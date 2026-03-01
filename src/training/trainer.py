"""Training pipeline for SpliceBind.

Multi-task training with comprehensive logging, early stopping,
and checkpointing.
"""

import os
import time
import logging
import json
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
from torch_geometric.loader import DataLoader as PyGDataLoader
from typing import Dict, List, Optional, Tuple
from pathlib import Path
from sklearn.metrics import roc_auc_score, average_precision_score

from ..modules.druggability_predictor import PocketGNN, MCDropoutPredictor
from ..utils.logging import TrainingLogger

logger = logging.getLogger(__name__)


class SpliceBindTrainer:
    """Trainer for SpliceBind druggability prediction model."""

    def __init__(self, model: PocketGNN, config: Dict, log_dir: str,
                 experiment_name: str = "splicebind"):
        self.model = model
        self.config = config
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.model = self.model.to(self.device)

        # Training config
        train_cfg = config.get("training", config.get("pretrain", {}))
        self.lr = train_cfg.get("lr", 1e-4)
        self.weight_decay = train_cfg.get("weight_decay", 1e-5)
        self.epochs = train_cfg.get("epochs", 100)
        self.patience = train_cfg.get("early_stopping_patience", 15)
        self.batch_size = train_cfg.get("batch_size", 64)

        # Optimizer
        self.optimizer = optim.AdamW(
            model.parameters(), lr=self.lr, weight_decay=self.weight_decay
        )

        # Scheduler
        self.scheduler = optim.lr_scheduler.CosineAnnealingLR(
            self.optimizer, T_max=self.epochs, eta_min=self.lr * 0.01
        )

        # Loss - use class-weighted BCE + focal loss for imbalanced data
        pos_weight = train_cfg.get("pos_weight", 10.0)  # upweight positives
        self.pos_weight = torch.tensor([pos_weight], dtype=torch.float32).to(self.device)
        self.criterion = nn.BCEWithLogitsLoss(pos_weight=self.pos_weight, reduction='none')
        self.focal_gamma = train_cfg.get("focal_gamma", 2.0)
        self.use_focal = train_cfg.get("use_focal", True)

        # Logging
        self.logger = TrainingLogger(log_dir, experiment_name)
        self.logger.log_hyperparams({
            "model": str(type(model).__name__),
            "lr": self.lr,
            "weight_decay": self.weight_decay,
            "epochs": self.epochs,
            "batch_size": self.batch_size,
            "patience": self.patience,
            "pos_weight": pos_weight,
            "focal_gamma": self.focal_gamma,
            "device": str(self.device),
        })

        # Checkpointing
        self.checkpoint_dir = os.path.join(log_dir, experiment_name, "checkpoints")
        os.makedirs(self.checkpoint_dir, exist_ok=True)

        # Early stopping
        self.best_val_loss = float("inf")
        self.best_val_auroc = 0.0
        self.patience_counter = 0

    def train(self, train_dataset, val_dataset, test_dataset=None) -> Dict:
        """Run full training loop.

        Returns dict with final metrics.
        """
        train_loader = PyGDataLoader(train_dataset, batch_size=self.batch_size,
                                      shuffle=True, drop_last=False)
        val_loader = PyGDataLoader(val_dataset, batch_size=self.batch_size,
                                    shuffle=False)

        self.logger.log_text("training", f"Starting training: {self.epochs} epochs, "
                            f"{len(train_dataset)} train / {len(val_dataset)} val samples")

        all_train_metrics = []
        all_val_metrics = []

        for epoch in range(1, self.epochs + 1):
            self.logger.epoch = epoch

            # Train epoch
            train_metrics = self._train_epoch(train_loader, epoch)
            all_train_metrics.append(train_metrics)

            # Validation
            val_metrics = self._validate(val_loader, epoch)
            all_val_metrics.append(val_metrics)

            # Log
            self.logger.log_epoch(epoch, train_metrics, val_metrics)

            # Learning rate
            current_lr = self.optimizer.param_groups[0]["lr"]
            self.logger.log_metrics({"learning_rate": current_lr}, phase="lr")
            self.scheduler.step()

            # Checkpointing
            if val_metrics.get("auroc", 0) > self.best_val_auroc:
                self.best_val_auroc = val_metrics["auroc"]
                self.best_val_loss = val_metrics["loss"]
                self.patience_counter = 0
                self._save_checkpoint(epoch, val_metrics, is_best=True)
                self.logger.log_text("checkpoint",
                    f"New best model at epoch {epoch}: AUROC={val_metrics['auroc']:.4f}")
            else:
                self.patience_counter += 1
                if epoch % 10 == 0:
                    self._save_checkpoint(epoch, val_metrics, is_best=False)

            # Early stopping
            if self.patience_counter >= self.patience:
                self.logger.log_text("early_stop",
                    f"Early stopping at epoch {epoch} (patience={self.patience})")
                break

        # Final evaluation
        self.logger.log_text("training", "Training complete. Loading best model...")
        self._load_best_checkpoint()

        final_results = {
            "best_val_auroc": self.best_val_auroc,
            "best_val_loss": self.best_val_loss,
            "total_epochs": epoch,
            "train_history": all_train_metrics,
            "val_history": all_val_metrics,
        }

        if test_dataset:
            test_loader = PyGDataLoader(test_dataset, batch_size=self.batch_size, shuffle=False)
            test_metrics = self._validate(test_loader, epoch, phase="test")
            final_results["test_metrics"] = test_metrics
            self.logger.log_metrics(test_metrics, phase="test")
            self.logger.log_text("test", f"Test AUROC: {test_metrics.get('auroc', 'N/A')}")

        # Save final results
        results_path = os.path.join(self.checkpoint_dir, "..", "final_results.json")
        with open(results_path, "w") as f:
            json.dump(_make_serializable(final_results), f, indent=2)

        self.logger.close()
        return final_results

    def _train_epoch(self, loader: DataLoader, epoch: int) -> Dict:
        """Run one training epoch."""
        self.model.train()
        total_loss = 0
        all_preds = []
        all_labels = []
        num_batches = 0

        for batch_idx, batch in enumerate(loader):
            batch = batch.to(self.device)
            self.optimizer.zero_grad()

            output = self.model(batch)
            logits = output["druggability_logits"].squeeze()
            pred = output["druggability_score"].squeeze()
            labels = batch.y.squeeze()

            # Ensure same shape
            if logits.dim() == 0:
                logits = logits.unsqueeze(0)
            if pred.dim() == 0:
                pred = pred.unsqueeze(0)
            if labels.dim() == 0:
                labels = labels.unsqueeze(0)

            # Weighted BCE loss (per-sample)
            bce_loss = self.criterion(logits, labels)

            # Apply focal weighting: (1-p_t)^gamma
            if self.use_focal:
                p_t = pred * labels + (1 - pred) * (1 - labels)
                focal_weight = (1 - p_t) ** self.focal_gamma
                loss = (focal_weight * bce_loss).mean()
            else:
                loss = bce_loss.mean()

            loss.backward()

            # Gradient clipping
            torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)

            self.optimizer.step()

            total_loss += loss.item()
            all_preds.extend(pred.detach().cpu().numpy().tolist())
            all_labels.extend(labels.detach().cpu().numpy().tolist())
            num_batches += 1

            # Step-level logging
            step = (epoch - 1) * len(loader) + batch_idx
            if batch_idx % 10 == 0:
                self.logger.log_metrics(
                    {"batch_loss": loss.item()},
                    step=step,
                    phase="train_step",
                )

        avg_loss = total_loss / max(num_batches, 1)
        metrics = {"loss": avg_loss}

        # Compute metrics
        if len(set(all_labels)) > 1:
            try:
                metrics["auroc"] = roc_auc_score(all_labels, all_preds)
                metrics["auprc"] = average_precision_score(all_labels, all_preds)
            except ValueError:
                metrics["auroc"] = 0.0
                metrics["auprc"] = 0.0
        else:
            metrics["auroc"] = 0.0
            metrics["auprc"] = 0.0

        # Accuracy
        preds_binary = [1.0 if p > 0.5 else 0.0 for p in all_preds]
        metrics["accuracy"] = sum(1 for p, l in zip(preds_binary, all_labels) if p == l) / max(len(all_labels), 1)

        return metrics

    @torch.no_grad()
    def _validate(self, loader: DataLoader, epoch: int, phase: str = "val") -> Dict:
        """Run validation/test evaluation."""
        self.model.eval()
        total_loss = 0
        all_preds = []
        all_labels = []
        all_embeddings = []
        num_batches = 0

        for batch in loader:
            batch = batch.to(self.device)
            output = self.model(batch)
            logits = output["druggability_logits"].squeeze()
            pred = output["druggability_score"].squeeze()
            labels = batch.y.squeeze()

            if logits.dim() == 0:
                logits = logits.unsqueeze(0)
            if pred.dim() == 0:
                pred = pred.unsqueeze(0)
            if labels.dim() == 0:
                labels = labels.unsqueeze(0)

            bce_loss = self.criterion(logits, labels)
            loss = bce_loss.mean()
            total_loss += loss.item()

            all_preds.extend(pred.cpu().numpy().tolist())
            all_labels.extend(labels.cpu().numpy().tolist())
            all_embeddings.append(output["pocket_embedding"].cpu().numpy())
            num_batches += 1

        avg_loss = total_loss / max(num_batches, 1)
        metrics = {"loss": avg_loss}

        if len(set(all_labels)) > 1:
            try:
                metrics["auroc"] = roc_auc_score(all_labels, all_preds)
                metrics["auprc"] = average_precision_score(all_labels, all_preds)
            except ValueError:
                metrics["auroc"] = 0.0
                metrics["auprc"] = 0.0
        else:
            metrics["auroc"] = 0.0
            metrics["auprc"] = 0.0

        preds_binary = [1.0 if p > 0.5 else 0.0 for p in all_preds]
        metrics["accuracy"] = sum(1 for p, l in zip(preds_binary, all_labels) if p == l) / max(len(all_labels), 1)
        metrics["num_samples"] = len(all_labels)

        return metrics

    def _save_checkpoint(self, epoch: int, metrics: Dict, is_best: bool = False):
        """Save model checkpoint."""
        checkpoint = {
            "epoch": epoch,
            "model_state_dict": self.model.state_dict(),
            "optimizer_state_dict": self.optimizer.state_dict(),
            "scheduler_state_dict": self.scheduler.state_dict(),
            "metrics": metrics,
            "best_val_auroc": self.best_val_auroc,
        }

        if is_best:
            path = os.path.join(self.checkpoint_dir, "best_model.pt")
        else:
            path = os.path.join(self.checkpoint_dir, f"checkpoint_epoch_{epoch}.pt")

        torch.save(checkpoint, path)

    def _load_best_checkpoint(self):
        """Load best model checkpoint."""
        path = os.path.join(self.checkpoint_dir, "best_model.pt")
        if os.path.exists(path):
            checkpoint = torch.load(path, map_location=self.device, weights_only=False)
            self.model.load_state_dict(checkpoint["model_state_dict"])
            logger.info(f"Loaded best model from epoch {checkpoint['epoch']}")

    def evaluate_with_uncertainty(self, dataset, num_mc_samples: int = 20) -> Dict:
        """Evaluate model with MC Dropout uncertainty estimation."""
        mc_predictor = MCDropoutPredictor(self.model, num_samples=num_mc_samples)
        loader = PyGDataLoader(dataset, batch_size=1, shuffle=False)

        results = []
        for data in loader:
            data = data.to(self.device)
            uc_result = mc_predictor.predict_with_uncertainty(data)
            uc_result["true_label"] = data.y.item()
            results.append(uc_result)

        # Aggregate
        predictions = [r["mean_prediction"] for r in results]
        uncertainties = [r["std_prediction"] for r in results]
        labels = [r["true_label"] for r in results]

        if isinstance(predictions[0], np.ndarray):
            predictions = [float(p.mean()) for p in predictions]
            uncertainties = [float(u.mean()) for u in uncertainties]

        metrics = {
            "mean_uncertainty": float(np.mean(uncertainties)),
            "std_uncertainty": float(np.std(uncertainties)),
        }

        if len(set(labels)) > 1:
            try:
                metrics["auroc"] = roc_auc_score(labels, predictions)
                metrics["auprc"] = average_precision_score(labels, predictions)
            except ValueError:
                pass

        # Calibration: higher uncertainty for incorrect predictions?
        correct = [1 if (p > 0.5) == (l > 0.5) else 0 for p, l in zip(predictions, labels)]
        if any(c == 0 for c in correct) and any(c == 1 for c in correct):
            correct_unc = [u for u, c in zip(uncertainties, correct) if c == 1]
            incorrect_unc = [u for u, c in zip(uncertainties, correct) if c == 0]
            metrics["mean_uncertainty_correct"] = float(np.mean(correct_unc))
            metrics["mean_uncertainty_incorrect"] = float(np.mean(incorrect_unc))
            metrics["uncertainty_calibrated"] = (
                np.mean(incorrect_unc) > np.mean(correct_unc)
            )

        return {
            "metrics": metrics,
            "per_sample": results,
        }


def _make_serializable(obj):
    """Make object JSON serializable."""
    if isinstance(obj, dict):
        return {k: _make_serializable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_make_serializable(v) for v in obj]
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (bool, np.bool_)):
        return bool(obj)
    return obj
