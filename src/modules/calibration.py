"""Uncertainty calibration for SpliceBind predictions.

Implements temperature scaling and deep ensemble aggregation
to improve prediction calibration.
"""

import logging
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from typing import Dict, List, Optional, Tuple
from torch_geometric.loader import DataLoader as PyGDataLoader
from sklearn.metrics import roc_auc_score, log_loss

logger = logging.getLogger(__name__)


class TemperatureScaling(nn.Module):
    """Post-hoc temperature scaling for calibration.

    Learns a single temperature parameter T such that
    calibrated_logits = logits / T, minimizing NLL on a validation set.
    """

    def __init__(self):
        super().__init__()
        self.temperature = nn.Parameter(torch.ones(1) * 1.5)

    def forward(self, logits: torch.Tensor) -> torch.Tensor:
        return logits / self.temperature

    def calibrated_probability(self, logits: torch.Tensor) -> torch.Tensor:
        return torch.sigmoid(self.forward(logits))


def fit_temperature_scaling(model, val_dataset, device: str = "cpu",
                             max_iter: int = 100, lr: float = 0.01) -> TemperatureScaling:
    """Fit temperature scaling on validation set.

    Args:
        model: Trained PocketGNN model
        val_dataset: Validation dataset
        device: torch device
        max_iter: Max LBFGS iterations
        lr: Learning rate

    Returns:
        Fitted TemperatureScaling module
    """
    model.eval()
    loader = PyGDataLoader(val_dataset, batch_size=len(val_dataset), shuffle=False)

    # Collect all logits and labels
    all_logits = []
    all_labels = []

    with torch.no_grad():
        for batch in loader:
            batch = batch.to(device)
            output = model(batch)
            logits = output["druggability_logits"].squeeze()
            labels = batch.y.squeeze()

            if logits.dim() == 0:
                logits = logits.unsqueeze(0)
            if labels.dim() == 0:
                labels = labels.unsqueeze(0)

            all_logits.append(logits.cpu())
            all_labels.append(labels.cpu())

    logits = torch.cat(all_logits)
    labels = torch.cat(all_labels)

    # Pre-calibration NLL
    pre_probs = torch.sigmoid(logits).numpy()
    pre_nll = _safe_log_loss(labels.numpy(), pre_probs)
    logger.info(f"Pre-calibration NLL: {pre_nll:.4f}")

    # Fit temperature
    temp_model = TemperatureScaling()
    criterion = nn.BCEWithLogitsLoss()

    optimizer = optim.LBFGS([temp_model.temperature], lr=lr, max_iter=max_iter)

    def closure():
        optimizer.zero_grad()
        scaled_logits = temp_model(logits)
        loss = criterion(scaled_logits, labels)
        loss.backward()
        return loss

    optimizer.step(closure)

    # Post-calibration NLL
    with torch.no_grad():
        post_probs = temp_model.calibrated_probability(logits).numpy()
    post_nll = _safe_log_loss(labels.numpy(), post_probs)

    logger.info(f"Post-calibration NLL: {post_nll:.4f}")
    logger.info(f"Learned temperature: {temp_model.temperature.item():.4f}")
    logger.info(f"NLL improvement: {pre_nll - post_nll:.4f}")

    return temp_model


def _safe_log_loss(y_true, y_pred, eps=1e-7):
    """Compute log loss safely."""
    y_pred = np.clip(y_pred, eps, 1 - eps)
    try:
        return log_loss(y_true, y_pred)
    except Exception:
        return float('inf')


class DeepEnsemble:
    """Deep ensemble for improved uncertainty estimation.

    Trains N models with different random seeds and aggregates predictions.
    """

    def __init__(self, models: List[nn.Module]):
        self.models = models

    @torch.no_grad()
    def predict(self, data, device: str = "cpu") -> Dict:
        """Predict with ensemble, returning mean and uncertainty."""
        all_probs = []

        data = data.to(device)

        for model in self.models:
            model.eval()
            model = model.to(device)
            output = model(data)
            prob = output["druggability_score"].squeeze().cpu().numpy()
            all_probs.append(prob)

        all_probs = np.array(all_probs)

        mean_pred = all_probs.mean(axis=0)
        std_pred = all_probs.std(axis=0)

        # Predictive entropy
        mean_clipped = np.clip(mean_pred, 1e-7, 1 - 1e-7)
        entropy = -(mean_clipped * np.log(mean_clipped) +
                   (1 - mean_clipped) * np.log(1 - mean_clipped))

        # Mutual information (epistemic uncertainty)
        individual_entropies = []
        for p in all_probs:
            p_clipped = np.clip(p, 1e-7, 1 - 1e-7)
            h = -(p_clipped * np.log(p_clipped) + (1 - p_clipped) * np.log(1 - p_clipped))
            individual_entropies.append(h)
        mean_entropy = np.mean(individual_entropies, axis=0)
        mutual_info = entropy - mean_entropy

        return {
            "mean_prediction": float(mean_pred) if mean_pred.ndim == 0 else mean_pred,
            "std_prediction": float(std_pred) if std_pred.ndim == 0 else std_pred,
            "predictive_entropy": float(entropy) if np.ndim(entropy) == 0 else entropy,
            "epistemic_uncertainty": float(mutual_info) if np.ndim(mutual_info) == 0 else mutual_info,
            "all_predictions": all_probs.tolist(),
        }


def evaluate_calibration(predictions: List[float], labels: List[float],
                          n_bins: int = 10) -> Dict:
    """Evaluate calibration using Expected Calibration Error (ECE).

    Args:
        predictions: Predicted probabilities
        labels: True binary labels
        n_bins: Number of bins for calibration

    Returns:
        Dict with calibration metrics
    """
    predictions = np.array(predictions)
    labels = np.array(labels)

    bin_boundaries = np.linspace(0, 1, n_bins + 1)
    ece = 0.0
    bin_stats = []

    for i in range(n_bins):
        low = bin_boundaries[i]
        high = bin_boundaries[i + 1]
        mask = (predictions >= low) & (predictions < high)

        if mask.sum() == 0:
            continue

        bin_preds = predictions[mask]
        bin_labels = labels[mask]

        avg_confidence = bin_preds.mean()
        avg_accuracy = bin_labels.mean()
        bin_size = mask.sum()

        ece += (bin_size / len(predictions)) * abs(avg_accuracy - avg_confidence)

        bin_stats.append({
            "bin_low": float(low),
            "bin_high": float(high),
            "avg_confidence": float(avg_confidence),
            "avg_accuracy": float(avg_accuracy),
            "count": int(bin_size),
        })

    # AUROC
    auroc = 0.0
    if len(set(labels)) > 1:
        try:
            auroc = roc_auc_score(labels, predictions)
        except ValueError:
            pass

    # Brier score
    brier = float(np.mean((predictions - labels) ** 2))

    return {
        "ece": float(ece),
        "brier_score": brier,
        "auroc": auroc,
        "n_bins": n_bins,
        "bin_stats": bin_stats,
    }
