"""Logging utilities for SpliceBind training and evaluation."""

import os
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional

import numpy as np


def setup_logger(name: str, log_dir: str, level: int = logging.INFO) -> logging.Logger:
    """Set up a logger with file and console handlers."""
    os.makedirs(log_dir, exist_ok=True)

    logger = logging.getLogger(name)
    logger.setLevel(level)

    if logger.handlers:
        return logger

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    fh = logging.FileHandler(os.path.join(log_dir, f"{name}_{timestamp}.log"))
    fh.setLevel(level)

    ch = logging.StreamHandler()
    ch.setLevel(level)

    formatter = logging.Formatter(
        "%(asctime)s | %(name)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger


class TrainingLogger:
    """Structured training logger that writes metrics to JSON lines + TensorBoard."""

    def __init__(self, log_dir: str, experiment_name: str, use_tensorboard: bool = True):
        self.log_dir = Path(log_dir) / experiment_name
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.metrics_file = self.log_dir / "metrics.jsonl"
        self.text_logger = setup_logger(experiment_name, str(self.log_dir))

        self.tb_writer = None
        if use_tensorboard:
            try:
                from torch.utils.tensorboard import SummaryWriter
                self.tb_writer = SummaryWriter(log_dir=str(self.log_dir / "tensorboard"))
            except ImportError:
                self.text_logger.warning("TensorBoard not available, skipping.")

        self.step = 0
        self.epoch = 0

    def log_metrics(self, metrics: Dict[str, float], step: Optional[int] = None,
                    phase: str = "train"):
        """Log metrics to JSONL file and TensorBoard."""
        if step is not None:
            self.step = step

        record = {
            "timestamp": datetime.now().isoformat(),
            "epoch": self.epoch,
            "step": self.step,
            "phase": phase,
            **metrics,
        }

        with open(self.metrics_file, "a") as f:
            f.write(json.dumps(record) + "\n")

        if self.tb_writer:
            for key, value in metrics.items():
                if isinstance(value, (int, float)):
                    self.tb_writer.add_scalar(f"{phase}/{key}", value, self.step)

        metric_str = " | ".join(f"{k}: {v:.4f}" if isinstance(v, float) else f"{k}: {v}"
                                for k, v in metrics.items())
        self.text_logger.info(f"[{phase}] epoch={self.epoch} step={self.step} | {metric_str}")

    def log_epoch(self, epoch: int, train_metrics: Dict, val_metrics: Dict):
        """Log full epoch summary."""
        self.epoch = epoch
        self.log_metrics(train_metrics, phase="train")
        self.log_metrics(val_metrics, phase="val")

    def log_text(self, tag: str, text: str):
        """Log text message."""
        self.text_logger.info(f"[{tag}] {text}")
        if self.tb_writer:
            self.tb_writer.add_text(tag, text, self.step)

    def log_hyperparams(self, hparams: Dict[str, Any]):
        """Log hyperparameters."""
        with open(self.log_dir / "hparams.json", "w") as f:
            json.dump(_serialize(hparams), f, indent=2)
        self.text_logger.info(f"Hyperparameters: {hparams}")

    def close(self):
        if self.tb_writer:
            self.tb_writer.close()


def _serialize(obj):
    """Make object JSON serializable."""
    if isinstance(obj, dict):
        return {k: _serialize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_serialize(v) for v in obj]
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, Path):
        return str(obj)
    return obj
