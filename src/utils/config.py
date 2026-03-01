"""Configuration management for SpliceBind."""

import os
from pathlib import Path
from omegaconf import OmegaConf, DictConfig


def load_config(config_path: str = None) -> DictConfig:
    """Load configuration from YAML file with defaults."""
    project_root = Path(__file__).parent.parent.parent
    default_path = project_root / "configs" / "default.yaml"

    cfg = OmegaConf.load(str(default_path))

    if config_path and os.path.exists(config_path):
        override = OmegaConf.load(config_path)
        cfg = OmegaConf.merge(cfg, override)

    # Resolve paths relative to project root
    for key in ["raw_dir", "processed_dir", "embeddings_dir", "structures_dir"]:
        path = cfg.data[key]
        if not os.path.isabs(path):
            cfg.data[key] = str(project_root / path)

    cfg.logging.log_dir = str(project_root / cfg.logging.log_dir)

    return cfg


def get_device(cfg: DictConfig):
    """Get torch device based on config."""
    import torch

    if cfg.project.device == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(cfg.project.device)
