"""Validation framework for SpliceBind.

Implements Tier 1 (hold-out kinase families), Tier 2 (literature case studies),
and Tier 3 (prospective predictions) validation.
"""

import logging
import json
import numpy as np
import pandas as pd
import torch
from typing import Dict, List, Optional, Tuple
from sklearn.metrics import (
    roc_auc_score, average_precision_score, precision_recall_curve,
    classification_report, confusion_matrix
)
from scipy.stats import spearmanr
from pathlib import Path

logger = logging.getLogger(__name__)


class Tier1Validator:
    """Hold-out kinase family validation.

    Tests generalization by evaluating on kinase families excluded from training.
    """

    def __init__(self, train_families: List[str], test_families: List[str]):
        self.train_families = train_families
        self.test_families = test_families

    def validate(self, predictions: List[Dict], labels: List[float],
                 metadata: List[Dict]) -> Dict:
        """Run Tier 1 validation.

        Args:
            predictions: list of {mean_prediction, std_prediction, ...}
            labels: true druggability labels
            metadata: list of {gene, family, ...}
        """
        results = {"overall": {}, "per_family": {}}

        # Overall metrics
        preds = [p["mean_prediction"] if isinstance(p, dict) else p for p in predictions]
        if isinstance(preds[0], (list, np.ndarray)):
            preds = [float(np.mean(p)) for p in preds]

        if len(set(labels)) > 1:
            results["overall"]["auroc"] = roc_auc_score(labels, preds)
            results["overall"]["auprc"] = average_precision_score(labels, preds)

        preds_binary = [1.0 if p > 0.5 else 0.0 for p in preds]
        results["overall"]["accuracy"] = sum(
            1 for p, l in zip(preds_binary, labels) if p == l
        ) / len(labels)
        results["overall"]["num_samples"] = len(labels)

        # Per-family metrics
        for family in set(m["family"] for m in metadata):
            family_idx = [i for i, m in enumerate(metadata) if m["family"] == family]
            if len(family_idx) < 2:
                continue

            f_preds = [preds[i] for i in family_idx]
            f_labels = [labels[i] for i in family_idx]

            f_results = {"num_samples": len(family_idx)}

            if len(set(f_labels)) > 1:
                try:
                    f_results["auroc"] = roc_auc_score(f_labels, f_preds)
                    f_results["auprc"] = average_precision_score(f_labels, f_preds)
                except ValueError:
                    f_results["auroc"] = None
                    f_results["auprc"] = None

            f_binary = [1.0 if p > 0.5 else 0.0 for p in f_preds]
            f_results["accuracy"] = sum(
                1 for p, l in zip(f_binary, f_labels) if p == l
            ) / len(f_labels)

            # Uncertainty analysis
            if isinstance(predictions[0], dict) and "std_prediction" in predictions[0]:
                uncertainties = [predictions[i].get("std_prediction", 0) for i in family_idx]
                if isinstance(uncertainties[0], (list, np.ndarray)):
                    uncertainties = [float(np.mean(u)) for u in uncertainties]
                f_results["mean_uncertainty"] = float(np.mean(uncertainties))

            results["per_family"][family] = f_results

        return results


class Tier2Validator:
    """Literature case study validation.

    Tests known splicing-drug interactions to validate model predictions.
    """

    # Known case studies from research plan
    CASE_STUDIES = [
        {
            "gene": "BCR-ABL",
            "splicing_event": "35-bp insertion in kinase domain",
            "drug": "Imatinib",
            "expected_effect": "resistance",
            "expected_delta_druggability": "negative",
            "reference": "Multiple studies",
            "description": "35INS variant inserts residues near ATP-binding pocket",
        },
        {
            "gene": "AR",
            "splicing_event": "Exon 4-8 skipping (AR-V7)",
            "drug": "Enzalutamide",
            "expected_effect": "resistance",
            "expected_delta_druggability": "negative",
            "reference": "Antonarakis et al. 2014",
            "description": "AR-V7 lacks ligand binding domain entirely",
        },
        {
            "gene": "EGFR",
            "splicing_event": "Exon 19 deletion",
            "drug": "Gefitinib",
            "expected_effect": "sensitivity",
            "expected_delta_druggability": "preserved",
            "reference": "Lynch et al. 2004",
            "description": "Exon 19 deletion preserves but alters ATP binding pocket",
        },
        {
            "gene": "PIK3CD",
            "splicing_event": "Exon 20 skipping (PIK3CD-S)",
            "drug": "Idelalisib",
            "expected_effect": "resistance",
            "expected_delta_druggability": "negative",
            "reference": "Durandy et al.",
            "description": "Short isoform lacks key binding residues",
        },
        {
            "gene": "BRAF",
            "splicing_event": "p61 splice variant (exon 4-8 skip)",
            "drug": "Vemurafenib",
            "expected_effect": "resistance",
            "expected_delta_druggability": "altered",
            "reference": "Poulikakos et al. 2011",
            "description": "Truncated BRAF promotes RAS-independent dimerization",
        },
        {
            "gene": "BCL2L1",
            "splicing_event": "BCL-XL to BCL-XS switching",
            "drug": "Navitoclax",
            "expected_effect": "changed_function",
            "expected_delta_druggability": "negative",
            "reference": "Boise et al. 1993",
            "description": "BCL-XS lacks BH1/BH2 domains, opposite function to BCL-XL",
        },
    ]

    def validate(self, predictions: Dict[str, Dict]) -> Dict:
        """Run Tier 2 validation against known case studies.

        Args:
            predictions: dict mapping gene -> {
                delta_druggability, pocket_preserved, confidence, ...
            }
        """
        results = []

        for case in self.CASE_STUDIES:
            gene = case["gene"]
            pred = predictions.get(gene, {})

            result = {
                "gene": gene,
                "drug": case["drug"],
                "expected_effect": case["expected_effect"],
                "expected_direction": case["expected_delta_druggability"],
                "has_prediction": bool(pred),
            }

            if pred:
                delta_drug = pred.get("delta_druggability", 0)
                result["predicted_delta_druggability"] = delta_drug
                result["predicted_pocket_preserved"] = pred.get("pocket_preserved", None)
                result["confidence"] = pred.get("confidence", None)

                # Check correctness
                expected = case["expected_delta_druggability"]
                if expected == "negative":
                    result["correct"] = delta_drug < 0
                elif expected == "preserved":
                    result["correct"] = abs(delta_drug) < 0.2
                elif expected == "altered":
                    result["correct"] = abs(delta_drug) > 0.1
                else:
                    result["correct"] = None
            else:
                result["correct"] = None
                result["note"] = "Gene not in prediction set"

            results.append(result)

        # Summary
        evaluated = [r for r in results if r["correct"] is not None]
        num_correct = sum(1 for r in evaluated if r["correct"])
        total = len(evaluated)

        summary = {
            "total_cases": len(self.CASE_STUDIES),
            "evaluated_cases": total,
            "correct_predictions": num_correct,
            "accuracy": num_correct / total if total > 0 else 0,
            "target_accuracy": 0.8,  # Success criterion: ≥80%
            "meets_criterion": (num_correct / total >= 0.8) if total > 0 else False,
            "per_case": results,
        }

        return summary


class Tier3Validator:
    """Prospective prediction validation.

    Generates and ranks novel predictions for experimental follow-up.
    """

    def generate_predictions(self, druggability_scores: Dict,
                              expression_data: pd.DataFrame,
                              cancer_types: List[str],
                              top_k: int = 50) -> pd.DataFrame:
        """Generate ranked prospective predictions.

        Args:
            druggability_scores: {gene: {isoform: score}}
            expression_data: DataFrame with cancer expression
            cancer_types: list of cancer types to evaluate

        Returns:
            DataFrame of top-K predictions ranked by therapeutic potential.
        """
        from ..modules.selectivity_calculator import SelectiveTherapeuticIndex

        sti_calculator = SelectiveTherapeuticIndex()

        predictions = []

        for gene, iso_scores in druggability_scores.items():
            for cancer in cancer_types:
                # Get expression for this gene/cancer
                gene_expr = expression_data[expression_data["gene"] == gene]
                if gene_expr.empty:
                    continue

                cancer_expr = gene_expr[
                    gene_expr.get("cancer_type", gene_expr.get("tissue", pd.Series())) == cancer
                ]

                tumor_expr = {}
                normal_expr = {}

                for _, row in cancer_expr.iterrows():
                    iso = row.get("isoform", "canonical")
                    if row.get("condition", "") == "tumor":
                        tumor_expr[iso] = row.get("tpm", 0)
                    elif row.get("condition", "") == "normal":
                        if iso not in normal_expr:
                            normal_expr[iso] = {}
                        normal_expr[iso]["primary"] = row.get("tpm", 0)

                sti_results = sti_calculator.compute_sti(iso_scores, tumor_expr, normal_expr)

                for iso_id, sti_data in sti_results.items():
                    predictions.append({
                        "gene": gene,
                        "isoform": iso_id,
                        "cancer_type": cancer,
                        **sti_data,
                    })

        df = pd.DataFrame(predictions)
        if not df.empty:
            df = df.sort_values("sti", ascending=False).head(top_k)

        return df


class ValidationReport:
    """Generate comprehensive validation report."""

    def __init__(self, output_dir: str):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def generate_report(self, tier1_results: Dict, tier2_results: Dict,
                        tier3_results: Optional[pd.DataFrame] = None,
                        training_history: Optional[Dict] = None) -> str:
        """Generate full validation report.

        Returns report as string and saves to file.
        """
        lines = []
        lines.append("=" * 70)
        lines.append("SpliceBind Validation Report")
        lines.append("=" * 70)
        lines.append("")

        # Tier 1
        lines.append("## Tier 1: Hold-out Kinase Family Validation")
        lines.append("-" * 50)
        overall = tier1_results.get("overall", {})
        lines.append(f"Overall AUROC: {overall.get('auroc', 'N/A')}")
        lines.append(f"Overall AUPRC: {overall.get('auprc', 'N/A')}")
        lines.append(f"Overall Accuracy: {overall.get('accuracy', 'N/A')}")
        lines.append(f"Total samples: {overall.get('num_samples', 'N/A')}")
        lines.append("")

        for family, metrics in tier1_results.get("per_family", {}).items():
            lines.append(f"  {family}: AUROC={metrics.get('auroc', 'N/A')}, "
                        f"Acc={metrics.get('accuracy', 'N/A'):.3f}, "
                        f"n={metrics.get('num_samples', '?')}")
        lines.append("")

        # Tier 2
        lines.append("## Tier 2: Literature Case Study Validation")
        lines.append("-" * 50)
        lines.append(f"Evaluated cases: {tier2_results.get('evaluated_cases', 0)}/{tier2_results.get('total_cases', 0)}")
        lines.append(f"Correct predictions: {tier2_results.get('correct_predictions', 0)}")
        lines.append(f"Accuracy: {tier2_results.get('accuracy', 0):.1%}")
        lines.append(f"Target: ≥80%")
        lines.append(f"Meets criterion: {tier2_results.get('meets_criterion', False)}")
        lines.append("")

        for case in tier2_results.get("per_case", []):
            status = "PASS" if case.get("correct") else ("FAIL" if case.get("correct") is False else "SKIP")
            lines.append(f"  [{status}] {case['gene']}/{case['drug']}: "
                        f"expected={case['expected_direction']}, "
                        f"predicted={case.get('predicted_delta_druggability', 'N/A')}")
        lines.append("")

        # Tier 3
        if tier3_results is not None and not tier3_results.empty:
            lines.append("## Tier 3: Prospective Predictions (Top 10)")
            lines.append("-" * 50)
            for _, row in tier3_results.head(10).iterrows():
                lines.append(f"  {row['gene']}/{row['isoform']} ({row.get('cancer_type', '?')}): "
                            f"STI={row.get('sti', 0):.4f}, "
                            f"Drug={row.get('druggability', 0):.3f}")
            lines.append("")

        # Training summary
        if training_history:
            lines.append("## Training Summary")
            lines.append("-" * 50)
            lines.append(f"Best val AUROC: {training_history.get('best_val_auroc', 'N/A')}")
            lines.append(f"Best val loss: {training_history.get('best_val_loss', 'N/A')}")
            lines.append(f"Total epochs: {training_history.get('total_epochs', 'N/A')}")
            if "test_metrics" in training_history:
                test = training_history["test_metrics"]
                lines.append(f"Test AUROC: {test.get('auroc', 'N/A')}")
                lines.append(f"Test AUPRC: {test.get('auprc', 'N/A')}")
                lines.append(f"Test Accuracy: {test.get('accuracy', 'N/A')}")

        report = "\n".join(lines)

        # Save
        report_path = self.output_dir / "validation_report.txt"
        with open(report_path, "w") as f:
            f.write(report)

        # Save JSON
        json_path = self.output_dir / "validation_results.json"
        with open(json_path, "w") as f:
            json.dump({
                "tier1": _make_serializable(tier1_results),
                "tier2": _make_serializable(tier2_results),
                "tier3": tier3_results.to_dict("records") if tier3_results is not None else [],
                "training": _make_serializable(training_history) if training_history else {},
            }, f, indent=2, default=str)

        logger.info(f"Validation report saved to {report_path}")
        return report


def _make_serializable(obj):
    """Make object JSON-serializable."""
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
    if isinstance(obj, pd.DataFrame):
        return obj.to_dict("records")
    return obj
