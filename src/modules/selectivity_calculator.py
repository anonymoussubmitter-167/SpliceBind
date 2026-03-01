"""Module 4: Selectivity Index Calculator

Integrates druggability predictions with tissue/cancer-specific
expression data to compute Selective Therapeutic Index (STI).
"""

import logging
import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


class SelectiveTherapeuticIndex:
    """Compute Selective Therapeutic Index (STI) for isoform targeting.

    STI = Druggability(cancer_isoform) × Expression(cancer_isoform, tumor)
          ─────────────────────────────────────────────────────────────────
          Σ [Druggability(isoform_i) × Expression(isoform_i, normal_tissue)]

    Higher STI = better therapeutic window for selective targeting.
    """

    def __init__(self, min_expression: float = 1.0, pseudocount: float = 0.1):
        self.min_expression = min_expression
        self.pseudocount = pseudocount

    def compute_sti(self, druggability_scores: Dict[str, float],
                    tumor_expression: Dict[str, float],
                    normal_expression: Dict[str, Dict[str, float]]) -> Dict:
        """Compute STI for a gene's isoforms in a specific cancer type.

        Args:
            druggability_scores: {isoform_id: druggability_score}
            tumor_expression: {isoform_id: TPM_in_tumor}
            normal_expression: {isoform_id: {tissue: TPM}}

        Returns:
            Dict with STI scores and components.
        """
        results = {}

        for iso_id in druggability_scores:
            drug_score = druggability_scores[iso_id]
            tumor_expr = tumor_expression.get(iso_id, 0)

            # Numerator: druggability × tumor expression
            numerator = drug_score * (tumor_expr + self.pseudocount)

            # Denominator: sum over all isoforms of (druggability × max normal expression)
            denominator = self.pseudocount
            for other_iso in druggability_scores:
                other_drug = druggability_scores[other_iso]
                if other_iso in normal_expression:
                    max_normal = max(normal_expression[other_iso].values()) if normal_expression[other_iso] else 0
                else:
                    max_normal = 0
                denominator += other_drug * (max_normal + self.pseudocount)

            sti = numerator / denominator if denominator > 0 else 0

            results[iso_id] = {
                "sti": float(sti),
                "druggability": float(drug_score),
                "tumor_expression": float(tumor_expr),
                "numerator": float(numerator),
                "denominator": float(denominator),
            }

        return results

    def compute_sti_matrix(self, gene_data: Dict, expression_data: pd.DataFrame,
                           cancer_types: List[str]) -> pd.DataFrame:
        """Compute STI matrix across genes and cancer types.

        Returns DataFrame with columns: gene, isoform, cancer_type, sti, ...
        """
        records = []

        for gene, gdata in gene_data.items():
            isoforms = list(gdata.get("druggability", {}).keys())
            if not isoforms:
                continue

            drug_scores = gdata["druggability"]

            for cancer in cancer_types:
                # Get tumor expression
                tumor_mask = (
                    (expression_data["gene"] == gene) &
                    (expression_data.get("cancer_type", expression_data.get("tissue", "")) == cancer) &
                    (expression_data.get("condition", "tumor") == "tumor")
                )
                tumor_expr = {}
                if tumor_mask.any():
                    for _, row in expression_data[tumor_mask].iterrows():
                        tumor_expr[row.get("isoform", "canonical")] = row.get("tpm", 0)

                # Get normal expression across tissues
                normal_expr = {}
                normal_mask = (
                    (expression_data["gene"] == gene) &
                    (expression_data.get("condition", "normal") == "normal")
                )
                if "tissue" in expression_data.columns:
                    for _, row in expression_data[normal_mask].iterrows():
                        iso = row.get("isoform", "canonical")
                        tissue = row.get("tissue", "unknown")
                        if iso not in normal_expr:
                            normal_expr[iso] = {}
                        normal_expr[iso][tissue] = row.get("tpm", 0)

                # Compute STI
                sti_results = self.compute_sti(drug_scores, tumor_expr, normal_expr)

                for iso_id, sti_data in sti_results.items():
                    records.append({
                        "gene": gene,
                        "isoform": iso_id,
                        "cancer_type": cancer,
                        **sti_data,
                    })

        return pd.DataFrame(records)

    def rank_targets(self, sti_df: pd.DataFrame, top_k: int = 50) -> pd.DataFrame:
        """Rank isoform targets by STI across cancer types.

        Returns top-K targets sorted by STI.
        """
        if sti_df.empty:
            return sti_df

        # Get best cancer type per isoform
        ranked = (sti_df
                  .sort_values("sti", ascending=False)
                  .drop_duplicates(subset=["gene", "isoform"], keep="first")
                  .head(top_k))

        return ranked


class DifferentialExpressionAnalyzer:
    """Analyze differential isoform expression between tumor and normal."""

    def __init__(self, min_fold_change: float = 1.5, min_expression: float = 1.0):
        self.min_fold_change = min_fold_change
        self.min_expression = min_expression

    def compute_differential_expression(self, expression_data: pd.DataFrame) -> pd.DataFrame:
        """Compute differential expression for each gene/isoform/cancer combination."""
        results = []

        if "condition" not in expression_data.columns:
            logger.warning("No 'condition' column in expression data, cannot compute DE")
            return pd.DataFrame()

        # Group by gene, isoform, cancer_type
        group_cols = ["gene", "isoform"]
        if "cancer_type" in expression_data.columns:
            group_cols.append("cancer_type")

        for name, group in expression_data.groupby(group_cols):
            if len(group_cols) == 3:
                gene, isoform, cancer = name
            else:
                gene, isoform = name
                cancer = "unknown"

            tumor_data = group[group["condition"] == "tumor"]
            normal_data = group[group["condition"] == "normal"]

            if tumor_data.empty or normal_data.empty:
                continue

            tumor_tpm = tumor_data["tpm"].mean()
            normal_tpm = normal_data["tpm"].mean()

            if normal_tpm < self.min_expression:
                fold_change = tumor_tpm / self.min_expression if tumor_tpm > 0 else 0
            else:
                fold_change = tumor_tpm / normal_tpm

            log2fc = np.log2(fold_change + 1e-10)

            results.append({
                "gene": gene,
                "isoform": isoform,
                "cancer_type": cancer,
                "tumor_tpm": tumor_tpm,
                "normal_tpm": normal_tpm,
                "fold_change": fold_change,
                "log2_fold_change": log2fc,
                "significant": abs(log2fc) > np.log2(self.min_fold_change),
                "direction": "up" if log2fc > 0 else "down",
            })

        return pd.DataFrame(results)

    def identify_cancer_specific_isoforms(self, de_results: pd.DataFrame,
                                          min_fc: float = 2.0) -> pd.DataFrame:
        """Find isoforms that are specifically upregulated in cancer."""
        if de_results.empty:
            return de_results

        cancer_specific = de_results[
            (de_results["fold_change"] >= min_fc) &
            (de_results["tumor_tpm"] >= self.min_expression) &
            (de_results["direction"] == "up")
        ].copy()

        cancer_specific["specificity_score"] = (
            cancer_specific["log2_fold_change"] *
            np.log2(cancer_specific["tumor_tpm"] + 1)
        )

        return cancer_specific.sort_values("specificity_score", ascending=False)
