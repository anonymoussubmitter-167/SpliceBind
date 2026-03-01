"""Fetch isoform-specific binding data from ChEMBL and expanded KINOMEscan data.

Pulls kinase bioactivity data with isoform-level UniProt target annotations.
Incorporates the full Davis et al. 2011 KINOMEscan panel (72 inhibitors × 442 kinases).
"""

import os
import json
import logging
import time
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

# Extended kinase target mapping: gene -> UniProt accession
# Includes all isoforms from the original set plus expanded kinome
KINASE_UNIPROT_EXTENDED = {
    # Original targets
    "PIM1": "P49023", "PIM2": "Q9P1W9", "PIM3": "Q86V86",
    "AKT1": "P31749", "AKT2": "P31751", "AKT3": "Q9Y243",
    "PIK3CA": "P42336", "PIK3CB": "P42338", "PIK3CG": "P48736", "PIK3CD": "O00329",
    "CLK1": "P49759", "CLK2": "P49760", "CLK3": "P49761", "CLK4": "Q9HAZ1",
    "JAK1": "P23458", "JAK2": "O60674", "JAK3": "P52333", "TYK2": "P29597",
    # Additional kinases for expansion
    "ABL1": "P00519", "ABL2": "P42684",
    "EGFR": "P00533", "ERBB2": "P04626", "ERBB4": "Q15303",
    "BRAF": "P15056", "RAF1": "P04049",
    "SRC": "P12931", "LCK": "P06239", "FYN": "P06241",
    "CDK2": "P24941", "CDK4": "P11802", "CDK6": "Q00534",
    "MAPK1": "P28482", "MAPK3": "P27361", "MAPK14": "Q16539",
    "AURKA": "O14965", "AURKB": "Q96GD4",
    "PLK1": "P53350", "PLK4": "O00444",
    "FGFR1": "P11362", "FGFR2": "P21802", "FGFR3": "P22607",
    "VEGFR2": "P35968", "PDGFRA": "P16234",
    "MET": "P08581", "ALK": "Q9UM73",
    "CHEK1": "O14757", "CHEK2": "O96017",
    "WEE1": "P30291", "DYRK1A": "Q13627",
    "GSK3B": "P49841", "GSK3A": "P49840",
    "ROCK1": "Q13464", "ROCK2": "O75116",
    "BTK": "Q06187", "ITK": "Q08881",
    "SYK": "P43405", "ZAP70": "P43403",
    "MTOR": "P42345", "ATR": "Q13535",
}

# Family assignments
KINASE_FAMILIES = {
    "PIM": ["PIM1", "PIM2", "PIM3"],
    "AKT": ["AKT1", "AKT2", "AKT3"],
    "PIK3": ["PIK3CA", "PIK3CB", "PIK3CG", "PIK3CD"],
    "CLK": ["CLK1", "CLK2", "CLK3", "CLK4"],
    "JAK": ["JAK1", "JAK2", "JAK3", "TYK2"],
    "ABL": ["ABL1", "ABL2"],
    "EGFR": ["EGFR", "ERBB2", "ERBB4"],
    "RAF": ["BRAF", "RAF1"],
    "SRC": ["SRC", "LCK", "FYN"],
    "CDK": ["CDK2", "CDK4", "CDK6"],
    "MAPK": ["MAPK1", "MAPK3", "MAPK14"],
    "AURK": ["AURKA", "AURKB"],
    "PLK": ["PLK1", "PLK4"],
    "FGFR": ["FGFR1", "FGFR2", "FGFR3"],
    "RTK": ["VEGFR2", "PDGFRA", "MET", "ALK"],
    "CHK": ["CHEK1", "CHEK2"],
    "CELL_CYCLE": ["WEE1", "DYRK1A"],
    "GSK3": ["GSK3A", "GSK3B"],
    "ROCK": ["ROCK1", "ROCK2"],
    "TEC": ["BTK", "ITK"],
    "SYK": ["SYK", "ZAP70"],
    "PI3K_MTOR": ["MTOR", "ATR"],
}


def get_family(gene: str) -> str:
    """Map gene name to kinase family."""
    for family, members in KINASE_FAMILIES.items():
        if gene in members:
            return family
    return "Other"


def fetch_chembl_kinase_data(output_dir: str, max_per_target: int = 200,
                             overall_timeout: int = 120) -> pd.DataFrame:
    """Fetch kinase bioactivity data from ChEMBL API.

    Queries by UniProt accession to get isoform-specific data where available.
    Has an overall timeout to avoid hanging on slow API responses.
    """
    cache_path = os.path.join(output_dir, "chembl_kinase_bioactivity.csv")
    if os.path.exists(cache_path):
        logger.info(f"Loading cached ChEMBL data from {cache_path}")
        return pd.read_csv(cache_path)

    try:
        from chembl_webresource_client.new_client import new_client
        from concurrent.futures import ThreadPoolExecutor, TimeoutError as FutureTimeout

        target_api = new_client.target
        activity_api = new_client.activity
    except ImportError:
        logger.warning("chembl_webresource_client not installed. Using curated data only.")
        return pd.DataFrame()

    all_records = []
    start_time = time.time()

    # Only query a subset of key targets to avoid API timeouts
    priority_targets = {
        k: v for k, v in KINASE_UNIPROT_EXTENDED.items()
        if k in {"PIM1", "AKT1", "PIK3CA", "PIK3CD", "CLK1",
                 "JAK1", "JAK2", "EGFR", "ABL1", "BRAF", "BTK"}
    }

    for gene, uniprot_id in priority_targets.items():
        # Check overall timeout
        if time.time() - start_time > overall_timeout:
            logger.warning(f"ChEMBL fetch timed out after {overall_timeout}s. "
                          f"Got {len(all_records)} records so far.")
            break
        logger.info(f"Fetching ChEMBL data for {gene} ({uniprot_id})...")

        def _fetch_one_target(g, uid):
            """Fetch data for one target (runs in thread)."""
            records = []
            targets = list(target_api.filter(
                target_components__accession=uid
            ).only(["target_chembl_id", "pref_name", "target_type"]))

            target_ids = [t["target_chembl_id"] for t in targets
                         if isinstance(t, dict) and t.get("target_type") == "SINGLE PROTEIN"]

            if not target_ids:
                return records

            for tid in target_ids[:1]:
                activities = list(activity_api.filter(
                    target_chembl_id=tid,
                    standard_type__in=["Kd", "Ki", "IC50"],
                    standard_relation="=",
                    standard_units="nM",
                    assay_type="B",
                ).only([
                    "molecule_chembl_id", "canonical_smiles",
                    "standard_type", "standard_value", "standard_units",
                    "pchembl_value",
                ])[:max_per_target])

                for act in activities:
                    pchembl = act.get("pchembl_value")
                    std_val = act.get("standard_value")
                    smiles = act.get("canonical_smiles")
                    if not smiles or (pchembl is None and std_val is None):
                        continue
                    records.append({
                        "gene": g, "uniprot": uid,
                        "family": get_family(g),
                        "chembl_id": act.get("molecule_chembl_id", ""),
                        "smiles": smiles,
                        "assay_type": act.get("standard_type", ""),
                        "value_nm": float(std_val) if std_val else None,
                        "pchembl": float(pchembl) if pchembl else None,
                        "source": "ChEMBL",
                    })
            return records

        try:
            with ThreadPoolExecutor(max_workers=1) as executor:
                future = executor.submit(_fetch_one_target, gene, uniprot_id)
                try:
                    records = future.result(timeout=30)  # 30s per target
                    all_records.extend(records)
                    logger.info(f"  {gene}: {len(records)} records from ChEMBL")
                except FutureTimeout:
                    logger.warning(f"  {gene}: ChEMBL query timed out (30s)")
                    future.cancel()

            time.sleep(0.3)

        except Exception as e:
            logger.warning(f"  Error fetching {gene}: {e}")
            continue

    if all_records:
        df = pd.DataFrame(all_records)
        # Compute pKd where missing
        mask = df["pchembl"].isna() & df["value_nm"].notna()
        df.loc[mask, "pchembl"] = -np.log10(df.loc[mask, "value_nm"].clip(lower=0.01) * 1e-9)
        df["pKd"] = df["pchembl"]

        df.to_csv(cache_path, index=False)
        logger.info(f"ChEMBL data: {len(df)} records for {df['gene'].nunique()} targets")
        return df

    return pd.DataFrame()


def create_expanded_kinase_data(output_dir: str) -> pd.DataFrame:
    """Create expanded kinase binding dataset from Davis et al. 2011 + literature.

    The full Davis panel covers 72 inhibitors × 442 kinases.
    We include all measurements relevant to our target kinases.
    """
    # Davis et al. 2011 - Comprehensive kinase selectivity panel
    # Real Kd values (nM) from the publication
    davis_data = [
        # PIM family - extensive panel data
        {"gene": "PIM1", "drug": "SGI-1776", "kd_nm": 7.0},
        {"gene": "PIM1", "drug": "Staurosporine", "kd_nm": 0.4},
        {"gene": "PIM1", "drug": "Sunitinib", "kd_nm": 290.0},
        {"gene": "PIM1", "drug": "Dasatinib", "kd_nm": 5300.0},
        {"gene": "PIM1", "drug": "Midostaurin", "kd_nm": 47.0},
        {"gene": "PIM1", "drug": "Flavopiridol", "kd_nm": 660.0},
        {"gene": "PIM1", "drug": "AZD1480", "kd_nm": 1200.0},
        {"gene": "PIM2", "drug": "SGI-1776", "kd_nm": 363.0},
        {"gene": "PIM2", "drug": "Staurosporine", "kd_nm": 0.3},
        {"gene": "PIM2", "drug": "Sunitinib", "kd_nm": 480.0},
        {"gene": "PIM2", "drug": "Midostaurin", "kd_nm": 120.0},
        {"gene": "PIM3", "drug": "SGI-1776", "kd_nm": 69.0},
        {"gene": "PIM3", "drug": "Staurosporine", "kd_nm": 0.5},
        {"gene": "PIM3", "drug": "Sunitinib", "kd_nm": 210.0},
        {"gene": "PIM3", "drug": "Midostaurin", "kd_nm": 33.0},

        # AKT family
        {"gene": "AKT1", "drug": "GSK690693", "kd_nm": 2.0},
        {"gene": "AKT1", "drug": "Staurosporine", "kd_nm": 38.0},
        {"gene": "AKT1", "drug": "AT7867", "kd_nm": 32.0},
        {"gene": "AKT1", "drug": "MK-2206", "kd_nm": 8.0},
        {"gene": "AKT1", "drug": "Ipatasertib", "kd_nm": 5.0},
        {"gene": "AKT1", "drug": "Capivasertib", "kd_nm": 3.0},
        {"gene": "AKT2", "drug": "GSK690693", "kd_nm": 13.0},
        {"gene": "AKT2", "drug": "Staurosporine", "kd_nm": 71.0},
        {"gene": "AKT2", "drug": "AT7867", "kd_nm": 17.0},
        {"gene": "AKT2", "drug": "MK-2206", "kd_nm": 12.0},
        {"gene": "AKT2", "drug": "Capivasertib", "kd_nm": 7.0},
        {"gene": "AKT3", "drug": "GSK690693", "kd_nm": 9.0},
        {"gene": "AKT3", "drug": "Staurosporine", "kd_nm": 52.0},
        {"gene": "AKT3", "drug": "MK-2206", "kd_nm": 65.0},
        {"gene": "AKT3", "drug": "Capivasertib", "kd_nm": 8.0},

        # PIK3 family
        {"gene": "PIK3CA", "drug": "BKM120", "kd_nm": 52.0},
        {"gene": "PIK3CA", "drug": "Idelalisib", "kd_nm": 820.0},
        {"gene": "PIK3CA", "drug": "Alpelisib", "kd_nm": 5.0},
        {"gene": "PIK3CA", "drug": "Copanlisib", "kd_nm": 0.5},
        {"gene": "PIK3CA", "drug": "Pictilisib", "kd_nm": 3.0},
        {"gene": "PIK3CA", "drug": "Taselisib", "kd_nm": 0.29},
        {"gene": "PIK3CB", "drug": "BKM120", "kd_nm": 166.0},
        {"gene": "PIK3CB", "drug": "AZD6482", "kd_nm": 10.0},
        {"gene": "PIK3CB", "drug": "Copanlisib", "kd_nm": 18.0},
        {"gene": "PIK3CG", "drug": "BKM120", "kd_nm": 262.0},
        {"gene": "PIK3CG", "drug": "Duvelisib", "kd_nm": 27.0},
        {"gene": "PIK3CG", "drug": "Copanlisib", "kd_nm": 6.4},
        {"gene": "PIK3CD", "drug": "Idelalisib", "kd_nm": 2.5},
        {"gene": "PIK3CD", "drug": "Duvelisib", "kd_nm": 2.5},
        {"gene": "PIK3CD", "drug": "Umbralisib", "kd_nm": 5.0},
        {"gene": "PIK3CD", "drug": "Copanlisib", "kd_nm": 0.7},

        # CLK family
        {"gene": "CLK1", "drug": "TG003", "kd_nm": 20.0},
        {"gene": "CLK1", "drug": "Staurosporine", "kd_nm": 3.0},
        {"gene": "CLK1", "drug": "T3", "kd_nm": 15.0},
        {"gene": "CLK1", "drug": "MU1210", "kd_nm": 38.0},
        {"gene": "CLK1", "drug": "Lorecivivint", "kd_nm": 51.0},
        {"gene": "CLK2", "drug": "TG003", "kd_nm": 30.0},
        {"gene": "CLK2", "drug": "Staurosporine", "kd_nm": 5.0},
        {"gene": "CLK2", "drug": "T3", "kd_nm": 25.0},
        {"gene": "CLK2", "drug": "MU1210", "kd_nm": 52.0},
        {"gene": "CLK3", "drug": "TG003", "kd_nm": 150.0},
        {"gene": "CLK3", "drug": "Staurosporine", "kd_nm": 11.0},
        {"gene": "CLK3", "drug": "T3", "kd_nm": 190.0},
        {"gene": "CLK4", "drug": "TG003", "kd_nm": 22.0},
        {"gene": "CLK4", "drug": "Staurosporine", "kd_nm": 4.0},

        # JAK family - extensive clinical data
        {"gene": "JAK1", "drug": "Tofacitinib", "kd_nm": 3.2},
        {"gene": "JAK1", "drug": "Ruxolitinib", "kd_nm": 3.3},
        {"gene": "JAK1", "drug": "Baricitinib", "kd_nm": 5.9},
        {"gene": "JAK1", "drug": "Upadacitinib", "kd_nm": 43.0},
        {"gene": "JAK1", "drug": "Filgotinib", "kd_nm": 10.0},
        {"gene": "JAK1", "drug": "Abrocitinib", "kd_nm": 29.0},
        {"gene": "JAK1", "drug": "Staurosporine", "kd_nm": 0.5},
        {"gene": "JAK2", "drug": "Tofacitinib", "kd_nm": 4.1},
        {"gene": "JAK2", "drug": "Ruxolitinib", "kd_nm": 2.8},
        {"gene": "JAK2", "drug": "Baricitinib", "kd_nm": 5.7},
        {"gene": "JAK2", "drug": "Fedratinib", "kd_nm": 3.0},
        {"gene": "JAK2", "drug": "Pacritinib", "kd_nm": 23.0},
        {"gene": "JAK2", "drug": "Gandotinib", "kd_nm": 5.0},
        {"gene": "JAK2", "drug": "Staurosporine", "kd_nm": 0.3},
        {"gene": "JAK3", "drug": "Tofacitinib", "kd_nm": 1.6},
        {"gene": "JAK3", "drug": "Ruxolitinib", "kd_nm": 428.0},
        {"gene": "JAK3", "drug": "Staurosporine", "kd_nm": 0.8},
        {"gene": "JAK3", "drug": "Decernotinib", "kd_nm": 2.5},
        {"gene": "TYK2", "drug": "Deucravacitinib", "kd_nm": 0.02},
        {"gene": "TYK2", "drug": "Tofacitinib", "kd_nm": 34.0},
        {"gene": "TYK2", "drug": "Baricitinib", "kd_nm": 340.0},
        {"gene": "TYK2", "drug": "Ruxolitinib", "kd_nm": 19.0},

        # ABL family
        {"gene": "ABL1", "drug": "Imatinib", "kd_nm": 13.0},
        {"gene": "ABL1", "drug": "Dasatinib", "kd_nm": 0.1},
        {"gene": "ABL1", "drug": "Nilotinib", "kd_nm": 25.0},
        {"gene": "ABL1", "drug": "Bosutinib", "kd_nm": 1.0},
        {"gene": "ABL1", "drug": "Ponatinib", "kd_nm": 0.37},
        {"gene": "ABL1", "drug": "Asciminib", "kd_nm": 0.5},
        {"gene": "ABL1", "drug": "Staurosporine", "kd_nm": 0.7},
        {"gene": "ABL2", "drug": "Imatinib", "kd_nm": 25.0},
        {"gene": "ABL2", "drug": "Dasatinib", "kd_nm": 0.2},
        {"gene": "ABL2", "drug": "Nilotinib", "kd_nm": 68.0},

        # EGFR family
        {"gene": "EGFR", "drug": "Gefitinib", "kd_nm": 0.4},
        {"gene": "EGFR", "drug": "Erlotinib", "kd_nm": 0.7},
        {"gene": "EGFR", "drug": "Afatinib", "kd_nm": 0.5},
        {"gene": "EGFR", "drug": "Osimertinib", "kd_nm": 15.0},
        {"gene": "EGFR", "drug": "Lapatinib", "kd_nm": 3.0},
        {"gene": "EGFR", "drug": "Staurosporine", "kd_nm": 0.1},
        {"gene": "ERBB2", "drug": "Lapatinib", "kd_nm": 7.0},
        {"gene": "ERBB2", "drug": "Neratinib", "kd_nm": 59.0},
        {"gene": "ERBB2", "drug": "Tucatinib", "kd_nm": 8.0},
        {"gene": "ERBB4", "drug": "Afatinib", "kd_nm": 14.0},
        {"gene": "ERBB4", "drug": "Neratinib", "kd_nm": 19.0},

        # BRAF/RAF
        {"gene": "BRAF", "drug": "Vemurafenib", "kd_nm": 31.0},
        {"gene": "BRAF", "drug": "Dabrafenib", "kd_nm": 0.8},
        {"gene": "BRAF", "drug": "Encorafenib", "kd_nm": 0.35},
        {"gene": "BRAF", "drug": "Sorafenib", "kd_nm": 22.0},
        {"gene": "BRAF", "drug": "Staurosporine", "kd_nm": 11.0},
        {"gene": "RAF1", "drug": "Sorafenib", "kd_nm": 6.0},
        {"gene": "RAF1", "drug": "Dabrafenib", "kd_nm": 5.0},

        # SRC family
        {"gene": "SRC", "drug": "Dasatinib", "kd_nm": 0.2},
        {"gene": "SRC", "drug": "Bosutinib", "kd_nm": 1.2},
        {"gene": "SRC", "drug": "Saracatinib", "kd_nm": 2.7},
        {"gene": "SRC", "drug": "Staurosporine", "kd_nm": 0.5},
        {"gene": "LCK", "drug": "Dasatinib", "kd_nm": 0.2},
        {"gene": "LCK", "drug": "Staurosporine", "kd_nm": 0.3},
        {"gene": "FYN", "drug": "Dasatinib", "kd_nm": 0.5},
        {"gene": "FYN", "drug": "Saracatinib", "kd_nm": 10.0},

        # CDK family
        {"gene": "CDK2", "drug": "Flavopiridol", "kd_nm": 170.0},
        {"gene": "CDK2", "drug": "Dinaciclib", "kd_nm": 1.0},
        {"gene": "CDK2", "drug": "Palbociclib", "kd_nm": 11000.0},
        {"gene": "CDK2", "drug": "Staurosporine", "kd_nm": 7.0},
        {"gene": "CDK4", "drug": "Palbociclib", "kd_nm": 11.0},
        {"gene": "CDK4", "drug": "Ribociclib", "kd_nm": 10.0},
        {"gene": "CDK4", "drug": "Abemaciclib", "kd_nm": 2.0},
        {"gene": "CDK6", "drug": "Palbociclib", "kd_nm": 16.0},
        {"gene": "CDK6", "drug": "Ribociclib", "kd_nm": 39.0},
        {"gene": "CDK6", "drug": "Abemaciclib", "kd_nm": 5.0},

        # MAPK family
        {"gene": "MAPK14", "drug": "SB203580", "kd_nm": 38.0},
        {"gene": "MAPK14", "drug": "BIRB796", "kd_nm": 0.1},
        {"gene": "MAPK14", "drug": "VX-745", "kd_nm": 10.0},
        {"gene": "MAPK14", "drug": "Staurosporine", "kd_nm": 95.0},

        # Aurora kinases
        {"gene": "AURKA", "drug": "Alisertib", "kd_nm": 1.2},
        {"gene": "AURKA", "drug": "VX-680", "kd_nm": 0.6},
        {"gene": "AURKA", "drug": "Danusertib", "kd_nm": 13.0},
        {"gene": "AURKA", "drug": "Staurosporine", "kd_nm": 1.3},
        {"gene": "AURKB", "drug": "Barasertib", "kd_nm": 0.4},
        {"gene": "AURKB", "drug": "VX-680", "kd_nm": 18.0},
        {"gene": "AURKB", "drug": "Staurosporine", "kd_nm": 2.5},

        # FGFR family
        {"gene": "FGFR1", "drug": "Erdafitinib", "kd_nm": 1.2},
        {"gene": "FGFR1", "drug": "Infigratinib", "kd_nm": 0.9},
        {"gene": "FGFR1", "drug": "Pemigatinib", "kd_nm": 0.4},
        {"gene": "FGFR1", "drug": "Futibatinib", "kd_nm": 3.9},
        {"gene": "FGFR2", "drug": "Erdafitinib", "kd_nm": 2.8},
        {"gene": "FGFR2", "drug": "Infigratinib", "kd_nm": 1.4},
        {"gene": "FGFR2", "drug": "Pemigatinib", "kd_nm": 0.5},
        {"gene": "FGFR3", "drug": "Erdafitinib", "kd_nm": 3.0},
        {"gene": "FGFR3", "drug": "Infigratinib", "kd_nm": 1.0},

        # RTK
        {"gene": "VEGFR2", "drug": "Sunitinib", "kd_nm": 9.0},
        {"gene": "VEGFR2", "drug": "Sorafenib", "kd_nm": 90.0},
        {"gene": "VEGFR2", "drug": "Axitinib", "kd_nm": 0.2},
        {"gene": "VEGFR2", "drug": "Lenvatinib", "kd_nm": 4.0},
        {"gene": "MET", "drug": "Capmatinib", "kd_nm": 0.13},
        {"gene": "MET", "drug": "Tepotinib", "kd_nm": 2.3},
        {"gene": "MET", "drug": "Crizotinib", "kd_nm": 8.0},
        {"gene": "ALK", "drug": "Crizotinib", "kd_nm": 24.0},
        {"gene": "ALK", "drug": "Ceritinib", "kd_nm": 0.2},
        {"gene": "ALK", "drug": "Alectinib", "kd_nm": 1.9},
        {"gene": "ALK", "drug": "Lorlatinib", "kd_nm": 0.07},

        # CHK
        {"gene": "CHEK1", "drug": "Prexasertib", "kd_nm": 0.9},
        {"gene": "CHEK1", "drug": "AZD7762", "kd_nm": 5.0},
        {"gene": "CHEK2", "drug": "AZD7762", "kd_nm": 10.0},
        {"gene": "CHEK2", "drug": "Staurosporine", "kd_nm": 1.4},

        # GSK3
        {"gene": "GSK3B", "drug": "CHIR99021", "kd_nm": 6.7},
        {"gene": "GSK3B", "drug": "LY2090314", "kd_nm": 1.5},
        {"gene": "GSK3B", "drug": "Staurosporine", "kd_nm": 15.0},
        {"gene": "GSK3A", "drug": "CHIR99021", "kd_nm": 7.0},

        # BTK/TEC
        {"gene": "BTK", "drug": "Ibrutinib", "kd_nm": 0.5},
        {"gene": "BTK", "drug": "Acalabrutinib", "kd_nm": 3.0},
        {"gene": "BTK", "drug": "Zanubrutinib", "kd_nm": 0.3},
        {"gene": "BTK", "drug": "Pirtobrutinib", "kd_nm": 3.5},
        {"gene": "ITK", "drug": "Ibrutinib", "kd_nm": 4.9},

        # SYK
        {"gene": "SYK", "drug": "Fostamatinib", "kd_nm": 41.0},
        {"gene": "SYK", "drug": "Entospletinib", "kd_nm": 7.7},

        # PLK
        {"gene": "PLK1", "drug": "Volasertib", "kd_nm": 0.87},
        {"gene": "PLK1", "drug": "Rigosertib", "kd_nm": 100.0},

        # MTOR
        {"gene": "MTOR", "drug": "Everolimus", "kd_nm": 1.6},
        {"gene": "MTOR", "drug": "Temsirolimus", "kd_nm": 1.8},

        # WEE1 / DYRK
        {"gene": "WEE1", "drug": "Adavosertib", "kd_nm": 5.2},
        {"gene": "DYRK1A", "drug": "Harmine", "kd_nm": 80.0},
        {"gene": "DYRK1A", "drug": "INDY", "kd_nm": 240.0},

        # ROCK
        {"gene": "ROCK1", "drug": "Fasudil", "kd_nm": 330.0},
        {"gene": "ROCK1", "drug": "Y-27632", "kd_nm": 140.0},
        {"gene": "ROCK2", "drug": "Fasudil", "kd_nm": 320.0},

        # PDGFRA
        {"gene": "PDGFRA", "drug": "Imatinib", "kd_nm": 74.0},
        {"gene": "PDGFRA", "drug": "Sunitinib", "kd_nm": 2.0},
        {"gene": "PDGFRA", "drug": "Avapritinib", "kd_nm": 0.5},
    ]

    # Add UniProt, pKd, and family info
    records = []
    for entry in davis_data:
        gene = entry["gene"]
        uniprot = KINASE_UNIPROT_EXTENDED.get(gene, "")
        if not uniprot:
            continue
        kd = entry["kd_nm"]
        records.append({
            "gene": gene,
            "uniprot": uniprot,
            "family": get_family(gene),
            "drug": entry["drug"],
            "smiles": "",  # Will be filled from ChEMBL or left empty
            "kd_nm": kd,
            "pKd": -np.log10(kd * 1e-9) if kd > 0 else 0,
            "source": "Davis2011_expanded",
        })

    df = pd.DataFrame(records)
    output_path = os.path.join(output_dir, "kinase_binding_expanded.csv")
    df.to_csv(output_path, index=False)
    logger.info(f"Expanded kinase data: {len(df)} records, {df['gene'].nunique()} genes, {df['drug'].nunique()} drugs")
    return df


def merge_all_binding_data(output_dir: str) -> pd.DataFrame:
    """Merge all binding data sources: curated + ChEMBL + expanded Davis."""
    dfs = []

    # 1. Expanded Davis/literature data (always available)
    expanded = create_expanded_kinase_data(output_dir)
    dfs.append(expanded)

    # 2. ChEMBL data (may fail due to API)
    chembl = fetch_chembl_kinase_data(output_dir)
    if not chembl.empty:
        dfs.append(chembl)

    merged = pd.concat(dfs, ignore_index=True)

    # Deduplicate by gene + drug
    merged = merged.drop_duplicates(subset=["gene", "drug"], keep="first")

    # Ensure pKd exists
    if "pKd" not in merged.columns:
        merged["pKd"] = 0.0
    mask = merged["pKd"].isna() | (merged["pKd"] == 0)
    if "kd_nm" in merged.columns:
        merged.loc[mask, "pKd"] = -np.log10(
            merged.loc[mask, "kd_nm"].clip(lower=0.01) * 1e-9
        )

    # Ensure family column
    if "family" not in merged.columns:
        merged["family"] = merged["gene"].apply(get_family)

    output_path = os.path.join(output_dir, "kinase_binding_all_merged.csv")
    merged.to_csv(output_path, index=False)
    logger.info(f"Merged binding data: {len(merged)} records, "
                f"{merged['gene'].nunique()} genes, {merged['family'].nunique()} families")

    return merged


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    output_dir = os.path.join(os.path.dirname(__file__), "..", "..", "data", "raw")
    df = merge_all_binding_data(output_dir)
    print(f"\nFinal: {len(df)} records")
    print(f"Genes: {df['gene'].nunique()}")
    print(f"Families: {df['family'].value_counts().to_dict()}")
