"""Fetch isoform sequences and annotations from UniProt.

Downloads protein sequences for all isoforms of target kinases,
identifies splicing events, and builds isoform-canonical alignments.
"""

import os
import json
import logging
import time
import requests
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from io import StringIO

logger = logging.getLogger(__name__)

UNIPROT_API = "https://rest.uniprot.org"


def fetch_uniprot_entry(uniprot_id: str, include_isoforms: bool = True) -> Dict:
    """Fetch full UniProt entry including isoform annotations."""
    url = f"{UNIPROT_API}/uniprotkb/{uniprot_id}.json"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        return resp.json()
    except Exception as e:
        logger.error(f"Failed to fetch UniProt entry {uniprot_id}: {e}")
        return {}


def fetch_isoform_sequences(uniprot_id: str) -> Dict[str, str]:
    """Fetch all isoform sequences for a UniProt entry.

    Returns dict mapping isoform_id -> sequence string.
    """
    # Fetch canonical sequence
    url = f"{UNIPROT_API}/uniprotkb/{uniprot_id}.fasta"
    sequences = {}

    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        for record in SeqIO.parse(StringIO(resp.text), "fasta"):
            sequences[f"{uniprot_id}-1"] = str(record.seq)
    except Exception as e:
        logger.error(f"Failed to fetch canonical sequence for {uniprot_id}: {e}")
        return sequences

    # Fetch isoform sequences
    url_iso = f"{UNIPROT_API}/uniprotkb/{uniprot_id}/isoforms"
    try:
        resp = requests.get(url_iso, timeout=30, headers={"Accept": "application/json"})
        if resp.status_code == 200:
            data = resp.json()
            for entry in data if isinstance(data, list) else [data]:
                if "sequence" in entry and "value" in entry["sequence"]:
                    iso_id = entry.get("isoformId", entry.get("accession", ""))
                    sequences[iso_id] = entry["sequence"]["value"]
    except Exception as e:
        logger.warning(f"Could not fetch isoforms for {uniprot_id}: {e}")

    # Alternative: try fetching FASTA with isoforms
    if len(sequences) <= 1:
        url_fasta = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta?include=yes"
        try:
            resp = requests.get(url_fasta, timeout=30)
            if resp.status_code == 200:
                for record in SeqIO.parse(StringIO(resp.text), "fasta"):
                    header = record.description
                    if "Isoform" in header or uniprot_id in record.id:
                        seq_id = record.id.split("|")[1] if "|" in record.id else record.id
                        sequences[seq_id] = str(record.seq)
        except Exception as e:
            logger.warning(f"Could not fetch isoform FASTA for {uniprot_id}: {e}")

    return sequences


def align_isoforms(canonical_seq: str, isoform_seq: str) -> Dict:
    """Align isoform sequence to canonical and identify differences.

    Returns dict with alignment info and splicing event details.
    """
    if canonical_seq == isoform_seq:
        return {"identical": True, "events": []}

    # Use pairwise alignment
    alignments = pairwise2.align.globalms(
        canonical_seq, isoform_seq,
        match=2, mismatch=-1, open=-10, extend=-0.5,
        one_alignment_only=True,
    )

    if not alignments:
        return {"identical": False, "events": [], "error": "alignment_failed"}

    aln = alignments[0]
    aln_can = aln.seqA
    aln_iso = aln.seqB

    events = []
    i_can = 0
    i_iso = 0

    in_gap_can = False
    in_gap_iso = False
    gap_start_can = -1
    gap_start_iso = -1

    for pos in range(len(aln_can)):
        c_can = aln_can[pos]
        c_iso = aln_iso[pos]

        if c_can == "-" and not in_gap_can:
            # Insertion in isoform (exon insertion)
            in_gap_can = True
            gap_start_can = i_can
            gap_start_iso = i_iso
        elif c_can != "-" and in_gap_can:
            in_gap_can = False
            events.append({
                "type": "insertion",
                "canonical_pos": gap_start_can,
                "isoform_pos_start": gap_start_iso,
                "isoform_pos_end": i_iso,
                "length": i_iso - gap_start_iso,
                "inserted_seq": isoform_seq[gap_start_iso:i_iso],
            })

        if c_iso == "-" and not in_gap_iso:
            # Deletion in isoform (exon skipping)
            in_gap_iso = True
            gap_start_can = i_can
            gap_start_iso = i_iso
        elif c_iso != "-" and in_gap_iso:
            in_gap_iso = False
            events.append({
                "type": "deletion",
                "canonical_pos_start": gap_start_can,
                "canonical_pos_end": i_can,
                "isoform_pos": gap_start_iso,
                "length": i_can - gap_start_can,
                "deleted_seq": canonical_seq[gap_start_can:i_can],
            })

        if c_can != "-":
            i_can += 1
        if c_iso != "-":
            i_iso += 1

    # Handle trailing gaps
    if in_gap_can:
        events.append({
            "type": "insertion",
            "canonical_pos": gap_start_can,
            "isoform_pos_start": gap_start_iso,
            "isoform_pos_end": i_iso,
            "length": i_iso - gap_start_iso,
        })
    if in_gap_iso:
        events.append({
            "type": "deletion",
            "canonical_pos_start": gap_start_can,
            "canonical_pos_end": i_can,
            "isoform_pos": gap_start_iso,
            "length": i_can - gap_start_can,
        })

    return {
        "identical": False,
        "canonical_length": len(canonical_seq),
        "isoform_length": len(isoform_seq),
        "num_events": len(events),
        "events": events,
        "sequence_identity": sum(1 for a, b in zip(aln_can, aln_iso) if a == b and a != "-") / max(len(canonical_seq), 1),
    }


def fetch_all_kinase_isoforms(kinase_targets: Dict, output_dir: str) -> Dict:
    """Fetch all isoform sequences for target kinases.

    Returns comprehensive isoform data structure.
    """
    os.makedirs(output_dir, exist_ok=True)
    cache_path = os.path.join(output_dir, "kinase_isoforms.json")

    if os.path.exists(cache_path):
        logger.info(f"Loading cached isoform data from {cache_path}")
        with open(cache_path) as f:
            return json.load(f)

    all_data = {}

    for family, members in kinase_targets.items():
        for gene, uniprot_id in members.items():
            logger.info(f"Fetching isoforms for {gene} ({uniprot_id})...")
            time.sleep(0.5)  # Rate limiting

            sequences = fetch_isoform_sequences(uniprot_id)

            # Identify canonical
            canonical_id = f"{uniprot_id}-1"
            canonical_seq = sequences.get(canonical_id, sequences.get(uniprot_id, ""))

            if not canonical_seq:
                # Try the first sequence as canonical
                if sequences:
                    canonical_id = list(sequences.keys())[0]
                    canonical_seq = sequences[canonical_id]

            gene_data = {
                "gene": gene,
                "family": family,
                "uniprot": uniprot_id,
                "canonical_id": canonical_id,
                "canonical_seq": canonical_seq,
                "canonical_length": len(canonical_seq),
                "isoforms": {},
            }

            # Analyze each isoform
            for iso_id, iso_seq in sequences.items():
                if iso_id == canonical_id:
                    continue
                alignment_info = align_isoforms(canonical_seq, iso_seq)
                gene_data["isoforms"][iso_id] = {
                    "sequence": iso_seq,
                    "length": len(iso_seq),
                    "alignment": alignment_info,
                }

            gene_data["num_isoforms"] = len(gene_data["isoforms"]) + 1  # +1 for canonical

            all_data[gene] = gene_data
            logger.info(f"  -> {gene_data['num_isoforms']} isoforms, "
                       f"canonical length: {len(canonical_seq)}")

            # Save sequences as FASTA
            fasta_path = os.path.join(output_dir, f"{gene}_isoforms.fasta")
            with open(fasta_path, "w") as f:
                f.write(f">{canonical_id} {gene} canonical\n{canonical_seq}\n")
                for iso_id, iso_info in gene_data["isoforms"].items():
                    f.write(f">{iso_id} {gene} isoform\n{iso_info['sequence']}\n")

    with open(cache_path, "w") as f:
        json.dump(all_data, f, indent=2)

    logger.info(f"Saved isoform data to {cache_path}")
    return all_data


def summarize_isoform_data(isoform_data: Dict) -> Dict:
    """Generate summary statistics of isoform data."""
    summary = {
        "total_genes": len(isoform_data),
        "total_isoforms": sum(d["num_isoforms"] for d in isoform_data.values()),
        "genes_with_multiple_isoforms": sum(1 for d in isoform_data.values() if d["num_isoforms"] > 1),
        "total_splicing_events": 0,
        "event_types": {"insertion": 0, "deletion": 0},
        "per_family": {},
    }

    for gene, data in isoform_data.items():
        family = data["family"]
        if family not in summary["per_family"]:
            summary["per_family"][family] = {"genes": 0, "isoforms": 0}
        summary["per_family"][family]["genes"] += 1
        summary["per_family"][family]["isoforms"] += data["num_isoforms"]

        for iso_id, iso_info in data["isoforms"].items():
            aln = iso_info.get("alignment", {})
            events = aln.get("events", [])
            summary["total_splicing_events"] += len(events)
            for evt in events:
                evt_type = evt.get("type", "unknown")
                summary["event_types"][evt_type] = summary["event_types"].get(evt_type, 0) + 1

    return summary


if __name__ == "__main__":
    from fetch_bindingdb import KINASE_TARGETS
    logging.basicConfig(level=logging.INFO)
    output_dir = os.path.join(os.path.dirname(__file__), "..", "..", "data", "raw")
    data = fetch_all_kinase_isoforms(KINASE_TARGETS, output_dir)
    summary = summarize_isoform_data(data)
    print(json.dumps(summary, indent=2))
