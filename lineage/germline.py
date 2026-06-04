"""
germline.py
Infers the unmutated germline ancestor for a clone.

Strategy (in order of preference):
  1. If the input data contains a 'germline_alignment_h' column (produced by
     IgBLAST / Change-O / IMGT tools), use the consensus of those sequences
     as the root. This is the most accurate approach and matches what R tools
     like Dowser/IgPhyML use.

  2. Fallback: find the cell(s) with the lowest total SHM count and use their
     sequence as a proxy for the germline. This is less accurate but works
     when no germline reference is available.

The returned germline dict is structured identically to a cell row so the
rest of the pipeline (distance, tree, visualize) can treat it uniformly.
"""

import pandas as pd
from collections import Counter


def infer_germline(clone_df: pd.DataFrame) -> dict:
    """
    Infer the germline ancestor for a clone.

    Returns a dict with keys matching the standard cell row format.
    """
    if "germline_alignment_h" in clone_df.columns:
        germline_seqs = clone_df["germline_alignment_h"].dropna().astype(str).tolist()
        germline_seqs = [s for s in germline_seqs if s and s not in ("nan", "")]
        if germline_seqs:
            return _germline_from_alignment(germline_seqs, clone_df)

    # Fallback: least-mutated cell
    return _germline_from_least_mutated(clone_df)


def _germline_from_alignment(germline_seqs: list, clone_df: pd.DataFrame) -> dict:
    """
    Build the germline root from the germline_alignment_h column.
    Takes the consensus across all cells in the clone (they should be
    near-identical since they share the same V gene).
    """
    consensus_seq = _consensus(germline_seqs)
    ref_row = clone_df.iloc[0]

    return {
        "cell_id":              "GERMLINE",
        "VDJ_sequence_H":       consensus_seq,
        "VDJ_sequence_L":       "",
        "VDJ_aa_sequence_H":    "",
        "VDJ_aa_sequence_L":    "",
        "mu_count_h":           0,
        "mu_count_l":           0,
        "cluster_annotated":    "Germline",
        "c_call":               "",
        "sample_id":            "germline",
        "v_call_h":             ref_row.get("v_call_h", ""),
        "j_call_h":             ref_row.get("j_call_h", ""),
        "v_call_l":             ref_row.get("v_call_l", ""),
        "j_call_l":             ref_row.get("j_call_l", ""),
        "mu_total":             0,
    }


def _germline_from_least_mutated(clone_df: pd.DataFrame) -> dict:
    """
    Fallback: use the least-mutated observed cell as a proxy for germline.
    This cell is included as a real node AND as the root — so it will appear
    once in the tree as the germline anchor.
    """
    clone_df = clone_df.copy()
    clone_df["mu_total"] = clone_df["mu_count_h"] + clone_df["mu_count_l"]
    min_mu = clone_df["mu_total"].min()
    candidates = clone_df[clone_df["mu_total"] == min_mu]

    if len(candidates) == 1:
        row = candidates.iloc[0]
        seq_h = row["VDJ_sequence_H"]
        aa_h  = row["VDJ_aa_sequence_H"]
    else:
        seq_h = _consensus([r["VDJ_sequence_H"] for _, r in candidates.iterrows()])
        aa_h  = _consensus([r["VDJ_aa_sequence_H"] for _, r in candidates.iterrows()])

    ref_row = candidates.iloc[0]
    return {
        "cell_id":              "GERMLINE",
        "VDJ_sequence_H":       seq_h,
        "VDJ_sequence_L":       str(ref_row.get("VDJ_sequence_L", "") or ""),
        "VDJ_aa_sequence_H":    aa_h,
        "VDJ_aa_sequence_L":    str(ref_row.get("VDJ_aa_sequence_L", "") or ""),
        "mu_count_h":           0,
        "mu_count_l":           0,
        "cluster_annotated":    "Germline",
        "c_call":               ref_row.get("c_call", ""),
        "sample_id":            "germline",
        "v_call_h":             ref_row.get("v_call_h", ""),
        "j_call_h":             ref_row.get("j_call_h", ""),
        "v_call_l":             ref_row.get("v_call_l", ""),
        "j_call_l":             ref_row.get("j_call_l", ""),
        "mu_total":             0,
    }


def _consensus(sequences: list) -> str:
    """
    Position-wise consensus from a list of sequences.
    Truncates to the shortest sequence length.
    At each position the most common character wins.
    """
    if not sequences:
        return ""
    min_len = min(len(s) for s in sequences)
    result = []
    for i in range(min_len):
        chars = [s[i] for s in sequences]
        result.append(Counter(chars).most_common(1)[0][0])
    return "".join(result)
