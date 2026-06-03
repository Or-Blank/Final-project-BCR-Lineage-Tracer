"""
germline.py
Infers the unmutated germline ancestor for a clone.

Strategy:
  1. Find the cell(s) with the lowest total SHM count (mu_count_h + mu_count_l).
     These are the closest observed sequences to the original naive B cell.
  2. If multiple cells tie at the minimum, take the consensus (most frequent
     nucleotide at each position).
  3. Return the germline sequence pair (H, L) as a named dict so the rest of
     the pipeline can use it as the tree root.

This avoids the need to download IMGT germline databases while still giving a
biologically meaningful root: the least-mutated observed ancestor is the
standard proxy used in many published lineage tools.
"""

import pandas as pd
from collections import Counter


def infer_germline(clone_df: pd.DataFrame) -> dict:
    """
    Infer the germline ancestor from the least-mutated cell(s) in the clone.

    Returns a dict with keys:
      cell_id          : "GERMLINE"
      VDJ_sequence_H   : heavy chain nucleotide sequence
      VDJ_sequence_L   : light chain nucleotide sequence
      VDJ_aa_sequence_H: heavy chain amino acid sequence
      VDJ_aa_sequence_L: light chain amino acid sequence
      mu_count_h       : 0  (by definition, the root has zero mutations)
      mu_count_l       : 0
      cluster_annotated: "Germline"
      c_call           : germline isotype (from least-mutated cell)
      sample_id        : "germline"
      v_call_h         : V gene of the root cell
      j_call_h         : J gene of the root cell
    """
    clone_df = clone_df.copy()
    clone_df["mu_total"] = clone_df["mu_count_h"] + clone_df["mu_count_l"]
    min_mu = clone_df["mu_total"].min()
    candidates = clone_df[clone_df["mu_total"] == min_mu]

    if len(candidates) == 1:
        row = candidates.iloc[0]
        seq_h = row["VDJ_sequence_H"]
        seq_l = row["VDJ_sequence_L"]
        aa_h  = row["VDJ_aa_sequence_H"]
        aa_l  = row["VDJ_aa_sequence_L"]
    else:
        # Build consensus from all tied candidates
        seq_h = _consensus([r["VDJ_sequence_H"] for _, r in candidates.iterrows()])
        seq_l = _consensus([r["VDJ_sequence_L"] for _, r in candidates.iterrows()])
        aa_h  = _consensus([r["VDJ_aa_sequence_H"] for _, r in candidates.iterrows()])
        aa_l  = _consensus([r["VDJ_aa_sequence_L"] for _, r in candidates.iterrows()])

    ref_row = candidates.iloc[0]
    return {
        "cell_id":           "GERMLINE",
        "VDJ_sequence_H":    seq_h,
        "VDJ_sequence_L":    seq_l,
        "VDJ_aa_sequence_H": aa_h,
        "VDJ_aa_sequence_L": aa_l,
        "mu_count_h":        0,
        "mu_count_l":        0,
        "cluster_annotated": "Germline",
        "c_call":            ref_row.get("c_call", ""),
        "sample_id":         "germline",
        "v_call_h":          ref_row.get("v_call_h", ""),
        "j_call_h":          ref_row.get("j_call_h", ""),
        "v_call_l":          ref_row.get("v_call_l", ""),
        "j_call_l":          ref_row.get("j_call_l", ""),
        "mu_total":          0,
    }


def _consensus(sequences: list) -> str:
    """
    Build a position-wise consensus string from a list of sequences.
    Sequences are truncated / padded to the length of the shortest one.
    At each position the most common character wins; ties go to the first.
    """
    if not sequences:
        return ""
    min_len = min(len(s) for s in sequences)
    result = []
    for i in range(min_len):
        chars = [s[i] for s in sequences]
        result.append(Counter(chars).most_common(1)[0][0])
    return "".join(result)
