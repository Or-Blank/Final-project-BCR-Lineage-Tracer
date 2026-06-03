"""
distance.py
Computes pairwise sequence distances between cells in a clone.

Distance metric:
  Combined Hamming distance on aligned nucleotide sequences (H + L).
  Sequences are zero-padded to equal length before comparison so that
  length differences also contribute to distance.

Returns a square numpy distance matrix plus an ordered list of cell IDs.
"""

import numpy as np
import pandas as pd


def hamming(seq_a: str, seq_b: str) -> int:
    """
    Hamming distance between two strings.
    Pads the shorter one with '-' so length differences count as mismatches.
    """
    len_a, len_b = len(seq_a), len(seq_b)
    if len_a < len_b:
        seq_a = seq_a + "-" * (len_b - len_a)
    elif len_b < len_a:
        seq_b = seq_b + "-" * (len_a - len_b)
    return sum(c1 != c2 for c1, c2 in zip(seq_a, seq_b))


def combined_distance(row_a: dict, row_b: dict) -> float:
    """
    Combined distance = Hamming(H nucleotide) + Hamming(L nucleotide).
    If L chain sequences are absent (heavy-only files), distance is H only.
    """
    d_h = hamming(str(row_a["VDJ_sequence_H"]), str(row_b["VDJ_sequence_H"]))
    seq_l_a = str(row_a.get("VDJ_sequence_L", "") or "")
    seq_l_b = str(row_b.get("VDJ_sequence_L", "") or "")
    d_l = hamming(seq_l_a, seq_l_b) if seq_l_a and seq_l_b else 0
    return float(d_h + d_l)


def build_distance_matrix(rows: list) -> tuple:
    """
    Build a symmetric pairwise distance matrix for a list of sequence dicts.

    Parameters
    ----------
    rows : list of dicts, each with VDJ_sequence_H and VDJ_sequence_L

    Returns
    -------
    matrix : np.ndarray shape (n, n)
    cell_ids : list of str  (same order as matrix rows/cols)
    """
    n = len(rows)
    matrix = np.zeros((n, n), dtype=float)
    cell_ids = [r["cell_id"] for r in rows]

    for i in range(n):
        for j in range(i + 1, n):
            d = combined_distance(rows[i], rows[j])
            matrix[i, j] = d
            matrix[j, i] = d

    return matrix, cell_ids
