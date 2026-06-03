"""
tests/test_distance.py
Unit tests for distance computation.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
from lineage.distance import hamming, combined_distance, build_distance_matrix


def test_hamming_identical():
    assert hamming("ATCG", "ATCG") == 0


def test_hamming_all_different():
    assert hamming("AAAA", "TTTT") == 4


def test_hamming_different_lengths():
    # Padding makes length difference count as mismatches
    assert hamming("AT", "ATCG") == 2


def test_combined_distance_zero():
    row = {"VDJ_sequence_H": "AAAA", "VDJ_sequence_L": "CCCC"}
    assert combined_distance(row, row) == 0.0


def test_combined_distance_positive():
    a = {"VDJ_sequence_H": "AAAA", "VDJ_sequence_L": "CCCC"}
    b = {"VDJ_sequence_H": "TTTT", "VDJ_sequence_L": "GGGG"}
    assert combined_distance(a, b) == 8.0


def test_build_distance_matrix_shape():
    rows = [
        {"cell_id": "A", "VDJ_sequence_H": "AAAA", "VDJ_sequence_L": "CC"},
        {"cell_id": "B", "VDJ_sequence_H": "AAAT", "VDJ_sequence_L": "CC"},
        {"cell_id": "C", "VDJ_sequence_H": "TTTT", "VDJ_sequence_L": "GG"},
    ]
    mat, ids = build_distance_matrix(rows)
    assert mat.shape == (3, 3)
    assert ids == ["A", "B", "C"]


def test_build_distance_matrix_symmetric():
    rows = [
        {"cell_id": "A", "VDJ_sequence_H": "AAAA", "VDJ_sequence_L": "CC"},
        {"cell_id": "B", "VDJ_sequence_H": "TTTT", "VDJ_sequence_L": "GG"},
    ]
    mat, _ = build_distance_matrix(rows)
    assert mat[0, 1] == mat[1, 0]


def test_build_distance_matrix_diagonal_zero():
    rows = [
        {"cell_id": "A", "VDJ_sequence_H": "AAAA", "VDJ_sequence_L": "CC"},
        {"cell_id": "B", "VDJ_sequence_H": "TTTT", "VDJ_sequence_L": "GG"},
    ]
    mat, _ = build_distance_matrix(rows)
    assert mat[0, 0] == 0.0
    assert mat[1, 1] == 0.0
