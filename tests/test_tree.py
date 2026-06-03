"""
tests/test_tree.py
Integration tests for the lineage tree builder.
Uses a small synthetic clone (5 cells) so tests run without real data files.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
from lineage.germline import infer_germline
from lineage.tree import build_tree, TreeNode


# ── Synthetic clone fixture ───────────────────────────────────────────────────

SYNTHETIC_CLONE = pd.DataFrame([
    {
        "cell_id": "cell_001", "clone_id": "TEST_CLONE", "clone_count": 5,
        "VDJ_sequence_H":    "ATCGATCGATCG",
        "VDJ_sequence_L":    "GCTAGCTAGCTA",
        "VDJ_aa_sequence_H": "IDID",
        "VDJ_aa_sequence_L": "ASAS",
        "v_call_h": "IGHV3-23*01", "d_call_h": "IGHD1-1*01",
        "j_call_h": "IGHJ4*01",    "v_call_l": "IGKV1-5*01", "j_call_l": "IGKJ1*01",
        "mu_count_h": 0,  "mu_count_l": 0,
        "cluster_annotated": "NBC", "c_call": "IGHM", "sample_id": "T1",
    },
    {
        "cell_id": "cell_002", "clone_id": "TEST_CLONE", "clone_count": 5,
        "VDJ_sequence_H":    "ATCGATCGTTCG",   # 1 mutation vs cell_001
        "VDJ_sequence_L":    "GCTAGCTAGCTA",
        "VDJ_aa_sequence_H": "IDIT",
        "VDJ_aa_sequence_L": "ASAS",
        "v_call_h": "IGHV3-23*01", "d_call_h": "IGHD1-1*01",
        "j_call_h": "IGHJ4*01",    "v_call_l": "IGKV1-5*01", "j_call_l": "IGKJ1*01",
        "mu_count_h": 1,  "mu_count_l": 0,
        "cluster_annotated": "GC_Cycling", "c_call": "IGHG1", "sample_id": "T1",
    },
    {
        "cell_id": "cell_003", "clone_id": "TEST_CLONE", "clone_count": 5,
        "VDJ_sequence_H":    "ATCGATCGTTCG",
        "VDJ_sequence_L":    "GCTAGCTAGGTA",   # 1 mutation vs cell_001
        "VDJ_aa_sequence_H": "IDIT",
        "VDJ_aa_sequence_L": "ASAG",
        "v_call_h": "IGHV3-23*01", "d_call_h": "IGHD1-1*01",
        "j_call_h": "IGHJ4*01",    "v_call_l": "IGKV1-5*01", "j_call_l": "IGKJ1*01",
        "mu_count_h": 1,  "mu_count_l": 1,
        "cluster_annotated": "MBC_Basal_CS", "c_call": "IGHG1", "sample_id": "T2",
    },
    {
        "cell_id": "cell_004", "clone_id": "TEST_CLONE", "clone_count": 5,
        "VDJ_sequence_H":    "ATCGATCGTTCG",
        "VDJ_sequence_L":    "GCTAGCTAGGTA",
        "VDJ_aa_sequence_H": "IDIT",
        "VDJ_aa_sequence_L": "ASAG",
        "v_call_h": "IGHV3-23*01", "d_call_h": "IGHD1-1*01",
        "j_call_h": "IGHJ4*01",    "v_call_l": "IGKV1-5*01", "j_call_l": "IGKJ1*01",
        "mu_count_h": 5,  "mu_count_l": 3,
        "cluster_annotated": "PC", "c_call": "IGHG2", "sample_id": "T2",
    },
    {
        "cell_id": "cell_005", "clone_id": "TEST_CLONE", "clone_count": 5,
        "VDJ_sequence_H":    "TTCGATCGTTCG",   # further mutated
        "VDJ_sequence_L":    "GCTAGCTAGGTA",
        "VDJ_aa_sequence_H": "FDIT",
        "VDJ_aa_sequence_L": "ASAG",
        "v_call_h": "IGHV3-23*01", "d_call_h": "IGHD1-1*01",
        "j_call_h": "IGHJ4*01",    "v_call_l": "IGKV1-5*01", "j_call_l": "IGKJ1*01",
        "mu_count_h": 8,  "mu_count_l": 3,
        "cluster_annotated": "PC", "c_call": "IGHG2", "sample_id": "T3",
    },
])


# ── Tests ─────────────────────────────────────────────────────────────────────

def test_infer_germline_picks_lowest_mu():
    germline = infer_germline(SYNTHETIC_CLONE)
    assert germline["cell_id"] == "GERMLINE"
    assert germline["mu_count_h"] == 0
    assert germline["mu_count_l"] == 0
    # Should match cell_001 (lowest mu)
    assert germline["VDJ_sequence_H"] == "ATCGATCGATCG"


def test_build_tree_returns_root():
    germline = infer_germline(SYNTHETIC_CLONE)
    root, all_nodes, all_mutations = build_tree(SYNTHETIC_CLONE, germline)
    assert isinstance(root, TreeNode)
    assert root.cell_data["cell_id"] == "GERMLINE"


def test_build_tree_has_children():
    germline = infer_germline(SYNTHETIC_CLONE)
    root, all_nodes, _ = build_tree(SYNTHETIC_CLONE, germline)
    assert len(root.children) > 0


def test_build_tree_all_nodes_present():
    germline = infer_germline(SYNTHETIC_CLONE)
    root, all_nodes, _ = build_tree(SYNTHETIC_CLONE, germline)
    node_ids = {n.node_id for n in all_nodes}
    # Cells with identical sequences are deduplicated intentionally.
    assert "GERMLINE" in node_ids
    assert len(node_ids) >= 3


def test_mutations_recorded():
    germline = infer_germline(SYNTHETIC_CLONE)
    _, _, all_mutations = build_tree(SYNTHETIC_CLONE, germline)
    assert len(all_mutations) > 0
    # Each mutation record should have required keys
    for rec in all_mutations:
        assert "parent_cell_id" in rec
        assert "child_cell_id"  in rec
        assert "chain"          in rec
        assert "position"       in rec
        assert "ref"            in rec
        assert "alt"            in rec
