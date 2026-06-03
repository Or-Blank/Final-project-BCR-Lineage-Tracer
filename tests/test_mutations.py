"""
tests/test_mutations.py
Unit tests for mutation detection logic.
"""
import pytest
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from lineage.mutations import diff_sequences, get_branch_mutations, mutations_to_label, Mutation


def test_diff_identical_sequences():
    muts = diff_sequences("ATCG", "ATCG", chain="H", level="nt")
    assert muts == []


def test_diff_single_mutation():
    muts = diff_sequences("ATCG", "ATGG", chain="H", level="nt")
    assert len(muts) == 1
    assert muts[0].position == 3
    assert muts[0].ref == "C"
    assert muts[0].alt == "G"
    assert muts[0].chain == "H"


def test_diff_multiple_mutations():
    muts = diff_sequences("AAAA", "TTTT", chain="L", level="aa")
    assert len(muts) == 4


def test_diff_different_lengths():
    # Shorter seq – only compare up to min length
    muts = diff_sequences("ATCG", "ATC", chain="H", level="nt")
    assert all(m.position <= 3 for m in muts)


def test_get_branch_mutations_returns_both_chains():
    parent = {
        "VDJ_sequence_H": "AAAA",
        "VDJ_sequence_L": "CCCC",
        "VDJ_aa_sequence_H": "MM",
        "VDJ_aa_sequence_L": "KK",
    }
    child = {
        "VDJ_sequence_H": "AAAT",     # 1 nt change in H
        "VDJ_sequence_L": "CCCC",     # no change in L
        "VDJ_aa_sequence_H": "MM",
        "VDJ_aa_sequence_L": "KK",
    }
    muts = get_branch_mutations(parent, child)
    chains = [m.chain for m in muts]
    assert "H" in chains


def test_mutations_to_label_empty():
    assert mutations_to_label([]) == ""


def test_mutations_to_label_silent_only():
    muts = [Mutation(chain="H", level="nt", position=1, ref="A", alt="T")]
    label = mutations_to_label(muts)
    assert "silent" in label


def test_mutations_to_label_aa():
    muts = [Mutation(chain="H", level="aa", position=56, ref="T", alt="A")]
    label = mutations_to_label(muts)
    assert "T56A" in label


def test_mutation_str():
    m = Mutation(chain="H", level="aa", position=56, ref="T", alt="A")
    assert str(m) == "H:T56A"
