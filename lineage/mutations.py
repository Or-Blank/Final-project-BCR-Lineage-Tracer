"""
mutations.py
Detects mutations between two sequences (parent → child).

For each branch in the lineage tree we report:
  - Nucleotide changes  : position, ref base, alt base  (e.g. A123G)
  - Amino acid changes  : position, ref aa,  alt aa      (e.g. T56A)
  - Whether the change is synonymous (silent) or non-synonymous

Chain is either "H" (heavy) or "L" (light).
"""

from dataclasses import dataclass


@dataclass
class Mutation:
    chain: str       # "H" or "L"
    level: str       # "nt" or "aa"
    position: int    # 1-based
    ref: str         # original character
    alt: str         # mutated character
    silent: bool = False   # only meaningful for aa level

    def __str__(self):
        return f"{self.chain}:{self.ref}{self.position}{self.alt}"


def diff_sequences(parent_seq: str, child_seq: str, chain: str, level: str) -> list:
    """
    Return a list of Mutation objects for positions that differ between
    parent_seq and child_seq.

    Parameters
    ----------
    parent_seq : str  (nucleotide or amino acid)
    child_seq  : str
    chain      : "H" or "L"
    level      : "nt" or "aa"
    """
    mutations = []
    min_len = min(len(parent_seq), len(child_seq))

    for i in range(min_len):
        ref = parent_seq[i]
        alt = child_seq[i]
        if ref != alt and ref != "-" and alt != "-":
            mutations.append(Mutation(
                chain=chain,
                level=level,
                position=i + 1,
                ref=ref,
                alt=alt,
            ))

    return mutations


def get_branch_mutations(parent: dict, child: dict) -> list:
    """
    Collect all mutations on the branch from parent to child node.
    Compares both nucleotide (nt) and amino acid (aa) sequences for H and L.

    Parameters
    ----------
    parent, child : dicts with keys VDJ_sequence_H/L, VDJ_aa_sequence_H/L

    Returns
    -------
    List of Mutation objects
    """
    mutations = []

    # Nucleotide level
    mutations += diff_sequences(
        parent["VDJ_sequence_H"], child["VDJ_sequence_H"], chain="H", level="nt"
    )
    mutations += diff_sequences(
        parent["VDJ_sequence_L"], child["VDJ_sequence_L"], chain="L", level="nt"
    )

    # Amino acid level
    mutations += diff_sequences(
        parent["VDJ_aa_sequence_H"], child["VDJ_aa_sequence_H"], chain="H", level="aa"
    )
    mutations += diff_sequences(
        parent["VDJ_aa_sequence_L"], child["VDJ_aa_sequence_L"], chain="L", level="aa"
    )

    return mutations


def mutations_to_label(mutations: list, max_show: int = 4) -> str:
    """
    Build a short human-readable label for a branch edge.
    Shows aa-level mutations only (more biologically meaningful for labelling).
    Falls back to 'silent' if only nt changes exist.
    """
    aa_muts = [m for m in mutations if m.level == "aa"]
    if not aa_muts:
        nt_count = len([m for m in mutations if m.level == "nt"])
        return f"{nt_count} silent" if nt_count else ""

    labels = [str(m) for m in aa_muts[:max_show]]
    if len(aa_muts) > max_show:
        labels.append(f"+{len(aa_muts) - max_show} more")
    return ", ".join(labels)


def mutations_to_records(parent_id: str, child_id: str, mutations: list) -> list:
    """
    Convert a list of Mutation objects to flat dicts for CSV export.
    """
    records = []
    for m in mutations:
        records.append({
            "parent_cell_id": parent_id,
            "child_cell_id":  child_id,
            "chain":          m.chain,
            "level":          m.level,
            "position":       m.position,
            "ref":            m.ref,
            "alt":            m.alt,
            "notation":       str(m),
        })
    return records
