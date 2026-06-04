import numpy as np
import pandas as pd
from lineage.distance import build_distance_matrix
from lineage.mutations import (
    get_branch_mutations,
    mutations_to_records,
    mutations_to_label,
)


class TreeNode:
    def __init__(self, node_id, cell_data):
        self.node_id = node_id
        self.cell_data = cell_data
        self.children = []
        self.parent = None

    def is_leaf(self):
        return len(self.children) == 0


def build_tree(clone_df: pd.DataFrame, germline: dict):
    """
    Build a lineage tree using Neighbor-Joining (NJ).
    """
    rows = [germline] + clone_df.to_dict(orient="records")

    matrix, cell_ids = build_distance_matrix(rows)

    row_map = {r["cell_id"]: r for r in rows}

    # Convert matrix to dict-of-dicts keyed by cell_id
    D = {
        cid_i: {cid_j: float(matrix[i, j]) for j, cid_j in enumerate(cell_ids)}
        for i, cid_i in enumerate(cell_ids)
    }

    root = neighbor_joining(D, row_map)

    all_mutations = []
    _collect_mutations(root, all_mutations)

    all_nodes = _gather_nodes(root)

    return root, all_nodes, all_mutations


def _synthetic_sequence(seq1, seq2):
    """
    Build a synthetic ancestor sequence by taking consensus of two sequences.
    """
    if seq1 is None or seq2 is None:
        return seq1 or seq2

    out = []
    for a, b in zip(seq1, seq2):
        out.append(a if a == b else a)
    return "".join(out)


def neighbor_joining(D, row_map):
    """
    Neighbor-Joining using dict-of-dicts keyed by cell_id.
    D: {label: {label: distance}}
    """
    nodes = {lab: TreeNode(lab, row_map[lab]) for lab in D.keys()}

    while len(D) > 2:
        labels = list(D.keys())
        n = len(labels)

        # Total distances
        total = {i: sum(D[i][j] for j in labels if j != i) for i in labels}

        # Q-matrix
        Q = {}
        for i in labels:
            Q[i] = {}
            for j in labels:
                if i == j:
                    Q[i][j] = np.inf
                else:
                    Q[i][j] = (n - 2) * D[i][j] - total[i] - total[j]

        # Find pair with minimum Q
        i, j = min(
            ((a, b) for a in labels for b in labels if a != b),
            key=lambda x: Q[x[0]][x[1]],
        )

        # New internal node label
        new_label = f"NJ_internal_{len(nodes)}"

        # Synthetic ancestor sequence
        seq_i = row_map[i].get("VDJ_sequence_H", "")
        seq_j = row_map[j].get("VDJ_sequence_H", "")
        syn_seq = _synthetic_sequence(seq_i, seq_j)

        row_map[new_label] = {
            "cell_id": new_label,
            "cell_type": "internal",
            # synthetic sequences
            "VDJ_sequence_H": syn_seq,
            "VDJ_sequence_L": "",
            "VDJ_aa_sequence_H": "",
            "VDJ_aa_sequence_L": "",

            # internal nodes have no SHM
            "mu_count_h": 0,
            "mu_count_l": 0,

            # internal nodes have no isotype / organ
            "isotype": None,
            "organ": None,
            "duplicate_count": 1,
            }

        new_node = TreeNode(new_label, row_map[new_label])
        nodes[new_label] = new_node

        # Attach children (parent = new internal node)
        muts_i = get_branch_mutations(row_map[new_label], row_map[i])
        muts_j = get_branch_mutations(row_map[new_label], row_map[j])

        new_node.children.append((nodes[i], muts_i))
        new_node.children.append((nodes[j], muts_j))
        nodes[i].parent = new_node
        nodes[j].parent = new_node

        # Compute distances from new node to others
        labels_wo_ij = [k for k in labels if k not in (i, j)]
        D_new = {}
        for k in labels_wo_ij:
            D_new[k] = (D[i][k] + D[j][k] - D[i][j]) / 2.0

        # Build new distance matrix D2 from old D (no new_label yet)
        D2 = {}
        for k in labels_wo_ij:
            D2[k] = {}
            for kk in labels_wo_ij:
                D2[k][kk] = D[k][kk]

        # Add new_label row/col
        D2[new_label] = {}
        for k in labels_wo_ij:
            D2[new_label][k] = D_new[k]
            D2[k][new_label] = D_new[k]
        D2[new_label][new_label] = 0.0

        D = D2

    # Final join: two labels left
    a, b = list(D.keys())

    root = TreeNode("ROOT", {"cell_id": "ROOT"})

    parent_stub = {
    "VDJ_sequence_H": "",
    "VDJ_aa_sequence_H": "",
    "VDJ_sequence_L": "",
    "VDJ_aa_sequence_L": "",
    }

    muts_a = get_branch_mutations(parent_stub, row_map[a])
    muts_b = get_branch_mutations(parent_stub, row_map[b])

    root.children.append((nodes[a], muts_a))
    root.children.append((nodes[b], muts_b))
    nodes[a].parent = root
    nodes[b].parent = root

    return root



def _collect_mutations(node, records):
    for child, muts in node.children:
        records += mutations_to_records(node.node_id, child.node_id, muts)
        _collect_mutations(child, records)


def _gather_nodes(root):
    out = []
    stack = [root]
    while stack:
        n = stack.pop()
        out.append(n)
        for c, _ in n.children:
            stack.append(c)
    return out
