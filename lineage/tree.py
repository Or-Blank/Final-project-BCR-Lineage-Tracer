import numpy as np
import pandas as pd
from itertools import combinations
from lineage.mutations import (
    get_branch_mutations,
    mutations_to_records,
)


class TreeNode:
    def __init__(self, node_id, cell_data):
        self.node_id = node_id
        self.cell_data = cell_data  # dict with sequences, metadata
        self.children = []          # list of (child_node, mutations)
        self.parent = None

    def is_leaf(self):
        return len(self.children) == 0


# ---------- PUBLIC ENTRY POINT ----------

def build_tree(clone_df: pd.DataFrame, germline: dict):
    """
    Build a BCR lineage tree MiXCR-style:
    - collapse identical sequences into clonotypes
    - build NJ tree on clonotypes + germline
    """
    # 1) collapse to clonotypes
    clonotypes = _collapse_to_clonotypes(clone_df)

    # 2) add germline as a pseudo-row
    rows = [germline] + clonotypes.to_dict(orient="records")

    # 3) build distance matrix (simple Hamming on heavy-chain nt)
    dist_matrix, labels = _build_distance_matrix(rows)

    # 4) map label -> row dict
    row_map = {r["cell_id"]: r for r in rows}

    # 5) run Neighbor-Joining
    root = _neighbor_joining(dist_matrix, labels, row_map)

    # 6) collect mutations per branch
    all_mutations = []
    _collect_mutations(root, all_mutations)

    # 7) gather all nodes
    all_nodes = _gather_nodes(root)

    return root, all_nodes, all_mutations


# ---------- STEP 1: CLONOTYPES ----------

def _collapse_to_clonotypes(clone_df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse identical sequences into clonotypes.
    Group by heavy+light nt + aa + isotype + organ.
    Sum duplicate_count, keep first SHM counts.
    """
    group_cols = [
        "VDJ_sequence_H",
        "VDJ_sequence_L",
        "VDJ_aa_sequence_H",
        "VDJ_aa_sequence_L",
        "isotype",
        "organ",
    ]

    agg_df = (
        clone_df
        .groupby(group_cols, as_index=False)
        .agg({
            "mu_count_h": "first",
            "mu_count_l": "first",
            "duplicate_count": "sum",
        })
    )

    # give each clonotype its own ID
    agg_df["cell_id"] = [
        f"clonotype_{i+1}" for i in range(len(agg_df))
    ]

    return agg_df


# ---------- STEP 2: DISTANCE MATRIX ----------

def _hamming(a: str, b: str) -> int:
    if a is None or b is None:
        return 0
    a = str(a)
    b = str(b)
    n = min(len(a), len(b))
    return sum(1 for i in range(n) if a[i] != b[i]) + abs(len(a) - len(b))


def _build_distance_matrix(rows):
    """
    rows: list of dicts with 'cell_id' and 'VDJ_sequence_H'
    returns: (numpy matrix, labels list)
    """
    labels = [r["cell_id"] for r in rows]
    n = len(labels)
    M = np.zeros((n, n), dtype=float)

    for i, j in combinations(range(n), 2):
        seq_i = rows[i].get("VDJ_sequence_H", "")
        seq_j = rows[j].get("VDJ_sequence_H", "")
        d = _hamming(seq_i, seq_j)
        M[i, j] = M[j, i] = d

    return M, labels


# ---------- STEP 3: NEIGHBOR-JOINING ----------

def _synthetic_sequence(seq1, seq2):
    """
    Simple consensus: if equal, keep; else keep seq1.
    """
    if seq1 is None and seq2 is None:
        return ""
    if seq1 is None:
        return str(seq2)
    if seq2 is None:
        return str(seq1)
    a = str(seq1)
    b = str(seq2)
    out = []
    for x, y in zip(a, b):
        out.append(x if x == y else x)
    return "".join(out)


def _neighbor_joining(D, labels, row_map):
    """
    Classic NJ on a numeric matrix D with labels.
    Returns root TreeNode.
    """
    # convert to dict-of-dicts keyed by label
    Ddict = {
        lab_i: {lab_j: float(D[i, j]) for j, lab_j in enumerate(labels)}
        for i, lab_i in enumerate(labels)
    }

    nodes = {lab: TreeNode(lab, row_map[lab]) for lab in labels}

    while len(Ddict) > 2:
        labs = list(Ddict.keys())
        n = len(labs)

        # total distances
        total = {i: sum(Ddict[i][j] for j in labs if j != i) for i in labs}

        # Q-matrix
        Q = {}
        for i in labs:
            Q[i] = {}
            for j in labs:
                if i == j:
                    Q[i][j] = float("inf")
                else:
                    Q[i][j] = (n - 2) * Ddict[i][j] - total[i] - total[j]

        # find min Q
        i, j = min(
            ((a, b) for a in labs for b in labs if a != b),
            key=lambda x: Q[x[0]][x[1]],
        )

        # new internal node
        new_label = f"NJ_internal_{len(nodes)}"

        # synthetic ancestor sequence (heavy only)
        seq_i = row_map[i].get("VDJ_sequence_H", "")
        seq_j = row_map[j].get("VDJ_sequence_H", "")
        syn_seq = _synthetic_sequence(seq_i, seq_j)

        row_map[new_label] = {
            "cell_id": new_label,
            "cell_type": "internal",
            "VDJ_sequence_H": syn_seq,
            "VDJ_sequence_L": "",
            "VDJ_aa_sequence_H": "",
            "VDJ_aa_sequence_L": "",
            "mu_count_h": 0,
            "mu_count_l": 0,
            "isotype": None,
            "organ": None,
            "duplicate_count": 1,
        }

        new_node = TreeNode(new_label, row_map[new_label])
        nodes[new_label] = new_node

        # attach children (parent = new internal node)
        muts_i = get_branch_mutations(row_map[new_label], row_map[i])
        muts_j = get_branch_mutations(row_map[new_label], row_map[j])

        new_node.children.append((nodes[i], muts_i))
        new_node.children.append((nodes[j], muts_j))
        nodes[i].parent = new_node
        nodes[j].parent = new_node

        # compute distances from new node to others
        labs_wo_ij = [k for k in labs if k not in (i, j)]
        D_new = {}
        for k in labs_wo_ij:
            D_new[k] = (Ddict[i][k] + Ddict[j][k] - Ddict[i][j]) / 2.0

        # rebuild Ddict without i,j, then add new_label
        D2 = {}
        for k in labs_wo_ij:
            D2[k] = {}
            for kk in labs_wo_ij:
                D2[k][kk] = Ddict[k][kk]

        D2[new_label] = {}
        for k in labs_wo_ij:
            D2[new_label][k] = D_new[k]
            D2[k][new_label] = D_new[k]
        D2[new_label][new_label] = 0.0

        Ddict = D2

    # final join
    a, b = list(Ddict.keys())

    root = TreeNode("ROOT", {"cell_id": "ROOT", "cell_type": "root"})

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


# ---------- STEP 4: MUTATION COLLECTION ----------

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
