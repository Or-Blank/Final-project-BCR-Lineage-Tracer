"""
tree.py
Builds a phylogenetic lineage tree for a BCR clone.

Algorithm: UPGMA (Unweighted Pair Group Method with Arithmetic mean).
  - Simple, deterministic, well-understood.
  - Works well on BCR data where the germline root is explicitly defined.

The tree is represented as a list of nodes and a list of edges.
Each node carries the original cell metadata.
Each edge carries the list of mutations that happened on that branch.
"""

import numpy as np
import pandas as pd
from lineage.distance import build_distance_matrix, combined_distance
from lineage.mutations import get_branch_mutations, mutations_to_label, mutations_to_records


class TreeNode:
    def __init__(self, node_id, cell_data: dict):
        self.node_id = node_id
        self.cell_data = cell_data          # original row as dict
        self.children = []                  # list of (TreeNode, mutations_list)
        self.parent = None

    @property
    def label(self):
        cid = self.cell_data.get("cell_id", self.node_id)
        ct  = self.cell_data.get("cluster_annotated", "")
        iso = self.cell_data.get("c_call", "")
        mu  = self.cell_data.get("mu_count_h", "")
        return f"{cid}\n{ct} | {iso}\nSHM-H:{mu}"

    def is_leaf(self):
        return len(self.children) == 0


def build_tree(clone_df: pd.DataFrame, germline: dict) -> tuple:
    """
    Build a lineage tree for the clone.

    Parameters
    ----------
    clone_df  : DataFrame of all cells in the clone
    germline  : dict representing the germline root (from germline.py)

    Returns
    -------
    root      : TreeNode (the germline node)
    all_nodes : list of all TreeNode objects
    all_mutations : list of flat mutation records (for CSV export)
    """
    # Build list of all sequence dicts: germline first, then cells
    rows = [germline] + clone_df.to_dict(orient="records")

    # Deduplicate by VDJ_sequence_H + VDJ_sequence_L
    seen = {}
    unique_rows = []
    for r in rows:
        key = str(r["VDJ_sequence_H"]) + "|" + str(r.get("VDJ_sequence_L", ""))
        if key not in seen:
            seen[key] = r
            unique_rows.append(r)

    matrix, cell_ids = build_distance_matrix(unique_rows)
    id_to_row = {r["cell_id"]: r for r in unique_rows}

    # Create one node per unique sequence
    nodes = {cid: TreeNode(cid, id_to_row[cid]) for cid in cell_ids}

    # UPGMA clustering
    root = _upgma(matrix, cell_ids, nodes, unique_rows)

    # Collect all mutation records for CSV
    all_mutations = []
    _collect_mutations(root, all_mutations)

    all_nodes = list(nodes.values())
    return root, all_nodes, all_mutations


def _upgma(matrix: np.ndarray, cell_ids: list, nodes: dict, rows: list) -> TreeNode:
    """
    UPGMA implementation. Returns the root TreeNode.
    The germline node (cell_id == 'GERMLINE') is pinned as the root.
    """
    n = len(cell_ids)
    active = list(range(n))
    id_map = {i: cell_ids[i] for i in range(n)}
    row_map = {i: rows[i] for i in range(n)}
    cluster_size = {i: 1 for i in range(n)}
    dist = matrix.copy()

    while len(active) > 1:
        # Find closest pair
        min_d = np.inf
        pair = (0, 1)
        for i in range(len(active)):
            for j in range(i + 1, len(active)):
                ii, jj = active[i], active[j]
                if dist[ii, jj] < min_d:
                    min_d = dist[ii, jj]
                    pair = (ii, jj)

        i, j = pair
        # The node with lower mu_total is the parent
        mu_i = row_map[i].get("mu_count_h", 0) + row_map[i].get("mu_count_l", 0)
        mu_j = row_map[j].get("mu_count_h", 0) + row_map[j].get("mu_count_l", 0)

        if mu_i <= mu_j:
            parent_idx, child_idx = i, j
        else:
            parent_idx, child_idx = j, i

        parent_id = id_map[parent_idx]
        child_id  = id_map[child_idx]

        parent_node = nodes[parent_id]
        child_node  = nodes[child_id]

        if child_node.parent is None:
            muts = get_branch_mutations(row_map[parent_idx], row_map[child_idx])
            parent_node.children.append((child_node, muts))
            child_node.parent = parent_node

        # Update distances (UPGMA average)
        sz_i = cluster_size[parent_idx]
        sz_j = cluster_size[child_idx]
        for k in active:
            if k == i or k == j:
                continue
            new_d = (sz_i * dist[parent_idx, k] + sz_j * dist[child_idx, k]) / (sz_i + sz_j)
            dist[parent_idx, k] = new_d
            dist[k, parent_idx] = new_d

        cluster_size[parent_idx] = sz_i + sz_j
        active.remove(child_idx)

    # The remaining active node is the root
    root_idx = active[0]
    root_id  = id_map[root_idx]
    root     = nodes[root_id]

    # Re-root toward germline if needed
    germline_node = nodes.get("GERMLINE")
    if germline_node and root_id != "GERMLINE":
        root = _reroot(root, germline_node)

    return root


def _reroot(current_root: TreeNode, target: TreeNode) -> TreeNode:
    """Make 'target' the new root by reversing edges along the path."""
    path = _find_path(current_root, target)
    if not path:
        return current_root

    for i in range(len(path) - 1):
        parent = path[i]
        child  = path[i + 1]
        parent.children = [(c, m) for c, m in parent.children if c.node_id != child.node_id]
        muts = get_branch_mutations(child.cell_data, parent.cell_data)
        child.children.append((parent, muts))
        parent.parent = child

    target.parent = None
    return target


def _find_path(root: TreeNode, target: TreeNode) -> list:
    """DFS to find path from root to target."""
    if root.node_id == target.node_id:
        return [root]
    for child, _ in root.children:
        path = _find_path(child, target)
        if path:
            return [root] + path
    return []


def _collect_mutations(node: TreeNode, records: list):
    """Recursively collect mutation records from all branches."""
    for child, muts in node.children:
        records += mutations_to_records(node.cell_data["cell_id"],
                                        child.cell_data["cell_id"],
                                        muts)
        _collect_mutations(child, records)


def get_mutation_label(parent: TreeNode, child_muts: list) -> str:
    return mutations_to_label(child_muts)
