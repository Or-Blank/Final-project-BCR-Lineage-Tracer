import numpy as np
import pandas as pd
from lineage.distance import build_distance_matrix
from lineage.mutations import get_branch_mutations, mutations_to_records, mutations_to_label


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
    id_to_row = {r["cell_id"]: r for r in rows}

    # Build NJ tree
    root = neighbor_joining(matrix, cell_ids, id_to_row)

    # Collect mutation records
    all_mutations = []
    _collect_mutations(root, all_mutations)

    # Return all nodes
    all_nodes = _gather_nodes(root)
    return root, all_nodes, all_mutations


def neighbor_joining(D, labels, row_map):
    """
    Classic Neighbor-Joining implementation.
    D: distance matrix
    labels: list of node IDs
    row_map: mapping from node ID to cell data
    """
    D = D.astype(float)
    nodes = {lab: TreeNode(lab, row_map[lab]) for lab in labels}

    while len(labels) > 2:
        n = len(labels)
        total_dist = {i: np.sum(D[i]) for i in range(n)}

        # Build Q-matrix
        Q = np.zeros_like(D)
        for i in range(n):
            for j in range(n):
                if i != j:
                    Q[i, j] = (n - 2) * D[i, j] - total_dist[i] - total_dist[j]

        # Find minimum Q
        i, j = np.unravel_index(np.argmin(Q), Q.shape)
        if i == j:
            break

        li, lj = labels[i], labels[j]

        # Create new internal node
        new_label = f"NJ_internal_{len(nodes)}"
        new_node = TreeNode(new_label, {"cell_id": new_label})
        nodes[new_label] = new_node

        # Attach children
        muts_i = get_branch_mutations(row_map[li], row_map[lj])
        muts_j = get_branch_mutations(row_map[lj], row_map[li])

        new_node.children.append((nodes[li], muts_i))
        new_node.children.append((nodes[j], muts_j))
        nodes[li].parent = new_node
        nodes[j].parent = new_node

        # Compute distances to new node
        new_row = []
        for k in range(n):
            if k != i and k != j:
                d = (D[i, k] + D[j, k] - D[i, j]) / 2
                new_row.append(d)

        # Build new matrix
        new_D = np.zeros((n - 1, n - 1))
        new_labels = [lab for idx, lab in enumerate(labels) if idx not in (i, j)]
        new_labels.append(new_label)

        # Fill matrix
        idx_map = {lab: idx for idx, lab in enumerate(new_labels)}

        for a in range(n - 1):
            for b in range(n - 1):
                if new_labels[a] == new_label:
                    if new_labels[b] == new_label:
                        new_D[a, b] = 0
                    else:
                        old_idx = labels.index(new_labels[b])
                        new_D[a, b] = new_row[old_idx if old_idx < max(i, j) else old_idx - 1]
                elif new_labels[b] == new_label:
                    old_idx = labels.index(new_labels[a])
                    new_D[a, b] = new_row[old_idx if old_idx < max(i, j) else old_idx - 1]
                else:
                    old_a = labels.index(new_labels[a])
                    old_b = labels.index(new_labels[b])
                    new_D[a, b] = D[old_a, old_b]

        D = new_D
        labels = new_labels

    # Final join
    a, b = labels
    root = TreeNode("ROOT", {"cell_id": "ROOT"})
    muts_a = get_branch_mutations(row_map[a], row_map[b])
    muts_b = get_branch_mutations(row_map[b], row_map[a])
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
    nodes = []
    stack = [root]
    while stack:
        n = stack.pop()
        nodes.append(n)
        for c, _ in n.children:
            stack.append(c)
    return nodes
