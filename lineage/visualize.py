"""
visualize.py
Renders the lineage tree as a PNG image using matplotlib.

Color scheme: fully dynamic — colors are auto-assigned from a qualitative
palette based on whatever cell types appear in the data. No hardcoded names.
The germline root is always drawn in grey.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from lineage.tree import TreeNode
from lineage.mutations import mutations_to_label

# ── Fixed colors ──────────────────────────────────────────────────────────────
GERMLINE_COLOR = "#888888"

# Qualitative palette — distinct, colorblind-friendly, enough for ~20 types
_PALETTE = [
    "#E63946", "#457B9D", "#2A9D8F", "#E9C46A", "#F4A261",
    "#264653", "#A8DADC", "#6A4C93", "#1D3557", "#80B918",
    "#F72585", "#4CC9F0", "#7209B7", "#3A0CA3", "#4361EE",
    "#B5179E", "#560BAD", "#480CA8", "#3F37C9", "#4895EF",
]

# Node circle radius (smaller than before)
NODE_RADIUS = 0.10


def _build_color_map(root: TreeNode) -> dict:
    seen = set()
    _collect_cell_types(root, seen)
    clean = sorted(
        [str(ct) for ct in seen if ct and str(ct) not in ("nan", "None", "")],
        key=str.lower
    )
    color_map = {"Germline": GERMLINE_COLOR}
    palette_cycle = _PALETTE * ((len(clean) // len(_PALETTE)) + 1)
    for i, ct in enumerate(clean):
        if ct == "Germline":
            continue
        color_map[ct] = palette_cycle[i]
    return color_map


def _get_color(cell_type, color_map: dict) -> str:
    ct = str(cell_type) if cell_type else ""
    return color_map.get(ct, GERMLINE_COLOR)


# ── Layout ────────────────────────────────────────────────────────────────────

def _assign_positions(node: TreeNode, depth: int, counter: list) -> dict:
    positions = {}
    if node.is_leaf():
        y = counter[0]
        counter[0] += 1
        positions[node.node_id] = (depth, y)
    else:
        child_ys = []
        for child, _ in node.children:
            child_pos = _assign_positions(child, depth + 1, counter)
            positions.update(child_pos)
            child_ys.append(child_pos[child.node_id][1])
        y = np.mean(child_ys)
        positions[node.node_id] = (depth, y)
    return positions


# ── Drawing ───────────────────────────────────────────────────────────────────

def _draw_tree(ax, node: TreeNode, positions: dict, color_map: dict):
    x, y = positions[node.node_id]

    for child, muts in node.children:
        cx, cy = positions[child.node_id]

        # L-shaped connector
        ax.plot([x, cx], [y, y], color="#CCCCCC", linewidth=1.0, zorder=1)
        ax.plot([cx, cx], [y, cy], color="#CCCCCC", linewidth=1.0, zorder=1)

        _draw_tree(ax, child, positions, color_map)

    # ── Node circle (small, crisp) ──
    cell_type = node.cell_data.get("cluster_annotated", "")
    color = _get_color(cell_type, color_map)
    circle = plt.Circle((x, y), NODE_RADIUS, color=color, zorder=3,
                         linewidth=0.8, edgecolor="white")
    ax.add_patch(circle)

    # ── Label to the RIGHT of the dot ──
    mu_h  = node.cell_data.get("mu_count_h", 0)
    mu_l  = node.cell_data.get("mu_count_l", 0)
    ct    = str(cell_type)[:14] if cell_type and str(cell_type) not in ("nan", "") else ""
    iso   = str(node.cell_data.get("c_call", ""))

    # Build the per-node mutation label (aa changes on the branch FROM parent)
    # These are stored on the child node for display purposes
    branch_label = node.cell_data.get("_branch_label", "")

    # Compose multi-line label
    lines = []
    if ct:
        lines.append(ct)
    if iso and iso not in ("nan", ""):
        lines.append(iso)
    lines.append(f"SHM: {mu_h} (H) / {mu_l} (L)")
    if branch_label:
        lines.append(branch_label)

    label_text = "\n".join(lines)
    ax.text(x + NODE_RADIUS + 0.06, y, label_text,
            fontsize=4.8, ha="left", va="center", color="#111111",
            linespacing=1.4)


def _attach_branch_labels(node: TreeNode):
    """
    Walk the tree and store the mutation label for each branch
    directly on the child node's cell_data under '_branch_label',
    so _draw_tree can access it when drawing the node.
    """
    for child, muts in node.children:
        label = mutations_to_label(muts, max_show=4)
        child.cell_data["_branch_label"] = label
        _attach_branch_labels(child)


def render_tree(root: TreeNode, output_path: str, clone_id: str = ""):
    """
    Render the lineage tree to a PNG file.

    Parameters
    ----------
    root        : TreeNode  (germline root)
    output_path : str       path for the output PNG
    clone_id    : str       used in the figure title
    """
    color_map = _build_color_map(root)

    # Attach branch mutation labels to each child node before drawing
    root.cell_data["_branch_label"] = ""   # germline has no incoming branch
    _attach_branch_labels(root)

    # Layout
    counter = [0]
    positions = _assign_positions(root, depth=0, counter=counter)
    n_leaves = counter[0]

    fig_w = max(14, (max(x for x, y in positions.values()) + 2) * 3.2)
    fig_h = max(6,  n_leaves * 0.62)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    _draw_tree(ax, root, positions, color_map)

    all_x = [x for x, y in positions.values()]
    all_y = [y for x, y in positions.values()]
    ax.set_xlim(min(all_x) - 0.5, max(all_x) + 2.8)
    ax.set_ylim(min(all_y) - 0.8, max(all_y) + 0.8)
    ax.axis("off")

    # ── Legend: cell type colors ──
    legend_patches = [
        mpatches.Patch(color=color, label=ct)
        for ct, color in sorted(color_map.items())
    ]

    # ── Notation legend (below cell type legend) ──
    notation_lines = [
        "─── Notation ───",
        "H:XnY  mutation in Heavy chain",
        "L:XnY  mutation in Light chain",
        "  X = original AA,  n = position",
        "  Y = new AA",
        "SHM  somatic hypermutation count",
    ]
    notation_text = "\n".join(notation_lines)

    cell_type_legend = ax.legend(
        handles=legend_patches,
        loc="upper right",
        fontsize=7,
        framealpha=0.92,
        title="Cell type",
        title_fontsize=7.5,
        borderpad=0.8,
    )
    ax.add_artist(cell_type_legend)

    # Add notation box below the legend using a text box
    ax.text(
        1.01, 0.02, notation_text,
        transform=ax.transAxes,
        fontsize=6.2,
        va="bottom", ha="left",
        family="monospace",
        bbox=dict(boxstyle="round,pad=0.5", fc="#F8F8F8", ec="#CCCCCC", lw=0.8),
    )

    title = f"BCR Lineage Tree — {clone_id}" if clone_id else "BCR Lineage Tree"
    ax.set_title(title, fontsize=11, pad=10)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Tree image saved → {output_path}")


def _collect_cell_types(node: TreeNode, seen: set):
    seen.add(node.cell_data.get("cluster_annotated", ""))
    for child, _ in node.children:
        _collect_cell_types(child, seen)
