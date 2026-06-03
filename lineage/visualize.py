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
    "#F72585", "#4CC9F0", "#7209B7", "#3A0CA3", "#4362EEF8",
    "#B5179E", "#560BAD", "#480CA8", "#3F37C9", "#4895EF",
]


def _build_color_map(root: TreeNode) -> dict:
    """
    Collect all unique cell types from the tree and assign each a color
    dynamically from the palette. Germline always gets grey.
    Returns a dict: cell_type_string → hex_color.
    """
    seen = set()
    _collect_cell_types(root, seen)

    # Clean: remove None / NaN, convert to str
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
        ax.plot([x, cx], [y, y], color="#CCCCCC", linewidth=1.2, zorder=1)
        ax.plot([cx, cx], [y, cy], color="#CCCCCC", linewidth=1.2, zorder=1)

        # Mutation label
        mut_label = mutations_to_label(muts, max_show=3)
        if mut_label:
            mid_y = (y + cy) / 2
            ax.text(cx + 0.05, mid_y, mut_label,
                    fontsize=5.5, color="#333333", va="center", ha="left",
                    style="italic",
                    bbox=dict(boxstyle="round,pad=0.1", fc="white", alpha=0.6, lw=0))

        _draw_tree(ax, child, positions, color_map)

    # Node circle
    cell_type = node.cell_data.get("cluster_annotated", "")
    color = _get_color(cell_type, color_map)
    circle = plt.Circle((x, y), 0.18, color=color, zorder=3,
                         linewidth=1, edgecolor="white")
    ax.add_patch(circle)

    # Node label: short cell barcode + cell type + SHM count
    cid  = str(node.cell_data.get("cell_id", ""))[-8:]   # last 8 chars only
    ct   = str(cell_type)[:12] if cell_type else ""
    mu_h = node.cell_data.get("mu_count_h", "")
    ax.text(x, y - 0.28, f"{cid}\n{ct}\nSHM:{mu_h}",
            fontsize=4.5, ha="center", va="top", color="#111111")


def render_tree(root: TreeNode, output_path: str, clone_id: str = ""):
    """
    Render the lineage tree to a PNG file.

    Parameters
    ----------
    root        : TreeNode  (germline root)
    output_path : str       path for the output PNG
    clone_id    : str       used in the figure title
    """
    # Build color map from whatever cell types are in THIS tree
    color_map = _build_color_map(root)

    # Layout
    counter = [0]
    positions = _assign_positions(root, depth=0, counter=counter)
    n_leaves = counter[0]

    fig_w = max(12, (max(x for x, y in positions.values()) + 2) * 2.5)
    fig_h = max(6,  n_leaves * 0.55)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    _draw_tree(ax, root, positions, color_map)

    all_x = [x for x, y in positions.values()]
    all_y = [y for x, y in positions.values()]
    ax.set_xlim(min(all_x) - 0.5, max(all_x) + 1.5)
    ax.set_ylim(min(all_y) - 0.8, max(all_y) + 0.8)
    ax.axis("off")

    # Legend — one patch per cell type actually present
    legend_patches = [
        mpatches.Patch(color=color, label=ct)
        for ct, color in sorted(color_map.items())
    ]
    ax.legend(handles=legend_patches, loc="upper right",
              fontsize=7, framealpha=0.9, title="Cell type", title_fontsize=7)

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
