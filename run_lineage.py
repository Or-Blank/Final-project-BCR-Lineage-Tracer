"""
run_lineage.py
Command-line entry point for BCR Lineage Tracer.

Usage:
    python run_lineage.py --input your_file.xlsx
    python run_lineage.py --input your_file.xlsx --clone PT7_sJRC_138_G1G2
    python run_lineage.py --input your_file.xlsx --all-clones --min-cells 5
"""

import argparse
import os
import sys
import pandas as pd
from pathlib import Path

from lineage.loader import load_file, list_clones, get_clone, get_largest_clone
from lineage.germline import infer_germline
from lineage.tree import build_tree
from lineage.visualize import render_tree


def run_clone(clone_df: pd.DataFrame, clone_id: str, outdir: str):
    """Run the full pipeline for a single clone and save outputs."""
    clone_outdir = Path(outdir) / clone_id.replace("/", "_")
    clone_outdir.mkdir(parents=True, exist_ok=True)

    print(f"\n── Clone: {clone_id}  ({len(clone_df)} cells) ──")

    # 1. Infer germline root
    print("  Inferring germline ancestor...")
    germline = infer_germline(clone_df)

    # 2. Build tree
    print("  Building lineage tree...")
    root, all_nodes, all_mutations = build_tree(clone_df, germline)

    # 3. Render tree image
    tree_path = str(clone_outdir / "lineage_tree.png")
    print("  Rendering tree...")
    render_tree(root, tree_path, clone_id=clone_id)

    # 4. Save mutation table
    if all_mutations:
        mut_df = pd.DataFrame(all_mutations)
        mut_path = clone_outdir / "mutations_per_branch.xlsx"
        mut_df.to_excel(str(mut_path), index=False)
        print(f"  Mutation table saved → {mut_path}")

    # 5. Save clone summary
    summary = {
        "clone_id":         clone_id,
        "total_cells":      len(clone_df),
        "unique_sequences": len(set(clone_df["VDJ_sequence_H"] + "|" + clone_df["VDJ_sequence_L"])),
        "mean_mu_h":        round(clone_df["mu_count_h"].mean(), 2),
        "mean_mu_l":        round(clone_df["mu_count_l"].mean(), 2),
        "isotypes":         ", ".join(sorted(clone_df["c_call"].dropna().unique())),
        "cell_types":       ", ".join(sorted(clone_df["cluster_annotated"].dropna().unique())),
        "timepoints":       ", ".join(sorted(clone_df["sample_id"].dropna().unique())),
        "total_mutations":  len(all_mutations),
    }
    summary_df = pd.DataFrame([summary])
    summary_path = clone_outdir / "clone_summary.xlsx"
    summary_df.to_excel(str(summary_path), index=False)
    print(f"  Clone summary saved → {summary_path}")

    return clone_outdir


def main():
    parser = argparse.ArgumentParser(
        description="BCR Lineage Tracer — reconstruct B cell clonal lineage trees"
    )
    parser.add_argument("--input",      required=True,  help="Path to input .xlsx / .csv / .tsv file")
    parser.add_argument("--clone",      default=None,   help="Clone ID to analyse (default: largest clone)")
    parser.add_argument("--all-clones", action="store_true", help="Run on all clones >= --min-cells")
    parser.add_argument("--min-cells",  type=int, default=5,  help="Minimum clone size (default: 5)")
    parser.add_argument("--outdir",     default="results",    help="Output directory (default: results/)")
    args = parser.parse_args()

    # Load data
    print(f"Loading {args.input} ...")
    df = load_file(args.input)
    print(f"  Loaded {len(df)} cells.")

    os.makedirs(args.outdir, exist_ok=True)

    if args.all_clones:
        clones = list_clones(df, min_cells=args.min_cells)
        print(f"\nFound {len(clones)} clones with >= {args.min_cells} cells.")
        for _, row in clones.iterrows():
            clone_df = get_clone(df, row["clone_id"])
            run_clone(clone_df, row["clone_id"], args.outdir)
    else:
        clone_id = args.clone or get_largest_clone(df)
        if args.clone is None:
            print(f"No --clone specified. Using largest clone: {clone_id}")
        clone_df = get_clone(df, clone_id)
        run_clone(clone_df, clone_id, args.outdir)

    print(f"\nDone. Results saved to: {args.outdir}/")


if __name__ == "__main__":
    main()
