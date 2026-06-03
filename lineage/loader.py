"""
loader.py
Reads the input file (.xlsx / .csv / .tsv) and returns a filtered
DataFrame ready for lineage analysis.
"""

import pandas as pd
from pathlib import Path

REQUIRED_COLUMNS = [
    "cell_id", "clone_id", "clone_count",
    "VDJ_sequence_H", "VDJ_sequence_L",
    "VDJ_aa_sequence_H", "VDJ_aa_sequence_L",
    "v_call_h", "j_call_h", "v_call_l", "j_call_l",
    "mu_count_h", "mu_count_l",
    "cluster_annotated", "c_call", "sample_id",
]

OPTIONAL_COLUMNS = ["d_call_h"]


def load_file(path: str) -> pd.DataFrame:
    """Load input file into a DataFrame. Supports .xlsx, .csv, .tsv."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    suffix = p.suffix.lower()
    if suffix == ".xlsx":
        df = pd.read_excel(path, engine="openpyxl")
    elif suffix == ".csv":
        df = pd.read_csv(path)
    elif suffix == ".tsv":
        df = pd.read_csv(path, sep="\t")
    else:
        raise ValueError(f"Unsupported file format: {suffix}. Use .xlsx, .csv, or .tsv")

    _validate_columns(df, path)
    df = _clean(df)
    return df


def _validate_columns(df: pd.DataFrame, path: str):
    """Raise an informative error if required columns are missing."""
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns in {path}:\n  " + "\n  ".join(missing)
        )


def _clean(df: pd.DataFrame) -> pd.DataFrame:
    """Drop rows with missing sequence data and normalise types."""
    # Only require heavy chain sequence — light chain may be absent (heavy-only files)
    df = df.dropna(subset=["VDJ_sequence_H"])
    df["VDJ_sequence_H"]    = df["VDJ_sequence_H"].fillna("").astype(str)
    df["VDJ_sequence_L"]    = df["VDJ_sequence_L"].fillna("").astype(str)
    df["VDJ_aa_sequence_H"] = df["VDJ_aa_sequence_H"].fillna("").astype(str)
    df["VDJ_aa_sequence_L"] = df["VDJ_aa_sequence_L"].fillna("").astype(str)
    df["mu_count_h"] = pd.to_numeric(df["mu_count_h"], errors="coerce").fillna(0).astype(int)
    df["mu_count_l"] = pd.to_numeric(df["mu_count_l"], errors="coerce").fillna(0).astype(int)
    df["clone_count"] = pd.to_numeric(df["clone_count"], errors="coerce").fillna(1).astype(int)
    df = df.reset_index(drop=True)
    return df


def list_clones(df: pd.DataFrame, min_cells: int = 2) -> pd.DataFrame:
    """
    Return a summary DataFrame of clones with at least min_cells cells,
    sorted by size descending.
    """
    summary = (
        df.groupby("clone_id")
        .agg(
            cell_count=("cell_id", "count"),
            isotypes=("c_call", lambda x: ", ".join(sorted(x.dropna().unique()))),
            cell_types=("cluster_annotated", lambda x: ", ".join(sorted(x.dropna().unique()))),
            timepoints=("sample_id", lambda x: ", ".join(sorted(x.dropna().unique()))),
            mean_mu_h=("mu_count_h", "mean"),
            mean_mu_l=("mu_count_l", "mean"),
        )
        .reset_index()
        .query("cell_count >= @min_cells")
        .sort_values("cell_count", ascending=False)
    )
    return summary


def get_clone(df: pd.DataFrame, clone_id: str) -> pd.DataFrame:
    """Return all cells belonging to a given clone_id."""
    sub = df[df["clone_id"] == clone_id].copy()
    if sub.empty:
        raise ValueError(f"Clone '{clone_id}' not found in data.")
    return sub.reset_index(drop=True)


def get_largest_clone(df: pd.DataFrame) -> str:
    """Return the clone_id of the clone with the most cells (minimum 2)."""
    counts = df.groupby("clone_id")["cell_id"].count()
    counts = counts[counts >= 2]
    if counts.empty:
        raise ValueError("No clones with 2 or more cells found in the data.")
    return counts.idxmax()
