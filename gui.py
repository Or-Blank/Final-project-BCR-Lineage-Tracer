"""
gui.py
Graphical User Interface for BCR Lineage Tracer.
Built with tkinter (included in standard Python).

Layout:
  ┌──────────────────────────────────────┐
  │  File Upload Panel                   │
  │  Clone Selector (dropdown)           │
  │  Settings (min cells, output dir)    │
  │  [Run] button                        │
  ├──────────────────────────────────────┤
  │  Log / progress panel                │
  ├──────────────────────────────────────┤
  │  Tree image preview (after run)      │
  └──────────────────────────────────────┘
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from pathlib import Path
import threading
import os
import sys

# Allow running gui.py from the project root
sys.path.insert(0, str(Path(__file__).parent))

from lineage.loader import load_file, list_clones, get_clone, get_largest_clone
from lineage.germline import infer_germline
from lineage.tree import build_tree
from lineage.visualize import render_tree
import pandas as pd


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("BCR Lineage Tracer")
        self.geometry("820x700")
        self.resizable(True, True)
        self.configure(bg="#F5F5F5")

        self.df = None
        self.clone_list = []
        self.last_tree_path = None

        self._build_ui()

    # ── UI Construction ───────────────────────────────────────────────────────

    def _build_ui(self):
        pad = dict(padx=12, pady=6)

        # ── Title ──
        title_frame = tk.Frame(self, bg="#1D3557")
        title_frame.pack(fill="x")
        tk.Label(title_frame, text="BCR Lineage Tracer",
                 font=("Helvetica", 16, "bold"), fg="white", bg="#1D3557",
                 pady=10).pack()
        tk.Label(title_frame,
                 text="B Cell Clonal Lineage Tree Maker",
                 font=("Helvetica", 10), fg="#A8DADC", bg="#1D3557",
                 pady=2).pack()

        # ── File upload ──
        file_frame = ttk.LabelFrame(self, text="1. Input File", padding=10)
        file_frame.pack(fill="x", **pad)

        self.file_var = tk.StringVar(value="No file selected")
        tk.Label(file_frame, textvariable=self.file_var,
                 fg="#555", font=("Helvetica", 9)).pack(side="left", expand=True)
        ttk.Button(file_frame, text="Browse...",
                   command=self._browse_file).pack(side="right")

        # ── Clone selector ──
        clone_frame = ttk.LabelFrame(self, text="2. Select Clone", padding=10)
        clone_frame.pack(fill="x", **pad)

        tk.Label(clone_frame, text="Clone ID:", font=("Helvetica", 9)).pack(side="left")
        self.clone_var = tk.StringVar(value="— load a file first —")
        self.clone_combo = ttk.Combobox(clone_frame, textvariable=self.clone_var,
                                        state="disabled", width=45)
        self.clone_combo.pack(side="left", padx=8)

        self.clone_info = tk.StringVar(value="")
        tk.Label(clone_frame, textvariable=self.clone_info,
                 font=("Helvetica", 8), fg="#666").pack(side="left", padx=4)
        self.clone_combo.bind("<<ComboboxSelected>>", self._on_clone_selected)

        # ── Settings ──
        settings_frame = ttk.LabelFrame(self, text="3. Settings", padding=10)
        settings_frame.pack(fill="x", **pad)

        tk.Label(settings_frame, text="Output directory:", font=("Helvetica", 9)).grid(
            row=0, column=0, sticky="w")
        self.outdir_var = tk.StringVar(value=str(Path.cwd() / "results"))
        tk.Entry(settings_frame, textvariable=self.outdir_var, width=40).grid(
            row=0, column=1, padx=6, sticky="w")
        ttk.Button(settings_frame, text="Browse",
                   command=self._browse_outdir).grid(row=0, column=2, padx=4)

        tk.Label(settings_frame, text="Run all clones with ≥",
                 font=("Helvetica", 9)).grid(row=1, column=0, sticky="w", pady=4)
        self.min_cells_var = tk.IntVar(value=5)
        tk.Spinbox(settings_frame, from_=2, to=500, textvariable=self.min_cells_var,
                   width=5).grid(row=1, column=1, sticky="w", padx=6)
        tk.Label(settings_frame, text="cells", font=("Helvetica", 9)).grid(
            row=1, column=1, sticky="w", padx=55)

        self.all_clones_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(settings_frame, text="Run all clones (not just selected)",
                        variable=self.all_clones_var).grid(
            row=2, column=0, columnspan=3, sticky="w")

        # ── Run button ──
        btn_frame = tk.Frame(self, bg="#F5F5F5")
        btn_frame.pack(fill="x", padx=12, pady=4)
        self.run_btn = ttk.Button(btn_frame, text="▶  Run Analysis",
                                  command=self._run, style="Accent.TButton")
        self.run_btn.pack(side="left")
        self.open_btn = ttk.Button(btn_frame, text="📂  Open Results Folder",
                                   command=self._open_results, state="disabled")
        self.open_btn.pack(side="left", padx=10)

        # ── Log panel ──
        log_frame = ttk.LabelFrame(self, text="Log", padding=6)
        log_frame.pack(fill="both", expand=False, **pad)
        self.log_text = tk.Text(log_frame, height=8, state="disabled",
                                font=("Courier", 8), bg="#1e1e1e", fg="#d4d4d4",
                                relief="flat")
        scrollbar = ttk.Scrollbar(log_frame, command=self.log_text.yview)
        self.log_text.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side="right", fill="y")
        self.log_text.pack(fill="both", expand=True)

        # ── Tree preview label ──
        self.preview_label = tk.Label(self, text="Tree preview will appear here after running.",
                                      bg="#F5F5F5", fg="#888", font=("Helvetica", 9))
        self.preview_label.pack(pady=4)

    # ── Actions ───────────────────────────────────────────────────────────────

    def _browse_file(self):
        path = filedialog.askopenfilename(
            title="Select input file",
            filetypes=[("Supported files", "*.xlsx *.csv *.tsv"), ("All files", "*.*")]
        )
        if not path:
            return
        self.file_var.set(path)
        self._log(f"Loading {Path(path).name} ...")
        try:
            self.df = load_file(path)
            self._log(f"  Loaded {len(self.df)} cells.")
            self._populate_clones()
        except Exception as e:
            messagebox.showerror("Load error", str(e))
            self._log(f"ERROR: {e}")

    def _populate_clones(self):
        if self.df is None:
            return
        clones_df = list_clones(self.df, min_cells=1)
        self.clone_list = clones_df["clone_id"].tolist()
        self.clone_combo["values"] = self.clone_list
        self.clone_combo.config(state="readonly")

        # Default: largest clone
        largest = get_largest_clone(self.df)
        self.clone_var.set(largest)
        self._on_clone_selected()

    def _on_clone_selected(self, event=None):
        if self.df is None:
            return
        cid = self.clone_var.get()
        try:
            sub = self.df[self.df["clone_id"] == cid]
            n   = len(sub)
            iso = ", ".join(sorted(sub["c_call"].dropna().unique()))
            cts = ", ".join(sorted(sub["cluster_annotated"].dropna().unique()))
            self.clone_info.set(f"{n} cells | {iso} | {cts}")
        except Exception:
            pass

    def _browse_outdir(self):
        d = filedialog.askdirectory(title="Select output directory")
        if d:
            self.outdir_var.set(d)

    def _run(self):
        if self.df is None:
            messagebox.showwarning("No file", "Please load an input file first.")
            return
        self.run_btn.config(state="disabled")
        self._log("─" * 50)
        thread = threading.Thread(target=self._run_thread, daemon=True)
        thread.start()

    def _run_thread(self):
        try:
            outdir = self.outdir_var.get()
            os.makedirs(outdir, exist_ok=True)

            if self.all_clones_var.get():
                clones_df = list_clones(self.df, min_cells=self.min_cells_var.get())
                self._log(f"Running {len(clones_df)} clones...")
                last_path = None
                for _, row in clones_df.iterrows():
                    clone_df = get_clone(self.df, row["clone_id"])
                    last_path = self._run_one(clone_df, row["clone_id"], outdir)
                if last_path:
                    self.last_tree_path = last_path
            else:
                clone_id = self.clone_var.get()
                clone_df = get_clone(self.df, clone_id)
                self.last_tree_path = self._run_one(clone_df, clone_id, outdir)

            self._log("✓ Analysis complete!")
            self.after(0, lambda: self.open_btn.config(state="normal"))
            if self.last_tree_path:
                self.after(0, lambda: self._show_preview(self.last_tree_path))
        except Exception as e:
            self._log(f"ERROR: {e}")
            import traceback
            self._log(traceback.format_exc())
        finally:
            self.after(0, lambda: self.run_btn.config(state="normal"))

    def _run_one(self, clone_df, clone_id, outdir) -> str:
        """Run pipeline for one clone and return path to tree image."""
        from pathlib import Path as P
        clone_outdir = P(outdir) / clone_id.replace("/", "_")
        clone_outdir.mkdir(parents=True, exist_ok=True)

        self._log(f"  Clone: {clone_id} ({len(clone_df)} cells)")
        germline = infer_germline(clone_df)
        root, _, all_mutations = build_tree(clone_df, germline)

        tree_path = str(clone_outdir / "lineage_tree.png")
        render_tree(root, tree_path, clone_id=clone_id)
        self._log(f"    Tree → {tree_path}")

        if all_mutations:
            mut_df = pd.DataFrame(all_mutations)
            mut_path = clone_outdir / "mutations_per_branch.xlsx"
            mut_df.to_excel(str(mut_path), index=False)
            self._log(f"    Mutations → {mut_path}")

        summary = {
            "clone_id":       clone_id,
            "total_cells":    len(clone_df),
            "mean_mu_h":      round(clone_df["mu_count_h"].mean(), 2),
            "mean_mu_l":      round(clone_df["mu_count_l"].mean(), 2),
            "isotypes":       ", ".join(sorted(clone_df["c_call"].dropna().unique())),
            "cell_types":     ", ".join(sorted(clone_df["cluster_annotated"].dropna().unique())),
            "timepoints":     ", ".join(sorted(clone_df["sample_id"].dropna().unique())),
            "total_mutations": len(all_mutations),
        }
        pd.DataFrame([summary]).to_excel(
            str(clone_outdir / "clone_summary.xlsx"), index=False)

        return tree_path

    def _show_preview(self, image_path: str):
        """Show a small preview of the tree image."""
        try:
            from PIL import Image, ImageTk
            img = Image.open(image_path)
            img.thumbnail((780, 200))
            photo = ImageTk.PhotoImage(img)
            self.preview_label.config(image=photo, text="")
            self.preview_label.image = photo  # keep reference
        except ImportError:
            self.preview_label.config(
                text=f"Tree saved to: {image_path}\n(Install Pillow for in-app preview: pip install Pillow)"
            )

    def _open_results(self):
        outdir = self.outdir_var.get()
        if os.path.exists(outdir):
            if sys.platform == "darwin":
                os.system(f"open '{outdir}'")
            elif sys.platform == "win32":
                os.startfile(outdir)
            else:
                os.system(f"xdg-open '{outdir}'")

    def _log(self, message: str):
        def _insert():
            self.log_text.config(state="normal")
            self.log_text.insert("end", message + "\n")
            self.log_text.see("end")
            self.log_text.config(state="disabled")
        self.after(0, _insert)


if __name__ == "__main__":
    app = App()
    app.mainloop()
