#!/usr/bin/env python3
"""Interactive ECD lineshape viewer — toggle .dat files, switch units, export."""

from __future__ import annotations

import ast
import itertools
import json
import os
import sys
import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)
from plot_class_functions import SpectrumData

STYLE_FILE      = os.path.join(_SCRIPT_DIR, "lineshapes.mplstyle")
FILE_PREFIX     = "lineshape-ECD-VH-internal-Qi-"
INTERACTIVE_DPI = 96

_COLORS     = ["#56B4E9", "#E69F00", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"]
_LINESTYLES = ["-", "--", "-.", ":", "-", "--", "-."]


def _get_smart_stem(fname: str) -> str:
    """Strip known prefix and .dat extension; return short label."""
    stem = os.path.splitext(fname)[0]
    if stem.startswith(FILE_PREFIX):
        stem = stem[len(FILE_PREFIX):]
    return stem


def _scan_folder(folder: str, color_it, ls_it) -> list:
    """Return records for all .dat files in folder, consuming colors/styles."""
    records = []
    try:
        fnames = sorted(f for f in os.listdir(folder) if f.endswith(".dat"))
    except FileNotFoundError:
        print(f"Warning: folder not found: {folder}", file=sys.stderr)
        return records
    for fname in fnames:
        path = os.path.join(folder, fname)
        stem = _get_smart_stem(fname)
        try:
            data = SpectrumData(path)
        except Exception as exc:
            print(f"Warning: could not load {path}: {exc}", file=sys.stderr)
            continue
        records.append({
            "path":      path,
            "stem":      stem,
            "data":      data,
            "color":     next(color_it),
            "linestyle": next(ls_it),
        })
    return records


# ---------------------------------------------------------------------------
# PlotState — pure data, no widgets
# ---------------------------------------------------------------------------

@dataclass
class PlotState:
    x_unit:        str             = "eV"
    title:         str             = ""
    xmin:          Optional[float] = None
    xmax:          Optional[float] = None
    ymin:          Optional[float] = None
    ymax:          Optional[float] = None
    # styles: [{"path": str, "enabled": bool}] — saved so config restores them
    styles:        list            = field(default_factory=list)
    style_preview: bool            = False
    style_export:  bool            = True
    # folders: [{"path": str, "enabled": bool}] — saved so config restores them
    folders:       list            = field(default_factory=list)
    # selected / labels keyed by *file path* to stay unique across folders
    selected:      list            = field(default_factory=list)
    labels:        dict            = field(default_factory=dict)


# ---------------------------------------------------------------------------
# PlotEngine — stateless renderer
# ---------------------------------------------------------------------------

class PlotEngine:
    def render(self, state: PlotState, records: list, fig,
               use_cycler: bool = False) -> None:
        fig.clear()
        ax = fig.add_subplot(111)

        selected_set = set(state.selected)
        for rec in records:
            if rec["path"] not in selected_set:
                continue
            data  = rec["data"]
            x     = data.energy if state.x_unit == "eV" else data.cm_inv
            label = state.labels.get(rec["path"], rec["stem"])
            if use_cycler:
                ax.plot(x, data.me, marker="None", label=label)
            else:
                ax.plot(x, data.me, marker="None",
                        label=label,
                        color=rec["color"],
                        linestyle=rec["linestyle"])

        ax.set_xlabel(r"Energy (eV)" if state.x_unit == "eV"
                      else r"$\tilde{\nu}\ (\mathrm{cm^{-1}})$")
        ax.set_ylabel(r"Lineshape (a.u.)")
        if state.title:
            ax.set_title(state.title)

        handles, _ = ax.get_legend_handles_labels()
        if handles:
            ax.legend()

        if state.xmin is not None and state.xmax is not None:
            ax.set_xlim(state.xmin, state.xmax)
        if state.ymin is not None and state.ymax is not None:
            ax.set_ylim(state.ymin, state.ymax)

        if hasattr(fig, "canvas"):
            fig.canvas.draw_idle()


# ---------------------------------------------------------------------------
# Config I/O
# ---------------------------------------------------------------------------

_CONFIG_FIELDS = (
    "x_unit", "title", "xmin", "xmax", "ymin", "ymax",
    "styles", "style_preview", "style_export",
    "folders", "selected", "labels",
)


def _export_config(state: PlotState, path: str) -> None:
    cfg = {f: getattr(state, f) for f in _CONFIG_FIELDS}
    Path(path).write_text(json.dumps(cfg, indent=2))


def _load_config(path: str) -> dict:
    return json.loads(Path(path).read_text())


# ---------------------------------------------------------------------------
# Script generation
# ---------------------------------------------------------------------------

def _generate_script(state: PlotState, records: list) -> str:
    """Generate a portable standalone .py that reproduces the current plot."""
    selected_set  = set(state.selected)
    selected_recs = [r for r in records if r["path"] in selected_set]

    x_expr      = "data.energy"      if state.x_unit == "eV" else "data.cm_inv"
    x_label_str = r"r'Energy (eV)'"  if state.x_unit == "eV" \
                  else r"r'$\tilde{\nu}\ (\mathrm{cm^{-1}})$'"

    lines = [
        "#!/usr/bin/env python3",
        '"""Snapshot — auto-generated by lineshape-dynamic-plotting.py',
        "",
        "Run from any directory.  If plot_class_functions.py is not found",
        "beside this script, pass its location with --scripts-dir.",
        '"""',
        "import sys, argparse",
        "from pathlib import Path",
        "",
        "_ap = argparse.ArgumentParser(description='Render ECD lineshape figure to PDF')",
        "_ap.add_argument('-sd', '--scripts-dir', default='.',",
        "    metavar='DIR',",
        "    help='directory that contains plot_class_functions.py (default: .)')",
        "_args = _ap.parse_args()",
        "",
        "import matplotlib",
        'matplotlib.use("pgf")',
        "import matplotlib.pyplot as plt",
        "",
        "# ---- Locate plot_class_functions.py ----",
        f"_FALLBACK = Path({repr(_SCRIPT_DIR)})",
        "_HERE = Path(__file__).parent",
        "_SCRIPTS_DIR = Path(_args.scripts_dir).resolve()",
        "for _d in [_HERE, _SCRIPTS_DIR, _FALLBACK]:",
        "    if (_d / 'plot_class_functions.py').exists():",
        "        sys.path.insert(0, str(_d))",
        "        break",
        "else:",
        "    print('Error: plot_class_functions.py not found.', file=sys.stderr)",
        "    print('Use: python snapshot.py --scripts-dir /path/to/scripts', file=sys.stderr)",
        "    sys.exit(1)",
        "",
        "from plot_class_functions import SpectrumData",
        "",
        "_active_styles = " + repr([s["path"] for s in state.styles if s.get("enabled", True)]),
        "if _active_styles:",
        "    plt.style.use(_active_styles)",
        "",
        "# ---- Data ----",
        "paths = [",
    ]
    for rec in selected_recs:
        lines.append(f"    {repr(rec['path'])},")
    lines += [
        "]",
        "labels = [",
    ]
    for rec in selected_recs:
        lines.append(f"    {repr(state.labels.get(rec['path'], rec['stem']))},")
    lines += [
        "]",
        "datasets = [SpectrumData(p) for p in paths]",
        "",
        "# ---- Plot ----",
        "fig, ax = plt.subplots()",
        "for data, label in zip(datasets, labels):",
        f"    x = {x_expr}",
        "    ax.plot(x, data.me, marker='None', label=label)",
        "",
        f"ax.set_xlabel({x_label_str})",
        "ax.set_ylabel(r'Lineshape (a.u.)')",
    ]
    if state.xmin is not None and state.xmax is not None:
        lines.append(f"ax.set_xlim({state.xmin!r}, {state.xmax!r})")
    if state.ymin is not None and state.ymax is not None:
        lines.append(f"ax.set_ylim({state.ymin!r}, {state.ymax!r})")
    if state.title:
        lines.append(f"ax.set_title({repr(state.title)})")
    lines += [
        "",
        "handles, labels_ = ax.get_legend_handles_labels()",
        "if handles:",
        "    ax.legend()",
        "",
        '_out = Path(__file__).with_suffix(".pdf")',
        'fig.savefig(_out, bbox_inches="tight")',
        'print(f"Saved {_out}")',
    ]

    src = "\n".join(lines) + "\n"
    ast.parse(src)   # validate before returning
    return src


# ---------------------------------------------------------------------------
# Application
# ---------------------------------------------------------------------------

class InteractivePlotApp:

    PANEL_WIDTH = 280

    def __init__(self, root, initial_folder: Optional[str] = None):
        import tkinter as tk
        from tkinter import ttk, filedialog, messagebox
        from matplotlib.backends.backend_tkagg import (
            FigureCanvasTkAgg, NavigationToolbar2Tk,
        )

        self._tk                   = tk
        self._ttk                  = ttk
        self._fd                   = filedialog
        self._mb                   = messagebox
        self._FigureCanvasTkAgg    = FigureCanvasTkAgg
        self._NavigationToolbar2Tk = NavigationToolbar2Tk

        self.root = root
        root.title("ECD Spectrum Viewer")
        root.geometry("1200x700")

        # Global iterators — colors assigned in the order folders are opened
        self._color_it = itertools.cycle(_COLORS)
        self._ls_it    = itertools.cycle(_LINESTYLES)

        # folder_entries: [{"path", "enabled" (BooleanVar), "records", "_frame"}]
        self.folder_entries: list = []
        # style_entries: [{"path": str, "enabled": BooleanVar}]
        self.style_entries: list = []

        self.state  = PlotState()
        self.engine = PlotEngine()
        self._suppress = False

        self._build_layout()

        # Auto-load the default style (same folder as this script)
        self._add_style(STYLE_FILE, warn=False)

        if initial_folder:
            self._add_folder(os.path.abspath(initial_folder))

    # ── Helpers ──────────────────────────────────────────────────────────────

    def _all_records(self) -> list:
        """Flat list of records from all *enabled* folders."""
        out = []
        for fe in self.folder_entries:
            if fe["enabled"].get():
                out.extend(fe["records"])
        return out

    def _scroll(self, direction: int):
        self._sb_canvas.yview_scroll(direction, "units")

    # ── Folder management ─────────────────────────────────────────────────────

    def _add_folder(self, folder: str):
        for fe in self.folder_entries:
            if fe["path"] == folder:
                self._mb.showinfo("Info", f"Folder already open:\n{folder}")
                return

        records = _scan_folder(folder, self._color_it, self._ls_it)
        if not records:
            self._mb.showwarning("Warning", f"No .dat files found in:\n{folder}")
            return

        tk = self._tk
        enabled_var = tk.BooleanVar(value=True)

        fe: dict = {
            "path":    folder,
            "enabled": enabled_var,
            "records": records,
            "_frame":  None,   # dataset sub-frame, set below
        }
        self.folder_entries.append(fe)

        # Add folder checkbox row to Folders panel
        self._append_folder_row(fe)
        # Add dataset rows to Datasets panel
        self._append_dataset_rows(fe)

        # Wire the toggle *after* the frame exists
        enabled_var.trace_add("write", lambda *_: self._on_folder_toggle())

        # Mark all as selected with default labels
        for rec in records:
            self.state.selected.append(rec["path"])
            self.state.labels[rec["path"]] = rec["stem"]

        self._on_change()

    def _on_folder_toggle(self):
        if self._suppress:
            return
        self._update_dataset_visibility()
        self._on_change()

    def _update_dataset_visibility(self):
        for fe in self.folder_entries:
            if fe["_frame"] is not None:
                if fe["enabled"].get():
                    fe["_frame"].pack(fill=self._tk.X, padx=2)
                else:
                    fe["_frame"].pack_forget()

    def _on_open_folder(self):
        folder = self._fd.askdirectory(title="Select folder with .dat files")
        if folder:
            self._add_folder(os.path.abspath(folder))

    # ── Style management ──────────────────────────────────────────────────────

    def _add_style(self, path: str, warn: bool = True):
        if any(se["path"] == path for se in self.style_entries):
            if warn:
                self._mb.showinfo("Info", f"Style already loaded:\n{path}")
            return
        if not Path(path).exists():
            if warn:
                self._mb.showwarning("Warning", f"Style file not found:\n{path}")
            return
        enabled_var = self._tk.BooleanVar(value=True)
        se = {"path": path, "enabled": enabled_var}
        self.style_entries.append(se)
        self._append_style_row(se)
        self._on_change()

    def _append_style_row(self, se: dict):
        row = self._ttk.Frame(self._style_list_frm)
        row.pack(fill=self._tk.X, pady=1)
        self._ttk.Checkbutton(row, variable=se["enabled"],
                               command=self._on_change).pack(side=self._tk.LEFT)
        short = Path(se["path"]).name
        self._ttk.Label(row, text=short, anchor="w").pack(
            side=self._tk.LEFT, fill=self._tk.X, expand=True)

    def _active_styles(self) -> list:
        return [se["path"] for se in self.style_entries if se["enabled"].get()]

    # ── Layout ───────────────────────────────────────────────────────────────

    def _build_layout(self):
        tk  = self._tk
        FigureCanvasTkAgg    = self._FigureCanvasTkAgg
        NavigationToolbar2Tk = self._NavigationToolbar2Tk

        self.root.rowconfigure(0, weight=1)
        self.root.columnconfigure(0, weight=1)

        main = tk.Frame(self.root)
        main.grid(row=0, column=0, sticky="nsew")
        main.rowconfigure(0, weight=1)
        main.columnconfigure(1, weight=1)

        # ── Left: scrollable sidebar ─────────────────────────────────────────
        left_outer = tk.Frame(main, width=self.PANEL_WIDTH, bg="#f0f0f0")
        left_outer.grid(row=0, column=0, sticky="ns", padx=(4, 2), pady=4)
        left_outer.pack_propagate(False)

        scrollbar = tk.Scrollbar(left_outer, orient=tk.VERTICAL)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        self._sb_canvas = tk.Canvas(left_outer, bg="#f0f0f0",
                                    highlightthickness=0,
                                    yscrollcommand=scrollbar.set)
        self._sb_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.config(command=self._sb_canvas.yview)

        self._inner = tk.Frame(self._sb_canvas, bg="#f0f0f0")
        self._sb_canvas.create_window((0, 0), window=self._inner, anchor="nw")
        self._inner.bind(
            "<Configure>",
            lambda e: self._sb_canvas.configure(
                scrollregion=self._sb_canvas.bbox("all")
            ),
        )

        # Scroll bindings: Button-4/5 for Linux/X11, MouseWheel for Win/Mac
        self.root.bind_all("<Button-4>",    lambda e: self._scroll(-1))
        self.root.bind_all("<Button-5>",    lambda e: self._scroll( 1))
        self.root.bind_all("<MouseWheel>",
                           lambda e: self._scroll(-1 * (e.delta // 120)))

        pad = {"padx": 4, "pady": 3}
        self._build_config_panel(pad)
        self._build_folders_panel(pad)
        self._build_datasets_panel(pad)
        self._build_axis_panel(pad)
        self._build_styles_panel(pad)
        self._build_style_export_panel(pad)

        # ── Right: toolbar + canvas ──────────────────────────────────────────
        right = tk.Frame(main)
        right.grid(row=0, column=1, sticky="nsew", padx=(2, 4), pady=4)
        right.rowconfigure(1, weight=1)
        right.columnconfigure(0, weight=1)

        self._fig = plt.figure()
        plt.rcParams["figure.dpi"] = INTERACTIVE_DPI

        self._canvas = FigureCanvasTkAgg(self._fig, master=right)
        toolbar = NavigationToolbar2Tk(self._canvas, right)
        toolbar.update()
        toolbar.grid(row=0, column=0, sticky="ew")
        self._canvas.get_tk_widget().grid(row=1, column=0, sticky="nsew")

    # ── Sidebar panels ────────────────────────────────────────────────────────

    def _build_config_panel(self, pad):
        tk  = self._tk
        ttk = self._ttk

        frm = ttk.LabelFrame(self._inner, text="Config")
        frm.pack(fill=tk.X, **pad)

        row = ttk.Frame(frm)
        row.pack(fill=tk.X, padx=2, pady=4)
        ttk.Button(row, text="Load config",
                   command=self._on_load_config).pack(
                       side=tk.LEFT, expand=True, fill=tk.X, padx=(0, 2))
        ttk.Button(row, text="Export config",
                   command=self._on_export_config).pack(
                       side=tk.LEFT, expand=True, fill=tk.X)

    def _build_folders_panel(self, pad):
        tk  = self._tk
        ttk = self._ttk

        frm = ttk.LabelFrame(self._inner, text="Folders")
        frm.pack(fill=tk.X, **pad)

        ttk.Button(frm, text="Open folder\u2026",
                   command=self._on_open_folder).pack(
                       fill=tk.X, padx=2, pady=(4, 2))

        # Rows added dynamically by _append_folder_row
        self._folder_list_frm = ttk.Frame(frm)
        self._folder_list_frm.pack(fill=tk.X, padx=2, pady=(0, 4))

    def _append_folder_row(self, fe: dict):
        ttk = self._ttk

        row = ttk.Frame(self._folder_list_frm)
        row.pack(fill=self._tk.X, pady=1)
        ttk.Checkbutton(row, variable=fe["enabled"]).pack(side=self._tk.LEFT)
        short = Path(fe["path"]).name or fe["path"]
        ttk.Label(row, text=short, anchor="w").pack(
            side=self._tk.LEFT, fill=self._tk.X, expand=True)

    def _build_datasets_panel(self, pad):
        tk  = self._tk
        ttk = self._ttk

        frm = ttk.LabelFrame(self._inner, text="Datasets")
        frm.pack(fill=tk.X, **pad)
        self._datasets_frm = frm

        # Per-folder sub-frames are packed here by _append_dataset_rows
        self._dataset_body = ttk.Frame(frm)
        self._dataset_body.pack(fill=tk.X, padx=0, pady=(2, 0))

        btn_row = ttk.Frame(frm)
        btn_row.pack(fill=tk.X, padx=2, pady=(4, 2))
        ttk.Button(btn_row, text="All",
                   command=self._select_all).pack(
                       side=tk.LEFT, expand=True, fill=tk.X, padx=(0, 2))
        ttk.Button(btn_row, text="None",
                   command=self._select_none).pack(
                       side=tk.LEFT, expand=True, fill=tk.X)

    def _append_dataset_rows(self, fe: dict):
        """Create a sub-frame for fe's records and pack dataset rows into it."""
        tk  = self._tk
        ttk = self._ttk

        sub = ttk.Frame(self._dataset_body)
        sub.pack(fill=tk.X, padx=2)
        fe["_frame"] = sub

        for rec in fe["records"]:
            row = ttk.Frame(sub)
            row.pack(fill=tk.X, pady=1)

            var = tk.BooleanVar(value=True)
            rec["checked"] = var
            ttk.Checkbutton(row, variable=var,
                            command=self._on_change).pack(side=tk.LEFT)

            lv = tk.StringVar(value=rec["stem"])
            rec["label_var"] = lv
            lv.trace_add("write", lambda *_: self._on_change())
            ttk.Entry(row, textvariable=lv, width=19).pack(
                side=tk.LEFT, fill=tk.X, expand=True, padx=(2, 0))

    def _build_axis_panel(self, pad):
        tk  = self._tk
        ttk = self._ttk

        frm = ttk.LabelFrame(self._inner, text="Axis")
        frm.pack(fill=tk.X, **pad)

        row_xu = ttk.Frame(frm)
        row_xu.pack(fill=tk.X, padx=2, pady=2)
        ttk.Label(row_xu, text="X:").pack(side=tk.LEFT)
        self._x_unit_var = tk.StringVar(value="eV")
        ttk.Radiobutton(row_xu, text="eV",
                        variable=self._x_unit_var, value="eV",
                        command=self._on_change).pack(side=tk.LEFT)
        ttk.Radiobutton(row_xu, text="cm\u207b\u00b9",
                        variable=self._x_unit_var, value="cm-1",
                        command=self._on_change).pack(side=tk.LEFT)

        def _rng_row(label, a_from, a_to):
            row = ttk.Frame(frm)
            row.pack(fill=tk.X, padx=2, pady=1)
            ttk.Label(row, text=label, width=8).pack(side=tk.LEFT)
            vf = tk.StringVar(value="")
            vt = tk.StringVar(value="")
            ttk.Entry(row, textvariable=vf, width=7).pack(side=tk.LEFT, padx=(0, 2))
            ttk.Label(row, text="\u2013").pack(side=tk.LEFT)
            ttk.Entry(row, textvariable=vt, width=7).pack(side=tk.LEFT, padx=(2, 0))
            setattr(self, a_from, vf)
            setattr(self, a_to,   vt)

        _rng_row("X range:", "_xmin_var", "_xmax_var")
        _rng_row("Y range:", "_ymin_var", "_ymax_var")

        row_t = ttk.Frame(frm)
        row_t.pack(fill=tk.X, padx=2, pady=2)
        ttk.Label(row_t, text="Title:", width=8).pack(side=tk.LEFT)
        self._title_var = tk.StringVar(value="")
        ttk.Entry(row_t, textvariable=self._title_var).pack(
            side=tk.LEFT, fill=tk.X, expand=True)

        ttk.Button(frm, text="Apply",
                   command=self._on_change).pack(pady=(2, 4))

    def _build_styles_panel(self, pad):
        tk  = self._tk
        ttk = self._ttk

        frm = ttk.LabelFrame(self._inner, text="Styles")
        frm.pack(fill=tk.X, **pad)

        ttk.Button(frm, text="Browse\u2026",
                   command=self._on_style_browse).pack(
                       fill=tk.X, padx=2, pady=(4, 2))

        self._style_list_frm = ttk.Frame(frm)
        self._style_list_frm.pack(fill=tk.X, padx=2, pady=(0, 4))

    def _build_style_export_panel(self, pad):
        tk  = self._tk
        ttk = self._ttk

        frm = ttk.LabelFrame(self._inner, text="Export")
        frm.pack(fill=tk.X, **pad)

        row_tog = ttk.Frame(frm)
        row_tog.pack(fill=tk.X, padx=2, pady=(4, 2))
        self._style_preview_var = tk.BooleanVar(value=False)
        self._style_export_var  = tk.BooleanVar(value=True)
        ttk.Checkbutton(row_tog, text="Preview style",
                        variable=self._style_preview_var,
                        command=self._on_change).pack(side=tk.LEFT, padx=2)
        ttk.Checkbutton(row_tog, text="Export style",
                        variable=self._style_export_var).pack(side=tk.LEFT, padx=2)

        ttk.Separator(frm, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=4, pady=2)
        ttk.Button(frm, text="Export PDF",
                   command=self._on_export_pdf).pack(fill=tk.X, padx=2, pady=1)
        ttk.Button(frm, text="Snapshot Script",
                   command=self._on_export_script).pack(fill=tk.X, padx=2, pady=1)
        ttk.Separator(frm, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=4, pady=2)
        ttk.Button(frm, text="Exit",
                   command=self._quit).pack(fill=tk.X, padx=2, pady=(1, 6))

    # ── Event handlers ────────────────────────────────────────────────────────

    def _on_change(self, *_):
        if self._suppress:
            return
        self._sync_state()
        self._render()

    def _select_all(self):
        self._suppress = True
        for rec in self._all_records():
            rec["checked"].set(True)
        self._suppress = False
        self._on_change()

    def _select_none(self):
        self._suppress = True
        for rec in self._all_records():
            rec["checked"].set(False)
        self._suppress = False
        self._on_change()

    def _on_style_browse(self):
        path = self._fd.askopenfilename(
            title="Select .mplstyle file",
            filetypes=[("Matplotlib style", "*.mplstyle"), ("All files", "*.*")],
        )
        if path:
            self._add_style(os.path.abspath(path))

    def _on_load_config(self):
        path = self._fd.askopenfilename(
            title="Load config",
            filetypes=[("JSON config", "*.json"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            cfg = _load_config(path)
        except Exception as exc:
            self._mb.showerror("Error", f"Could not read config:\n{exc}")
            return
        self._suppress = True
        self._apply_config(cfg)
        self._suppress = False
        self._update_dataset_visibility()
        self._sync_state()
        self._render()

    def _apply_config(self, cfg: dict):
        """Restore folders and all widget state from a config dict."""
        # Open any folders listed in the config that aren't already loaded
        for fc in cfg.get("folders", []):
            fpath   = fc.get("path", "")
            enabled = fc.get("enabled", True)
            if not fpath:
                continue
            existing = next(
                (fe for fe in self.folder_entries if fe["path"] == fpath), None
            )
            if existing is None:
                if os.path.isdir(fpath):
                    records = _scan_folder(fpath, self._color_it, self._ls_it)
                    if records:
                        tk = self._tk
                        ev = tk.BooleanVar(value=enabled)
                        fe = {"path": fpath, "enabled": ev,
                              "records": records, "_frame": None}
                        self.folder_entries.append(fe)
                        self._append_folder_row(fe)
                        self._append_dataset_rows(fe)
                        # Trace wired after frame exists
                        ev.trace_add("write", lambda *_: self._on_folder_toggle())
                        for rec in records:
                            self.state.selected.append(rec["path"])
                            self.state.labels[rec["path"]] = rec["stem"]
            else:
                existing["enabled"].set(enabled)

        # Axis / title
        self._x_unit_var.set(cfg.get("x_unit", "eV"))
        self._title_var.set(cfg.get("title", ""))

        def _s(v):
            return "" if v is None else str(v)

        self._xmin_var.set(_s(cfg.get("xmin")))
        self._xmax_var.set(_s(cfg.get("xmax")))
        self._ymin_var.set(_s(cfg.get("ymin")))
        self._ymax_var.set(_s(cfg.get("ymax")))

        # Styles
        for sc in cfg.get("styles", []):
            spath    = sc.get("path", "")
            senabled = sc.get("enabled", True)
            if not spath:
                continue
            existing = next((se for se in self.style_entries if se["path"] == spath), None)
            if existing is None and Path(spath).exists():
                self._add_style(spath, warn=False)
                existing = next((se for se in self.style_entries if se["path"] == spath), None)
            if existing is not None:
                existing["enabled"].set(senabled)
        self._style_preview_var.set(cfg.get("style_preview", False))
        self._style_export_var.set(cfg.get("style_export", True))

        # Dataset selection and labels
        selected = set(cfg.get("selected", []))
        labels   = cfg.get("labels", {})
        for rec in self._all_records():
            if selected:
                rec["checked"].set(rec["path"] in selected)
            if rec["path"] in labels:
                rec["label_var"].set(labels[rec["path"]])

    def _on_export_config(self):
        self._sync_state()
        path = self._fd.asksaveasfilename(
            title="Save config",
            defaultextension=".json",
            filetypes=[("JSON config", "*.json"), ("All files", "*.*")],
        )
        if not path:
            return
        _export_config(self.state, path)
        self._mb.showinfo("Export", f"Config saved: {path}")

    def _on_export_pdf(self):
        path = self._fd.asksaveasfilename(
            defaultextension=".pdf",
            filetypes=[("PDF files", "*.pdf")],
            initialfile="spectrum.pdf",
        )
        if not path:
            return
        self._sync_state()
        active = self._active_styles()
        if self.state.style_export and active:
            with plt.style.context(active):
                fig2 = plt.figure()
                self.engine.render(self.state, self._all_records(), fig2,
                                   use_cycler=True)
                fig2.savefig(path, bbox_inches="tight")
        else:
            fig2 = plt.figure()
            self.engine.render(self.state, self._all_records(), fig2)
            fig2.savefig(path, bbox_inches="tight")
        plt.close(fig2)
        print(f"Saved PDF: {path}")
        self._mb.showinfo("Export", f"Saved: {path}")

    def _on_export_script(self):
        self._sync_state()
        path = self._fd.asksaveasfilename(
            defaultextension=".py",
            filetypes=[("Python scripts", "*.py")],
            initialfile="snapshot.py",
        )
        if not path:
            return
        try:
            src = _generate_script(self.state, self._all_records())
        except SyntaxError as exc:
            self._mb.showerror("Error", f"Generated script has syntax error:\n{exc}")
            return
        Path(path).write_text(src)
        print(f"Snapshot script saved: {path}")
        self._mb.showinfo("Export", f"Script saved: {path}")

    # ── State sync ────────────────────────────────────────────────────────────

    def _sync_state(self):
        s = self.state
        s.x_unit = self._x_unit_var.get()
        s.title  = self._title_var.get()

        def _pf(sv):
            try:
                return float(sv.get())
            except ValueError:
                return None

        s.xmin = _pf(self._xmin_var)
        s.xmax = _pf(self._xmax_var)
        s.ymin = _pf(self._ymin_var)
        s.ymax = _pf(self._ymax_var)

        s.styles        = [{"path": se["path"], "enabled": se["enabled"].get()}
                           for se in self.style_entries]
        s.style_preview = self._style_preview_var.get()
        s.style_export  = self._style_export_var.get()

        s.folders = [
            {"path": fe["path"], "enabled": fe["enabled"].get()}
            for fe in self.folder_entries
        ]

        all_recs  = self._all_records()
        s.selected = [r["path"] for r in all_recs if r["checked"].get()]
        s.labels   = {r["path"]: r["label_var"].get() for r in all_recs}

    # ── Rendering ─────────────────────────────────────────────────────────────

    def _render(self):
        active = self._active_styles()
        recs   = self._all_records()
        if self.state.style_preview and active:
            plt.style.use(["default"] + active)
        else:
            plt.style.use("default")
        plt.rcParams["figure.dpi"] = INTERACTIVE_DPI
        self.engine.render(self.state, recs, self._fig)

    def _quit(self):
        plt.close("all")
        self.root.destroy()
        sys.exit(0)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Interactive ECD lineshape viewer",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "folder", nargs="?", default=None,
        help="Initial folder containing .dat lineshape files (optional)",
    )
    args = parser.parse_args()

    import tkinter as tk
    root = tk.Tk()
    initial = os.path.abspath(args.folder) if args.folder else None
    app = InteractivePlotApp(root, initial)
    root.protocol("WM_DELETE_WINDOW", app._quit)
    root.mainloop()


if __name__ == "__main__":
    main()
