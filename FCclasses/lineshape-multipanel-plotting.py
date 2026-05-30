#!/usr/bin/env python3
"""Interactive ECD lineshape viewer — multi-panel figures."""

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
    stem = os.path.splitext(fname)[0]
    if stem.startswith(FILE_PREFIX):
        stem = stem[len(FILE_PREFIX):]
    return stem


def _scan_folder(folder: str, color_it, ls_it) -> list:
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
# State
# ---------------------------------------------------------------------------

@dataclass
class MultiPlotState:
    rows:          int  = 1
    cols:          int  = 2
    sharex:        bool = False
    sharey:        bool = False
    sup_xlabel:    bool = False
    sup_ylabel:    bool = False
    panels:        list = field(default_factory=list)   # list of dict, one per panel
    styles:        list = field(default_factory=list)
    style_preview: bool = False
    style_export:  bool = True
    folders:       list = field(default_factory=list)

    def ensure_panels(self) -> None:
        n = self.rows * self.cols
        while len(self.panels) < n:
            self.panels.append({})
        if len(self.panels) > n:
            self.panels = self.panels[:n]


def _panel_from_dict(d: dict):
    return {
        "x_unit":      d.get("x_unit", "eV"),
        "title":       d.get("title", ""),
        "xmin":        d.get("xmin"),
        "xmax":        d.get("xmax"),
        "ymin":        d.get("ymin"),
        "ymax":        d.get("ymax"),
        "show_xlabel": d.get("show_xlabel", True),
        "show_ylabel": d.get("show_ylabel", True),
        "selected":    list(d.get("selected", [])),
        "labels":      dict(d.get("labels", {})),
    }


# ---------------------------------------------------------------------------
# Engine
# ---------------------------------------------------------------------------

class MultiPanelEngine:
    def render(self, state: MultiPlotState, records: list, fig,
               use_cycler: bool = False) -> None:
        fig.clear()
        state.ensure_panels()

        raw = fig.subplots(state.rows, state.cols,
                           sharex=state.sharex, sharey=state.sharey)
        if state.rows == 1 and state.cols == 1:
            axes_flat = [raw]
        elif state.rows == 1 or state.cols == 1:
            axes_flat = list(raw)
        else:
            axes_flat = [ax for row in raw for ax in row]

        for ax, panel_dict in zip(axes_flat, state.panels):
            panel = _panel_from_dict(panel_dict)

            selected_set = set(panel["selected"])
            for rec in records:
                if rec["path"] not in selected_set:
                    continue
                data  = rec["data"]
                x     = data.energy if panel["x_unit"] == "eV" else data.cm_inv
                label = panel["labels"].get(rec["path"], rec["stem"])
                if use_cycler:
                    ax.plot(x, data.me, marker="None", label=label)
                else:
                    ax.plot(x, data.me, marker="None",
                            label=label,
                            color=rec["color"],
                            linestyle=rec["linestyle"])

            if panel["show_xlabel"] and not state.sup_xlabel:
                ax.set_xlabel(r"Energy (eV)" if panel["x_unit"] == "eV"
                              else r"$\tilde{\nu}\ (\mathrm{cm^{-1}})$")
            if panel["show_ylabel"] and not state.sup_ylabel:
                ax.set_ylabel(r"Lineshape (a.u.)")
            if panel["title"]:
                ax.set_title(panel["title"])

            handles, _ = ax.get_legend_handles_labels()
            if handles:
                ax.legend()

            if panel["xmin"] is not None and panel["xmax"] is not None:
                ax.set_xlim(panel["xmin"], panel["xmax"])
            if panel["ymin"] is not None and panel["ymax"] is not None:
                ax.set_ylim(panel["ymin"], panel["ymax"])

        if state.sup_xlabel:
            first_unit = state.panels[0].get("x_unit", "eV") if state.panels else "eV"
            fig.supxlabel(r"Energy (eV)" if first_unit == "eV"
                          else r"$\tilde{\nu}\ (\mathrm{cm^{-1}})$")
        if state.sup_ylabel:
            fig.supylabel(r"Lineshape (a.u.)")

        fig.tight_layout()
        if hasattr(fig, "canvas"):
            fig.canvas.draw_idle()


# ---------------------------------------------------------------------------
# Config I/O
# ---------------------------------------------------------------------------

_CONFIG_FIELDS = (
    "rows", "cols", "sharex", "sharey", "sup_xlabel", "sup_ylabel",
    "panels", "styles", "style_preview", "style_export", "folders",
)


def _export_config(state: MultiPlotState, path: str) -> None:
    cfg = {f: getattr(state, f) for f in _CONFIG_FIELDS}
    Path(path).write_text(json.dumps(cfg, indent=2))


def _load_config_file(path: str) -> dict:
    return json.loads(Path(path).read_text())


# ---------------------------------------------------------------------------
# Snapshot script generation
# ---------------------------------------------------------------------------

def _generate_script(state: MultiPlotState, records: list) -> str:
    state.ensure_panels()
    rows, cols = state.rows, state.cols
    n = rows * cols

    all_paths: set[str] = set()
    for d in state.panels:
        all_paths.update(d.get("selected", []))
    path_to_rec = {r["path"]: r for r in records if r["path"] in all_paths}

    sx = repr(state.sharex)
    sy = repr(state.sharey)
    if n == 1:
        axes_setup = [
            f"fig, _ax = plt.subplots({rows}, {cols}, sharex={sx}, sharey={sy})",
            "axes_flat = [_ax]",
        ]
    elif rows == 1 or cols == 1:
        axes_setup = [
            f"fig, _axes = plt.subplots({rows}, {cols}, sharex={sx}, sharey={sy})",
            "axes_flat = list(_axes)",
        ]
    else:
        axes_setup = [
            f"fig, _axes = plt.subplots({rows}, {cols}, sharex={sx}, sharey={sy})",
            "axes_flat = [ax for row in _axes for ax in row]",
        ]

    lines = [
        "#!/usr/bin/env python3",
        '"""Snapshot — auto-generated by lineshape-multipanel-plotting.py"""',
        "import sys, argparse",
        "from pathlib import Path",
        "",
        "_ap = argparse.ArgumentParser()",
        "_ap.add_argument('-sd', '--scripts-dir', default='.',",
        "    metavar='DIR', help='directory containing plot_class_functions.py')",
        "_args = _ap.parse_args()",
        "",
        "import matplotlib",
        'matplotlib.use("pgf")',
        "import matplotlib.pyplot as plt",
        "",
        f"_FALLBACK = Path({repr(_SCRIPT_DIR)})",
        "_HERE = Path(__file__).parent",
        "_SCRIPTS_DIR = Path(_args.scripts_dir).resolve()",
        "for _d in [_HERE, _SCRIPTS_DIR, _FALLBACK]:",
        "    if (_d / 'plot_class_functions.py').exists():",
        "        sys.path.insert(0, str(_d))",
        "        break",
        "else:",
        "    print('Error: plot_class_functions.py not found.', file=sys.stderr)",
        "    sys.exit(1)",
        "from plot_class_functions import SpectrumData",
        "",
        "_active_styles = " + repr([s["path"] for s in state.styles if s.get("enabled", True)]),
        "if _active_styles:",
        "    plt.style.use(_active_styles)",
        "",
    ] + axes_setup + [""]

    if path_to_rec:
        lines += ["_ds = {}"]
        for path in path_to_rec:
            lines.append(f"_ds[{repr(path)}] = SpectrumData({repr(path)})")
        lines.append("")

    for i, panel_dict in enumerate(state.panels):
        panel   = _panel_from_dict(panel_dict)
        x_attr  = "energy" if panel["x_unit"] == "eV" else "cm_inv"
        x_label = r"r'Energy (eV)'" if panel["x_unit"] == "eV" \
                  else r"r'$\tilde{\nu}\ (\mathrm{cm^{-1}})$'"
        lines.append(f"ax = axes_flat[{i}]")
        for path in panel["selected"]:
            if path in path_to_rec:
                label = panel["labels"].get(path, path_to_rec[path]["stem"])
                lines.append(
                    f"ax.plot(_ds[{repr(path)}].{x_attr}, "
                    f"_ds[{repr(path)}].me, marker='None', label={repr(label)})"
                )
        if panel["show_xlabel"] and not state.sup_xlabel:
            lines.append(f"ax.set_xlabel({x_label})")
        if panel["show_ylabel"] and not state.sup_ylabel:
            lines.append("ax.set_ylabel(r'Lineshape (a.u.)')")
        if panel["title"]:
            lines.append(f"ax.set_title({repr(panel['title'])})")
        if panel["xmin"] is not None and panel["xmax"] is not None:
            lines.append(f"ax.set_xlim({panel['xmin']!r}, {panel['xmax']!r})")
        if panel["ymin"] is not None and panel["ymax"] is not None:
            lines.append(f"ax.set_ylim({panel['ymin']!r}, {panel['ymax']!r})")
        lines += [
            "_h, _ = ax.get_legend_handles_labels()",
            "if _h: ax.legend()",
            "",
        ]

    if state.sup_xlabel:
        first_unit = state.panels[0].get("x_unit", "eV") if state.panels else "eV"
        sup_x_text = r"r'Energy (eV)'" if first_unit == "eV" \
                     else r"r'$\tilde{\nu}\ (\mathrm{cm^{-1}})$'"
        lines.append(f"fig.supxlabel({sup_x_text})")
    if state.sup_ylabel:
        lines.append("fig.supylabel(r'Lineshape (a.u.)')")

    lines += [
        "fig.tight_layout()",
        '_out = Path(__file__).with_suffix(".pdf")',
        'fig.savefig(_out, bbox_inches="tight")',
        'print(f"Saved {_out}")',
    ]

    src = "\n".join(lines) + "\n"
    ast.parse(src)
    return src


# ---------------------------------------------------------------------------
# Application
# ---------------------------------------------------------------------------

class InteractiveMultiPanelApp:

    PANEL_WIDTH = 290

    def __init__(self, root, initial_folder: Optional[str] = None):
        import tkinter as tk
        from tkinter import ttk, filedialog, messagebox
        from matplotlib.backends.backend_tkagg import (
            FigureCanvasTkAgg, NavigationToolbar2Tk,
        )

        self._tk  = tk
        self._ttk = ttk
        self._fd  = filedialog
        self._mb  = messagebox
        self._FigureCanvasTkAgg    = FigureCanvasTkAgg
        self._NavigationToolbar2Tk = NavigationToolbar2Tk

        self.root = root
        root.title("ECD Multi-Panel Viewer")
        root.geometry("1300x750")

        self._color_it = itertools.cycle(_COLORS)
        self._ls_it    = itertools.cycle(_LINESTYLES)

        self.folder_entries: list = []
        self.style_entries:  list = []

        self.state          = MultiPlotState()
        self.state.ensure_panels()
        self.engine         = MultiPanelEngine()
        self._suppress      = False
        self._active_panel  = 0

        self._build_layout()
        self._add_style(STYLE_FILE, warn=False)

        if initial_folder:
            self._add_folder(os.path.abspath(initial_folder))

    # ── Helpers ──────────────────────────────────────────────────────────────

    def _all_records(self) -> list:
        out = []
        for fe in self.folder_entries:
            if fe["enabled"].get():
                out.extend(fe["records"])
        return out

    def _scroll(self, direction: int):
        self._sb_canvas.yview_scroll(direction, "units")

    @staticmethod
    def _pf(sv) -> Optional[float]:
        try:
            return float(sv.get())
        except ValueError:
            return None

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
        fe: dict = {"path": folder, "enabled": enabled_var,
                    "records": records, "_frame": None}
        self.folder_entries.append(fe)
        self._append_folder_row(fe)
        self._append_dataset_rows(fe)
        enabled_var.trace_add("write", lambda *_: self._on_folder_toggle())

        # Select new records in active panel by default
        for rec in records:
            rec["checked"].set(True)

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
        self._ttk.Label(row, text=Path(se["path"]).name, anchor="w").pack(
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

        self.root.bind_all("<Button-4>",    lambda e: self._scroll(-1))
        self.root.bind_all("<Button-5>",    lambda e: self._scroll( 1))
        self.root.bind_all("<MouseWheel>",
                           lambda e: self._scroll(-1 * (e.delta // 120)))

        pad = {"padx": 4, "pady": 3}
        self._build_config_panel(pad)
        self._build_folders_panel(pad)
        self._build_grid_panel(pad)
        self._build_panel_selector(pad)
        self._build_datasets_panel(pad)
        self._build_axis_panel(pad)
        self._build_styles_panel(pad)
        self._build_export_panel(pad)

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
        ttk.Button(frm, text="Open folder…",
                   command=self._on_open_folder).pack(
                       fill=tk.X, padx=2, pady=(4, 2))
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

    def _build_grid_panel(self, pad):
        tk  = self._tk
        ttk = self._ttk
        frm = ttk.LabelFrame(self._inner, text="Grid")
        frm.pack(fill=tk.X, **pad)
        row = ttk.Frame(frm)
        row.pack(fill=tk.X, padx=2, pady=4)
        ttk.Label(row, text="Rows:").pack(side=tk.LEFT)
        self._rows_var = tk.IntVar(value=self.state.rows)
        ttk.Spinbox(row, from_=1, to=6, width=4,
                    textvariable=self._rows_var).pack(side=tk.LEFT, padx=(2, 8))
        ttk.Label(row, text="Cols:").pack(side=tk.LEFT)
        self._cols_var = tk.IntVar(value=self.state.cols)
        ttk.Spinbox(row, from_=1, to=6, width=4,
                    textvariable=self._cols_var).pack(side=tk.LEFT, padx=(2, 0))
        row2 = ttk.Frame(frm)
        row2.pack(fill=tk.X, padx=2, pady=(0, 2))
        self._sharex_var = tk.BooleanVar(value=self.state.sharex)
        self._sharey_var = tk.BooleanVar(value=self.state.sharey)
        ttk.Checkbutton(row2, text="Share X",
                        variable=self._sharex_var).pack(side=tk.LEFT, padx=(0, 8))
        ttk.Checkbutton(row2, text="Share Y",
                        variable=self._sharey_var).pack(side=tk.LEFT)

        row3 = ttk.Frame(frm)
        row3.pack(fill=tk.X, padx=2, pady=(0, 2))
        self._sup_xlabel_var = tk.BooleanVar(value=self.state.sup_xlabel)
        self._sup_ylabel_var = tk.BooleanVar(value=self.state.sup_ylabel)
        ttk.Checkbutton(row3, text="Sup X label",
                        variable=self._sup_xlabel_var,
                        command=self._on_sup_label_change).pack(side=tk.LEFT, padx=(0, 4))
        ttk.Checkbutton(row3, text="Sup Y label",
                        variable=self._sup_ylabel_var,
                        command=self._on_sup_label_change).pack(side=tk.LEFT)

        ttk.Button(frm, text="Apply grid",
                   command=self._on_apply_grid).pack(fill=tk.X, padx=2, pady=(0, 4))

    def _build_panel_selector(self, pad):
        tk  = self._tk
        ttk = self._ttk
        self._panel_selector_frm = ttk.LabelFrame(self._inner, text="Active panel")
        self._panel_selector_frm.pack(fill=tk.X, **pad)
        self._panel_var         = tk.IntVar(value=0)
        self._panel_radio_frame = ttk.Frame(self._panel_selector_frm)
        self._panel_radio_frame.pack(fill=tk.X, padx=2, pady=(2, 4))
        self._rebuild_panel_radios()

    def _rebuild_panel_radios(self):
        ttk = self._ttk
        for w in self._panel_radio_frame.winfo_children():
            w.destroy()
        n = self.state.rows * self.state.cols
        for i in range(n):
            ttk.Radiobutton(
                self._panel_radio_frame,
                text=f"P{i + 1}",
                variable=self._panel_var,
                value=i,
                command=lambda idx=i: self._on_panel_select(idx),
            ).pack(side=self._tk.LEFT, padx=2)

    def _build_datasets_panel(self, pad):
        tk  = self._tk
        ttk = self._ttk
        frm = ttk.LabelFrame(self._inner, text="Datasets (active panel)")
        frm.pack(fill=tk.X, **pad)
        self._datasets_frm = frm

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
        tk  = self._tk
        ttk = self._ttk
        sub = ttk.Frame(self._dataset_body)
        sub.pack(fill=tk.X, padx=2)
        fe["_frame"] = sub

        for rec in fe["records"]:
            row = ttk.Frame(sub)
            row.pack(fill=tk.X, pady=1)

            var = tk.BooleanVar(value=False)
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
        frm = ttk.LabelFrame(self._inner, text="Axis (active panel)")
        frm.pack(fill=tk.X, **pad)

        row_lbl = ttk.Frame(frm)
        row_lbl.pack(fill=tk.X, padx=2, pady=(4, 2))
        self._show_xlabel_var = tk.BooleanVar(value=True)
        self._show_ylabel_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(row_lbl, text="X label",
                        variable=self._show_xlabel_var,
                        command=self._on_change).pack(side=tk.LEFT, padx=(0, 8))
        ttk.Checkbutton(row_lbl, text="Y label",
                        variable=self._show_ylabel_var,
                        command=self._on_change).pack(side=tk.LEFT)

        row_xu = ttk.Frame(frm)
        row_xu.pack(fill=tk.X, padx=2, pady=2)
        ttk.Label(row_xu, text="X unit:").pack(side=tk.LEFT)
        self._x_unit_var = tk.StringVar(value="eV")
        ttk.Radiobutton(row_xu, text="eV",
                        variable=self._x_unit_var, value="eV",
                        command=self._on_change).pack(side=tk.LEFT)
        ttk.Radiobutton(row_xu, text="cm⁻¹",
                        variable=self._x_unit_var, value="cm-1",
                        command=self._on_change).pack(side=tk.LEFT)

        def _rng_row(label, a_from, a_to):
            row = ttk.Frame(frm)
            row.pack(fill=tk.X, padx=2, pady=1)
            ttk.Label(row, text=label, width=8).pack(side=tk.LEFT)
            vf = tk.StringVar(value="")
            vt = tk.StringVar(value="")
            ttk.Entry(row, textvariable=vf, width=7).pack(side=tk.LEFT, padx=(0, 2))
            ttk.Label(row, text="–").pack(side=tk.LEFT)
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

        ttk.Button(frm, text="Apply", command=self._on_change).pack(pady=(2, 4))

    def _build_styles_panel(self, pad):
        tk  = self._tk
        ttk = self._ttk
        frm = ttk.LabelFrame(self._inner, text="Styles")
        frm.pack(fill=tk.X, **pad)
        ttk.Button(frm, text="Browse…",
                   command=self._on_style_browse).pack(
                       fill=tk.X, padx=2, pady=(4, 2))
        self._style_list_frm = ttk.Frame(frm)
        self._style_list_frm.pack(fill=tk.X, padx=2, pady=(0, 4))

    def _build_export_panel(self, pad):
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

    # ── Per-panel state ───────────────────────────────────────────────────────

    def _save_panel_to_state(self, idx: int):
        self.state.ensure_panels()
        if idx < 0 or idx >= len(self.state.panels):
            return
        all_recs = self._all_records()
        panel = {
            "x_unit":      self._x_unit_var.get(),
            "title":       self._title_var.get(),
            "xmin":        self._pf(self._xmin_var),
            "xmax":        self._pf(self._xmax_var),
            "ymin":        self._pf(self._ymin_var),
            "ymax":        self._pf(self._ymax_var),
            "show_xlabel": self._show_xlabel_var.get(),
            "show_ylabel": self._show_ylabel_var.get(),
            "selected":    [r["path"] for r in all_recs if r["checked"].get()],
            "labels":      {r["path"]: r["label_var"].get() for r in all_recs},
        }
        self.state.panels[idx] = panel

        if self.state.sharex:
            for i, p in enumerate(self.state.panels):
                if i != idx:
                    p["xmin"] = panel["xmin"]
                    p["xmax"] = panel["xmax"]
        if self.state.sharey:
            for i, p in enumerate(self.state.panels):
                if i != idx:
                    p["ymin"] = panel["ymin"]
                    p["ymax"] = panel["ymax"]

    def _load_panel_from_state(self, idx: int):
        self.state.ensure_panels()
        if idx < 0 or idx >= len(self.state.panels):
            return
        d        = self.state.panels[idx]
        selected = set(d.get("selected", []))
        labels   = d.get("labels", {})
        self._suppress = True
        for rec in self._all_records():
            rec["checked"].set(rec["path"] in selected)
            rec["label_var"].set(labels.get(rec["path"], rec["stem"]))
        self._x_unit_var.set(d.get("x_unit", "eV"))
        self._title_var.set(d.get("title", ""))
        self._show_xlabel_var.set(d.get("show_xlabel", True))
        self._show_ylabel_var.set(d.get("show_ylabel", True))
        def _s(v): return "" if v is None else str(v)
        self._xmin_var.set(_s(d.get("xmin")))
        self._xmax_var.set(_s(d.get("xmax")))
        self._ymin_var.set(_s(d.get("ymin")))
        self._ymax_var.set(_s(d.get("ymax")))
        self._suppress = False

    # ── Event handlers ────────────────────────────────────────────────────────

    def _on_change(self, *_):
        if self._suppress:
            return
        self._sync_state()
        self._render()

    def _on_panel_select(self, idx: int):
        if self._suppress:
            return
        self._save_panel_to_state(self._active_panel)
        self._active_panel = idx
        self._load_panel_from_state(idx)
        self._on_change()

    def _on_sup_label_change(self):
        self.state.sup_xlabel = self._sup_xlabel_var.get()
        self.state.sup_ylabel = self._sup_ylabel_var.get()
        self._on_change()

    def _on_apply_grid(self):
        self._save_panel_to_state(self._active_panel)
        try:
            rows = max(1, int(self._rows_var.get()))
            cols = max(1, int(self._cols_var.get()))
        except (ValueError, self._tk.TclError):
            return
        self.state.rows      = rows
        self.state.cols      = cols
        self.state.sharex    = self._sharex_var.get()
        self.state.sharey    = self._sharey_var.get()
        self.state.sup_xlabel = self._sup_xlabel_var.get()
        self.state.sup_ylabel = self._sup_ylabel_var.get()
        self.state.ensure_panels()
        n = rows * cols
        self._active_panel = min(self._active_panel, n - 1)
        self._panel_var.set(self._active_panel)
        self._rebuild_panel_radios()
        self._load_panel_from_state(self._active_panel)
        self._on_change()

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
            cfg = _load_config_file(path)
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
                        ev.trace_add("write", lambda *_: self._on_folder_toggle())
            else:
                existing["enabled"].set(enabled)

        rows = cfg.get("rows", 1)
        cols = cfg.get("cols", 2)
        self._rows_var.set(rows)
        self._cols_var.set(cols)
        self.state.rows   = rows
        self.state.cols   = cols
        self.state.sharex     = cfg.get("sharex", False)
        self.state.sharey     = cfg.get("sharey", False)
        self.state.sup_xlabel = cfg.get("sup_xlabel", False)
        self.state.sup_ylabel = cfg.get("sup_ylabel", False)
        self._sharex_var.set(self.state.sharex)
        self._sharey_var.set(self.state.sharey)
        self._sup_xlabel_var.set(self.state.sup_xlabel)
        self._sup_ylabel_var.set(self.state.sup_ylabel)
        self.state.panels = cfg.get("panels", [])
        self.state.ensure_panels()
        self._active_panel = min(self._active_panel, rows * cols - 1)
        self._panel_var.set(self._active_panel)
        self._rebuild_panel_radios()

        for sc in cfg.get("styles", []):
            spath    = sc.get("path", "")
            senabled = sc.get("enabled", True)
            if not spath:
                continue
            existing = next(
                (se for se in self.style_entries if se["path"] == spath), None
            )
            if existing is None and Path(spath).exists():
                self._add_style(spath, warn=False)
                existing = next(
                    (se for se in self.style_entries if se["path"] == spath), None
                )
            if existing is not None:
                existing["enabled"].set(senabled)
        self._style_preview_var.set(cfg.get("style_preview", False))
        self._style_export_var.set(cfg.get("style_export", True))

        self._load_panel_from_state(self._active_panel)

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
            initialfile="spectrum-multipanel.pdf",
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
            initialfile="snapshot_multipanel.py",
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
        self._save_panel_to_state(self._active_panel)
        s.styles        = [{"path": se["path"], "enabled": se["enabled"].get()}
                           for se in self.style_entries]
        s.style_preview = self._style_preview_var.get()
        s.style_export  = self._style_export_var.get()
        s.folders       = [{"path": fe["path"], "enabled": fe["enabled"].get()}
                           for fe in self.folder_entries]

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
        description="Interactive ECD multi-panel lineshape viewer",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "folder", nargs="?", default=None,
        help="Initial folder containing .dat lineshape files (optional)",
    )
    parser.add_argument(
        "-c", "--config", default=None, metavar="FILE",
        help="JSON config file to load on startup",
    )
    args = parser.parse_args()

    import tkinter as tk
    root = tk.Tk()
    initial = os.path.abspath(args.folder) if args.folder else None
    app = InteractiveMultiPanelApp(root, initial)

    if args.config:
        try:
            cfg = _load_config_file(os.path.abspath(args.config))
        except Exception as exc:
            print(f"Error loading config {args.config!r}: {exc}", file=sys.stderr)
            sys.exit(1)
        app._suppress = True
        app._apply_config(cfg)
        app._suppress = False
        app._update_dataset_visibility()
        app._sync_state()
        app._render()

    root.protocol("WM_DELETE_WINDOW", app._quit)
    root.mainloop()


if __name__ == "__main__":
    main()
