#!/usr/bin/env python3
"""
castep_extract.py — Extract key data from CASTEP .castep output files.

Usage:
    castep_extract.py [FILE ...] [DIR ...]     # scan dirs recursively for *.castep
    castep_extract.py -f table  FILE [...]
    castep_extract.py -f csv    FILE [...]
    castep_extract.py -f json   FILE [...]     (default)
    castep_extract.py -o out.json FILE [...]
"""

import re
import sys
import json
import csv as csv_mod
import argparse
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _float(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def _int(s):
    try:
        return int(s)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Extraction
# ---------------------------------------------------------------------------

def extract(path: Path) -> dict:
    text = path.read_text(encoding="utf-8", errors="replace")
    d: dict = {"file": str(path), "name": path.stem}

    # --- Job metadata ---
    m = re.search(r"Job started on host\s+(\S+)", text)
    d["host"] = m.group(1) if m else None

    m = re.search(r"at\s+\w+\s+(\w+\s+\d+\s+\d+:\d+:\d+\s+\d+)", text)
    d["job_date"] = m.group(1).strip() if m else None

    m = re.search(r"CASTEP version\s+(\S+)", text)
    d["castep_version"] = m.group(1) if m else None

    # --- Calculation setup ---
    m = re.search(r"type of calculation\s*:\s*(.+)", text)
    d["calc_type"] = m.group(1).strip() if m else None

    m = re.search(r"using functional\s*:\s*(.+)", text)
    d["functional"] = m.group(1).strip() if m else None

    m = re.search(r"DFT\+D: Semi-empirical dispersion correction\s*:\s*(\w+)", text)
    d["dispersion"] = (m.group(1).strip() == "on") if m else False

    m = re.search(r"SEDC with\s*:\s*(.+)", text)
    d["dispersion_scheme"] = m.group(1).strip() if m else None

    m = re.search(r"plane wave basis set cut-off\s*:\s*([\d.]+)\s*eV", text)
    d["cutoff_eV"] = _float(m.group(1)) if m else None

    m = re.search(r"Calculation parallelised over\s+(\d+)\s+processes", text)
    d["n_processes"] = _int(m.group(1)) if m else None

    # --- System ---
    m = re.search(r"number of +electrons\s*:\s*([\d.]+)", text)
    d["n_electrons"] = _float(m.group(1)) if m else None

    m = re.search(r"number of bands\s*:\s*(\d+)", text)
    d["n_bands"] = _int(m.group(1)) if m else None

    m = re.search(r"Total number of ions in cell\s*=\s*(\d+)", text)
    d["n_ions"] = _int(m.group(1)) if m else None

    m = re.search(r"Total number of species in cell\s*=\s*(\d+)", text)
    d["n_species"] = _int(m.group(1)) if m else None

    # --- Cell parameters (last occurrence = optimized geometry if geomopt) ---
    cell_matches = re.findall(
        r"a\s*=\s*([\d.]+)\s+alpha\s*=\s*([\d.]+)\s*\n"
        r"\s*b\s*=\s*([\d.]+)\s+beta\s*=\s*([\d.]+)\s*\n"
        r"\s*c\s*=\s*([\d.]+)\s+gamma\s*=\s*([\d.]+)",
        text,
    )
    if cell_matches:
        a, alpha, b, beta, c, gamma = cell_matches[-1]
        d["cell_a_A"]      = _float(a)
        d["cell_b_A"]      = _float(b)
        d["cell_c_A"]      = _float(c)
        d["cell_alpha_deg"] = _float(alpha)
        d["cell_beta_deg"]  = _float(beta)
        d["cell_gamma_deg"] = _float(gamma)
    else:
        d.update(dict.fromkeys(
            ["cell_a_A","cell_b_A","cell_c_A","cell_alpha_deg","cell_beta_deg","cell_gamma_deg"]
        ))

    cell_vols = re.findall(r"Current cell volume\s*=\s*([\d.]+)\s+A\*\*3", text)
    d["cell_volume_A3"] = _float(cell_vols[-1]) if cell_vols else None

    # --- Energies (last occurrence = final step) ---
    final_energies = re.findall(r"^Final energy\s*=\s*([-\d.E+]+)\s+eV", text, re.MULTILINE)
    d["final_energy_eV"] = _float(final_energies[-1]) if final_energies else None

    sedc_corrs = re.findall(r"Total Energy Correction\s*:\s*([-\d.E+]+)\s+eV", text)
    d["sedc_correction_eV"] = _float(sedc_corrs[-1]) if sedc_corrs else None

    disp_energies = re.findall(
        r"Dispersion corrected final energy\*\s*=\s*([-\d.E+]+)\s+eV", text
    )
    d["dispersion_corrected_energy_eV"] = _float(disp_energies[-1]) if disp_energies else None

    # --- SCF cycles in the last SCF run ---
    # Delimited by lines of dashes ending "<-- SCF"
    scf_delims = [m.start() for m in re.finditer(r"-{40,} <-- SCF", text)]
    if len(scf_delims) >= 2:
        block = text[scf_delims[-2] : scf_delims[-1]]
        d["final_scf_cycles"] = len(re.findall(r"^\s+\d+\s+-[\d.E+]+", block, re.MULTILINE))
    else:
        d["final_scf_cycles"] = None

    # --- GeomOpt / LBFGS ---
    lbfgs_iters = re.findall(
        r"LBFGS: finished iteration\s+(\d+)\s+with enthalpy=\s*([-\d.E+]+)\s+eV", text
    )
    if lbfgs_iters:
        d["geomopt_n_iterations"] = _int(lbfgs_iters[-1][0])
        d["geomopt_final_enthalpy_eV"] = _float(lbfgs_iters[-1][1])
    else:
        d["geomopt_n_iterations"] = None
        d["geomopt_final_enthalpy_eV"] = None

    # Last LBFGS convergence table — search the block after the last finished-iteration line
    last_iter_pos = None
    for m in re.finditer(r"LBFGS: finished iteration\s+\d+", text):
        last_iter_pos = m.end()

    if last_iter_pos is not None:
        snippet = text[last_iter_pos : last_iter_pos + 800]
        # Each row has a fixed name; |F|max / |dR|max contain literal | so match explicitly
        patterns = [
            ("dE_ion_eV",  r"\|  dE/ion\s+\|\s+([\d.E+-]+)\s+\|\s+([\d.E+-]+)\s+\|[^|]+\|\s+(Yes|No)\s+\| <-- LBFGS"),
            ("Fmax_eVA",   r"\|  \|F\|max\s+\|\s+([\d.E+-]+)\s+\|\s+([\d.E+-]+)\s+\|[^|]+\|\s+(Yes|No)\s+\| <-- LBFGS"),
            ("dRmax_A",    r"\|  \|dR\|max\s+\|\s+([\d.E+-]+)\s+\|\s+([\d.E+-]+)\s+\|[^|]+\|\s+(Yes|No)\s+\| <-- LBFGS"),
        ]
        ok_flags = []
        for col, pat in patterns:
            m = re.search(pat, snippet)
            d[f"geomopt_{col}"]           = _float(m.group(1)) if m else None
            d[f"geomopt_{col}_tol"]       = _float(m.group(2)) if m else None
            d[f"geomopt_{col}_converged"] = (m.group(3) == "Yes") if m else None
            if m:
                ok_flags.append(m.group(3) == "Yes")
        d["geomopt_converged"] = all(ok_flags) if ok_flags else None
    else:
        for k in [
            "geomopt_dE_ion_eV","geomopt_dE_ion_eV_tol","geomopt_dE_ion_eV_converged",
            "geomopt_Fmax_eVA","geomopt_Fmax_eVA_tol","geomopt_Fmax_eVA_converged",
            "geomopt_dRmax_A","geomopt_dRmax_A_tol","geomopt_dRmax_A_converged",
            "geomopt_converged",
        ]:
            d[k] = None

    # --- Timing ---
    m = re.search(r"Total time\s*=\s*([\d.]+)\s*s", text)
    d["total_time_s"] = _float(m.group(1)) if m else None

    m = re.search(r"Calculation time\s*=\s*([\d.]+)\s*s", text)
    d["calc_time_s"] = _float(m.group(1)) if m else None

    m = re.search(r"Peak Memory Use\s*=\s*(\d+)\s*kB", text)
    d["peak_memory_MB"] = _int(m.group(1)) // 1024 if m else None

    m = re.search(r"Overall parallel efficiency rating:\s*\w+\s*\((\d+)%\)", text)
    d["parallel_efficiency_pct"] = _int(m.group(1)) if m else None

    # --- Job status ---
    # "Total time" line is only written on clean termination; absent = killed/crashed
    d["job_completed"] = bool(re.search(r"Total time\s*=\s*[\d.]+\s*s", text))

    # SCF non-convergence warning (CASTEP writes this when max SCF cycles hit)
    scf_warn = re.search(
        r"Warning: max\. number of SCF cycles performed, but system has not converged",
        text, re.IGNORECASE,
    )
    d["scf_converged"] = not bool(scf_warn)

    return d


# ---------------------------------------------------------------------------
# Formatters
# ---------------------------------------------------------------------------

def _yn(val):
    if val is None:
        return "N/A"
    return "Yes" if val else "No"


def format_table(results: list[dict]) -> str:
    lines = []
    sep = "=" * 60

    for d in results:
        lines.append(sep)
        lines.append(f"File : {d['file']}")
        lines.append(f"Name : {d['name']}")
        lines.append(f"Date : {d.get('job_date') or 'N/A'}   Host: {d.get('host') or 'N/A'}")
        lines.append(f"CASTEP version : {d.get('castep_version') or 'N/A'}")
        completed_str = "OK" if d.get("job_completed") else "*** INCOMPLETE / KILLED ***"
        scf_str = "OK" if d.get("scf_converged") else "*** SCF DID NOT CONVERGE ***"
        lines.append(f"Job status     : {completed_str}   SCF: {scf_str}")
        lines.append("")

        lines.append("-- System --")
        lines.append(f"  Calc type    : {d.get('calc_type') or 'N/A'}")
        lines.append(f"  Ions         : {d.get('n_ions')}   Species: {d.get('n_species')}   Electrons: {d.get('n_electrons')}")
        lines.append(f"  Cell (Å/°)   : a={d.get('cell_a_A')}  b={d.get('cell_b_A')}  c={d.get('cell_c_A')}")
        lines.append(f"               : α={d.get('cell_alpha_deg')}  β={d.get('cell_beta_deg')}  γ={d.get('cell_gamma_deg')}")
        lines.append(f"  Volume       : {d.get('cell_volume_A3')} Å³")
        lines.append("")

        lines.append("-- Setup --")
        lines.append(f"  Functional   : {d.get('functional') or 'N/A'}")
        disp_str = f"{d.get('dispersion_scheme') or 'N/A'}" if d.get("dispersion") else "off"
        lines.append(f"  Dispersion   : {disp_str}")
        lines.append(f"  Cutoff       : {d.get('cutoff_eV')} eV")
        lines.append(f"  Processes    : {d.get('n_processes')}")
        lines.append("")

        lines.append("-- Energies --")
        lines.append(f"  Final energy          : {d.get('final_energy_eV')} eV")
        if d.get("dispersion"):
            lines.append(f"  SEDC correction       : {d.get('sedc_correction_eV')} eV")
            lines.append(f"  Dispersion-corrected  : {d.get('dispersion_corrected_energy_eV')} eV")
        lines.append(f"  Final SCF cycles      : {d.get('final_scf_cycles')}")
        lines.append("")

        if d.get("geomopt_n_iterations") is not None:
            lines.append("-- GeomOpt --")
            lines.append(f"  Converged    : {_yn(d.get('geomopt_converged'))}   Iterations: {d.get('geomopt_n_iterations')}")
            lines.append(f"  Enthalpy     : {d.get('geomopt_final_enthalpy_eV')} eV")
            c = d.get
            lines.append(
                f"  dE/ion  : {c('geomopt_dE_ion_eV'):.4E} eV  "
                f"(tol {c('geomopt_dE_ion_eV_tol'):.1E})  {_yn(c('geomopt_dE_ion_eV_converged'))}"
                if c("geomopt_dE_ion_eV") is not None else "  dE/ion  : N/A"
            )
            lines.append(
                f"  |F|max   : {c('geomopt_Fmax_eVA'):.4E} eV/Å  "
                f"(tol {c('geomopt_Fmax_eVA_tol'):.1E})  {_yn(c('geomopt_Fmax_eVA_converged'))}"
                if c("geomopt_Fmax_eVA") is not None else "  |F|max  : N/A"
            )
            lines.append(
                f"  |dR|max  : {c('geomopt_dRmax_A'):.4E} Å    "
                f"(tol {c('geomopt_dRmax_A_tol'):.1E})  {_yn(c('geomopt_dRmax_A_converged'))}"
                if c("geomopt_dRmax_A") is not None else "  |dR|max : N/A"
            )
            lines.append("")

        lines.append("-- Timing --")
        lines.append(f"  Total / Calc : {d.get('total_time_s')} s / {d.get('calc_time_s')} s")
        lines.append(f"  Peak memory  : {d.get('peak_memory_MB')} MB")
        lines.append(f"  Parallel eff : {d.get('parallel_efficiency_pct')} %")
        lines.append("")

    lines.append(sep)
    return "\n".join(lines)


def format_json(results: list[dict]) -> str:
    return json.dumps(results if len(results) > 1 else results[0], indent=2)


def format_csv(results: list[dict]) -> str:
    if not results:
        return ""
    import io
    buf = io.StringIO()
    fields = list(results[0].keys())
    writer = csv_mod.DictWriter(buf, fieldnames=fields, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(results)
    return buf.getvalue()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def collect_files(paths: list[str]) -> list[Path]:
    files = []
    for p in paths:
        p = Path(p)
        if p.is_dir():
            files.extend(sorted(p.rglob("*.castep")))
        elif p.suffix == ".castep":
            files.append(p)
        else:
            print(f"Warning: skipping {p} (not a .castep file or directory)", file=sys.stderr)
    return files


def main():
    parser = argparse.ArgumentParser(
        description="Extract key data from CASTEP .castep output files."
    )
    parser.add_argument(
        "paths", nargs="+", metavar="FILE_OR_DIR",
        help=".castep files or directories (scanned recursively for *.castep)"
    )
    parser.add_argument(
        "-f", "--format", choices=["json", "table", "csv"], default="json",
        help="Output format (default: json)"
    )
    parser.add_argument(
        "-o", "--output", metavar="FILE",
        help="Write output to FILE instead of stdout"
    )
    args = parser.parse_args()

    files = collect_files(args.paths)
    if not files:
        sys.exit("No .castep files found.")

    results = []
    for f in files:
        try:
            results.append(extract(f))
        except Exception as e:
            print(f"Error processing {f}: {e}", file=sys.stderr)

    if args.format == "json":
        output = format_json(results)
    elif args.format == "table":
        output = format_table(results)
    else:
        output = format_csv(results)

    if args.output:
        Path(args.output).write_text(output, encoding="utf-8")
        print(f"Written to {args.output}")
    else:
        print(output)


if __name__ == "__main__":
    main()
