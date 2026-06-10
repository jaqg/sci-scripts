#!/usr/bin/env python3
"""
gamess_extract.py — Extract key data from GAMESS .log output files.

Usage:
    gamess_extract.py file.log                # print summary + write file.xyz
    gamess_extract.py file.log -o out.xyz     # custom XYZ output path
    gamess_extract.py file.log --no-xyz       # skip XYZ file

Extracts:
  - Job metadata: SCF type, run type, functional, basis, charge, multiplicity
  - Geometry: input coordinates (converted Bohr → Å), written as .xyz
  - Energies: final SCF energy, total energy, energy components
  - SCF convergence: iterations, density converged flag
  - Timing: CPU time, wall clock time
  - Job status: completed / failed / killed
"""

import re
import sys
import argparse
from pathlib import Path


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BOHR_TO_ANGSTROM = 0.529177210903

ATOMIC_SYMBOLS: dict[int, str] = {
    1: "H",   2: "He",  3: "Li",  4: "Be",  5: "B",   6: "C",
    7: "N",   8: "O",   9: "F",  10: "Ne", 11: "Na", 12: "Mg",
   13: "Al", 14: "Si", 15: "P",  16: "S",  17: "Cl", 18: "Ar",
   19: "K",  20: "Ca", 21: "Sc", 22: "Ti", 23: "V",  24: "Cr",
   25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu", 30: "Zn",
   31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br", 36: "Kr",
   37: "Rb", 38: "Sr", 39: "Y",  40: "Zr", 41: "Nb", 42: "Mo",
   43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd",
   49: "In", 50: "Sn", 51: "Sb", 52: "Te", 53: "I",  54: "Xe",
   55: "Cs", 56: "Ba", 57: "La", 58: "Ce", 59: "Pr", 60: "Nd",
   61: "Pm", 62: "Sm", 63: "Eu", 64: "Gd", 65: "Tb", 66: "Dy",
   67: "Ho", 68: "Er", 69: "Tm", 70: "Yb", 71: "Lu", 72: "Hf",
   73: "Ta", 74: "W",  75: "Re", 76: "Os", 77: "Ir", 78: "Pt",
   79: "Au", 80: "Hg", 81: "Tl", 82: "Pb", 83: "Bi", 84: "Po",
   85: "At", 86: "Rn", 87: "Fr", 88: "Ra", 89: "Ac", 90: "Th",
   91: "Pa", 92: "U",  93: "Np", 94: "Pu", 95: "Am", 96: "Cm",
   97: "Bk", 98: "Cf", 99: "Es", 100: "Fm", 101: "Md", 102: "No",
  103: "Lr", 104: "Rf", 105: "Db", 106: "Sg", 107: "Bh", 108: "Hs",
  109: "Mt", 110: "Ds", 111: "Rg", 112: "Cn", 113: "Nh", 114: "Fl",
  115: "Mc", 116: "Lv", 117: "Ts", 118: "Og",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _float(s: str) -> float | None:
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def _int(s: str) -> int | None:
    try:
        return int(s)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Extraction
# ---------------------------------------------------------------------------

def extract(path: Path) -> dict:
    """Parse a GAMESS .log file and return a dictionary of extracted data."""
    text = path.read_text(encoding="utf-8", errors="replace")
    d: dict = {"file": str(path), "name": path.stem}

    # -- Run metadata (printed after input parsing) --
    m = re.search(r"SCFTYP=(\S+)\s+RUNTYP=(\S+)\s+EXETYP=(\S+)", text)
    if m:
        d["scftyp"] = m.group(1)
        d["runtyp"] = m.group(2)
        d["exetyp"] = m.group(3)
    else:
        d["scftyp"] = d["runtyp"] = d["exetyp"] = None

    # DFT functional (absent for pure wavefunction methods)
    m = re.search(r"DFTTYP=(\S+)", text)
    d["dfttyp"] = m.group(1) if m and m.group(1) != "NONE" else None

    # Charge and multiplicity
    m = re.search(r"MULT\s*=\s*(\d+)\s+ICHARG\s*=\s*(-?\d+)", text)
    if m:
        d["mult"] = _int(m.group(1))
        d["icharg"] = _int(m.group(2))
    else:
        d["mult"] = d["icharg"] = None

    # Basis set
    m = re.search(r"GBASIS=(\S+)", text)
    d["gbasis"] = m.group(1) if m else None
    m = re.search(r"IGAUSS=\s*(\d+)", text)
    d["igauss"] = _int(m.group(1)) if m else None

    # System info
    m = re.search(r"TOTAL NUMBER OF ATOMS\s*=\s*(\d+)", text)
    d["n_atoms"] = _int(m.group(1)) if m else None

    m = re.search(r"NUMBER OF ELECTRONS\s*=\s*(\d+)", text)
    d["n_electrons"] = _int(m.group(1)) if m else None

    m = re.search(r"CHARGE OF MOLECULE\s*=\s*(-?\d+)", text)
    d["molecule_charge"] = _int(m.group(1)) if m else None

    m = re.search(r"SPIN MULTIPLICITY\s*=\s*(\d+)", text)
    d["spin_mult"] = _int(m.group(1)) if m else None

    m = re.search(r"THE POINT GROUP OF THE MOLECULE IS (\S+)", text)
    d["point_group"] = m.group(1) if m else None

    # -- Geometry: input coordinates (always printed in Bohr) --
    atoms = _parse_coordinates_bohr(text)
    d["atoms"] = atoms  # list of (symbol, atomic_num, x_bohr, y_bohr, z_bohr)

    # -- SCF iterations --
    # RHF/DFT iteration table has header: ITER EX DEM     TOTAL ENERGY ...
    # GVB iteration table has header: ITER EX     TOTAL ENERGY       E CHANGE        SQCDF       ORB. GRAD
    # Count iterations by matching lines starting with a number after the header
    d["scf_iterations"] = _count_scf_iterations(text)
    d["scf_converged"] = bool(re.search(r"DENSITY CONVERGED", text))

    # -- Final SCF energy (printed right after convergence) --
    m = re.search(
        r"FINAL\s+(?:R-B3LYP|R-HF|U-HF|ROHF|GVB|MCSCF|CASSCF|R\-?\w*)\s+ENERGY IS\s+([-\d.E+]+)\s+AFTER\s+(\d+)\s+ITERATIONS",
        text,
    )
    if m:
        d["final_scf_energy"] = _float(m.group(1))
        d["final_scf_iterations"] = _int(m.group(2))
    else:
        d["final_scf_energy"] = None
        d["final_scf_iterations"] = None

    # -- Total energy from energy components block (most reliable final energy) --
    m = re.search(r"TOTAL ENERGY =\s+([-\d.E+]+)", text)
    d["total_energy"] = _float(m.group(1)) if m else None

    # -- Energy components --
    m = re.search(r"ONE ELECTRON ENERGY =\s+([-\d.E+]+)", text)
    d["one_electron_energy"] = _float(m.group(1)) if m else None

    m = re.search(r"TWO ELECTRON ENERGY =\s+([-\d.E+]+)", text)
    d["two_electron_energy"] = _float(m.group(1)) if m else None

    m = re.search(r"NUCLEAR REPULSION ENERGY =\s+([-\d.E+]+)", text)
    d["nuclear_repulsion_energy"] = _float(m.group(1)) if m else None

    # -- Timing: last occurrence of CPU / wall clock --
    cpu_matches = re.findall(
        r"CPU\s+\d+:\s+STEP CPU TIME=\s+([\d.]+)\s+TOTAL CPU TIME=\s+([\d.]+)\s+\(\s*([\d.]+)\s+MIN\)",
        text,
    )
    if cpu_matches:
        d["last_step_cpu_s"] = _float(cpu_matches[-1][0])
        d["total_cpu_s"] = _float(cpu_matches[-1][1])
        d["total_cpu_min"] = _float(cpu_matches[-1][2])
    else:
        d["last_step_cpu_s"] = d["total_cpu_s"] = d["total_cpu_min"] = None

    wall_matches = re.findall(
        r"TOTAL WALL CLOCK TIME=\s+([\d.]+)\s+SECONDS,\s+CPU UTILIZATION IS\s+([\d.]+)%",
        text,
    )
    if wall_matches:
        d["total_wall_s"] = _float(wall_matches[-1][0])
        d["cpu_utilization_pct"] = _float(wall_matches[-1][1])
    else:
        d["total_wall_s"] = d["cpu_utilization_pct"] = None

    # -- Job status --
    d["job_completed"] = bool(re.search(r"ddikick\.x:\s+exited gracefully", text))
    d["job_error"] = bool(re.search(
        r"ddikick\.x:\s+(?:application process \d+ quit unexpectedly|Execution terminated due to error)",
        text,
    ))

    # -- GAMESS version --
    m = re.search(r"GAMESS VERSION =\s*(\S[^*\n]+?)\s*\*?\s*$", text, re.MULTILINE)
    d["gamess_version"] = m.group(1).strip() if m else None

    return d


# ---------------------------------------------------------------------------
# Internal parsers
# ---------------------------------------------------------------------------

def _parse_coordinates_bohr(text: str) -> list[tuple[str, int, float, float, float]]:
    """Parse the COORDINATES (BOHR) table.
    Returns list of (symbol, atomic_num, x_bohr, y_bohr, z_bohr).
    """
    # Find the coordinates block header
    m = re.search(r"ATOM\s+ATOMIC\s+COORDINATES\s+\(BOHR\)", text)
    if not m:
        return []

    # Start searching after the header + one more line (the column labels)
    start = m.end()
    # Skip the label line
    label_end = text.index("\n", start) + 1

    atoms = []
    # Pattern for a coordinate line:  C           6.0    -3.0155434730        0.0004440856       -0.6505608480
    coord_re = re.compile(
        r"^\s*(\S+)\s+(\d+\.?\d*)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"
    )
    for line in text[label_end:].splitlines():
        line = line.strip()
        if not line or line.startswith("---") or line.startswith("INTERNUCLEAR"):
            break
        cm = coord_re.match(line)
        if cm:
            symbol_raw = cm.group(1)
            atomic_num = _int(float(cm.group(2)))  # handle "6.0" → 6
            x = float(cm.group(3))
            y = float(cm.group(4))
            z = float(cm.group(5))
            # Normalize symbol (e.g., "C" stays "C"; handle edge cases)
            if atomic_num and atomic_num in ATOMIC_SYMBOLS:
                symbol = ATOMIC_SYMBOLS[atomic_num]
            else:
                symbol = symbol_raw
            atoms.append((symbol, atomic_num, x, y, z))

    return atoms


def _count_scf_iterations(text: str) -> int | None:
    """Count SCF iterations by finding the iteration table and counting data rows."""
    # Find the iteration table header
    m = re.search(r"ITER\s+EX(?:\s+DEM)?\s+TOTAL ENERGY", text)
    if not m:
        return None

    # Find the end of the iteration block: "DENSITY CONVERGED" or "FINAL ... ENERGY"
    end_pattern = re.compile(r"DENSITY CONVERGED|FINAL\s+\w+(?:\s+\w+)*\s+ENERGY IS")

    # Start after the header line
    start = text.index("\n", m.start()) + 1
    count = 0
    for line in text[start:].splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if end_pattern.search(stripped):
            break
        # Iteration lines start with a number (the iteration count)
        if re.match(r"^\s*\d+\s", line):
            count += 1

    return count if count > 0 else None


# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------

def _yn(val: bool | None) -> str:
    if val is None:
        return "N/A"
    return "Yes" if val else "No"


def _fmt(val, fmt_spec: str = ".6f") -> str:
    if val is None:
        return "N/A"
    return f"{val:{fmt_spec}}"


def format_summary(d: dict) -> str:
    """Format extracted data as a human-readable summary block."""
    lines = []
    sep = "─" * 62

    lines.append(sep)
    lines.append(f"  GAMESS Extract  —  {d['name']}")
    lines.append(f"  File: {d['file']}")
    if d.get("gamess_version"):
        lines.append(f"  Version: {d['gamess_version']}")
    lines.append("")

    # Job status
    if d.get("job_completed"):
        status = "✓ COMPLETED"
    elif d.get("job_error"):
        status = "✗ ERROR / KILLED"
    else:
        status = "? UNKNOWN (check log)"
    scf_status = "✓ converged" if d.get("scf_converged") else "✗ NOT CONVERGED"
    lines.append(f"  Status:  {status}    SCF: {scf_status}")
    lines.append("")

    # Run setup
    lines.append("  ── Setup ──")
    scftyp = d.get("scftyp") or "N/A"
    runtyp = d.get("runtyp") or "N/A"
    dft = d.get("dfttyp")
    method = f"{scftyp}" + (f"/{dft}" if dft else "")
    lines.append(f"  Method:       {method}")
    lines.append(f"  Run type:     {runtyp}")
    gbasis = d.get("gbasis") or "N/A"
    igauss = d.get("igauss")
    basis_str = gbasis
    if igauss is not None and igauss > 0:
        basis_str += f"-{igauss}G"
    lines.append(f"  Basis:        {basis_str}")
    lines.append(f"  Charge: {d.get('icharg')}    Mult: {d.get('mult')}")
    lines.append(f"  Atoms: {d.get('n_atoms')}    Electrons: {d.get('n_electrons')}")
    lines.append(f"  Point group:  {d.get('point_group') or 'N/A'}")
    lines.append("")

    # Energies
    lines.append("  ── Energies (Hartree) ──")
    lines.append(f"  Final SCF energy:        {_fmt(d.get('final_scf_energy'), '.10f')}")
    lines.append(f"  Total energy:            {_fmt(d.get('total_energy'), '.10f')}")
    lines.append(f"    1-electron:            {_fmt(d.get('one_electron_energy'), '.10f')}")
    lines.append(f"    2-electron:            {_fmt(d.get('two_electron_energy'), '.10f')}")
    lines.append(f"    Nuclear repulsion:     {_fmt(d.get('nuclear_repulsion_energy'), '.10f')}")
    scf_it = d.get('final_scf_iterations') or d.get('scf_iterations')
    lines.append(f"  SCF iterations:  {scf_it or 'N/A'}")

    # Check energy consistency: total vs sum of components
    e1 = d.get("one_electron_energy")
    e2 = d.get("two_electron_energy")
    enuc = d.get("nuclear_repulsion_energy")
    etot = d.get("total_energy")
    if all(v is not None for v in (e1, e2, enuc, etot)):
        check = e1 + e2 + enuc
        diff = abs(check - etot)
        if diff > 1e-8:
            lines.append(f"  ⚠ Energy sum check: {e1} + {e2} + {enuc} = {check:.10f}  (diff = {diff:.2e})")

    lines.append("")

    # Geometry summary
    atoms = d.get("atoms", [])
    if atoms:
        lines.append(f"  ── Geometry ({len(atoms)} atoms) ──")
        lines.append(f"  {'Atom':<4s} {'Z':>3s}  {'X (Å)':>12s}  {'Y (Å)':>12s}  {'Z (Å)':>12s}")
        lines.append(f"  {'─'*4:<4s} {'─'*3:>3s}  {'─'*12:>12s}  {'─'*12:>12s}  {'─'*12:>12s}")
        for sym, z, x, y, z_bohr in atoms:
            xa = x * BOHR_TO_ANGSTROM
            ya = y * BOHR_TO_ANGSTROM
            za = z_bohr * BOHR_TO_ANGSTROM
            lines.append(f"  {sym:<4s} {z:>3d}  {xa:12.6f}  {ya:12.6f}  {za:12.6f}")
        lines.append("")
    else:
        lines.append("  ── Geometry: not found ──")
        lines.append("")

    # Timing
    lines.append("  ── Timing ──")
    lines.append(f"  Total CPU time:  {_fmt(d.get('total_cpu_s'), '.1f')} s  ({_fmt(d.get('total_cpu_min'), '.1f')} min)")
    lines.append(f"  Wall clock time: {_fmt(d.get('total_wall_s'), '.1f')} s")
    lines.append(f"  CPU utilization: {_fmt(d.get('cpu_utilization_pct'), '.1f')} %")
    lines.append("")

    lines.append(sep)
    return "\n".join(lines)


def write_xyz(d: dict, outpath: Path):
    """Write geometry as .xyz file (coordinates converted to Ångströms)."""
    atoms = d.get("atoms", [])
    if not atoms:
        print("Warning: no geometry to write XYZ.", file=sys.stderr)
        return

    with open(outpath, "w") as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"Geometry from {d['name']}  (GAMESS, units: Å)\n")
        for sym, _z, x, y, z_bohr in atoms:
            xa = x * BOHR_TO_ANGSTROM
            ya = y * BOHR_TO_ANGSTROM
            za = z_bohr * BOHR_TO_ANGSTROM
            f.write(f"{sym:<3s}  {xa:12.6f}  {ya:12.6f}  {za:12.6f}\n")

    print(f"\nXYZ written to: {outpath}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Extract key data from a GAMESS .log output file.",
        epilog="Prints a summary to stdout and writes an .xyz geometry file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "logfile", metavar="file.log",
        help="GAMESS log file to read",
    )
    parser.add_argument(
        "-o", "--output", metavar="file.xyz", default=None,
        help="Output .xyz file path (default: <file.log>.xyz)",
    )
    parser.add_argument(
        "--no-xyz", action="store_true",
        help="Skip writing the .xyz file",
    )
    args = parser.parse_args()

    logpath = Path(args.logfile)
    if not logpath.exists():
        sys.exit(f"Error: file not found: {args.logfile}")

    data = extract(logpath)

    print(format_summary(data))

    if not args.no_xyz:
        xyzpath = Path(args.output) if args.output else logpath.with_suffix(".xyz")
        write_xyz(data, xyzpath)


if __name__ == "__main__":
    main()
