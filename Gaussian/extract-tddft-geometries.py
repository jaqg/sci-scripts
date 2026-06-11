#!/usr/bin/env python3
"""
Extract optimized excited-state geometries from Gaussian TDDFT log files.

Usage:
    python3 extract-tddft-geometries.py file1.log [file2.log ...]

For each optimized geometry found, writes a .xyz file labeled by state
(e.g., mol_S1.xyz, mol_T3.xyz). Handles:
  - Single .log with one Opt=TD job
  - Multiple .log files (one per root) passed as arguments
  - Multi-step .log files (--Link1-- chaining)
  - Ground-state optimizations (labeled "GS")

State detection: parses "Excited State N: Singlet/Triplet-..." followed by
"This state for optimization" to identify the optimized root.
"""

import sys
import re
import argparse
import os

ATOMIC_SYMBOLS = {
    1: 'H',   6: 'C',   7: 'N',   8: 'O',   9: 'F',
   14: 'Si', 15: 'P',  16: 'S',  17: 'Cl', 35: 'Br', 53: 'I',
}

# ---------------------------------------------------------------------------
# Geometry parsing (reused from extract-optimized-geom-gaussian.py)
# ---------------------------------------------------------------------------

def parse_orientation_block(lines, start):
    """Parse one orientation table starting at the header line.
    Returns list of (symbol, x, y, z) or None if malformed."""
    i = start + 5
    atoms = []
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('---'):
            break
        parts = line.split()
        if len(parts) < 6:
            break
        try:
            atomic_num = int(parts[1])
            x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
        except (ValueError, IndexError):
            break
        symbol = ATOMIC_SYMBOLS.get(atomic_num, str(atomic_num))
        atoms.append((symbol, x, y, z))
        i += 1
    return atoms if atoms else None


# ---------------------------------------------------------------------------
# Log-file splitting (multi-step --Link1-- jobs)
# ---------------------------------------------------------------------------

def split_into_steps(lines):
    """Split a list of lines into sub-lists, one per Gaussian job step,
    delimited by 'Link1:  Proceeding to internal job step'.
    Returns list of (step_lines, step_number) tuples."""
    steps = []
    current_start = 0

    for i, line in enumerate(lines):
        if 'Proceeding to internal job step number' in line:
            if i > current_start:
                steps.append(lines[current_start:i])
            current_start = i

    # Last step
    if current_start < len(lines):
        steps.append(lines[current_start:])

    # If no Link1 boundaries found, whole file is one step
    if not steps:
        steps = [lines]

    return steps


# ---------------------------------------------------------------------------
# TDDFT state detection
# ---------------------------------------------------------------------------

def find_tddft_state(lines):
    """Look for 'Excited State N: <mult>-<sym>' followed within ~5 lines by
    'This state for optimization'. Returns (state_num, multiplicity_str) or
    (None, None) if not found."""
    for i, line in enumerate(lines):
        m = re.search(r'Excited State\s+(\d+):\s+(\S+)', line)
        if not m:
            continue
        state_num = int(m.group(1))
        mult_full = m.group(2)  # e.g. "Singlet-A" or "Triplet-A"

        # Check next ~6 lines for "This state for optimization"
        for j in range(i + 1, min(i + 8, len(lines))):
            if 'This state for optimization' in lines[j]:
                # Extract multiplicity word
                mult_word = mult_full.split('-')[0]  # "Singlet" or "Triplet"
                return state_num, mult_word

    return None, None


# ---------------------------------------------------------------------------
# Extract geometry from one optimization step
# ---------------------------------------------------------------------------

def extract_step_geometry(step_lines):
    """Within a single optimization step, find the final geometry after
    'Optimization completed' or 'Stationary point found'.
    Returns (atoms, state_num, multiplicity) or (None, None, None)."""

    # 1. Check this step actually ran an optimization
    has_opt = any('Optimization completed' in l or
                  'Stationary point found' in l
                  for l in step_lines)

    if not has_opt:
        return None, None, None

    # 2. Detect TDDFT state being optimized (if any)
    state_num, multiplicity = find_tddft_state(step_lines)

    # 3. Find the FINAL geometry: search for orientation blocks,
    #    take the LAST one (after last convergence marker is the true final)
    atoms = None
    for i, line in enumerate(step_lines):
        if re.search(r'(Standard|Input) orientation:', line):
            parsed = parse_orientation_block(step_lines, i)
            if parsed:
                atoms = parsed

    return atoms, state_num, multiplicity


# ---------------------------------------------------------------------------
# XYZ output
# ---------------------------------------------------------------------------

def write_xyz(outpath, atoms, comment, state_label, multiplicity):
    """Write an .xyz file."""
    xyz_comment = f"{comment} | State: {state_label} | Multiplicity: {multiplicity}"
    with open(outpath, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"{xyz_comment}\n")
        for sym, x, y, z in atoms:
            f.write(f"{sym:<3s}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")
    print(f"  -> {outpath}")


MULT_TO_PREFIX = {
    'Singlet': 'S', 'Triplet': 'T', 'Quintet': 'Q',
    'Septet': 'Sep', 'Nonet': 'N',
}


def build_state_label(state_num, multiplicity):
    """Build a short label like 'S1', 'T3', or 'GS'."""
    if state_num is None:
        return 'GS', 'Ground'
    prefix = MULT_TO_PREFIX.get(multiplicity, '?')
    return f"{prefix}{state_num}", multiplicity


# ---------------------------------------------------------------------------
# Disambiguate labels when multiple steps produce same label
# ---------------------------------------------------------------------------

def disambiguate_labels(results):
    """If two results have the same label (e.g. same state optimized twice),
    append _stepN to distinguish them."""
    seen = {}
    for r in results:
        lbl = r['label']
        if lbl in seen:
            # Both get disambiguated
            if not seen[lbl]['disambiguated']:
                seen[lbl]['results'][0]['label'] += '_step1'
                seen[lbl]['disambiguated'] = True
            step_tag = f"_step{r['step_num']}"
            r['label'] = lbl + step_tag
        else:
            seen[lbl] = {'results': [r], 'disambiguated': False}


# ---------------------------------------------------------------------------
# Main processing: one log file -> list of (label, atoms, multiplicity, step)
# ---------------------------------------------------------------------------

def process_log(logfile):
    """Process a single Gaussian .log file.
    Returns list of dicts: {label, atoms, multiplicity, step_num}."""
    with open(logfile) as f:
        lines = f.readlines()

    steps = split_into_steps(lines)
    results = []

    for step_idx, step_lines in enumerate(steps, 1):
        atoms, state_num, multiplicity = extract_step_geometry(step_lines)
        if atoms is None:
            continue

        label_long, _ = build_state_label(state_num, multiplicity)
        results.append({
            'label': label_long,
            'atoms': atoms,
            'multiplicity': multiplicity or 'Unknown',
            'step_num': step_idx,
        })

    # Disambiguate if multiple steps produced same label
    if len(results) > 1:
        disambiguate_labels(results)

    return results


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Extract optimized excited-state geometries from Gaussian "
                    "TDDFT log files.",
        epilog="Examples:\n"
               "  %(prog)s mol_S1.log\n"
               "  %(prog)s mol_S1.log mol_S2.log mol_S3.log\n"
               "  %(prog)s --outdir ./geoms *.log",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("logfiles", nargs='+', metavar="file.log",
                        help="One or more Gaussian .log files")
    parser.add_argument("-o", "--outdir", default=".", metavar="DIR",
                        help="Output directory for .xyz files (default: .)")
    parser.add_argument("--prefix", default=None, metavar="STR",
                        help="Prefix for output filenames (default: derived from log filename)")
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="Suppress printing of geometry blocks to stdout")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    total_geoms = 0
    multi_file = len(args.logfiles) > 1

    for logfile in args.logfiles:
        if not os.path.isfile(logfile):
            print(f"WARNING: file not found, skipping: {logfile}", file=sys.stderr)
            continue

        print(f"Processing: {logfile}")
        results = process_log(logfile)

        if not results:
            print(f"  No optimized geometry found.")
            continue

        log_basename = os.path.splitext(os.path.basename(logfile))[0]
        if args.prefix:
            # When prefix is given for multiple files, include log basename
            # to avoid collisions (e.g. prefix_mol1_S1.xyz)
            basename = args.prefix if not multi_file else f"{args.prefix}_{log_basename}"
        else:
            basename = log_basename

        for r in results:
            outname = f"{basename}_{r['label']}.xyz"
            outpath = os.path.join(args.outdir, outname)
            write_xyz(outpath, r['atoms'], f"TDDFT geometry from {logfile}",
                      r['label'], r['multiplicity'])
            total_geoms += 1

        # Print Gaussian-ready block to stdout (unless --quiet)
        if not args.quiet:
            print("")
            for r in results:
                print(f"# State: {r['label']}  ({r['multiplicity']})  from {logfile}")
                for sym, x, y, z in r['atoms']:
                    print(f" {sym:<3s}  {x:12.6f}  {y:12.6f}  {z:12.6f}")
                print("")

    print(f"Done. {total_geoms} geometry file(s) written to {os.path.abspath(args.outdir)}")


if __name__ == '__main__':
    main()
