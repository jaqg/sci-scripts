#!/usr/bin/env python3
"""
Extract the final optimized geometry from a Gaussian log file.

Usage:
    python3 extract_geom.py <file.log>

Outputs the geometry block ready to paste into a new Gaussian .com file,
and also writes <file.log>.xyz
"""

import sys
import re
import argparse

ATOMIC_SYMBOLS = {
    1: 'H',  6: 'C',  7: 'N',  8: 'O',  9: 'F',
   14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 35: 'Br', 53: 'I',
}

def parse_orientation_block(lines, start):
    """Parse one orientation table starting at the header line.
    Returns list of (symbol, x, y, z) or None if malformed."""
    # skip separator, two header lines, separator → data starts at start+5
    i = start + 5
    atoms = []
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('---'):
            break
        parts = line.split()
        if len(parts) < 6:
            break
        atomic_num = int(parts[1])
        x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
        symbol = ATOMIC_SYMBOLS.get(atomic_num, str(atomic_num))
        atoms.append((symbol, x, y, z))
        i += 1
    return atoms if atoms else None


def extract_last_geometry(logfile):
    with open(logfile) as f:
        lines = f.readlines()

    # Collect all orientation blocks (prefer Input orientation when nosymm)
    last_atoms = None
    for i, line in enumerate(lines):
        if re.search(r'(Standard|Input) orientation:', line):
            atoms = parse_orientation_block(lines, i)
            if atoms:
                last_atoms = atoms

    return last_atoms


def main():
    parser = argparse.ArgumentParser(
        description="Extract the final optimized geometry from a Gaussian log file.",
        epilog="Outputs the geometry block ready to paste into a new Gaussian .com file,\n"
               "and also writes <file.log>.xyz",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("logfile", metavar="file.log",
                        help="Gaussian log file to read")
    parser.add_argument("-o", "--output", metavar="file.xyz", default=None,
                        help="output .xyz file (default: <file.log>.xyz)")
    args = parser.parse_args()

    logfile = args.logfile
    atoms = extract_last_geometry(logfile)

    if atoms is None:
        print("ERROR: no orientation block found in", logfile)
        sys.exit(1)

    # Print Gaussian-ready block
    print("# Optimized geometry extracted from:", logfile)
    for sym, x, y, z in atoms:
        print(f"{sym:<3s}  {x:12.6f}  {y:12.6f}  {z:12.6f}")

    # Write .xyz file
    xyzfile = args.output if args.output else logfile + ".xyz"
    with open(xyzfile, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"Optimized geometry from {logfile}\n")
        for sym, x, y, z in atoms:
            f.write(f"{sym:<3s}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")
    print(f"\nAlso written to: {xyzfile}")


if __name__ == '__main__':
    main()
