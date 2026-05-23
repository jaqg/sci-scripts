#!/usr/bin/env python3
"""Remove large CASTEP binary output files (.castep_bin, .cst_esp) from a tree.

Safe after calculations are complete: energies live in .castep, geometries in
.xsd/.geom. These binaries are only needed to restart interrupted runs.

Usage:
    python castep-clean-largefiles.py [ROOT] [--skip-dir DIR] [--delete]

    ROOT        directory to scan (default: current directory)
    --skip-dir  directory name to exclude entirely (default: bin)
    --delete    actually delete; without this flag runs in dry-run mode
"""

import argparse
import os
from pathlib import Path

BINARY_EXTENSIONS = {'.castep_bin', '.cst_esp'}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('root', nargs='?', default='.',
                        help='Root directory to scan (default: .)')
    parser.add_argument('--skip-dir', default='bin',
                        help='Directory name to skip recursively (default: bin)')
    parser.add_argument('--delete', action='store_true',
                        help='Delete files; without this flag runs dry-run only')
    args = parser.parse_args()

    root = Path(args.root).resolve()
    total_size = 0
    found = []

    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d != args.skip_dir]
        for fname in filenames:
            p = Path(dirpath) / fname
            if p.suffix in BINARY_EXTENSIONS:
                size = p.stat().st_size
                found.append((p, size))
                total_size += size

    if not found:
        print('No CASTEP binary files found.')
        return

    mode = 'DELETE' if args.delete else 'DRY-RUN'
    print(f'[{mode}]  root: {root}  skip-dir: {args.skip_dir}')
    print(f'{"File":<80} {"Size":>10}')
    print('-' * 92)
    for p, size in found:
        rel = p.relative_to(root)
        print(f'{str(rel):<80} {size / 1024**3:>9.2f} GB')
        if args.delete:
            p.unlink()

    action = 'Deleted' if args.delete else 'Would delete'
    print('-' * 92)
    print(f'{action} {len(found)} files, {total_size / 1024**3:.2f} GB total')
    if not args.delete:
        print('Re-run with --delete to remove.')


if __name__ == '__main__':
    main()
