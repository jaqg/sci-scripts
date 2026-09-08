#!/usr/bin/env python3
"""Scan fchk file pairs for orientation-frame mismatches (report mode).

Wraps check_frame_alignment() from gen-d2num-findiff-dipfiles.py, which aborts
on the first mismatch (pipeline guard). This script instead reports the verdict
for every pair and exits 1 if any mismatch was found.

Auto mode (default): scan ROOT dirs recursively for S0-at-S1-geometry files
(*S0atS1*.fchk, *S0FC*.fchk) and pair each with its S1 twin in the same
directory (name substitution: S0atS1 -> S1, S0FC -> S1, trailing _R dropped).

Pair mode: pass explicit fchk paths in pairs (A B [C D ...]).

Usage:
    check-fchk-frames.py                       # auto: ~/R1R2-*
    check-fchk-frames.py -r ~/R1R2-internal    # auto: one root (repeatable)
    check-fchk-frames.py S1.fchk S0atS1.fchk   # explicit pair
"""

import argparse
import glob
import importlib.util
import os
import sys

# Import the guard from the sibling script (dashes in filename -> importlib)
_HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    'gen_d2num', os.path.join(_HERE, 'gen-d2num-findiff-dipfiles.py'))
_g = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_g)
check_frame_alignment = _g.check_frame_alignment


def s1_twin(path):
    """Name of the S1 counterpart of an S0-at-S1-geometry fchk, same dir."""
    stem = os.path.basename(path)
    for pat in ('S0atS1', 'S0FC'):
        if pat in stem:
            return os.path.join(os.path.dirname(path),
                                stem.replace(pat, 'S1').replace('_R', ''))
    return None


def find_pairs(roots):
    """Yield (s1, s0) fchk pairs from recursive scan of roots."""
    for root in roots:
        for s0 in sorted(glob.glob(os.path.join(root, '**', '*S0atS1*.fchk'),
                                   recursive=True)
                         + glob.glob(os.path.join(root, '**', '*S0FC*.fchk'),
                                     recursive=True)):
            s1 = s1_twin(s0)
            if s1 and os.path.exists(s1):
                yield s1, s0
            else:
                print(f'?? no S1 twin found for {s0}')


def check_pair(s1, s0):
    """Return None if aligned, else the guard's message."""
    try:
        check_frame_alignment(s0, s1, labels=(os.path.basename(s0),
                                              os.path.basename(s1)))
        return None
    except SystemExit as e:
        return str(e)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-r', '--root', action='append', default=[],
                    help='Root dir to scan in auto mode (repeatable)')
    ap.add_argument('files', nargs='*',
                    help='Explicit fchk pairs: A B [C D ...] (disables auto mode)')
    args = ap.parse_args()

    if args.files:
        if len(args.files) % 2 != 0:
            sys.exit('Error: explicit mode needs pairs of fchk paths')
        pairs = [(args.files[i], args.files[i + 1])
                 for i in range(0, len(args.files), 2)]
    else:
        roots = args.root or glob.glob(os.path.expanduser('~/R1R2-*'))
        if not roots:
            sys.exit('Error: no roots. Use -r ROOT or pass explicit pairs.')
        pairs = list(find_pairs(roots))

    n_bad = 0
    for s1, s0 in pairs:
        msg = check_pair(s1, s0)
        if msg is None:
            print(f'ALIGNED     {os.path.basename(s0)}  (vs {os.path.basename(s1)})')
        else:
            n_bad += 1
            print(f'MISALIGNED  {os.path.basename(s0)}  (vs {os.path.basename(s1)})')
            for line in msg.splitlines():
                print(f'    {line}')

    print(f'\n{len(pairs)} pair(s) checked, {n_bad} mismatched')
    sys.exit(1 if n_bad else 0)


if __name__ == '__main__':
    main()
