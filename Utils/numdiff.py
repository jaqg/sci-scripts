#!/usr/bin/env python3
"""
diff two files, suppressing lines that differ only by numerical noise.

Usage: numdiff.py [-r RTOL] [-a ATOL] file1 file2
  -r  relative tolerance (default: 1e-3)
  -a  absolute tolerance (default: 0, not used unless set)
"""
import sys
import re
import argparse


_NUM = re.compile(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][+-]?\d+)?')


def extract_floats(line):
    return [float(x.replace('d', 'e').replace('D', 'E')) for x in _NUM.findall(line)]


def only_numerical_difference(l1, l2, rtol, atol):
    f1, f2 = extract_floats(l1), extract_floats(l2)
    if not f1 and not f2:
        return False
    if len(f1) != len(f2):
        return False
    for a, b in zip(f1, f2):
        if abs(a - b) > atol + rtol * max(abs(a), abs(b)):
            return False
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('file1')
    parser.add_argument('file2')
    parser.add_argument('-r', '--rtol', type=float, default=1e-3)
    parser.add_argument('-a', '--atol', type=float, default=1e-12)
    args = parser.parse_args()

    lines1 = open(args.file1).readlines()
    lines2 = open(args.file2).readlines()

    found = False
    for i, (l1, l2) in enumerate(zip(lines1, lines2), 1):
        if l1 == l2:
            continue
        if only_numerical_difference(l1, l2, args.rtol, args.atol):
            continue
        print(f"{i}c{i}")
        print(f"< {l1.rstrip()}")
        print(f"---")
        print(f"> {l2.rstrip()}")
        found = True

    for i, l in enumerate(lines1[len(lines2):], len(lines2) + 1):
        print(f"{i}d")
        print(f"< {l.rstrip()}")
        found = True

    for i, l in enumerate(lines2[len(lines1):], len(lines1) + 1):
        print(f"{i}a")
        print(f"> {l.rstrip()}")
        found = True

    sys.exit(1 if found else 0)


if __name__ == '__main__':
    main()
