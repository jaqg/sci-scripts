#!/usr/bin/env python3
"""
Sync a directory between local machine and a remote SSH host via rsync.

Usage:
    folder-sync.py -H HOST -ld LOCAL_DIR -hd HOST_DIR (--up | --down) [--dry-run]

Required:
    -H  / --host        SSH host (e.g. cluster, or user@hostname)
    -ld / --local-dir   Local directory
    -hd / --host-dir    Remote directory
    --up / --upload     Sync local → remote
    --down / --download Sync remote → local

Options:
    -u  / --user        Remote username (default: user already in --host, or SSH config)
    -e  / --exclude     Exclude pattern (can be used multiple times)
                        e.g. -e '*.txt' -e 'dir/'
    -of / --only-files  Only sync files matching pattern (can be used multiple times)
                        e.g. -of '*.txt' -of 'data/*.csv'
                        (equivalent to excluding everything else)
    --dry-run           Show what would be transferred without doing it
    --delete            Delete files on the destination that are absent on the source
    --skip-empty-dirs   Skip empty directories (passes --prune-empty-dirs to rsync)
"""

import argparse
import os
import subprocess
import sys


def sync(src, dst, *, dry_run=False, delete=False, exclude=None, only_files=None, skip_empty_dirs=False):
    cmd = ["rsync", "-avz", "--progress"]
    if dry_run:
        cmd.append("--dry-run")
    if delete:
        cmd.append("--delete")
    if skip_empty_dirs:
        cmd.append("--prune-empty-dirs")
    # First-match-wins: specific excludes beat the broad include/exclude below
    for pattern in (exclude or []):
        cmd += ["--exclude", pattern]
    if only_files:
        cmd += ["--include", "*/"]  # keep recursing into directories
        for pattern in only_files:
            cmd += ["--include", pattern]
        cmd += ["--exclude", "*"]   # drop everything not included
    # Trailing slash on src means "contents of dir", not "dir itself"
    cmd += [src.rstrip("/") + "/", dst.rstrip("/") + "/"]
    print("  $", " ".join(cmd))
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit(result.returncode)


def main():
    parser = argparse.ArgumentParser(
        description="Sync a directory between local machine and a remote SSH host.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("-H", "--host", required=True, metavar="HOST",
                        help="SSH host (e.g. cluster, or user@hostname)")
    parser.add_argument("-u", "--user", default=None, metavar="USER",
                        help="Remote username (overrides any user in --host)")
    parser.add_argument("-ld", "--local-dir", required=True, metavar="LOCAL_DIR",
                        help="Local directory")
    parser.add_argument("-hd", "--host-dir", required=True, metavar="HOST_DIR",
                        help="Remote directory")

    direction = parser.add_mutually_exclusive_group(required=True)
    direction.add_argument("--up", "--upload", dest="upload", action="store_true",
                           help="Upload: local → remote")
    direction.add_argument("--down", "--download", dest="upload", action="store_false",
                           help="Download: remote → local")

    parser.add_argument("-e", "--exclude", action="append", default=[], metavar="PATTERN",
                        help="Exclude pattern, e.g. '*.txt' or 'dir/' (repeatable)")
    parser.add_argument("-of", "--only-files", action="append", nargs="+", default=[], metavar="PATTERN",
                        help="Only sync files matching PATTERN, e.g. '-of a.sh b.sh' or "
                             "'-of *.inp' (repeatable)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be transferred without doing it")
    parser.add_argument("--delete", action="store_true",
                        help="Delete destination files absent from the source")
    parser.add_argument("--skip-empty-dirs", action="store_true",
                        help="Skip empty directories (passes --prune-empty-dirs to rsync)")
    args = parser.parse_args()
    # --only-files uses append+nargs='+': flatten [[a,b],[c]] → [a,b,c]
    args.only_files = [p for group in args.only_files for p in group]

    # Build remote target: if -u given, prepend user@ (stripping any existing user@ from --host)
    host = args.host.split("@")[-1] if args.user else args.host
    if args.user:
        host = f"{args.user}@{host}"

    # If --host-dir was written with ~ but the shell expanded it to the local home,
    # convert it back so rsync expands it on the remote instead
    local_home = os.path.expanduser("~")
    host_dir = args.host_dir
    if host_dir.startswith(local_home):
        host_dir = "~" + host_dir[len(local_home):]

    remote = f"{host}:{host_dir}"

    if args.upload:
        print(f"Uploading  {args.local_dir}/ → {remote}/")
        sync(args.local_dir, remote, dry_run=args.dry_run, delete=args.delete, exclude=args.exclude, only_files=args.only_files, skip_empty_dirs=args.skip_empty_dirs)
    else:
        print(f"Downloading {remote}/ → {args.local_dir}/")
        sync(remote, args.local_dir, dry_run=args.dry_run, delete=args.delete, exclude=args.exclude, only_files=args.only_files, skip_empty_dirs=args.skip_empty_dirs)


if __name__ == "__main__":
    main()
