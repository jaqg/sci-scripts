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
    -u  / --user    Remote username (default: user already in --host, or SSH config)
    --dry-run       Show what would be transferred without doing it
    --delete        Delete files on the destination that are absent on the source
"""

import argparse
import os
import subprocess
import sys


def sync(src, dst, *, dry_run=False, delete=False):
    cmd = ["rsync", "-avz", "--progress"]
    if dry_run:
        cmd.append("--dry-run")
    if delete:
        cmd.append("--delete")
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

    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be transferred without doing it")
    parser.add_argument("--delete", action="store_true",
                        help="Delete destination files absent from the source")
    args = parser.parse_args()

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
        sync(args.local_dir, remote, dry_run=args.dry_run, delete=args.delete)
    else:
        print(f"Downloading {remote}/ → {args.local_dir}/")
        sync(remote, args.local_dir, dry_run=args.dry_run, delete=args.delete)


if __name__ == "__main__":
    main()
