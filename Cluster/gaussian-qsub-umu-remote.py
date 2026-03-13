#!/usr/bin/env python3
"""
Submit Gaussian calculations from a local machine to a remote PBS cluster.

Copies the .com files listed in FILE_LIST and the submission script to the
remote host, then launches gaussian-qsub-umu.sh there.

To retrieve results afterwards use folder-sync.py --down.

Usage:
    gaussian-qsub-umu-remote.py -H HOST [options]

Required:
    -H / --host     SSH host to submit to (e.g. cluster, or user@host)

Options:
    -d / --dir      Remote working directory
                    (default: $HOME/remote-calculations-tmp)
    -f / --files    Local file list (default: FILES)
    -N              PBS job name            (default: numders-fw)
    -n              Cluster node            (default: nv81)
    -p              Cores per node / ppn    (default: 32)
    -c              Cores per calculation   (default: 8)
    -m              Max parallel calcs      (default: 4)
"""

import argparse
import os
import subprocess
import sys

# gaussian-qsub-umu.sh lives next to this script
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_QSUB_SCRIPT = os.path.join(_SCRIPT_DIR, "gaussian-qsub-umu.sh")
_REMOTE_SCRIPT_NAME = "gaussian-qsub-umu.sh"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(cmd, *, input=None):
    """Print and execute a command; exit on failure."""
    print("  $", " ".join(cmd))
    result = subprocess.run(cmd, input=input, text=True)
    if result.returncode != 0:
        print(f"Error: command exited with code {result.returncode}", file=sys.stderr)
        sys.exit(result.returncode)


def read_file_list(path):
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip()]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Submit Gaussian calculations from local machine to a remote PBS cluster.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("-H", "--host", required=True, metavar="HOST",
                        help="SSH host (e.g. cluster, or user@hostname)")
    parser.add_argument("-d", "--dir", default="$HOME/remote-calculations-tmp",
                        metavar="REMOTE_DIR",
                        help="Remote working directory (default: $HOME/remote-calculations-tmp)")
    parser.add_argument("-f", "--files", default="FILES", metavar="FILE_LIST",
                        help="Local file list (default: FILES)")
    parser.add_argument("-N", dest="job_name", default="numders-fw", metavar="JOB_NAME",
                        help="PBS job name (default: numders-fw)")
    parser.add_argument("-n", dest="node", default="nv81", metavar="NODE",
                        help="Cluster node (default: nv81)")
    parser.add_argument("-p", dest="ppn", default=32, type=int, metavar="PPN",
                        help="Cores per node (default: 32)")
    parser.add_argument("-c", dest="cores_per_calc", default=8, type=int, metavar="CORES",
                        help="Cores per Gaussian calculation (default: 8)")
    parser.add_argument("-m", dest="max_parallel", default=4, type=int, metavar="MAX",
                        help="Max parallel calculations (default: 4)")
    args = parser.parse_args()

    # --- Validate local files ---
    if not os.path.isfile(args.files):
        print(f"Error: file list '{args.files}' not found.", file=sys.stderr)
        sys.exit(1)

    com_files = read_file_list(args.files)
    if not com_files:
        print(f"Error: '{args.files}' is empty.", file=sys.stderr)
        sys.exit(1)

    missing = [f for f in com_files if not os.path.isfile(f)]
    if missing:
        print("Error: the following input files are missing locally:", file=sys.stderr)
        for f in missing:
            print(f"  {f}", file=sys.stderr)
        sys.exit(1)

    if not os.path.isfile(_QSUB_SCRIPT):
        print(f"Error: submission script not found at {_QSUB_SCRIPT}", file=sys.stderr)
        sys.exit(1)

    # --- Step 1: create remote directory ---
    print(f"\n[1/3] Creating remote directory {args.host}:{args.dir}")
    run(["ssh", args.host, f"mkdir -p {args.dir}"])

    # --- Step 2: copy files ---
    # Upload the .com files listed in FILES and the submission script.
    # FILES is rewritten on the remote with basenames only, since all files
    # land flat in remote_dir regardless of their local paths.
    print(f"\n[2/3] Copying files to {args.host}:{args.dir}/")
    run(["rsync", "-avz", "--progress"] + com_files + [_QSUB_SCRIPT] + [f"{args.host}:{args.dir}/"])

    remote_file_list = os.path.basename(args.files)
    basenames = "\n".join(os.path.basename(f) for f in com_files) + "\n"
    run(["ssh", args.host, f"cat > {args.dir}/{remote_file_list}"], input=basenames)

    # --- Step 3: submit ---
    print(f"\n[3/3] Submitting job on {args.host}")
    remote_script = f"{args.dir}/{_REMOTE_SCRIPT_NAME}"
    submit_cmd = (
        f"cd {args.dir} && bash {remote_script}"
        f" -N {args.job_name}"
        f" -n {args.node}"
        f" -p {args.ppn}"
        f" -f {remote_file_list}"
        f" -c {args.cores_per_calc}"
        f" -m {args.max_parallel}"
    )
    run(["ssh", args.host, submit_cmd])

    print(f"\nDone. To retrieve results:\n"
          f"  folder-sync.py -H {args.host} -hd {args.dir} -ld <local_output_dir> --down")


if __name__ == "__main__":
    main()
