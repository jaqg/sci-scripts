#!/usr/bin/env python3
"""
Monitor PBS job progress on an HPC cluster.

Usage:
    qstat-info.py <job_id> [<job_id> ...]      query specific job IDs
    qstat-info.py -u [USER]                    query all jobs for USER (default: $USER)
    qstat-info.py -u --host scc               run qstat/cat via SSH on host 'scc'
    qstat-info.py -u --watch [SECONDS]        auto-refresh every N seconds (default: 30)
"""

import argparse
import os
import shlex
import subprocess
import sys
import time


# ---------------------------------------------------------------------------
# Command execution (local or via SSH)
# ---------------------------------------------------------------------------

def make_runner(host):
    """Return a run_cmd(args_list) callable that executes locally or via SSH."""
    def run_cmd(args_list):
        if host:
            cmd = ["ssh", host, shlex.join(args_list)]
        else:
            cmd = args_list
        result = subprocess.run(cmd, capture_output=True, text=True)
        return result.stdout  # empty string on error; callers handle missing data
    return run_cmd


# ---------------------------------------------------------------------------
# qstat parsers
# ---------------------------------------------------------------------------

def get_job_ids_for_user(user, run_cmd):
    output = run_cmd(["qstat", "-u", user])
    ids = []
    for i, line in enumerate(output.splitlines()):
        if i < 5:
            continue
        parts = line.split()
        if parts:
            ids.append(parts[0])
    return ids


_FIELDS_OF_INTEREST = {
    "Job_Name", "resources_used.walltime", "job_state", "exec_host", "Output_Path"
}


def _join_continuation_lines(text):
    """PBS qstat -f continues long values on lines starting with a literal tab."""
    joined = []
    for line in text.splitlines():
        if line.startswith("\t") and joined:
            joined[-1] = joined[-1].rstrip() + line[1:].rstrip()
        else:
            joined.append(line.rstrip())
    return joined


def parse_job_details(job_id, run_cmd):
    """Return a dict of selected qstat -f fields, or empty dict if not found."""
    raw = run_cmd(["qstat", "-f", job_id])
    fields = {}
    for line in _join_continuation_lines(raw):
        stripped = line.strip()
        if " = " in stripped:
            key, _, value = stripped.partition(" = ")
            key = key.strip()
            if key in _FIELDS_OF_INTEREST:
                fields[key] = value.strip()
    return fields


def _output_dir(fields):
    """Derive the calculation directory from the Output_Path field."""
    raw = fields.get("Output_Path", "")
    # Strip "hostname:" prefix (e.g. "scc_ser_1:/path/..." → "/path/...")
    path = raw.split(":", 1)[1] if ":" in raw else raw
    return os.path.dirname(path)


# ---------------------------------------------------------------------------
# progress.log / FILES parsers
# ---------------------------------------------------------------------------

def read_files_list(output_dir, run_cmd):
    content = run_cmd(["cat", f"{output_dir}/FILES"])
    return [l for l in content.splitlines() if l.strip()]


def parse_progress_log(output_dir, job_id, run_cmd):
    content = run_cmd(["cat", f"{output_dir}/progress_{job_id}.log"])
    running, completed, failed = [], [], []
    for line in content.splitlines():
        if line.startswith("Running: "):
            running.append(line[len("Running: "):].strip())
        elif line.startswith("Completed: "):
            completed.append(line[len("Completed: "):].strip())
        elif line.startswith("Failed: "):
            failed.append(line[len("Failed: "):].strip())
    return {"running": running, "completed": completed, "failed": failed}


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

_WIDTH = 62
_DIVIDER = "+" + "-" * _WIDTH + "+"


def _progress_bar(n_completed, n_total, n_failed, bar_width=28):
    if n_total == 0:
        return "[no files]"
    frac = n_completed / n_total
    filled = int(frac * bar_width)
    if filled < bar_width:
        bar = "=" * filled + ">" + " " * (bar_width - filled - 1)
    else:
        bar = "=" * bar_width
    pct = int(frac * 100)
    s = f"[{bar}] {n_completed}/{n_total} ({pct}%)"
    if n_failed:
        s += f"  |  {n_failed} failed"
    return s


def render_job_block(job_id, run_cmd):
    lines = [_DIVIDER, f"| Job: {job_id:<{_WIDTH - 6}}|", _DIVIDER]

    fields = parse_job_details(job_id, run_cmd)
    if not fields:
        lines.append(f"  No information found for job {job_id}")
        return "\n".join(lines)

    state    = fields.get("job_state", "?")
    walltime = fields.get("resources_used.walltime", "N/A")
    host     = fields.get("exec_host", "N/A")
    name     = fields.get("Job_Name", "N/A")
    out_dir  = _output_dir(fields)

    lines += [
        f"  Name:      {name}",
        f"  State:     {state}",
        f"  Walltime:  {walltime}",
        f"  Host:      {host}",
        f"  Directory: {out_dir}",
    ]

    if not out_dir:
        lines.append("  (Could not determine output directory)")
        return "\n".join(lines)

    files_list = read_files_list(out_dir, run_cmd)
    if not files_list:
        lines.append(f"\n  [FILES not found in {out_dir}]")
        return "\n".join(lines)

    n_total = len(files_list)
    lines.append(f"\n  FILES ({n_total} calculations):")
    for f in files_list:
        lines.append(f"    {f}")

    progress = parse_progress_log(out_dir, job_id, run_cmd)
    completed = progress["completed"]
    failed    = progress["failed"]
    running   = progress["running"]

    if not (completed or failed or running):
        lines.append(f"\n  [progress_{job_id}.log not found or empty]")
        return "\n".join(lines)

    lines.append(f"\n  Progress: {_progress_bar(len(completed), n_total, len(failed))}")

    if running:
        lines.append(f"\n  Running ({len(running)}):")
        for f in running:
            lines.append(f"    {f}")

    if completed:
        lines.append(f"\n  Completed ({len(completed)}):")
        for f in completed:
            lines.append(f"    {f}")

    if failed:
        lines.append(f"\n  Failed ({len(failed)}):")
        for f in failed:
            lines.append(f"    {f}")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _get_job_ids(args, run_cmd):
    if args.user is not None:
        user = args.user or os.environ.get("USER", "")
        return get_job_ids_for_user(user, run_cmd)
    return list(args.job_ids)


def main():
    parser = argparse.ArgumentParser(
        description="Monitor PBS job progress on an HPC cluster.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "job_ids", nargs="*",
        help="Job IDs to query (e.g. 152689.scc_ser_1)",
    )
    parser.add_argument(
        "-u", "--user", nargs="?", const="", metavar="USER",
        help="Query all jobs for USER (omit value to use $USER)",
    )
    parser.add_argument(
        "--host", metavar="HOST",
        help="Run qstat/cat via SSH on HOST (e.g. scc)",
    )
    parser.add_argument(
        "--watch", nargs="?", const=30, type=int, metavar="SECONDS",
        help="Auto-refresh every N seconds (default: 30)",
    )
    args = parser.parse_args()

    if args.user is None and not args.job_ids:
        args.user = ""  # default to $USER

    run_cmd = make_runner(args.host)

    def run_once():
        job_ids = _get_job_ids(args, run_cmd)
        if not job_ids:
            user = args.user or os.environ.get("USER", "")
            print(f"No jobs found for user '{user}'.")
            return
        for job_id in job_ids:
            print(render_job_block(job_id, run_cmd))
            print()

    if args.watch is not None:
        interval = args.watch
        try:
            while True:
                print("\033[2J\033[H", end="", flush=True)
                ts = time.strftime("%Y-%m-%d %H:%M:%S")
                print(f"qstat-info  [{ts}]  refreshing every {interval}s  (Ctrl+C to quit)\n")
                run_once()
                time.sleep(interval)
        except KeyboardInterrupt:
            print("\nStopped.")
    else:
        run_once()


if __name__ == "__main__":
    main()
