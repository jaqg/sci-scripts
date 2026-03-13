#!/usr/bin/env python3
"""
Query PBS node status on an HPC cluster.

Usage:
    pbsnodes-info.py                     show all nodes
    pbsnodes-info.py -f                  show only free nodes
    pbsnodes-info.py --host HOST         run pbsnodes via SSH on HOST
    pbsnodes-info.py --watch [SECONDS]   auto-refresh every N seconds (default: 30)
"""

import argparse
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
        return result.stdout
    return run_cmd


# ---------------------------------------------------------------------------
# pbsnodes parser
# ---------------------------------------------------------------------------

def parse_pbsnodes(run_cmd):
    """Parse pbsnodes -a output into a list of node dicts."""
    raw = run_cmd(["pbsnodes", "-a"])
    nodes = []
    current = None

    for line in raw.splitlines():
        # Node name lines have no leading whitespace
        if line and not line[0].isspace():
            if current is not None:
                nodes.append(current)
            current = {"name": line.strip(), "state": "?", "np": 0, "jobs": []}
        elif current is not None and " = " in line:
            key, _, value = line.strip().partition(" = ")
            key   = key.strip()
            value = value.strip()
            if key == "state":
                current["state"] = value
            elif key == "np":
                try:
                    current["np"] = int(value)
                except ValueError:
                    pass
            elif key == "jobs":
                current["jobs"] = [j.strip() for j in value.split(",") if j.strip()]

    if current is not None:
        nodes.append(current)
    return nodes


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

_WIDTH   = 62
_DIVIDER = "+" + "-" * _WIDTH + "+"


def _state_tag(state):
    if state == "free":
        return "FREE"
    if "job" in state:
        return "BUSY"
    return state.upper()[:6]


def render(nodes, free_only=False):
    free_nodes  = [n for n in nodes if n["state"] == "free"]
    busy_nodes  = [n for n in nodes if "job" in n["state"]]
    down_nodes  = [n for n in nodes if n["state"] in ("down", "offline")]
    total_cores = sum(n["np"] for n in nodes)
    used_cores  = sum(len(n["jobs"]) for n in nodes)

    lines = [
        _DIVIDER,
        f"| {'PBS Node Status':<{_WIDTH - 1}}|",
        _DIVIDER,
        f"  Nodes : {len(free_nodes)} free  {len(busy_nodes)} busy  {len(down_nodes)} down"
        f"  (total: {len(nodes)})",
        f"  Cores : {total_cores - used_cores} free  {used_cores} used  (total: {total_cores})",
        "",
    ]

    display = free_nodes if free_only else nodes
    if not display:
        lines.append("  (no nodes to display)")
        return "\n".join(lines)

    name_w = max(len(n["name"]) for n in display)

    for node in display:
        n_used = len(node["jobs"])
        n_free = node["np"] - n_used
        tag    = _state_tag(node["state"])
        lines.append(
            f"  {node['name']:<{name_w}}  [{tag:<6}]  {n_free:>3}/{node['np']:<3} cores free"
        )
        if node["jobs"]:
            job_ids = sorted(set(j.split("/")[1] if "/" in j else j for j in node["jobs"]))
            lines.append(f"  {'':<{name_w}}           jobs: {', '.join(job_ids)}")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Query PBS node status on an HPC cluster.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "-f", "--free", action="store_true",
        help="Show only free nodes",
    )
    parser.add_argument(
        "--host", metavar="HOST",
        help="Run pbsnodes via SSH on HOST (e.g. scc)",
    )
    parser.add_argument(
        "--watch", nargs="?", const=30, type=int, metavar="SECONDS",
        help="Auto-refresh every N seconds (default: 30)",
    )
    args = parser.parse_args()

    run_cmd = make_runner(args.host)

    def run_once():
        nodes = parse_pbsnodes(run_cmd)
        if not nodes:
            print("No node information returned.")
            return
        print(render(nodes, free_only=args.free))

    if args.watch is not None:
        interval = args.watch
        try:
            while True:
                print("\033[2J\033[H", end="", flush=True)
                ts = time.strftime("%Y-%m-%d %H:%M:%S")
                print(f"pbsnodes-info  [{ts}]  refreshing every {interval}s  (Ctrl+C to quit)\n")
                run_once()
                time.sleep(interval)
        except KeyboardInterrupt:
            print("\nStopped.")
    else:
        run_once()


if __name__ == "__main__":
    main()
