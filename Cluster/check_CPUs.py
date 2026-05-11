#!/usr/bin/env python3
"""
Check CPU availability on a local or remote Linux machine.

Reads /proc/cpuinfo and /proc/stat for topology and live usage; shows
load averages and top processes by CPU. Useful for deciding how many
--workers to pass to a pipeline before launching it.

Usage:
    check_CPUs.py                     local machine
    check_CPUs.py --host HOST         query HOST via SSH (e.g. qtcovi)
    check_CPUs.py --top N             show top N processes (default: 10)
    check_CPUs.py --watch [SECONDS]   auto-refresh every N seconds (default: 10)
"""

import argparse
import subprocess
import sys
import time

_WIDTH   = 66
_DIVIDER = "+" + "-" * _WIDTH + "+"

# Single shell command: collect everything in one round-trip.
# Two /proc/stat reads 1 second apart for accurate CPU usage %.
_COLLECT = (
    "cat /proc/cpuinfo; echo '===SEP==='; "
    "cat /proc/loadavg; echo '===SEP==='; "
    "cat /proc/stat; echo '===SEP==='; "
    "sleep 1; cat /proc/stat; echo '===SEP==='; "
    "ps aux --sort=-%cpu --no-header"
)


def collect(host, top_n):
    script = _COLLECT.replace(
        "ps aux --sort=-%cpu --no-header",
        f"ps aux --sort=-%cpu --no-header | head -{top_n}",
    )
    cmd = ["ssh", host, script] if host else ["bash", "-c", script]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    if result.returncode != 0 and not result.stdout:
        print(f"ERROR: {result.stderr.strip()}", file=sys.stderr)
        sys.exit(1)
    return result.stdout


def parse_cpuinfo(text):
    physical_ids = set()
    core_ids = {}
    cur_phys = None
    logical  = 0

    for line in text.splitlines():
        line = line.strip()
        if line.startswith("processor"):
            logical += 1
        elif line.startswith("physical id"):
            cur_phys = line.split(":")[1].strip()
            physical_ids.add(cur_phys)
            core_ids.setdefault(cur_phys, set())
        elif line.startswith("core id") and cur_phys is not None:
            core_ids[cur_phys].add(line.split(":")[1].strip())

    physical_cores = sum(len(v) for v in core_ids.values()) if core_ids else logical
    return {
        "logical":        logical,
        "physical_cores": physical_cores,
        "sockets":        max(len(physical_ids), 1),
        "ht":             logical > physical_cores,
    }


def parse_loadavg(text):
    parts = text.strip().split()
    return float(parts[0]), float(parts[1]), float(parts[2])


def parse_stat(text):
    for line in text.splitlines():
        if line.startswith("cpu "):
            return list(map(int, line.split()[1:]))
    return None


def cpu_usage_pct(stat1, stat2):
    total = sum(stat2) - sum(stat1)
    idle  = stat2[3]   - stat1[3]
    return 100.0 * (1.0 - idle / total) if total else 0.0


def parse_ps(text, n):
    rows = []
    for line in text.splitlines()[:n]:
        parts = line.split(None, 10)
        if len(parts) < 11:
            continue
        rows.append((parts[0], parts[1], parts[2], parts[3], parts[10].strip()[:45]))
    return rows


def suggest_workers(physical_cores, ht, load1):
    # Physical cores are the practical limit for CPU-bound work.
    # Subtract current load, apply 90% safety margin.
    headroom  = max(0.0, physical_cores - load1)
    suggested = max(1, int(headroom * 0.9))
    return min(suggested, physical_cores)


def render(cpuinfo, load1, load5, load15, cpu_pct, ps_rows):
    logical  = cpuinfo["logical"]
    physical = cpuinfo["physical_cores"]
    ht       = cpuinfo["ht"]
    sockets  = cpuinfo["sockets"]
    suggested = suggest_workers(physical, ht, load1)

    lines = [
        _DIVIDER,
        f"| {'CPU Availability':<{_WIDTH - 1}}|",
        _DIVIDER,
        f"  Logical CPUs   : {logical}",
        f"  Physical cores : {physical}"
        + (f"  ({sockets} socket{'s' if sockets > 1 else ''}, HT enabled)" if ht
           else f"  ({sockets} socket{'s' if sockets > 1 else ''})"),
        "",
        f"  Load avg (1/5/15m) : {load1:.2f}  {load5:.2f}  {load15:.2f}",
        f"  CPU usage (1s)     : {cpu_pct:.1f}%",
        "",
        f"  Suggested --workers : {suggested}  (physical cores - load, 90% margin)",
        "",
        f"  {'USER':<12} {'PID':>7} {'%CPU':>6} {'%MEM':>6}  COMMAND",
        f"  {'-'*12} {'-'*7} {'-'*6} {'-'*6}  {'-'*35}",
    ]
    for user, pid, cpu, mem, cmd in ps_rows:
        lines.append(f"  {user:<12} {pid:>7} {cpu:>6} {mem:>6}  {cmd}")

    return "\n".join(lines)


def run_once(host, top_n):
    raw   = collect(host, top_n)
    parts = raw.split("===SEP===\n")
    if len(parts) < 5:
        print("ERROR: unexpected output — check SSH connection or /proc availability.")
        sys.exit(1)

    cpuinfo_txt, loadavg_txt, stat1_txt, stat2_txt, ps_txt = parts[:5]

    cpuinfo = parse_cpuinfo(cpuinfo_txt)
    load1, load5, load15 = parse_loadavg(loadavg_txt)
    stat1   = parse_stat(stat1_txt)
    stat2   = parse_stat(stat2_txt)
    cpu_pct = cpu_usage_pct(stat1, stat2) if stat1 and stat2 else 0.0
    ps_rows = parse_ps(ps_txt, top_n)

    print(render(cpuinfo, load1, load5, load15, cpu_pct, ps_rows))


def main():
    parser = argparse.ArgumentParser(
        description="Check CPU availability on a local or remote Linux machine.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--host", metavar="HOST",
                        help="Query HOST via SSH (e.g. qtcovi)")
    parser.add_argument("--top", type=int, default=10, metavar="N",
                        help="Show top N processes by CPU (default: 10)")
    parser.add_argument("--watch", nargs="?", const=10, type=int, metavar="SECONDS",
                        help="Auto-refresh every N seconds (default: 10)")
    args = parser.parse_args()

    if args.watch is not None:
        interval = args.watch
        try:
            while True:
                print("\033[2J\033[H", end="", flush=True)
                ts = time.strftime("%Y-%m-%d %H:%M:%S")
                print(f"check_CPUs  [{ts}]  refreshing every {interval}s  (Ctrl+C to quit)\n")
                run_once(args.host, args.top)
                time.sleep(interval)
        except KeyboardInterrupt:
            print("\nStopped.")
    else:
        run_once(args.host, args.top)


if __name__ == "__main__":
    main()
