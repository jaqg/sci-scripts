#!/usr/bin/env python3
"""
Check hardware specifications for running local AI models.

Reports CPU, RAM, disk, and GPU (NVIDIA/AMD/Intel) in one pass.
Includes a model-fit table based on detected VRAM or RAM (CPU-only).

Usage:
    check_specifications.py                 local machine
    check_specifications.py --host HOST     query HOST via SSH (e.g. qtcovi)
"""

import argparse
import subprocess
import sys

_WIDTH   = 70
_DIVIDER = "+" + "-" * _WIDTH + "+"

_COLLECT = (
    # CPU model
    "grep -m1 'model name' /proc/cpuinfo 2>/dev/null || echo ''; echo '===SEP==='; "
    # logical CPUs, physical cores
    "cat /proc/cpuinfo; echo '===SEP==='; "
    # RAM
    "cat /proc/meminfo; echo '===SEP==='; "
    # Disk (root)
    "df -h / 2>/dev/null; echo '===SEP==='; "
    # NVIDIA GPU
    "nvidia-smi --query-gpu=name,memory.total,memory.free,utilization.gpu "
    "--format=csv,noheader 2>/dev/null || echo 'NO_NVIDIA'; echo '===SEP==='; "
    # AMD GPU
    "rocm-smi --showmeminfo vram --csv 2>/dev/null || echo 'NO_AMD'; echo '===SEP==='; "
    # Intel GPU / generic fallback
    "lspci 2>/dev/null | grep -iE 'vga|3d|display' || echo 'NO_LSPCI'"
)


def collect(host):
    cmd = ["ssh", host, _COLLECT] if host else ["bash", "-c", _COLLECT]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    if result.returncode != 0 and not result.stdout:
        print(f"ERROR: {result.stderr.strip()}", file=sys.stderr)
        sys.exit(1)
    return result.stdout


def parse_cpu(cpumodel_txt, cpuinfo_txt):
    model = cpumodel_txt.strip()
    if ":" in model:
        model = model.split(":", 1)[1].strip()

    physical_ids = set()
    core_ids = {}
    cur_phys = None
    logical = 0
    for line in cpuinfo_txt.splitlines():
        line = line.strip()
        if line.startswith("processor"):
            logical += 1
        elif line.startswith("physical id"):
            cur_phys = line.split(":")[1].strip()
            physical_ids.add(cur_phys)
            core_ids.setdefault(cur_phys, set())
        elif line.startswith("core id") and cur_phys is not None:
            core_ids[cur_phys].add(line.split(":")[1].strip())

    physical = sum(len(v) for v in core_ids.values()) if core_ids else logical
    return {"model": model or "unknown", "logical": logical, "physical": physical}


def parse_ram(meminfo_txt):
    total_kb = free_kb = avail_kb = 0
    for line in meminfo_txt.splitlines():
        if line.startswith("MemTotal:"):
            total_kb = int(line.split()[1])
        elif line.startswith("MemFree:"):
            free_kb = int(line.split()[1])
        elif line.startswith("MemAvailable:"):
            avail_kb = int(line.split()[1])
    return {
        "total_gb": total_kb / 1024**2,
        "free_gb":  free_kb  / 1024**2,
        "avail_gb": avail_kb / 1024**2,
    }


def parse_disk(df_txt):
    lines = [l for l in df_txt.splitlines() if l and not l.startswith("Filesystem")]
    if not lines:
        return None
    parts = lines[0].split()
    if len(parts) < 6:
        return None
    return {"total": parts[1], "used": parts[2], "avail": parts[3], "pct": parts[4]}


def parse_nvidia(nvidia_txt):
    if nvidia_txt.strip() == "NO_NVIDIA" or not nvidia_txt.strip():
        return []
    gpus = []
    for line in nvidia_txt.strip().splitlines():
        parts = [p.strip() for p in line.split(",")]
        if len(parts) < 4:
            continue
        name, total, free, util = parts[:4]
        def to_mb(s):
            s = s.replace("MiB", "").replace("MB", "").strip()
            try:
                return int(float(s))
            except ValueError:
                return 0
        gpus.append({
            "name": name,
            "vram_total_mb": to_mb(total),
            "vram_free_mb":  to_mb(free),
            "util_pct":      util,
        })
    return gpus


def parse_amd(rocm_txt):
    if rocm_txt.strip() == "NO_AMD" or not rocm_txt.strip():
        return []
    gpus = []
    for line in rocm_txt.strip().splitlines():
        if "VRAM" in line.upper() and "Total" in line:
            parts = line.split(",")
            if len(parts) >= 3:
                name = parts[0].strip()
                try:
                    total_mb = int(parts[1].strip())
                    free_mb  = int(parts[2].strip())
                    gpus.append({"name": name, "vram_total_mb": total_mb, "vram_free_mb": free_mb, "util_pct": "N/A"})
                except ValueError:
                    pass
    return gpus


def parse_lspci(lspci_txt):
    if lspci_txt.strip() == "NO_LSPCI":
        return []
    return [l.strip() for l in lspci_txt.strip().splitlines() if l.strip()]


def model_fit_table(vram_gb=None, ram_gb=None):
    rows = [
        ("3B–7B  Q4",  4,   "Llama 3.2 3B, Mistral 7B, Qwen2.5 7B"),
        ("7B–13B Q4",  8,   "Llama 3.1 8B, Mistral 7B, CodeLlama 13B"),
        ("13B–34B Q4", 16,  "CodeLlama 34B, Mixtral 8x7B (partial)"),
        ("70B Q4",     24,  "Llama 3 70B, Qwen2.5 72B"),
        ("70B Q8",     48,  "High-quality 70B inference"),
        ("405B Q4",    96,  "Llama 3.1 405B (multi-GPU)"),
    ]
    budget = vram_gb if vram_gb else ram_gb
    label  = "VRAM" if vram_gb else "RAM (CPU-only — slow)"
    fits   = [(tier, req, ex) for tier, req, ex in rows if budget >= req] if budget else []

    lines = [
        "",
        f"  Model fit  ({label}: {budget:.0f} GB available):",
        f"  {'Tier':<14} {'Min GB':>6}  Examples",
        f"  {'-'*14} {'-'*6}  {'-'*40}",
    ]
    if fits:
        for tier, req, ex in fits:
            marker = " <-- max fit" if (tier, req, ex) == fits[-1] else ""
            lines.append(f"  {tier:<14} {req:>6}  {ex}{marker}")
    else:
        lines.append("  No standard model fits detected memory.")
    return "\n".join(lines)


def render(cpu, ram, disk, nvidia_gpus, amd_gpus, lspci_lines):
    lines = [
        _DIVIDER,
        f"| {'Hardware Specifications':<{_WIDTH - 1}}|",
        _DIVIDER,
        "",
        "  CPU",
        f"    Model   : {cpu['model']}",
        f"    Logical : {cpu['logical']}  /  Physical cores: {cpu['physical']}",
        "",
        "  RAM",
        f"    Total     : {ram['total_gb']:.1f} GB",
        f"    Available : {ram['avail_gb']:.1f} GB  (free: {ram['free_gb']:.1f} GB)",
    ]

    if disk:
        lines += [
            "",
            "  Disk  (/)",
            f"    Total : {disk['total']}  Used: {disk['used']} ({disk['pct']})  Free: {disk['avail']}",
        ]

    all_gpus = nvidia_gpus + amd_gpus
    if all_gpus:
        lines += ["", "  GPU"]
        for g in all_gpus:
            vram_total = g["vram_total_mb"] / 1024
            vram_free  = g["vram_free_mb"]  / 1024
            lines.append(f"    {g['name']}")
            lines.append(f"      VRAM : {vram_total:.1f} GB total  /  {vram_free:.1f} GB free")
            lines.append(f"      Util : {g['util_pct']}")
        best_vram = max(g["vram_free_mb"] for g in all_gpus) / 1024
        lines.append(model_fit_table(vram_gb=best_vram))
    else:
        if lspci_lines:
            lines += ["", "  GPU (no NVIDIA/AMD driver — lspci only)"]
            for l in lspci_lines:
                lines.append(f"    {l}")
        else:
            lines += ["", "  GPU : none detected"]
        lines.append(model_fit_table(ram_gb=ram["avail_gb"]))

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Check hardware specifications for running local AI models.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--host", metavar="HOST",
                        help="Query HOST via SSH (e.g. qtcovi)")
    args = parser.parse_args()

    raw   = collect(args.host)
    parts = raw.split("===SEP===\n")
    if len(parts) < 7:
        print("ERROR: unexpected output — check SSH connection or command availability.")
        sys.exit(1)

    cpumodel_txt, cpuinfo_txt, meminfo_txt, df_txt, nvidia_txt, amd_txt, lspci_txt = parts[:7]

    cpu        = parse_cpu(cpumodel_txt, cpuinfo_txt)
    ram        = parse_ram(meminfo_txt)
    disk       = parse_disk(df_txt)
    nvidia_gpus = parse_nvidia(nvidia_txt)
    amd_gpus   = parse_amd(amd_txt)
    lspci_lines = parse_lspci(lspci_txt)

    print(render(cpu, ram, disk, nvidia_gpus, amd_gpus, lspci_lines))


if __name__ == "__main__":
    main()
