#!/usr/bin/env python3
"""
Estimate progress of a running Gaussian calculation from its .log file.

Usage:
    python3 gaussian-progress.py job.log            # one-shot snapshot
    python3 gaussian-progress.py job.log --watch    # real-time monitor (2 s refresh)
    python3 gaussian-progress.py job.log -w -n 5    # refresh every 5 seconds
    python3 gaussian-progress.py job.log --expected-steps 10

Supports: SCF, Opt, Opt=TD, Freq, and multi-step (--Link1--) jobs.
"""

import sys
import re
import os
import time
import argparse
from datetime import timedelta


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def tail_lines(filepath, n=500):
    """Read the last n lines of a file (efficient for large files)."""
    try:
        with open(filepath) as f:
            f.seek(0, 2)  # EOF
            size = f.tell()
            # Start reading from estimated position
            pos = max(0, size - 8192)
            while pos > 0:
                f.seek(pos)
                chunk = f.read(size - pos)
                if chunk.count('\n') >= n:
                    break
                pos = max(0, pos - 8192)
            f.seek(pos)
            lines = f.read().split('\n')
            return lines[-n:]
    except (IOError, OSError):
        return []


def find_all_matches(lines, pattern, group=0):
    """Return all regex matches across lines."""
    results = []
    for i, line in enumerate(lines):
        m = re.search(pattern, line)
        if m:
            results.append((i, m.group(group) if group else m))
    return results


# ---------------------------------------------------------------------------
# Phase detection and progress parsing
# ---------------------------------------------------------------------------

def parse_progress(filepath):
    """Parse a Gaussian .log file and return a progress dict.
    Returns None if file not readable or empty."""
    try:
        with open(filepath) as f:
            content = f.read()
    except (IOError, OSError):
        return None

    if not content.strip():
        return None

    lines = content.split('\n')
    n_total = len(lines)

    info = {
        'filename': os.path.basename(filepath),
        'total_lines': n_total,
        'phase': 'unknown',
        'phase_detail': '',
        'step_current': 0,
        'step_max': 0,
        'scf_cycle': 0,
        'scf_energy': None,
        'scf_converged': False,
        'max_force': None,
        'rms_force': None,
        'max_disp': None,
        'rms_disp': None,
        'forces_converged': False,
        'disps_converged': False,
        'td_in_progress': False,
        'td_nstates': None,
        'opt_converged': False,
        'freq_started': False,
        'job_finished': False,
        'elapsed_wall': None,
        'percent': 0.0,
    }

    # --- Detect job completion ---
    # "Normal termination" appears at the end of each Link1 step too;
    # only true completion is when it appears AFTER the last Link1 (or no Link1).
    for i, line in enumerate(lines):
        if 'Normal termination of Gaussian' in line:
            # Check if followed by another Link1 step
            has_next_step = any('Proceeding to internal job step' in lines[j]
                               for j in range(i, min(i+20, len(lines))))
            if not has_next_step:
                info['job_finished'] = True
                info['phase'] = 'complete'
                info['phase_detail'] = 'Job finished successfully'
                info['percent'] = 100.0

    # --- Detect multi-step boundaries ---
    link1_boundaries = []
    for i, line in enumerate(lines):
        if 'Proceeding to internal job step number' in line:
            link1_boundaries.append(i)
    current_step_idx = len(link1_boundaries)  # 0 = step 1, 1 = step 2, etc.

    # --- If in step 2+, freq is likely started ---
    if link1_boundaries:
        info['freq_started'] = True
    
    # Check if step 2 route contains 'Freq'
    for i in link1_boundaries:
        for j in range(i, min(i+10, len(lines))):
            if re.search(r'[Ff]req', lines[j]):
                info['freq_started'] = True
                break

    # --- Scan for key markers (from end, most relevant first) ---

    # Find last "Optimization completed"
    opt_done = None
    for i in range(len(lines) - 1, -1, -1):
        if 'Optimization completed' in lines[i]:
            opt_done = i
            break
    if opt_done is not None:
        info['opt_converged'] = True

    # Find last "Stationary point found"
    if any('Stationary point found' in l for l in lines):
        info['opt_converged'] = True

    # Find last "Berny optimization" (marks an optimization step)
    last_berny = None
    for i in range(len(lines) - 1, -1, -1):
        if 'Berny optimization' in lines[i]:
            last_berny = i
            break

    # Find last "Step number X out of a maximum of N" 
    # (scope to the optimization phase, not post-opt freq internal steps)
    last_step_match = None
    if link1_boundaries:
        # Multi-step: find the last job step that contains an optimization,
        # then get its step number
        step_ranges = [0] + link1_boundaries + [len(lines)]
        for s in range(len(step_ranges)-1, 0, -1):
            seg_start = step_ranges[s-1]
            seg_end = step_ranges[s]
            seg_lines = lines[seg_start:seg_end]
            if any('Optimization completed' in l for l in seg_lines):
                for line in reversed(seg_lines):
                    m = re.search(r'Step number\s+(\d+)\s+out of a maximum of\s+(\d+)', line)
                    if m:
                        last_step_match = (int(m.group(1)), int(m.group(2)))
                        break
                if last_step_match:
                    break
    
    if last_step_match is None:
        # Single-step or fallback: find the step number near last Berny
        if last_berny is not None:
            for i in range(last_berny, min(last_berny + 30, len(lines))):
                m = re.search(r'Step number\s+(\d+)\s+out of a maximum of\s+(\d+)', lines[i])
                if m:
                    last_step_match = (int(m.group(1)), int(m.group(2)))
                    break
    
    # Final fallback: search entire file
    if last_step_match is None:
        for line in reversed(lines):
            m = re.search(r'Step number\s+(\d+)\s+out of a maximum of\s+(\d+)', line)
            if m:
                last_step_match = (int(m.group(1)), int(m.group(2)))
                break

    # Find last SCF Done
    last_scf = None
    for line in reversed(lines):
        m = re.search(r'SCF Done:\s+E\([^)]+\)\s*=\s*(\S+)\s+A\.U\.\s+after\s+(\d+)\s+cycles', line)
        if m:
            last_scf = (float(m.group(1)), int(m.group(2)))
            break

    # Find last convergence check
    for line in reversed(lines):
        mf = re.search(r'Maximum Force\s+(\S+)\s+\S+\s+(YES|NO)', line)
        rf = re.search(r'RMS\s+Force\s+(\S+)\s+\S+\s+(YES|NO)', line)
        md = re.search(r'Maximum Displacement\s+(\S+)\s+\S+\s+(YES|NO)', line)
        rd = re.search(r'RMS\s+Displacement\s+(\S+)\s+\S+\s+(YES|NO)', line)
        if mf:
            info['max_force'] = float(mf.group(1))
            info['forces_converged'] = (mf.group(2) == 'YES')
            break  # Use the last set

    # Actually, need to find the LAST full convergence block.
    # Scan for the block pattern: Maximum Force / RMS Force / Max Disp / RMS Disp
    for i in range(len(lines) - 1, -1, -1):
        if 'Maximum Force' in lines[i] and 'RMS' in lines[i] and 'Displacement' in lines[i]:
            # This line contains all four? No, they're separate lines.
            pass
    # Let me do a better scan
    conv_block = {'max_force': None, 'rms_force': None, 'max_disp': None, 'rms_disp': None}
    for i in range(len(lines) - 10, -1, -1):
        if i < 0:
            break
        mf = re.search(r'Maximum Force\s+(\S+)\s+\S+\s+(YES|NO)', lines[i])
        if mf:
            # Found a convergence block; read forward
            conv_block['max_force'] = float(mf.group(1))
            info['forces_converged'] = (mf.group(2) == 'YES')
            for j in range(i+1, min(i+5, len(lines))):
                rf = re.search(r'RMS\s+Force\s+(\S+)\s+\S+\s+(YES|NO)', lines[j])
                md = re.search(r'Maximum Displacement\s+(\S+)\s+\S+\s+(YES|NO)', lines[j])
                rd = re.search(r'RMS\s+Displacement\s+(\S+)\s+\S+\s+(YES|NO)', lines[j])
                if rf:
                    conv_block['rms_force'] = float(rf.group(1))
                if md:
                    conv_block['max_disp'] = float(md.group(1))
                    info['disps_converged'] = (md.group(2) == 'YES')
                if rd:
                    conv_block['rms_disp'] = float(rd.group(1))
            break

    info['max_force'] = conv_block['max_force'] or info['max_force']
    info['rms_force'] = conv_block['rms_force']

    # Check if TDDFT is in progress
    for i in range(max(0, len(lines)-50), len(lines)):
        if 'Excitation energies and oscillator strengths' in lines[i]:
            info['td_in_progress'] = True
            # Check if TDDFT output follows (excited states printed)
            for j in range(i, min(i+200, len(lines))):
                m = re.search(r'Excited State\s+\d+:', lines[j])
                if m:
                    info['td_in_progress'] = True
                    break
            break

    # Try to detect NStates from TD(NStates=...) in route
    for line in lines[:min(200, len(lines))]:
        m = re.search(r'[Tt][Dd]\s*\(\s*[Nn][Ss][Tt][Aa][Tt][Ee][Ss]\s*=\s*(\d+)', line)
        if m:
            info['td_nstates'] = int(m.group(1))
            break
    # Also check Link1 routes
    for i in link1_boundaries:
        for j in range(i, min(i+20, len(lines))):
            m = re.search(r'[Tt][Dd]\s*\(\s*[Nn][Ss][Tt][Aa][Tt][Ee][Ss]\s*=\s*(\d+)', lines[j])
            if m:
                info['td_nstates'] = int(m.group(1))
                break

    # --- Populate parsed values ---
    if last_scf:
        info['scf_energy'] = last_scf[0]
        info['scf_cycle'] = last_scf[1]

    if last_step_match:
        info['step_current'] = last_step_match[0]
        info['step_max'] = last_step_match[1]

    # SCF convergence: last SCF followed by "Convergence on wavefunction"
    for i in range(len(lines) - 1, -1, -1):
        if 'Convergence on wavefunction' in lines[i]:
            info['scf_converged'] = True
            break

    # Elapsed wall time (last occurrence)
    for line in reversed(lines):
        m = re.search(r'Elapsed time:\s+(.*)', line)
        if m:
            info['elapsed_wall'] = m.group(1).strip()
            break

    # --- Determine phase ---
    if info['job_finished']:
        pass  # keep 'complete'
    elif info['opt_converged'] and info['freq_started']:
        info['phase'] = 'frequency'
        info['phase_detail'] = 'Frequency calculation'
    elif info['opt_converged']:
        info['phase'] = 'opt_converged'
        info['phase_detail'] = 'Geometry optimization converged'
    elif last_step_match:
        info['phase'] = 'optimization'
        step_info = f'step {info["step_current"]}'
        if info['scf_converged'] and info['td_in_progress']:
            info['phase_detail'] = f'TDDFT ({step_info})'
        elif info['scf_converged']:
            info['phase_detail'] = f'Gradient calc ({step_info})'
        else:
            info['phase_detail'] = f'SCF ({step_info}, cycle {info["scf_cycle"]})'
    elif last_berny is not None:
        info['phase'] = 'optimization'
        info['phase_detail'] = 'Optimization initializing...'
    elif last_scf:
        info['phase'] = 'initial_scf'
        info['phase_detail'] = f'Initial SCF (cycle {info["scf_cycle"]})'
    elif any('SCF' in l for l in lines):
        info['phase'] = 'initial_scf'
        info['phase_detail'] = 'Initial SCF starting...'
    else:
        info['phase'] = 'startup'
        info['phase_detail'] = 'Calculation starting...'

    return info


# ---------------------------------------------------------------------------
# Progress percentage estimation
# ---------------------------------------------------------------------------

def estimate_percent(info, expected_steps=8,
                     weight_initial_scf=15.0,
                     weight_opt=70.0,
                     weight_freq=15.0):
    """Estimate completion percentage from parsed info."""

    if info['job_finished']:
        return 100.0

    w1, w2, w3 = weight_initial_scf, weight_opt, weight_freq
    total_w = w1 + w2 + w3

    if info['phase'] in ('complete',):
        return 100.0

    if info['phase'] == 'frequency':
        # In freq step: estimate based on SCF within freq
        # Freq SCF is typically faster than opt SCF
        if info['scf_cycle'] > 0:
            freq_scf_pct = min(info['scf_cycle'] / 12.0, 1.0)
        else:
            freq_scf_pct = 0.1
        return w1 + w2 + freq_scf_pct * w3

    if info['phase'] in ('optimization', 'opt_converged'):
        if info['opt_converged'] and not info['freq_started']:
            return w1 + w2

        # Estimate based on steps
        current = info['step_current']
        max_steps = info['step_max'] if info['step_max'] > 0 else expected_steps

        # Use expected_steps as the expected convergence point
        # (max_steps is the Gaussian hardcoded limit, usually 106 or 151)
        effective_max = expected_steps

        # Sub-step progress: based on SCF vs TDDFT vs gradient phases
        sub_pct = 0.0
        if info['scf_converged'] and info['td_in_progress']:
            sub_pct = 0.7  # 70% through this step (SCF done, doing TDDFT)
        elif info['scf_converged']:
            sub_pct = 0.85  # 85% through (SCF+TD done, doing gradient)
        elif info['scf_cycle'] > 0:
            # SCF part: assume 15 SCF cycles per step
            sub_pct = 0.5 * min(info['scf_cycle'] / 15.0, 1.0)
        else:
            sub_pct = 0.05

        # Blend: completed steps + fraction of current step
        if current > 0:
            step_progress = ((current - 1) + sub_pct) / effective_max
        else:
            step_progress = 0.0

        # Cap at 1.0 to prevent >100% when current > expected_steps
        step_progress = min(step_progress, 1.0)

        return w1 + step_progress * w2

    if info['phase'] == 'initial_scf':
        if info['scf_cycle'] > 0:
            scf_frac = min(info['scf_cycle'] / 20.0, 1.0)
        else:
            scf_frac = 0.05
        return scf_frac * w1

    if info['phase'] == 'startup':
        return 1.0

    return 0.0


# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------

def format_progress(info, expected_steps=8,
                    weight_initial_scf=15.0, weight_opt=70.0, weight_freq=15.0):
    """Format a one-line progress summary."""
    pct = estimate_percent(info, expected_steps,
                           weight_initial_scf, weight_opt, weight_freq)

    bar_width = 30

    phase_icon = {
        'startup': '⚙',
        'initial_scf': '🔬',
        'optimization': '📐',
        'opt_converged': '✅',
        'frequency': '🎵',
        'complete': '🏁',
    }.get(info['phase'], '❓')

    elapsed = info.get('elapsed_wall', '?') or '?'

    # Build extra detail
    # --- Build phase pipeline ---
    scf_status = '✓' if (info['phase'] in ('optimization','opt_converged','frequency','complete')
                          or (info['phase'] == 'initial_scf' and info['scf_converged'])) \
                 else ('▸' if info['phase'] == 'initial_scf' else '·')
    opt_status = '✓' if info['phase'] in ('opt_converged','frequency','complete') \
                 else ('▸' if info['phase'] == 'optimization' else '·')
    freq_status = '✓' if info['phase'] == 'complete' \
                  else ('▸' if info['phase'] == 'frequency' else '·')

    pipeline = f"[SCF{scf_status}]→[Opt{opt_status}]→[Freq{freq_status}]"

    # --- Build detail line ---
    detail_parts = []
    if info['step_current'] > 0:
        detail_parts.append(f"step {info['step_current']}")
    if info['scf_cycle'] > 0 and info['phase'] == 'initial_scf':
        detail_parts.append(f"c{info['scf_cycle']}")
    elif info['scf_cycle'] > 0:
        detail_parts.append(f"SCF c{info['scf_cycle']}")
    if info['max_force'] is not None:
        detail_parts.append(f"MaxF={info['max_force']:.6f}")
    if info['forces_converged']:
        detail_parts.append("F✓")
    elif info['max_force'] is not None:
        detail_parts.append("F✗")
    if info['td_in_progress']:
        detail_parts.append("TDDFT")
    if info.get('td_nstates'):
        detail_parts.append(f"NSt={info['td_nstates']}")

    detail = ' | '.join(detail_parts) if detail_parts else ''

    # Cap bar display at 100%
    pct_display = min(pct, 100.0)
    filled = min(int(pct_display / 100.0 * bar_width), bar_width)
    bar = '█' * filled + '░' * (bar_width - filled)

    return (f"{phase_icon} [{bar}] {pct_display:5.1f}%  "
            f"{pipeline:28s}  "
            f"{info['phase_detail']:42s}  "
            f"⌛ {elapsed}  "
            f"{'(' + detail + ')' if detail else ''}")


def format_detail_block(info):
    """Return multi-line detailed status."""
    lines = []
    lines.append(f"File:       {info['filename']}")
    lines.append(f"Phase:      {info['phase']} ({info['phase_detail']})")
    lines.append(f"Lines:      {info['total_lines']}")
    lines.append(f"Elapsed:    {info['elapsed_wall'] or 'unknown'}")
    if info['scf_energy']:
        lines.append(f"SCF Energy: {info['scf_energy']:.8f}")
    lines.append(f"SCF cycles: {info['scf_cycle']} (converged: {info['scf_converged']})")
    if info['step_current'] > 0:
        lines.append(f"Opt step:   {info['step_current']} / {info['step_max']}")
    if info['max_force'] is not None:
        lines.append(f"Max Force:  {info['max_force']:.6f}  (conv: {info['forces_converged']})")
    if info['rms_force'] is not None:
        lines.append(f"RMS Force:  {info['rms_force']:.6f}")
    lines.append(f"TDDFT:      {'active' if info['td_in_progress'] else 'idle'}" +
                 (f" (NStates={info['td_nstates']})" if info['td_nstates'] else ''))
    lines.append(f"Opt conv:   {info['opt_converged']}")
    lines.append(f"Freq start: {info['freq_started']}")
    lines.append(f"Finished:   {info['job_finished']}")
    return '\n'.join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Estimate progress of a running Gaussian calculation.",
        epilog="Examples:\n"
               "  %(prog)s job.log              # one-shot\n"
               "  %(prog)s job.log -w           # watch mode (2 s refresh)\n"
               "  %(prog)s job.log -w -n 5      # 5 s refresh\n"
               "  %(prog)s job.log -e 10        # expect 10 opt steps",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("logfile", help="Gaussian .log file to monitor")
    parser.add_argument("-w", "--watch", action="store_true",
                        help="Real-time monitor mode (refresh every N seconds)")
    parser.add_argument("-n", "--interval", type=float, default=2.0,
                        help="Refresh interval in seconds (default: 2)")
    parser.add_argument("-e", "--expected-steps", type=int, default=8,
                        help="Expected number of optimization steps (default: 8)")
    parser.add_argument("--detail", action="store_true",
                        help="Show detailed status block instead of one-line summary")
    parser.add_argument("--weights", default="15,70,15",
                        help="Phase weight percentages: initial_scf,opt,freq (default: 15,70,15)")
    args = parser.parse_args()

    # Parse weights
    try:
        w_parts = [float(x) for x in args.weights.split(',')]
        if len(w_parts) != 3:
            raise ValueError
        w_scf, w_opt, w_freq = w_parts
    except (ValueError, TypeError):
        print("ERROR: --weights must be three comma-separated numbers (e.g. 15,70,15)",
              file=sys.stderr)
        sys.exit(1)

    def show_progress():
        info = parse_progress(args.logfile)
        if info is None:
            print(f"ERROR: cannot read file: {args.logfile}", file=sys.stderr)
            return False
        if args.detail:
            print(format_detail_block(info))
            pct = estimate_percent(info, args.expected_steps,
                                   w_scf, w_opt, w_freq)
            print(f"\nEstimated progress: {pct:.1f}%")
        else:
            print(format_progress(info, args.expected_steps,
                                  w_scf, w_opt, w_freq))
        return True

    if not args.watch:
        ok = show_progress()
        if not ok:
            sys.exit(1)
    else:
        # Watch mode: refresh display in-place
        last_mtime = 0
        try:
            while True:
                # Only update if file changed (or first iteration)
                try:
                    current_mtime = os.path.getmtime(args.logfile)
                except OSError:
                    current_mtime = 0

                if current_mtime != last_mtime:
                    # Clear screen and move cursor home
                    sys.stdout.write('\033[2J\033[H')
                    show_progress()
                    sys.stdout.write('\n(Press Ctrl+C to stop. Refresh: '
                                     f'{args.interval}s)\n')
                    sys.stdout.flush()
                    last_mtime = current_mtime

                time.sleep(args.interval)
        except KeyboardInterrupt:
            print("\nStopped.")


if __name__ == '__main__':
    main()
