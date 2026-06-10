# Visualization

Scripts for molecular visualization and quantum chemistry data exploration.

## Orbital Visualizer

Interactive 3D molecular orbital viewer — equivalent to VMD's Orbital panel. Browse canonical and localized orbitals from GAMESS calculations.

### Quick start

```bash
# One-time setup
bash setup.sh

# Launch GUI
./orbital-visualizer calculation.log

# CLI mode (headless, saves PNG)
python orbital-visualizer.py --cli --orbital 40 --output homo.png calculation.log
```

### Usage

| Action | How |
|--------|-----|
| Open file | File → Open (`Ctrl+O`) or CLI argument |
| Browse orbitals | Click in gallery (Occupied / Virtual / Localized tabs) |
| Change isovalue | Slider at bottom (0.01–0.20) |
| Change grid resolution | Grid slider (0.05–0.40 Å, coarser = faster) |
| Jump to HOMO/LUMO | HOMO / LUMO buttons |
| Export image | File → Export Image (`Ctrl+E`) |
| CLI render | `--cli --orbital N --output file.png` |

### Progressive refinement

Coarse grid renders first (instant) → auto-refines to medium → fine. Clicking another orbital cancels pending refinement.

### CPU usage

Defaults to `(cores − 2)` threads. Override:
```bash
NUMBA_NUM_THREADS=4 ./orbital-visualizer calculation.log
```

### Supported formats

| Format | Status |
|--------|--------|
| GAMESS `.log` | ✅ Canonical + Boys/Pipek-Mezey/Edmiston-Ruedenberg localized |
| Other (ORCA, Gaussian, etc.) | 🔜 Adapter pattern ready, not yet implemented |

### Requirements

- Python 3.10+
- Dependencies (installed by `setup.sh`): `cclib`, `numba`, `vispy`, `scikit-image`, `PyQt6`

## Files

| File | Purpose |
|------|---------|
| `orbital-visualizer.py` | Main application |
| `orbital-visualizer` | Bash launcher (portable across clone locations) |
| `setup.sh` | One-time venv + dependency install |
| `CONTEXT.md` | Design decisions and glossary (not distributed) |
