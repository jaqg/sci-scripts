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
| Split viewport | Click **Split ▸** to compare two orbitals side by side |
| Multiple molecules | File → Open or Ctrl+T to open additional molecules in new tabs |
| Switch active viewport | Click viewport label (top bar) — blue border = active, gallery click goes there |
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
- Python packages (installed by `setup.sh`): `cclib`, `numba`, `vispy`, `scikit-image`, `PyQt6`
- **System OpenGL drivers** — Mesa or proprietary GPU drivers

| Distro | Install command |
|--------|----------------|
| Ubuntu / Debian | `sudo apt install libgl1-mesa-glx libegl1-mesa mesa-utils` |
| Fedora / RHEL | `sudo dnf install mesa-libGL mesa-libEGL glx-utils` |
| Arch | `sudo pacman -S mesa libglvnd` |
| openSUSE | `sudo zypper install Mesa-libGL1 Mesa-libEGL1` |

### Platform notes

- **X11**: Fully supported (GLX)
- **Wayland (GNOME, KDE)**: Supported via Qt6's EGL backend. Ensure `libegl1-mesa` is installed.
  The PyQt6 wheel ships its own `libqwayland.so` — no extra Qt packages needed.
- **WSL2**: Works with `wslg` (built-in X11/Wayland). Install Mesa if missing.
- **macOS**: Untested. PyQt6 + vispy work on macOS but scene camera may need tuning.
- **Windows**: Untested. Should work with PyQt6 wheel (ships ANGLE/OpenGL ES).
- **Headless / SSH**: Use `--cli` mode (no GUI, renders to PNG).
- **Troubleshooting**: If you get a black window or OpenGL errors, try:
  ```bash
  QT_QPA_PLATFORM=xcb ./orbital-visualizer file.log   # force X11 on Wayland
  NUMBA_NUM_THREADS=1 ./orbital-visualizer file.log   # limit to 1 thread
  ```

## Files

| File | Purpose |
|------|---------|
| `orbital-visualizer.py` | Main application |
| `orbital-visualizer` | Bash launcher (portable across clone locations) |
| `setup.sh` | One-time venv + dependency install |
| `CONTEXT.md` | Design decisions and glossary (not distributed) |
