#!/usr/bin/env bash
# Set up Python virtual environment for Orbital Visualizer.
# Run once after cloning the repo.
# Usage: bash setup.sh

set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# --- System dependencies check ---
echo "Checking system OpenGL libraries ..."
MISSING=""

# Try to load OpenGL
if ! python3 -c "import ctypes; ctypes.CDLL('libGL.so.1')" 2>/dev/null && \
   ! python3 -c "import ctypes; ctypes.CDLL('libOpenGL.so.0')" 2>/dev/null; then
    MISSING="OpenGL"
fi

# Check GL dispatch library (libglvnd)
if ! python3 -c "import ctypes; ctypes.CDLL('libEGL.so.1')" 2>/dev/null; then
    MISSING="${MISSING:+$MISSING, }EGL"
fi

if [ -n "$MISSING" ]; then
    echo ""
    echo "WARNING: Missing system libraries: $MISSING"
    echo "The visualizer needs GPU drivers and OpenGL libraries."
    echo ""
    echo "Install on Ubuntu/Debian:"
    echo "  sudo apt install libgl1-mesa-glx libegl1-mesa libglu1-mesa mesa-utils"
    echo ""
    echo "Install on Fedora/RHEL:"
    echo "  sudo dnf install mesa-libGL mesa-libEGL mesa-libGLU glx-utils"
    echo ""
    echo "Install on Arch:"
    echo "  sudo pacman -S mesa libglvnd glu"
    echo ""
    echo "Install on openSUSE:"
    echo "  sudo zypper install Mesa-libGL1 Mesa-libEGL1 Mesa-libGLU1"
    echo ""
    echo "Continuing with Python setup anyway (will fail at runtime if missing)."
    echo ""
fi

echo "Creating virtual environment in $SCRIPT_DIR/venv ..."
python3 -m venv venv

echo "Installing dependencies ..."
source venv/bin/activate
pip install cclib numba vispy scikit-image PyQt6

echo ""
echo "Done. Run the visualizer with:"
echo "  orbital-viewer [file.log]"
echo "  or: python orbital-visualizer.py [file.log]"
