#!/usr/bin/env bash
# Set up Python virtual environment for Orbital Visualizer.
# Run once after cloning the repo.
# Usage: bash setup.sh

set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# --- System dependencies check ---
echo "Checking system libraries ..."
MISSING=""
MISSING_PKG=""

# Try to load OpenGL
if ! python3 -c "import ctypes; ctypes.CDLL('libGL.so.1')" 2>/dev/null && \
   ! python3 -c "import ctypes; ctypes.CDLL('libOpenGL.so.0')" 2>/dev/null; then
    MISSING="${MISSING:+$MISSING, }OpenGL"
    MISSING_PKG="${MISSING_PKG} libgl1-mesa-glx"
fi

# Check GL dispatch library (libglvnd)
if ! python3 -c "import ctypes; ctypes.CDLL('libEGL.so.1')" 2>/dev/null; then
    MISSING="${MISSING:+$MISSING, }EGL"
    MISSING_PKG="${MISSING_PKG} libegl1-mesa"
fi

# Check XCB cursor (needed by Qt 6.5+ on X11)
if ! ldconfig -p 2>/dev/null | grep -q libxcb-cursor.so; then
    if ! python3 -c "import ctypes; ctypes.CDLL('libxcb-cursor.so.0')" 2>/dev/null; then
        MISSING="${MISSING:+$MISSING, }xcb-cursor"
        MISSING_PKG="${MISSING_PKG} libxcb-cursor0"
    fi
fi

if [ -n "$MISSING" ]; then
    echo ""
    echo "WARNING: Missing system libraries: $MISSING"
    echo "The visualizer needs GPU drivers, OpenGL, and XCB cursor libraries."
    echo ""
    echo "Install on Ubuntu/Debian:"
    echo "  sudo apt install${MISSING_PKG}"
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
