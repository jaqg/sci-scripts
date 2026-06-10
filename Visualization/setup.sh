#!/usr/bin/env bash
# Set up Python virtual environment for Orbital Visualizer.
# Run once after cloning the repo.
# Usage: bash setup.sh

set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "Creating virtual environment in $SCRIPT_DIR/venv ..."
python3 -m venv venv

echo "Installing dependencies ..."
source venv/bin/activate
pip install cclib numba vispy scikit-image PyQt6

echo ""
echo "Done. Run the visualizer with:"
echo "  orbital-viewer [file.log]"
echo "  or: python orbital-visualizer.py [file.log]"
