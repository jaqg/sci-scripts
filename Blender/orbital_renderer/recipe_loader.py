# -*- coding: utf-8 -*-
#
#  Blender Orbital Renderer — recipe_loader.py
#  Reads render_recipe.json and Gaussian cube files.
#
#  Cube file reader adapted from Molecular Blender:
#    (c) Shane Parker, Joshua Szekely — GPLv3
#
#  License: GPLv3
#

"""Load render recipes and volumetric cube files."""

import json
import os
import numpy as np


# ---------------------------------------------------------------------------
# Unit conversion
# ---------------------------------------------------------------------------

BOHR_TO_ANG = 1.8897259885789 ** -1  # 1 Bohr = 0.529177... Å


# ---------------------------------------------------------------------------
# Cube file reader
# ---------------------------------------------------------------------------

def load_cube_file(cube_path: str):
    """
    Read a Gaussian cube file and return grid data + metadata.

    Parameters
    ----------
    cube_path : str
        Absolute path to .cube file.

    Returns
    -------
    dict with keys:
        grid : np.ndarray, shape (nx, ny, nz)
            Volumetric data (float32).
        atoms : list of dict
            Each atom: {symbol, atomic_number, position: [x, y, z] (Å)}.
        origin : np.ndarray, shape (3,)
            Grid origin in Å.
        axes : np.ndarray, shape (3, 3)
            Grid axis vectors (step between adjacent points) in Å.
        npoints : list[int], length 3
            Number of grid points along each axis [nx, ny, nz].

    Raises
    ------
    FileNotFoundError
        If cube_path does not exist.
    ValueError
        If the cube file has sheared/oblique axes.
    """
    if not os.path.isfile(cube_path):
        raise FileNotFoundError(f"Cube file not found: {cube_path}")

    atoms = []
    npoints = [0, 0, 0]
    axes = np.zeros((3, 3), dtype=np.float32)

    with open(cube_path, "r") as f:
        # Lines 1-2: comments
        f.readline()
        f.readline()

        # Line 3: natoms, origin_x, origin_y, origin_z (Bohr)
        parts = f.readline().split()
        natoms = int(parts[0])
        has_dset_ids = natoms < 0
        natoms = abs(natoms)
        origin = np.array([float(x) for x in parts[1:4]], dtype=np.float32) * BOHR_TO_ANG

        # Lines 4-6: grid axes
        for i in range(3):
            parts = f.readline().split()
            npoints[i] = int(parts[0])
            axis_vec = np.array([float(x) for x in parts[1:4]], dtype=np.float32)
            if npoints[i] > 0:
                # Positive npoints → Bohr units
                axis_vec *= BOHR_TO_ANG
            else:
                # Negative npoints → Ångström units (use absolute value)
                npoints[i] = abs(npoints[i])
            axes[i, :] = axis_vec

        # Check for oblique axes (not supported)
        overlap = np.dot(axes, axes.T)
        for (a, b) in [(0, 1), (0, 2), (1, 2)]:
            if abs(overlap[a, b]) > 1e-6:
                raise ValueError(f"Cube file has sheared/oblique axes: {cube_path}")

        # Read atoms: atomic_number, charge, x, y, z (Bohr)
        atomic_numbers = []
        for _ in range(natoms):
            parts = f.readline().split()
            anum = int(float(parts[0]))
            pos = np.array([float(x) for x in parts[2:5]], dtype=np.float32) * BOHR_TO_ANG
            symbol = _atomic_number_to_symbol(anum)
            atoms.append({
                "symbol": symbol,
                "atomic_number": anum,
                "position": pos.tolist(),
            })
            atomic_numbers.append(anum)

        # Optional dataset IDs
        if has_dset_ids:
            f.readline()  # Skip dataset ID line

        # Volumetric data: z-y-x inner-loop order, 6 floats per line
        total_points = npoints[0] * npoints[1] * npoints[2]
        data = np.zeros(total_points, dtype=np.float32)
        idx = 0
        while idx < total_points:
            line = f.readline()
            values = [float(x) for x in line.split()]
            data[idx: idx + len(values)] = values
            idx += len(values)

    # Reshape to (nx, ny, nz) — following Molecular Blender's convention.
    # Cube files store data in z-y-x order; C-order reshape is the
    # convention used by Molecular Blender and produces correct results.
    grid = data.reshape(npoints).copy()

    return {
        "grid": grid,
        "atoms": atoms,
        "origin": origin,
        "axes": axes,
        "npoints": npoints,
    }


# ---------------------------------------------------------------------------
# Atomic number ↔ symbol
# ---------------------------------------------------------------------------

# Lookup: atomic number → symbol (subset covering common elements)
_ATOMIC_NUMBER_TO_SYMBOL = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F",
    10: "Ne", 11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl",
    18: "Ar", 19: "K", 20: "Ca", 21: "Sc", 22: "Ti", 23: "V", 24: "Cr", 25: "Mn",
    26: "Fe", 27: "Co", 28: "Ni", 29: "Cu", 30: "Zn", 31: "Ga", 32: "Ge", 33: "As",
    34: "Se", 35: "Br", 36: "Kr", 37: "Rb", 38: "Sr", 39: "Y", 40: "Zr", 41: "Nb",
    42: "Mo", 43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In",
    50: "Sn", 51: "Sb", 52: "Te", 53: "I", 54: "Xe", 55: "Cs", 56: "Ba",
    57: "La", 58: "Ce", 59: "Pr", 60: "Nd", 61: "Pm", 62: "Sm", 63: "Eu",
    64: "Gd", 65: "Tb", 66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb",
    71: "Lu", 72: "Hf", 73: "Ta", 74: "W", 75: "Re", 76: "Os", 77: "Ir",
    78: "Pt", 79: "Au", 80: "Hg", 81: "Tl", 82: "Pb", 83: "Bi", 84: "Po",
    85: "At", 86: "Rn", 87: "Fr", 88: "Ra", 89: "Ac", 90: "Th", 91: "Pa",
    92: "U", 93: "Np", 94: "Pu", 95: "Am", 96: "Cm", 97: "Bk", 98: "Cf",
    99: "Es", 100: "Fm", 101: "Md", 102: "No", 103: "Lr", 104: "Rf",
    105: "Db", 106: "Sg", 107: "Bh", 108: "Hs", 109: "Mt",
}


def _atomic_number_to_symbol(anum: int) -> str:
    """Convert atomic number to element symbol."""
    return _ATOMIC_NUMBER_TO_SYMBOL.get(anum, "Xx")


# ---------------------------------------------------------------------------
# Recipe loader
# ---------------------------------------------------------------------------

def load_recipe(recipe_path: str):
    """
    Load a render recipe JSON file.

    Parameters
    ----------
    recipe_path : str
        Path to render_recipe.json.

    Returns
    -------
    dict
        Parsed recipe with all paths resolved to absolute.

    Raises
    ------
    FileNotFoundError
        If recipe_path does not exist.
    ValueError
        If recipe JSON is malformed.
    """
    if not os.path.isfile(recipe_path):
        raise FileNotFoundError(f"Recipe file not found: {recipe_path}")

    recipe_dir = os.path.dirname(os.path.abspath(recipe_path))

    with open(recipe_path, "r") as f:
        recipe = json.load(f)

    # Validate version
    version = recipe.get("version", 1)
    if version != 1:
        raise ValueError(f"Unsupported recipe version: {version}")

    # Resolve cube file paths to absolute
    for orbital in recipe.get("orbitals", []):
        cube_rel = orbital.get("cube_file", "")
        cube_abs = os.path.normpath(os.path.join(recipe_dir, cube_rel))
        orbital["_cube_path"] = cube_abs

        # Apply defaults for missing fields
        orbital.setdefault("isovalue", 0.05)
        orbital.setdefault("grid_spacing", 0.08)
        orbital.setdefault("wtype", "canonical")
        orbital.setdefault("label", orbital.get("wtype", "canonical"))

    # Resolve source_log path
    if "source_log" in recipe:
        recipe["_source_log_abs"] = os.path.normpath(
            os.path.join(recipe_dir, recipe["source_log"])
        )

    return recipe
