#!/usr/bin/env python3
"""
Orbital Visualizer — interactive 3D molecular orbital viewer.
Equivalent to VMD's Graphical Representations → Orbital panel.
Usage: python orbital-visualizer.py [calculation.log]
"""

import sys
import re
import math
import argparse
from collections import namedtuple
from pathlib import Path

import os
import numpy as np

# Cap CPU threads: leave 2 cores free by default. Set NUMBA_NUM_THREADS to override.
if 'NUMBA_NUM_THREADS' not in os.environ:
    cpu_count = os.cpu_count() or 4
    os.environ['NUMBA_NUM_THREADS'] = str(max(1, cpu_count - 2))

from numba import njit, prange

# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

Atom = namedtuple('Atom', ['idx', 'symbol', 'atomic_number', 'x', 'y', 'z'])
Shell = namedtuple('Shell', ['atom_idx', 'ang_mom', 'primitives'])  # primitives = [(exp, coeff), ...]


class BasisFunction:
    """One Cartesian basis function within a shell."""
    __slots__ = ('atom_idx', 'lx', 'ly', 'lz', 'angular_norm', 'shell_idx')
    
    def __init__(self, atom_idx, lx, ly, lz, shell_idx):
        self.atom_idx = atom_idx
        self.lx = lx
        self.ly = ly
        self.lz = lz
        self.shell_idx = shell_idx
        # Angular normalization: sqrt(lx! ly! lz! / ((2lx)! (2ly)! (2lz)!))
        self.angular_norm = self._angular_norm()
    
    def _angular_norm(self):
        lx, ly, lz = self.lx, self.ly, self.lz
        num = math.factorial(lx) * math.factorial(ly) * math.factorial(lz)
        den = math.factorial(2*lx) * math.factorial(2*ly) * math.factorial(2*lz)
        return math.sqrt(num / den)


class BasisSet:
    """Collection of shells and basis functions."""
    
    # Cartesian expansion per angular momentum
    CARTESIAN_MAP = {
        'S': [(0, 0, 0)],                                                    # 1
        'P': [(1, 0, 0), (0, 1, 0), (0, 0, 1)],                             # 3
        'D': [(2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1)],  # 6
        'F': [(3, 0, 0), (0, 3, 0), (0, 0, 3),                              # 10
              (2, 1, 0), (2, 0, 1), (1, 2, 0), (0, 2, 1), (1, 0, 2), (0, 1, 2),
              (1, 1, 1)],
        'G': [(4, 0, 0), (0, 4, 0), (0, 0, 4),                              # 15
              (3, 1, 0), (3, 0, 1), (1, 3, 0), (0, 3, 1), (1, 0, 3), (0, 1, 3),
              (2, 2, 0), (2, 0, 2), (0, 2, 2),
              (2, 1, 1), (1, 2, 1), (1, 1, 2)],
    }
    
    def __init__(self, shells, atoms):
        self.shells = shells          # list of Shell
        self.atoms = atoms            # list of Atom
        self.nbasis = 0
        self.basis_functions = []     # list of BasisFunction
        self._expand()
    
    def _expand(self):
        """Expand shells into Cartesian basis functions."""
        self.basis_functions = []
        for sidx, shell in enumerate(self.shells):
            lmn_list = self.CARTESIAN_MAP.get(shell.ang_mom, [])
            if not lmn_list:
                raise ValueError(f"Unsupported angular momentum: {shell.ang_mom}")
            for lx, ly, lz in lmn_list:
                bf = BasisFunction(shell.atom_idx, lx, ly, lz, sidx)
                self.basis_functions.append(bf)
        self.nbasis = len(self.basis_functions)
    
    def get_lmn_array(self):
        """Return (nbasis, 3) array of (lx, ly, lz)."""
        arr = np.zeros((self.nbasis, 3), dtype=np.int32)
        for i, bf in enumerate(self.basis_functions):
            arr[i, 0] = bf.lx
            arr[i, 1] = bf.ly
            arr[i, 2] = bf.lz
        return arr
    
    def get_atom_idx_array(self):
        """Return (nbasis,) array of atom indices."""
        return np.array([bf.atom_idx for bf in self.basis_functions], dtype=np.int32)
    
    def get_shell_idx_array(self):
        """Return (nbasis,) array of shell indices."""
        return np.array([bf.shell_idx for bf in self.basis_functions], dtype=np.int32)
    
    def get_angular_norm_array(self):
        """Return (nbasis,) array of angular normalization factors."""
        return np.array([bf.angular_norm for bf in self.basis_functions], dtype=np.float64)
    
    def flatten_primitives(self):
        """Return flat arrays for primitive data.

        Returns:
            shell_prim_start: (nshells+1,) int32
            shell_L: (nshells,) int32
            shell_cutoff_r2: (nshells,) float64 — max cutoff distance² across primitives in shell
            prim_exp: (nprims,) float64
            prim_coeff: (nprims,) float64
            prim_radial_norm: (nprims,) float64
        """
        CUTOFF_THRESHOLD = 1e-12
        nprims_total = sum(len(shell.primitives) for shell in self.shells)
        nshells = len(self.shells)
        shell_prim_start = np.zeros(nshells + 1, dtype=np.int32)
        shell_L = np.zeros(nshells, dtype=np.int32)
        shell_cutoff_r2 = np.zeros(nshells, dtype=np.float64)
        prim_exp = np.zeros(nprims_total, dtype=np.float64)
        prim_coeff = np.zeros(nprims_total, dtype=np.float64)
        prim_radial_norm = np.zeros(nprims_total, dtype=np.float64)

        idx = 0
        for sidx, shell in enumerate(self.shells):
            shell_prim_start[sidx] = idx
            lmn_list = self.CARTESIAN_MAP.get(shell.ang_mom, [(0, 0, 0)])
            L = sum(lmn_list[0])
            shell_L[sidx] = L

            # Spatial cutoff: find most diffuse primitive (smallest alpha)
            alpha_min = min(alpha for alpha, _ in shell.primitives)
            cutoff_r = math.sqrt(-math.log(CUTOFF_THRESHOLD) / alpha_min)
            shell_cutoff_r2[sidx] = cutoff_r * cutoff_r

            for alpha, coeff in shell.primitives:
                prim_exp[idx] = alpha
                prim_coeff[idx] = coeff
                radial = (2.0 * alpha / math.pi) ** 0.75
                if L > 0:
                    radial *= math.sqrt((8.0 * alpha) ** L)
                prim_radial_norm[idx] = radial
                idx += 1

        shell_prim_start[-1] = idx
        return shell_prim_start, shell_L, shell_cutoff_r2, prim_exp, prim_coeff, prim_radial_norm


class Wavefunction:
    """Holds MO coefficients and metadata."""
    
    def __init__(self, coefficients, energies, labels=None):
        """
        coefficients: (nmo, nbasis) array — each row is an MO
        energies: (nmo,) array — orbital energies in Hartree
        labels: list of str or None — 'canonical', 'boys', 'pipek-mezey', etc.
        """
        self.coefficients = coefficients  # (nmo, nbasis)
        self.energies = energies          # (nmo,)
        self.nmo = coefficients.shape[0]
        self.nbasis = coefficients.shape[1]
        if labels is None:
            labels = ['orbital'] * self.nmo
        self.labels = labels
    
    def get_mo(self, idx):
        """Return coefficients for MO at index idx."""
        return self.coefficients[idx]


# ---------------------------------------------------------------------------
# GAMESS log parsing
# ---------------------------------------------------------------------------

# Cartesian angular momentum labels as printed by GAMESS
CART_LABELS_GAMESS = ['S', 'X', 'Y', 'Z',
                       'XX', 'YY', 'ZZ', 'XY', 'XZ', 'YZ',
                       'XXX', 'YYY', 'ZZZ', 'XXY', 'XXZ', 'YYX', 'YYZ', 'ZZX', 'ZZY', 'XYZ',
                       'XXXX', 'YYYY', 'ZZZZ', 'XXXY', 'XXXZ', 'XYYY', 'YYYZ', 'XZZZ', 'YZZZ',
                       'XXYY', 'XXZZ', 'YYZZ', 'XXYZ', 'YYXZ', 'ZZXY']

CART_LABELS_TO_LMN = {
    'S': (0, 0, 0),
    'X': (1, 0, 0), 'Y': (0, 1, 0), 'Z': (0, 0, 1),
    'XX': (2, 0, 0), 'YY': (0, 2, 0), 'ZZ': (0, 0, 2),
    'XY': (1, 1, 0), 'XZ': (1, 0, 1), 'YZ': (0, 1, 1),
    'XXX': (3, 0, 0), 'YYY': (0, 3, 0), 'ZZZ': (0, 0, 3),
    'XXY': (2, 1, 0), 'XXZ': (2, 0, 1), 'YYX': (1, 2, 0),
    'YYZ': (0, 2, 1), 'ZZX': (1, 0, 2), 'ZZY': (0, 1, 2),
    'XYZ': (1, 1, 1),
    'XXXX': (4, 0, 0), 'YYYY': (0, 4, 0), 'ZZZZ': (0, 0, 4),
    'XXXY': (3, 1, 0), 'XXXZ': (3, 0, 1), 'XYYY': (1, 3, 0),
    'YYYZ': (0, 3, 1), 'XZZZ': (1, 0, 3), 'YZZZ': (0, 1, 3),
    'XXYY': (2, 2, 0), 'XXZZ': (2, 0, 2), 'YYZZ': (0, 2, 2),
    'XXYZ': (2, 1, 1), 'YYXZ': (1, 2, 1), 'ZZXY': (1, 1, 2),
}


def parse_gamess_log(filepath):
    """Parse a GAMESS .log file and return Molecule, BasisSet, canonical and localized wavefunctions.
    
    Returns:
        atoms: list of Atom
        basis_set: BasisSet
        canon_wfn: Wavefunction (canonical MOs)
        local_wfn: Wavefunction or None (localized MOs)
        homo_idx: int (0-based index of HOMO)
    """
    import cclib
    from periodictable import elements
    
    data = cclib.io.ccread(str(filepath))
    
    # --- Atoms / Molecule ---
    symbols = [elements[z].symbol for z in data.atomnos]
    coords = data.atomcoords[0]  # (natom, 3) in Angstrom
    atoms = []
    for i in range(data.natom):
        atoms.append(Atom(i, symbols[i], int(data.atomnos[i]),
                          float(coords[i, 0]), float(coords[i, 1]), float(coords[i, 2])))
    
    # --- Basis Set from cclib's gbasis ---
    shells = []
    for atom_idx, atom_shells in enumerate(data.gbasis):
        for ang_mom, primitives in atom_shells:
            # Convert primitives from list of tuples to list of (exp, coeff) floats
            prims = [(float(e), float(c)) for e, c in primitives]
            shells.append(Shell(atom_idx, ang_mom, prims))
    
    basis_set = BasisSet(shells, atoms)
    
    # Verify basis set size
    nbasis_expected = data.nbasis
    if basis_set.nbasis != nbasis_expected:
        raise ValueError(f"Basis function count mismatch: expanded {basis_set.nbasis}, "
                         f"expected {nbasis_expected}")
    
    # --- Canonical MOs from cclib ---
    # mocoeffs is stored as list of arrays; for single-point calcs, mocoeffs[0] is (nmo, nbasis)
    mocoeffs_raw = data.mocoeffs
    if isinstance(mocoeffs_raw, list):
        mocoeffs_arr = np.array(mocoeffs_raw[0], dtype=np.float64)
    else:
        mocoeffs_arr = np.array(mocoeffs_raw, dtype=np.float64)
    
    moenergies_raw = data.moenergies
    if isinstance(moenergies_raw, list):
        moenergies_arr = np.array(moenergies_raw[0], dtype=np.float64)
    else:
        moenergies_arr = np.array(moenergies_raw, dtype=np.float64)
    
    homo_idx = data.homos[0]  # 0-based (cclib convention for single ref)
    canon_labels = ['canonical'] * mocoeffs_arr.shape[0]
    canon_wfn = Wavefunction(mocoeffs_arr, moenergies_arr, canon_labels)
    
    # --- Localized MOs from raw text ---
    local_wfn = _parse_localized_orbitals(filepath, data.nbasis, homo_idx)
    
    return atoms, basis_set, canon_wfn, local_wfn, homo_idx


def _parse_localized_orbitals(filepath, nbasis, homo_idx):
    """Parse localized orbitals from raw GAMESS log text.
    
    Returns Wavefunction or None if no localized orbitals found.
    """
    with open(filepath, 'r') as f:
        text = f.read()
    
    # Try Boys first, then Pipek-Mezey, then Edmiston-Ruedenberg
    markers = [
        'THE BOYS LOCALIZED ORBITALS ARE',
        'THE PIPEK-MEZEY POPULATION LOCALIZED ORBITALS ARE',
        'EDMISTON-RUEDENBERG ENERGY LOCALIZED ORBITALS',
    ]
    
    found_marker = None
    marker_pos = -1
    
    for marker in markers:
        pos = text.find(marker)
        if pos != -1:
            found_marker = marker
            marker_pos = pos
            break
    
    if found_marker is None:
        return None
    
    # Determine label
    if 'BOYS' in found_marker:
        label = 'boys'
    elif 'PIPEK' in found_marker:
        label = 'pipek-mezey'
    elif 'EDMISTON' in found_marker:
        label = 'edmiston-ruedenberg'
    else:
        label = 'localized'
    
    # Determine number of localized orbitals from the localization summary
    # Look for "HAS   XX ORBITALS" before the marker
    local_section = text[:marker_pos + 500]
    nlocal_match = re.search(r'HAS\s+(\d+)\s+ORBITALS', local_section)
    if nlocal_match:
        nlocal = int(nlocal_match.group(1))
    else:
        # Default: number of occupied orbitals
        nlocal = homo_idx + 1
    
    # Parse the coefficient blocks
    # Format: header line with MO numbers "   1          2          3          4          5"
    # Then 520 rows of "  atom_idx  symbol  shell_idx  label  coeff1  coeff2  coeff3  coeff4  coeff5"
    # No energy line
    lines = text[marker_pos:].split('\n')
    
    coefficients = np.zeros((nlocal, nbasis), dtype=np.float64)
    
    # Skip the marker line and find first header
    line_idx = 0
    # Skip marker line
    line_idx += 1
    
    blocks_read = 0
    
    while line_idx < len(lines) and blocks_read * 5 < nlocal:
        line = lines[line_idx].strip()
        
        # Look for header line: just numbers "1          2          3          4          5"
        if re.match(r'^\s*\d+(\s+\d+)+$', line) and not re.search(r'[A-Za-z]', line):
            # This is a header line - parse the MO indices
            mo_nums = [int(x) for x in line.split()]
            ncols = len(mo_nums)
            line_idx += 1
            
            # Now read nbasis lines of coefficients (skip blank lines)
            bf_idx = 0
            while bf_idx < nbasis and line_idx < len(lines):
                data_line = lines[line_idx].strip()
                line_idx += 1
                
                if not data_line:
                    continue
                
                # Parse: "idx  symbol  shell  label  c1  c2  c3  c4  c5"
                parts = data_line.split()
                # First 4 parts are labels, rest are coefficients
                if len(parts) < 4 + ncols:
                    continue
                
                coeff_strs = parts[4:4 + ncols]
                for col, coeff_str in enumerate(coeff_strs):
                    mo_global_idx = mo_nums[col] - 1  # 0-based
                    if mo_global_idx < nlocal:
                        try:
                            coefficients[mo_global_idx, bf_idx] = float(coeff_str)
                        except ValueError:
                            pass
                bf_idx += 1
            
            blocks_read += 1
        else:
            line_idx += 1
    
    labels_list = [label] * nlocal
    # No energies for localized orbitals
    energies_arr = np.zeros(nlocal, dtype=np.float64)
    
    return Wavefunction(coefficients, energies_arr, labels_list)


# ---------------------------------------------------------------------------
# CPK colors and covalent radii
# ---------------------------------------------------------------------------

# CPK atom colors (RGB, 0-1)
CPK_COLORS = {
    1:  (1.0, 1.0, 1.0),    # H - white
    2:  (0.85, 1.0, 1.0),   # He
    3:  (0.8, 0.5, 1.0),    # Li
    4:  (0.76, 1.0, 0.0),   # Be
    5:  (1.0, 0.71, 0.71),  # B
    6:  (0.35, 0.35, 0.35), # C - dark gray
    7:  (0.14, 0.14, 0.82), # N - blue
    8:  (1.0, 0.05, 0.05),  # O - red
    9:  (0.5, 0.7, 0.3),    # F
    10: (0.85, 1.0, 1.0),   # Ne
    11: (0.67, 0.36, 0.95), # Na
    12: (0.54, 1.0, 0.0),   # Mg
    13: (0.75, 0.65, 0.65), # Al
    14: (0.5, 0.6, 0.6),    # Si
    15: (1.0, 0.5, 0.0),    # P
    16: (1.0, 1.0, 0.19),   # S - yellow
    17: (0.12, 0.94, 0.12), # Cl - green
    35: (0.65, 0.16, 0.16), # Br
    53: (0.58, 0.0, 0.58),  # I
}

# Covalent radii in Angstrom (for bond detection)
COVALENT_RADII = {
    1: 0.31,  2: 0.28,
    3: 1.28,  4: 0.96,  5: 0.84,  6: 0.76,  7: 0.71,  8: 0.66,  9: 0.57, 10: 0.58,
    11: 1.66, 12: 1.41, 13: 1.21, 14: 1.11, 15: 1.07, 16: 1.05, 17: 1.02, 18: 1.06,
    35: 1.20, 53: 1.39,
}

BOND_CUTOFF_FACTOR = 1.2


# ---------------------------------------------------------------------------
# Grid and basis function evaluation (numba kernel)
# ---------------------------------------------------------------------------

@njit(parallel=True)
def _eval_mo_kernel(grid_points, atom_centers, basis_atom_idx, basis_shell_idx,
                    basis_lx, basis_ly, basis_lz, angular_norms,
                    shell_prim_start, shell_L, shell_cutoff_r2,
                    prim_exp, prim_coeff, prim_radial_norm,
                    mo_coeffs, out_values):
    """Evaluate one MO on all grid points with spatial cutoff."""
    nbasis = basis_atom_idx.shape[0]
    npoints = grid_points.shape[0]

    for p_idx in prange(npoints):
        x = grid_points[p_idx, 0]
        y = grid_points[p_idx, 1]
        z = grid_points[p_idx, 2]
        total = 0.0

        for bf_idx in range(nbasis):
            aidx = basis_atom_idx[bf_idx]
            ax = atom_centers[aidx, 0]
            ay = atom_centers[aidx, 1]
            az = atom_centers[aidx, 2]
            rx = x - ax
            ry = y - ay
            rz = z - az
            r2 = rx * rx + ry * ry + rz * rz

            # Spatial cutoff: skip if grid point is beyond shell's cutoff radius
            sidx = basis_shell_idx[bf_idx]
            if r2 > shell_cutoff_r2[sidx]:
                continue

            p_start = shell_prim_start[sidx]
            p_end = shell_prim_start[sidx + 1]

            bf_val = 0.0
            for p in range(p_start, p_end):
                alpha = prim_exp[p]
                exp_part = math.exp(-alpha * r2)
                if exp_part < 1e-300:
                    continue
                coeff = prim_coeff[p]
                rad_norm = prim_radial_norm[p]

                ang_part = 1.0
                lx = basis_lx[bf_idx]
                ly = basis_ly[bf_idx]
                lz = basis_lz[bf_idx]
                if lx > 0:
                    ang_part *= rx ** lx
                if ly > 0:
                    ang_part *= ry ** ly
                if lz > 0:
                    ang_part *= rz ** lz

                full_norm = rad_norm * angular_norms[bf_idx]
                bf_val += coeff * full_norm * ang_part * exp_part

            total += mo_coeffs[bf_idx] * bf_val

        out_values[p_idx] = total


@njit(parallel=True)
def _eval_basis_kernel(grid_points, atom_centers, basis_atom_idx, basis_shell_idx,
                       basis_lx, basis_ly, basis_lz, angular_norms,
                       shell_prim_start, shell_L, shell_cutoff_r2,
                       prim_exp, prim_coeff, prim_radial_norm,
                       out_basis):
    """Evaluate all basis functions on all grid points. Result: (npoints, nbasis)."""
    nbasis = basis_atom_idx.shape[0]
    npoints = grid_points.shape[0]

    for p_idx in prange(npoints):
        x = grid_points[p_idx, 0]
        y = grid_points[p_idx, 1]
        z = grid_points[p_idx, 2]

        for bf_idx in range(nbasis):
            aidx = basis_atom_idx[bf_idx]
            ax = atom_centers[aidx, 0]
            ay = atom_centers[aidx, 1]
            az = atom_centers[aidx, 2]
            rx = x - ax
            ry = y - ay
            rz = z - az
            r2 = rx * rx + ry * ry + rz * rz

            sidx = basis_shell_idx[bf_idx]
            if r2 > shell_cutoff_r2[sidx]:
                out_basis[p_idx, bf_idx] = 0.0
                continue

            p_start = shell_prim_start[sidx]
            p_end = shell_prim_start[sidx + 1]

            bf_val = 0.0
            for p in range(p_start, p_end):
                alpha = prim_exp[p]
                exp_part = math.exp(-alpha * r2)
                if exp_part < 1e-300:
                    continue
                coeff = prim_coeff[p]
                rad_norm = prim_radial_norm[p]

                ang_part = 1.0
                lx = basis_lx[bf_idx]
                ly = basis_ly[bf_idx]
                lz = basis_lz[bf_idx]
                if lx > 0:
                    ang_part *= rx ** lx
                if ly > 0:
                    ang_part *= ry ** ly
                if lz > 0:
                    ang_part *= rz ** lz

                full_norm = rad_norm * angular_norms[bf_idx]
                bf_val += coeff * full_norm * ang_part * exp_part

            out_basis[p_idx, bf_idx] = bf_val


@njit
def _project_mo_kernel(basis_values, mo_coeffs, out_values):
    """Fast MO projection: out = basis_values @ mo_coeffs."""
    npoints = basis_values.shape[0]
    nbasis = basis_values.shape[1]
    for p_idx in range(npoints):
        total = 0.0
        for bf_idx in range(nbasis):
            total += basis_values[p_idx, bf_idx] * mo_coeffs[bf_idx]
        out_values[p_idx] = total


# Basis function cache: shared across orbital switches at same grid spacing
_basis_cache = {}  # key: spacing (rounded), value: (grid_points, basis_values, origin, spacing, shape)


def _prepare_kernel_data(atoms, basis_set):
    """Precompute flat arrays for the numba kernels. Cached by basis_set identity."""
    natom = len(atoms)
    atom_centers = np.zeros((natom, 3), dtype=np.float64)
    for i, atom in enumerate(atoms):
        atom_centers[i, 0] = atom.x
        atom_centers[i, 1] = atom.y
        atom_centers[i, 2] = atom.z

    basis_atom_idx = basis_set.get_atom_idx_array()
    basis_shell_idx = basis_set.get_shell_idx_array()
    lmn = basis_set.get_lmn_array()
    basis_lx = lmn[:, 0]
    basis_ly = lmn[:, 1]
    basis_lz = lmn[:, 2]
    angular_norms = basis_set.get_angular_norm_array()

    shell_prim_start, shell_L, shell_cutoff_r2, prim_exp, prim_coeff, prim_radial_norm = \
        basis_set.flatten_primitives()

    return (atom_centers, basis_atom_idx, basis_shell_idx,
            basis_lx, basis_ly, basis_lz, angular_norms,
            shell_prim_start, shell_L, shell_cutoff_r2,
            prim_exp, prim_coeff, prim_radial_norm)


def _build_grid(atoms, grid_spacing, padding=4.0):
    """Build a 3D grid around the molecule."""
    coords = np.array([[a.x, a.y, a.z] for a in atoms], dtype=np.float64)
    xyz_min = coords.min(axis=0) - padding
    xyz_max = coords.max(axis=0) + padding

    nx = max(2, int(np.ceil((xyz_max[0] - xyz_min[0]) / grid_spacing)) + 1)
    ny = max(2, int(np.ceil((xyz_max[1] - xyz_min[1]) / grid_spacing)) + 1)
    nz = max(2, int(np.ceil((xyz_max[2] - xyz_min[2]) / grid_spacing)) + 1)

    x = np.linspace(xyz_min[0], xyz_max[0], nx, dtype=np.float64)
    y = np.linspace(xyz_min[1], xyz_max[1], ny, dtype=np.float64)
    z = np.linspace(xyz_min[2], xyz_max[2], nz, dtype=np.float64)

    XX, YY, ZZ = np.meshgrid(x, y, z, indexing='ij')
    grid_points = np.column_stack((XX.ravel(), YY.ravel(), ZZ.ravel()))

    return grid_points, xyz_min, grid_spacing, (nx, ny, nz)


def eval_mo_on_grid(atoms, basis_set, mo_coeffs, grid_spacing=0.1, padding=4.0,
                    use_cache=True):
    """Evaluate a molecular orbital on a 3D grid.

    Uses spatial cutoff for efficiency. Caches basis function values at each
    grid spacing so that switching orbitals only requires a fast dot product.
    """
    cache_key = round(grid_spacing, 2)

    # Check cache for basis function values
    if use_cache and cache_key in _basis_cache:
        cached_grid, cached_basis, cached_origin, cached_spacing, cached_shape = _basis_cache[cache_key]
        mo_coeffs_arr = np.asarray(mo_coeffs, dtype=np.float64).ravel()
        out_values = np.zeros(len(cached_grid), dtype=np.float64)
        _project_mo_kernel(cached_basis, mo_coeffs_arr, out_values)
        return out_values.reshape(cached_shape), cached_origin, cached_spacing

    # Build grid
    grid_points, origin, spacing, shape = _build_grid(atoms, grid_spacing, padding)

    # Prepare kernel data
    kdata = _prepare_kernel_data(atoms, basis_set)
    (atom_centers, basis_atom_idx, basis_shell_idx,
     basis_lx, basis_ly, basis_lz, angular_norms,
     shell_prim_start, shell_L, shell_cutoff_r2,
     prim_exp, prim_coeff, prim_radial_norm) = kdata

    mo_coeffs_arr = np.asarray(mo_coeffs, dtype=np.float64).ravel()
    out_values = np.zeros(len(grid_points), dtype=np.float64)

    # If caching is enabled, compute full basis matrix first
    if use_cache and grid_spacing >= 0.15:  # only cache coarse/medium grids (saves RAM)
        nbasis = basis_set.nbasis
        out_basis = np.zeros((len(grid_points), nbasis), dtype=np.float64)
        _eval_basis_kernel(
            grid_points, atom_centers, basis_atom_idx, basis_shell_idx,
            basis_lx, basis_ly, basis_lz, angular_norms,
            shell_prim_start, shell_L, shell_cutoff_r2,
            prim_exp, prim_coeff, prim_radial_norm,
            out_basis)
        # Cache for future orbital switches
        _basis_cache[cache_key] = (grid_points, out_basis, origin, spacing, shape)
        # Project to MO
        _project_mo_kernel(out_basis, mo_coeffs_arr, out_values)
    else:
        # Direct evaluation (no caching for fine grids — saves RAM)
        _eval_mo_kernel(
            grid_points, atom_centers, basis_atom_idx, basis_shell_idx,
            basis_lx, basis_ly, basis_lz, angular_norms,
            shell_prim_start, shell_L, shell_cutoff_r2,
            prim_exp, prim_coeff, prim_radial_norm,
            mo_coeffs_arr, out_values)

    return out_values.reshape(shape), origin, spacing


def clear_basis_cache():
    """Clear the global basis function cache (e.g., when loading a new file)."""
    _basis_cache.clear()


# ---------------------------------------------------------------------------
# Isosurface extraction
# ---------------------------------------------------------------------------

def extract_isosurface(grid_values, isovalue, origin, spacing):
    """Extract isosurface mesh using marching cubes.
    
    Parameters
    ----------
    grid_values : (nx, ny, nz) array
    isovalue : float
    origin : (3,) float — grid origin
    spacing : float — grid spacing
    
    Returns
    -------
    vertices : (N, 3) array or None
    faces : (M, 3) array or None
    normals : (N, 3) array or None
    """
    from skimage import measure
    
    try:
        verts, faces, normals, values = measure.marching_cubes(
            grid_values, level=isovalue, spacing=(spacing, spacing, spacing)
        )
        # Shift vertices to world coordinates
        verts = verts + origin
        return verts, faces, normals
    except (ValueError, RuntimeError):
        # No surface at this isovalue
        return None, None, None


# ---------------------------------------------------------------------------
# Bond detection
# ---------------------------------------------------------------------------

def detect_bonds(atoms, cutoff_factor=BOND_CUTOFF_FACTOR):
    """Detect bonds between atoms based on covalent radii.
    
    Returns list of (atom_idx1, atom_idx2) pairs.
    """
    bonds = []
    for i in range(len(atoms)):
        ri = COVALENT_RADII.get(atoms[i].atomic_number, 0.7)
        for j in range(i + 1, len(atoms)):
            rj = COVALENT_RADII.get(atoms[j].atomic_number, 0.7)
            cutoff = (ri + rj) * cutoff_factor
            
            dx = atoms[i].x - atoms[j].x
            dy = atoms[i].y - atoms[j].y
            dz = atoms[i].z - atoms[j].z
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)
            
            if dist < cutoff:
                bonds.append((i, j))
    
    return bonds


# ---------------------------------------------------------------------------
# vispy rendering
# ---------------------------------------------------------------------------

def _get_atom_color(atomic_number):
    """Return (r, g, b, alpha) for an atom."""
    rgb = CPK_COLORS.get(atomic_number, (0.7, 0.7, 0.7))
    return (rgb[0], rgb[1], rgb[2], 1.0)


def create_sphere_mesh(center, radius=0.3, color=(0.7, 0.7, 0.7, 1.0), 
                       rows=16, cols=16):
    """Create a sphere mesh as (vertices, faces, colors)."""
    from vispy.geometry import create_sphere
    
    mesh_data = create_sphere(radius=radius, rows=rows, cols=cols)
    verts = mesh_data.get_vertices() + center
    faces = mesh_data.get_faces()
    
    nv = len(verts)
    colors_arr = np.tile(np.array(color), (nv, 1))
    
    return verts, faces, colors_arr


def create_cylinder_mesh(p1, p2, radius=0.1, color=(0.7, 0.7, 0.7, 1.0),
                          rows=8, cols=8):
    """Create a cylinder mesh between two points."""
    from vispy.geometry import create_cylinder
    
    # Direction and length
    direction = np.array(p2) - np.array(p1)
    length = np.linalg.norm(direction)
    if length < 1e-6:
        return np.zeros((0, 3)), np.zeros((0, 3), dtype=np.int32), np.zeros((0, 4))
    
    direction = direction / length
    
    # Create cylinder along Z axis
    mesh_data = create_cylinder(rows=rows, cols=cols, radius=[radius, radius], length=length)
    verts = mesh_data.get_vertices()
    faces = mesh_data.get_faces()
    
    # Rotate from Z to direction
    z_axis = np.array([0.0, 0.0, 1.0])
    if np.allclose(direction, z_axis):
        rot_matrix = np.eye(3)
    elif np.allclose(direction, -z_axis):
        rot_matrix = np.diag([1.0, 1.0, -1.0])
    else:
        v = np.cross(z_axis, direction)
        s = np.linalg.norm(v)
        c = np.dot(z_axis, direction)
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        rot_matrix = np.eye(3) + vx + vx @ vx * ((1 - c) / (s * s))
    
    verts = verts @ rot_matrix.T + np.array(p1)
    
    nv = len(verts)
    colors_arr = np.tile(np.array(color), (nv, 1))
    
    return verts, faces, colors_arr


class OrbitalCanvas:
    """Manages the 3D scene with vispy, embedded in PyQt."""
    
    def __init__(self, parent=None):
        from vispy import scene
        self.canvas = scene.SceneCanvas(keys='interactive', size=(800, 600),
                                        bgcolor='black', show=False,
                                        parent=parent)
        self.view = self.canvas.central_widget.add_view()
        self.view.camera = scene.TurntableCamera(fov=60, distance=30)
        
        self.orbital_mesh_positive = None
        self.orbital_mesh_negative = None
        self.atom_markers = []
        self.bond_markers = []
        self._atoms_ref = None
        self._bonds_ref = None
        self._pos_color = (1.0, 0.2, 0.2, 0.6)
        self._neg_color = (0.2, 0.2, 1.0, 0.6)
    
    @property
    def native_widget(self):
        return self.canvas.native
    
    def clear_orbital(self):
        if self.orbital_mesh_positive is not None:
            self.orbital_mesh_positive.parent = None
            self.orbital_mesh_positive = None
        if self.orbital_mesh_negative is not None:
            self.orbital_mesh_negative.parent = None
            self.orbital_mesh_negative = None
    
    def clear_all(self):
        self.clear_orbital()
        for marker in self.atom_markers:
            marker.parent = None
        self.atom_markers = []
        for marker in self.bond_markers:
            marker.parent = None
        self.bond_markers = []
    
    def add_atoms_and_bonds(self, atoms, bonds):
        from vispy import scene
        self._atoms_ref = atoms
        self._bonds_ref = bonds
        for i, atom in enumerate(atoms):
            color = _get_atom_color(atom.atomic_number)
            radius = 0.2 if atom.atomic_number == 1 else 0.35 if atom.atomic_number == 6 else 0.3
            verts, faces, colors = create_sphere_mesh(
                (atom.x, atom.y, atom.z), radius=radius, color=color)
            mesh = scene.visuals.Mesh(vertices=verts, faces=faces,
                                      vertex_colors=colors, shading='smooth',
                                      parent=self.view.scene)
            self.atom_markers.append(mesh)
        for i, j in bonds:
            p1 = (atoms[i].x, atoms[i].y, atoms[i].z)
            p2 = (atoms[j].x, atoms[j].y, atoms[j].z)
            verts, faces, colors = create_cylinder_mesh(p1, p2, radius=0.12, color=(0.5, 0.5, 0.5, 1.0))
            if len(verts) > 0:
                mesh = scene.visuals.Mesh(vertices=verts, faces=faces,
                                          vertex_colors=colors, shading='smooth',
                                          parent=self.view.scene)
                self.bond_markers.append(mesh)
    
    def set_orbital_surface(self, verts_pos, faces_pos, verts_neg, faces_neg):
        from vispy import scene
        self.clear_orbital()
        if verts_pos is not None and len(verts_pos) > 0:
            nv = len(verts_pos)
            colors_arr = np.tile(np.array(self._pos_color), (nv, 1))
            self.orbital_mesh_positive = scene.visuals.Mesh(
                vertices=verts_pos, faces=faces_pos,
                vertex_colors=colors_arr, shading='smooth',
                parent=self.view.scene)
            self.orbital_mesh_positive.set_gl_state('translucent', depth_test=True, cull_face=False)
        if verts_neg is not None and len(verts_neg) > 0:
            nv = len(verts_neg)
            colors_arr = np.tile(np.array(self._neg_color), (nv, 1))
            self.orbital_mesh_negative = scene.visuals.Mesh(
                vertices=verts_neg, faces=faces_neg,
                vertex_colors=colors_arr, shading='smooth',
                parent=self.view.scene)
            self.orbital_mesh_negative.set_gl_state('translucent', depth_test=True, cull_face=False)
        self.canvas.update()
    
    def set_camera_center(self, atoms=None):
        if atoms is None:
            atoms = self._atoms_ref
        if atoms is not None:
            center = np.array([[a.x, a.y, a.z] for a in atoms]).mean(axis=0)
            self.view.camera.center = center
    
    def screenshot(self, filename='orbital.png'):
        img = self.canvas.render()
        from vispy.io import imsave
        imsave(filename, img)


# ---------------------------------------------------------------------------
# Computation worker (QThread)
# ---------------------------------------------------------------------------

from PyQt6.QtCore import QThread, pyqtSignal


class GridWorker(QThread):
    """Background thread for MO grid computation."""
    finished = pyqtSignal(object, object, float, int)  # grid_values, origin, spacing, mo_idx
    
    def __init__(self, atoms, basis_set, mo_coeffs, grid_spacing, mo_idx, parent=None):
        super().__init__(parent)
        self.atoms = atoms
        self.basis_set = basis_set
        self.mo_coeffs = np.asarray(mo_coeffs, dtype=np.float64)
        self.grid_spacing = grid_spacing
        self.mo_idx = mo_idx
        self._cancelled = False
    
    def cancel(self):
        self._cancelled = True
    
    def run(self):
        if self._cancelled:
            return
        grid_values, origin, spacing = eval_mo_on_grid(
            self.atoms, self.basis_set, self.mo_coeffs,
            grid_spacing=self.grid_spacing
        )
        if not self._cancelled:
            self.finished.emit(grid_values, origin, spacing, self.mo_idx)


# ---------------------------------------------------------------------------
# Main GUI window
# ---------------------------------------------------------------------------

from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QSplitter, QListWidget, QListWidgetItem, QTabWidget,
    QLabel, QSlider, QPushButton, QFileDialog, QStatusBar,
    QMessageBox, QApplication
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QAction


class OrbitalViewer(QMainWindow):
    """Main application window."""
    REFINEMENT_STEPS = [0.35, 0.20, 0.10]
    
    def __init__(self, logpath=None):
        super().__init__()
        self.setWindowTitle("Orbital Visualizer")
        self.resize(1200, 800)
        
        # Data
        self.atoms = None
        self.basis_set = None
        self.canon_wfn = None
        self.local_wfn = None
        self.homo_idx = 0
        self.active_wfn = None
        self.active_mo_idx = -1
        self._current_grid_values = None
        self._current_origin = None
        self._current_spacing = None
        self._grid_worker = None
        self._refine_step = 0
        self._generation = 0  # increment on each new orbital request
        self._isovalue = 0.05
        self._grid_spacing = 0.35  # coarse by default
        
        self._build_ui()
        self._build_menu()
        
        if logpath is not None:
            self._load_file(Path(logpath))
    
    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QVBoxLayout(central)
        main_layout.setContentsMargins(4, 4, 4, 4)
        
        # --- Splitter: gallery | 3D view ---
        splitter = QSplitter(Qt.Orientation.Horizontal)
        
        # Left: gallery
        self.gallery_tabs = QTabWidget()
        self.occupied_list = QListWidget()
        self.occupied_list.currentRowChanged.connect(self._on_occupied_selected)
        self.gallery_tabs.addTab(self.occupied_list, "Occupied")
        self.virtual_list = QListWidget()
        self.virtual_list.currentRowChanged.connect(self._on_virtual_selected)
        self.gallery_tabs.addTab(self.virtual_list, "Virtual")
        self.localized_list = QListWidget()
        self.localized_list.currentRowChanged.connect(self._on_localized_selected)
        splitter.addWidget(self.gallery_tabs)
        
        # Right: canvas + controls
        right_side = QWidget()
        right_layout = QVBoxLayout(right_side)
        right_layout.setContentsMargins(0, 0, 0, 0)
        self.orbital_canvas = OrbitalCanvas()
        right_layout.addWidget(self.orbital_canvas.native_widget, 1)
        
        # Control bar
        ctrl = QWidget()
        ctrl_layout = QHBoxLayout(ctrl)
        ctrl_layout.setContentsMargins(4, 2, 4, 2)
        ctrl_layout.addWidget(QLabel("Isovalue:"))
        self.isovalue_slider = QSlider(Qt.Orientation.Horizontal)
        self.isovalue_slider.setRange(10, 200)
        self.isovalue_slider.setValue(50)
        self.isovalue_slider.valueChanged.connect(self._on_isovalue_changed)
        ctrl_layout.addWidget(self.isovalue_slider)
        self.isovalue_label = QLabel("0.050")
        self.isovalue_label.setFixedWidth(40)
        ctrl_layout.addWidget(self.isovalue_label)
        ctrl_layout.addSpacing(20)
        ctrl_layout.addWidget(QLabel("Grid:"))
        self.grid_slider = QSlider(Qt.Orientation.Horizontal)
        self.grid_slider.setRange(5, 40)
        self.grid_slider.setValue(35)
        self.grid_slider.valueChanged.connect(self._on_grid_changed)
        ctrl_layout.addWidget(self.grid_slider)
        self.grid_label = QLabel("0.35 Å")
        self.grid_label.setFixedWidth(50)
        ctrl_layout.addWidget(self.grid_label)

        # Grid quality preset buttons
        self.btn_coarse = QPushButton("Coarse")
        self.btn_coarse.clicked.connect(lambda: self._set_grid_preset(0.35, 35))
        ctrl_layout.addWidget(self.btn_coarse)
        self.btn_medium = QPushButton("Medium")
        self.btn_medium.clicked.connect(lambda: self._set_grid_preset(0.20, 20))
        ctrl_layout.addWidget(self.btn_medium)
        self.btn_fine = QPushButton("Fine")
        self.btn_fine.clicked.connect(lambda: self._set_grid_preset(0.10, 10))
        ctrl_layout.addWidget(self.btn_fine)

        ctrl_layout.addStretch()
        self.homo_btn = QPushButton("HOMO")
        self.homo_btn.clicked.connect(self._go_to_homo)
        ctrl_layout.addWidget(self.homo_btn)
        self.lumo_btn = QPushButton("LUMO")
        self.lumo_btn.clicked.connect(self._go_to_lumo)
        ctrl_layout.addWidget(self.lumo_btn)
        right_layout.addWidget(ctrl)
        
        splitter.addWidget(right_side)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([250, 950])
        main_layout.addWidget(splitter)
        
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready — Open a GAMESS .log file (Ctrl+O)")
    
    def _build_menu(self):
        mb = self.menuBar()
        fm = mb.addMenu("&File")
        a = QAction("&Open...", self); a.setShortcut("Ctrl+O"); a.triggered.connect(self._file_open); fm.addAction(a)
        fm.addSeparator()
        a = QAction("&Export Image...", self); a.setShortcut("Ctrl+E"); a.triggered.connect(self._file_export); fm.addAction(a)
        fm.addSeparator()
        a = QAction("&Quit", self); a.setShortcut("Ctrl+Q"); a.triggered.connect(self.close); fm.addAction(a)
    
    def _load_file(self, logpath):
        self.status_bar.showMessage(f"Loading {logpath.name} ...")
        QApplication.processEvents()
        try:
            atoms, basis_set, canon_wfn, local_wfn, homo_idx = parse_gamess_log(logpath)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to parse file:\n{e}")
            self.status_bar.showMessage("Error loading file")
            return

        # Clear basis cache — new molecule, new basis set
        clear_basis_cache()

        self.atoms = atoms
        self.basis_set = basis_set
        self.canon_wfn = canon_wfn
        self.local_wfn = local_wfn
        self.homo_idx = homo_idx
        self.setWindowTitle(f"Orbital Visualizer — {logpath.name}")
        
        bonds = detect_bonds(atoms)
        self.orbital_canvas.add_atoms_and_bonds(atoms, bonds)
        self.orbital_canvas.set_camera_center(atoms)
        self._populate_gallery()
        self.status_bar.showMessage(
            f"Loaded: {len(atoms)} atoms, {basis_set.nbasis} bf, "
            f"{canon_wfn.nmo} MOs, HOMO={homo_idx + 1}")
    
    def _populate_gallery(self):
        wfn = self.canon_wfn
        if wfn is None:
            return
        self.occupied_list.clear()
        for i in range(self.homo_idx + 1):
            e = wfn.energies[i]
            item = QListWidgetItem(f"MO {i+1}  ({e:+.3f} Eh)")
            item.setData(Qt.ItemDataRole.UserRole, i)
            self.occupied_list.addItem(item)
        self.virtual_list.clear()
        for i in range(self.homo_idx + 1, wfn.nmo):
            e = wfn.energies[i]
            item = QListWidgetItem(f"MO {i+1}  ({e:+.3f} Eh)")
            item.setData(Qt.ItemDataRole.UserRole, i)
            self.virtual_list.addItem(item)
        # Remove localized tab if it exists (from previous load)
        for t in range(self.gallery_tabs.count()):
            if self.gallery_tabs.tabText(t) == "Localized":
                self.gallery_tabs.removeTab(t)
                break
        if self.local_wfn is not None:
            self.localized_list.clear()
            for i in range(self.local_wfn.nmo):
                item = QListWidgetItem(f"Loc {i+1}")
                item.setData(Qt.ItemDataRole.UserRole, i)
                self.localized_list.addItem(item)
            self.gallery_tabs.addTab(self.localized_list, "Localized")
    
    # -- Orbital selection --
    
    def _select_orbital(self, mo_idx, wtype='canonical'):
        if self.atoms is None:
            return
        self._cancel_current_computation()
        self._generation += 1
        gen = self._generation
        wfn = self.canon_wfn if wtype == 'canonical' else self.local_wfn
        if wfn is None or mo_idx < 0 or mo_idx >= wfn.nmo:
            return
        self.active_wfn = wfn
        self.active_mo_idx = mo_idx
        self._refine_step = 0
        self.status_bar.showMessage(f"Computing MO {mo_idx+1} ({wfn.labels[mo_idx]}) ...")
        # Start coarse first
        self._start_grid_computation(mo_idx, self.REFINEMENT_STEPS[0], gen)
    
    def _start_grid_computation(self, mo_idx, spacing, generation):
        mo_coeffs = self.active_wfn.get_mo(mo_idx)
        worker = GridWorker(self.atoms, self.basis_set, mo_coeffs, spacing, mo_idx)
        worker.finished.connect(
            lambda gv, o, s, mi: self._on_grid_done(gv, o, s, mi, generation))
        worker.finished.connect(lambda *a: self._cleanup_worker(worker))
        worker.start()
        self._grid_worker = worker
    
    def _cleanup_worker(self, worker):
        """Called when worker finishes — clean up if it's the current one."""
        if self._grid_worker is worker:
            worker.wait(500)  # ensure fully done
            if worker is self._grid_worker:
                self._grid_worker = None
    
    def _on_grid_done(self, grid_values, origin, spacing, mo_idx, generation):
        if generation != self._generation or mo_idx != self.active_mo_idx:
            return
        self._current_grid_values = grid_values
        self._current_origin = origin
        self._current_spacing = spacing
        self._update_surface()
        self._refine_step += 1
        if self._refine_step < len(self.REFINEMENT_STEPS):
            next_sp = max(self.REFINEMENT_STEPS[self._refine_step], self._grid_spacing)
            if next_sp < spacing:
                self.status_bar.showMessage(f"Refining MO {mo_idx+1} ({spacing:.2f}→{next_sp:.2f} Å) ...")
                self._start_grid_computation(mo_idx, next_sp, generation)
                return
        self.status_bar.showMessage(f"MO {mo_idx+1} ready ({spacing:.2f} Å grid)")
    
    def _update_surface(self):
        if self._current_grid_values is None:
            return
        vp, fp, _ = extract_isosurface(self._current_grid_values, +self._isovalue,
                                        self._current_origin, self._current_spacing)
        vn, fn, _ = extract_isosurface(self._current_grid_values, -self._isovalue,
                                        self._current_origin, self._current_spacing)
        self.orbital_canvas.set_orbital_surface(vp, fp, vn, fn)
    
    def _cancel_current_computation(self):
        if self._grid_worker is not None:
            if self._grid_worker.isRunning():
                self._grid_worker.cancel()
                # Wait up to 3 seconds for the thread to finish.
                # Don't destroy the QThread while numba kernel is still running.
                if not self._grid_worker.wait(3000):
                    # Still running — detach, don't null the reference.
                    # The generation counter will reject its stale result.
                    return
            # Only null the reference if thread is done
            self._grid_worker = None
    
    # -- Slots --
    
    def _on_occupied_selected(self, row):
        if row < 0: return
        mo_idx = self.occupied_list.item(row).data(Qt.ItemDataRole.UserRole)
        self._select_orbital(mo_idx, 'canonical')
    
    def _on_virtual_selected(self, row):
        if row < 0: return
        mo_idx = self.virtual_list.item(row).data(Qt.ItemDataRole.UserRole)
        self._select_orbital(mo_idx, 'canonical')
    
    def _on_localized_selected(self, row):
        if row < 0: return
        mo_idx = self.localized_list.item(row).data(Qt.ItemDataRole.UserRole)
        self._select_orbital(mo_idx, 'localized')
    
    def _on_isovalue_changed(self, value):
        self._isovalue = value / 1000.0
        self.isovalue_label.setText(f"{self._isovalue:.3f}")
        if self._current_grid_values is not None:
            self._update_surface()
    
    def _on_grid_changed(self, value):
        self._grid_spacing = value / 100.0
        self.grid_label.setText(f"{self._grid_spacing:.2f} Å")
        if self.active_mo_idx >= 0:
            self._cancel_current_computation()
            self._generation += 1
            gen = self._generation
            self._refine_step = 0
            self._start_grid_computation(self.active_mo_idx, self._grid_spacing, gen)

    def _set_grid_preset(self, spacing, slider_value):
        """Set grid to a preset quality (Coarse/Medium/Fine)."""
        self.grid_slider.setValue(slider_value)
        # Slider valueChanged will trigger _on_grid_changed
    
    def _go_to_homo(self):
        self._select_orbital(self.homo_idx, 'canonical')
        self.gallery_tabs.setCurrentIndex(0)
        self.occupied_list.setCurrentRow(self.homo_idx)
    
    def _go_to_lumo(self):
        lumo_idx = self.homo_idx + 1
        self._select_orbital(lumo_idx, 'canonical')
        self.gallery_tabs.setCurrentIndex(1)
        self.virtual_list.setCurrentRow(0)
    
    def closeEvent(self, event):
        """Safely shut down worker threads before closing."""
        if self._grid_worker is not None:
            if self._grid_worker.isRunning():
                self._grid_worker.cancel()
                self._grid_worker.wait(3000)
            self._grid_worker = None
        self.orbital_canvas.canvas.close()
        super().closeEvent(event)
    
    def _file_open(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open GAMESS Log", "", "GAMESS Log Files (*.log *.out);;All Files (*)")
        if path:
            self._cancel_current_computation()
            self.orbital_canvas.clear_all()
            self._current_grid_values = None
            self._load_file(Path(path))
    
    def _file_export(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Image", "orbital.png", "PNG Images (*.png);;All Files (*)")
        if path:
            self._export_with_overlay(path)
            self.status_bar.showMessage(f"Saved: {path}")

    def _export_with_overlay(self, path):
        """Render viewport and overlay orbital info, then save to PNG."""
        from PyQt6.QtGui import QImage, QPainter, QColor, QFont

        # Get vispy render
        img_array = self.orbital_canvas.canvas.render()
        h, w = img_array.shape[:2]

        # Convert to QImage (RGBA8888)
        img_8bit = (np.clip(img_array, 0, 1) * 255).astype(np.uint8)
        qimg = QImage(img_8bit.data, w, h, w * 4, QImage.Format.Format_RGBA8888).copy()

        # Build overlay text
        lines = []
        if self.active_mo_idx >= 0 and self.active_wfn is not None:
            label = self.active_wfn.labels[self.active_mo_idx]
            lines.append(f"Orbital {self.active_mo_idx + 1} ({label})")
            if self.active_wfn is self.canon_wfn and self.active_mo_idx < len(self.active_wfn.energies):
                e = self.active_wfn.energies[self.active_mo_idx]
                lines.append(f"Energy: {e:+.6f} Eh  ({e * 27.2114:+.2f} eV)")
            lines.append(f"Isovalue: ±{self._isovalue:.3f}")
        if self._current_spacing is not None:
            lines.append(f"Grid: {self._current_spacing:.2f} Å")

        if not lines:
            qimg.save(path, 'PNG')
            return

        # Draw overlay text
        painter = QPainter(qimg)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        font = QFont("monospace", 12)
        font.setBold(True)
        painter.setFont(font)

        # Semi-transparent background for readability
        padding = 8
        line_height = 20
        box_height = len(lines) * line_height + padding * 2
        box_width = 0
        metrics = painter.fontMetrics()
        for line in lines:
            box_width = max(box_width, metrics.horizontalAdvance(line))
        box_width += padding * 2

        # Draw background box (top-left)
        painter.fillRect(0, 0, box_width, box_height, QColor(0, 0, 0, 180))
        painter.setPen(QColor(255, 255, 255, 255))

        y = padding + metrics.ascent()
        for line in lines:
            painter.drawText(padding, y, line)
            y += line_height

        painter.end()
        qimg.save(path, 'PNG')


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def cli_main():
    """Command-line mode for headless rendering."""
    parser = argparse.ArgumentParser(description='Orbital Visualizer')
    parser.add_argument('logfile', nargs='?', default='spval.log',
                        help='GAMESS .log file to visualize')
    parser.add_argument('--orbital', type=int, default=None,
                        help='Orbital index to render (1-based)')
    parser.add_argument('--isovalue', type=float, default=0.05,
                        help='Isosurface value (default: 0.05)')
    parser.add_argument('--grid', type=float, default=0.15,
                        help='Grid spacing in Angstrom (default: 0.15)')
    parser.add_argument('--output', type=str, default='orbital.png',
                        help='Output image filename')
    parser.add_argument('--type', type=str, default='canonical',
                        choices=['canonical', 'localized'],
                        help='Orbital type (default: canonical)')
    parser.add_argument('--cli', action='store_true', default=True,
                        help='Run in CLI mode (default when --output given)')
    args = parser.parse_args()
    
    logpath = Path(args.logfile)
    if not logpath.exists():
        print(f"Error: file not found: {args.logfile}")
        sys.exit(1)
    
    print(f"Parsing {args.logfile} ...")
    atoms, basis_set, canon_wfn, local_wfn, homo_idx = parse_gamess_log(logpath)
    
    print(f"  Atoms: {len(atoms)}")
    print(f"  Basis functions: {basis_set.nbasis}")
    print(f"  Canonical MOs: {canon_wfn.nmo}")
    print(f"  HOMO index (0-based): {homo_idx}")
    if local_wfn:
        print(f"  Localized MOs: {local_wfn.nmo}")
    
    if args.type == 'localized' and local_wfn is not None:
        wfn = local_wfn
    else:
        wfn = canon_wfn
    
    if args.orbital is not None:
        mo_idx = args.orbital - 1
        if mo_idx < 0 or mo_idx >= wfn.nmo:
            print(f"Error: orbital index out of range (1-{wfn.nmo})")
            sys.exit(1)
    else:
        mo_idx = homo_idx
    
    print(f"\nEvaluating orbital {mo_idx + 1} ({wfn.labels[mo_idx]}) "
          f"on {args.grid:.2f} Å grid ...")
    
    mo_coeffs = wfn.get_mo(mo_idx)
    grid_values, origin, spacing = eval_mo_on_grid(
        atoms, basis_set, mo_coeffs, grid_spacing=args.grid
    )
    
    print(f"  Grid shape: {grid_values.shape}")
    print(f"  Value range: [{grid_values.min():.6f}, {grid_values.max():.6f}]")
    
    verts_pos, faces_pos, _ = extract_isosurface(grid_values, +args.isovalue, origin, spacing)
    verts_neg, faces_neg, _ = extract_isosurface(grid_values, -args.isovalue, origin, spacing)
    
    if verts_pos is not None:
        print(f"  Positive lobe: {len(verts_pos)} vertices, {len(faces_pos)} faces")
    if verts_neg is not None:
        print(f"  Negative lobe: {len(verts_neg)} vertices, {len(faces_neg)} faces")
    
    print(f"\nRendering ...")
    canvas = OrbitalCanvas()
    bonds = detect_bonds(atoms)
    canvas.add_atoms_and_bonds(atoms, bonds)
    canvas.set_orbital_surface(verts_pos, faces_pos, verts_neg, faces_neg)
    canvas.set_camera_center(atoms)
    
    print(f"  Saving to {args.output} ...")
    from vispy import app
    canvas.canvas.show()
    app.process_events()
    canvas.screenshot(args.output)
    canvas.canvas.close()
    print(f"Done. Output: {args.output}")


def main():
    """Launch the GUI or CLI depending on arguments."""
    parser = argparse.ArgumentParser(description='Orbital Visualizer')
    parser.add_argument('logfile', nargs='?', default=None,
                        help='GAMESS .log file to visualize')
    parser.add_argument('--cli', action='store_true', default=False,
                        help='Run in command-line mode (render to PNG)')
    args, remaining = parser.parse_known_args()
    
    if args.cli:
        sys.argv = [sys.argv[0]] + remaining
        cli_main()
        return
    
    # GUI mode: initialize vispy with PyQt6 backend before creating canvas
    from vispy.app import use_app
    use_app('pyqt6')
    
    from PyQt6.QtWidgets import QApplication
    app = QApplication(sys.argv)
    app.setApplicationName("Orbital Visualizer")
    
    viewer = OrbitalViewer(logpath=args.logfile)
    viewer.show()
    
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
