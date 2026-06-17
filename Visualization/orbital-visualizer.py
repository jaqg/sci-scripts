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
# FChk helper: read a data section from Gaussian formatted checkpoint text
# ---------------------------------------------------------------------------

def _read_fchk_array(text, field_name, dtype):
    """Read a data section from Gaussian fchk text (scalar or array).

    Returns numpy array of dtype (np.int32 or np.float64).
    """
    import re
    escaped = re.escape(field_name)

    # Try array field: "Field Name         I/R   N=   NNN"
    pattern_arr = escaped + r'\s+[IRC]\s+N=\s+(\d+)'
    match = re.search(pattern_arr, text)

    if match:
        nvals = int(match.group(1))
        data_start = text.find('\n', match.end()) + 1
        if data_start == 0:
            raise ValueError(f"Unexpected end of file after header '{field_name}'")

        # Find next header (line starting with letter, containing I/R/C type tag)
        next_header = float('inf')
        for m in re.finditer(r'\n([A-Za-z])', text[data_start:]):
            pos = data_start + m.start() + 1
            snippet = text[pos:pos + 80]
            if re.match(r'[A-Za-z].{20,}?\s+[IRC]\s', snippet):
                next_header = pos
                break

        data_block = text[data_start:next_header]

        if dtype == int:
            values = [int(x) for x in data_block.split()]
        else:
            values = [float(x.replace('D', 'E')) for x in data_block.split()]

        if len(values) != nvals:
            raise ValueError(
                f"Field '{field_name}': expected {nvals} values, got {len(values)}"
            )

        dtype_np = np.float64 if dtype == float else np.int32
        return np.array(values, dtype=dtype_np)

    # Scalar field: "Field Name         I/R            value"
    pattern_scalar = escaped + r'\s+[IRC]\s+(-?\d+\.?\d*(?:[EeDd][+-]?\d+)?)'
    match = re.search(pattern_scalar, text)
    if not match:
        raise ValueError(f"Field '{field_name}' not found in fchk file")

    val_str = match.group(1).replace('D', 'E')
    if dtype == int:
        return np.array([int(float(val_str))], dtype=np.int32)
    else:
        return np.array([float(val_str)], dtype=np.float64)


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
# Gaussian formatted checkpoint (.fchk) parsing
# ---------------------------------------------------------------------------

# Gaussian fchk angular momentum codes → our string labels
# Negative → spherical (pure); positive/non-negative → Cartesian
# 0=S (same either way), 1=P (Cartesian), -2=D (spherical), 2=D (Cartesian), etc.
_FCHK_ANG_MOM_LABEL = {0: 'S', 1: 'P', 2: 'D', 3: 'F', 4: 'G'}

# Number of Cartesian components per angular momentum
_FCHK_CART_SIZE = {'S': 1, 'P': 3, 'D': 6, 'F': 10, 'G': 15}
# Number of spherical components per angular momentum
_FCHK_SPH_SIZE = {'S': 1, 'P': 3, 'D': 5, 'F': 7, 'G': 9}


def parse_gaussian_fchk(filepath):
    """Parse a Gaussian .fchk file and return atoms, BasisSet, Wavefunction.

    Returns:
        atoms: list of Atom
        basis_set: BasisSet (Cartesian expansion for evaluation)
        canon_wfn: Wavefunction (MO coefficients transformed to Cartesian basis)
        local_wfn: None (localized orbitals not available in .fchk)
        homo_idx: int (0-based index of HOMO)
        mocoeffs_sph: (nmo, nsph) array — original spherical MO coeffs for NTO calc
        T_sph_to_cart: (nsph, ncart) array — transformation matrix
    """
    import cclib

    data = cclib.io.ccread(str(filepath))

    # --- Atoms / Molecule ---
    symbols = _symbols_from_atomnos(data.atomnos)
    coords = data.atomcoords[0]  # (natom, 3) in Angstrom (cclib converts from Bohr)
    atoms = []
    for i in range(data.natom):
        atoms.append(Atom(i, symbols[i], int(data.atomnos[i]),
                          float(coords[i, 0]), float(coords[i, 1]), float(coords[i, 2])))

    # --- Parse shells with spherical/Cartesian info ---
    with open(filepath, 'r') as f:
        fchk_text = f.read()

    shells, nbasis_fchk, is_spherical = _parse_fchk_basis_shells(fchk_text)

    # Build Cartesian BasisSet (all shells expanded to Cartesian for evaluation)
    basis_set = BasisSet(shells, atoms)

    # Build spherical function count to compute the transformation
    nsph, ncart = _count_basis_sizes(shells, is_spherical)

    if nsph != nbasis_fchk:
        raise ValueError(
            f"Spherical basis count mismatch: counted {nsph}, expected {nbasis_fchk}"
        )

    # --- Canonical MOs from cclib ---
    if isinstance(data.mocoeffs, list):
        mocoeffs_sph = np.array(data.mocoeffs[0], dtype=np.float64)
    else:
        mocoeffs_sph = np.array(data.mocoeffs, dtype=np.float64)

    if isinstance(data.moenergies, list):
        moenergies_arr = np.array(data.moenergies[0], dtype=np.float64)
    else:
        moenergies_arr = np.array(data.moenergies, dtype=np.float64)

    # Build spherical→Cartesian transformation matrix and transform MO coefficients
    T = _build_sph_to_cart_transform(shells, is_spherical)
    # T: (nsph, ncart), C_sph: (nmo, nsph) → C_cart: (nmo, ncart)
    mocoeffs_cart = np.dot(mocoeffs_sph, T)

    homo_idx = data.homos[0]
    canon_labels = ['canonical'] * mocoeffs_cart.shape[0]
    canon_wfn = Wavefunction(mocoeffs_cart, moenergies_arr, canon_labels)

    return atoms, basis_set, canon_wfn, None, homo_idx, mocoeffs_sph, T


def _parse_fchk_basis_shells(fchk_text):
    """Parse basis set shells from Gaussian fchk text.

    In Gaussian fchk:
      - Negative shell type → pure/spherical harmonics (fewer functions)
      - Positive shell type → Cartesian functions
      - 0 → S (same either way)

    Returns:
        shells: list of Shell namedtuples — one per contracted shell
        nbasis_expected: int — number of basis functions according to fchk header
        is_spherical: list of bool — True if shell uses spherical harmonics
    """
    shell_types = _read_fchk_array(fchk_text, 'Shell types', int)
    nprims_per_shell = _read_fchk_array(fchk_text, 'Number of primitives per shell', int)
    shell_to_atom = _read_fchk_array(fchk_text, 'Shell to atom map', int)
    prim_exps = _read_fchk_array(fchk_text, 'Primitive exponents', float)
    prim_coeffs = _read_fchk_array(fchk_text, 'Contraction coefficients', float)

    nshells = len(shell_types)
    nbasis_expected = int(_read_fchk_array(fchk_text, 'Number of basis functions', int)[0])

    shells = []
    is_spherical = []
    prim_offset = 0
    for ishell in range(nshells):
        ang_code = int(shell_types[ishell])
        # abs(code) = angular momentum L. sign indicates spherical vs Cartesian.
        # negative → spherical (pure), non-negative → Cartesian
        # Exception: 0 → S (always spherical since only 1 component)
        ang_mom_abs = abs(ang_code)
        ang_mom_label = _FCHK_ANG_MOM_LABEL.get(ang_mom_abs)
        if ang_mom_label is None:
            raise ValueError(f"Unsupported angular momentum L={ang_mom_abs} in shell {ishell}")

        use_spherical = ang_code < 0 or ang_code == 0
        is_spherical.append(use_spherical)

        nprim = int(nprims_per_shell[ishell])
        prims = []
        for iprim in range(nprim):
            exp = float(prim_exps[prim_offset + iprim])
            coeff = float(prim_coeffs[prim_offset + iprim])
            prims.append((exp, coeff))
        prim_offset += nprim

        atom_idx = int(shell_to_atom[ishell]) - 1  # 1-based → 0-based
        shells.append(Shell(atom_idx, ang_mom_label, prims))

    return shells, nbasis_expected, is_spherical


def _count_basis_sizes(shells, is_spherical):
    """Count total number of MO-coefficient (spherical side) and Cartesian basis functions."""
    nsph = 0
    ncart = 0
    for shell, sph in zip(shells, is_spherical):
        lbl = shell.ang_mom
        if sph and lbl in _FCHK_SPH_SIZE:
            # Spherical shells: fewer functions in MO basis
            nsph += _FCHK_SPH_SIZE[lbl]
        else:
            # Cartesian (or S/P): same number in both
            nsph += _FCHK_CART_SIZE.get(lbl, 1)
        ncart += _FCHK_CART_SIZE.get(lbl, 1)
    return nsph, ncart


def _build_sph_to_cart_transform(shells, is_spherical):
    """Build the spherical→Cartesian transformation matrix for the full basis.

    For each shell, constructs the (n_sph, n_cart) block that maps Cartesian
    basis functions to the spherical functions used by Gaussian. The MO coefficients
    from cclib are in the spherical basis; multiplying C_sph @ T gives C_cart.

    Returns:
        T: (nsph_total, ncart_total) numpy array
    """
    nsph_total, ncart_total = _count_basis_sizes(shells, is_spherical)
    T = np.zeros((nsph_total, ncart_total), dtype=np.float64)

    sph_offset = 0
    cart_offset = 0

    for shell, sph in zip(shells, is_spherical):
        lbl = shell.ang_mom
        n_cart_shell = _FCHK_CART_SIZE.get(lbl, 1)
        if sph and lbl in ('D', 'F'):
            T_block = _make_shell_transform(lbl)
            n_sph_shell = T_block.shape[0]
            T[sph_offset:sph_offset + n_sph_shell,
              cart_offset:cart_offset + n_cart_shell] = T_block
            sph_offset += n_sph_shell
        elif sph and lbl == 'S':
            T[sph_offset, cart_offset] = 1.0
            sph_offset += 1
        elif sph and lbl == 'P':
            # P shells in fchk are always Cartesian (3 spherical = 3 Cartesian)
            T[sph_offset:sph_offset + 3, cart_offset:cart_offset + 3] = np.eye(3)
            sph_offset += 3
        else:
            # Cartesian shell → Cartesian target: identity block
            T[sph_offset:sph_offset + n_cart_shell,
              cart_offset:cart_offset + n_cart_shell] = np.eye(n_cart_shell)
            sph_offset += n_cart_shell
        cart_offset += n_cart_shell

    return T


def _make_shell_transform(lbl):
    """Return the spherical→Cartesian transformation block for one shell.

    Rows: spherical functions (HORTON/Gaussian order)
    Cols: Cartesian functions (our internal order from CARTESIAN_MAP)

    Coefficients taken from HORTON's normalized transformation matrices,
    permuted to match our Cartesian ordering.
    """
    if lbl == 'D':
        # Spherical order: C₂₀, C₂₁, S₂₁, C₂₂, S₂₂
        # Our Cartesian order: xx(0), yy(1), zz(2), xy(3), xz(4), yz(5)
        # From HORTON (normalized):
        #   C₂₀ = -0.5·xx -0.5·yy + zz
        #   C₂₁ = xz
        #   S₂₁ = yz
        #   C₂₂ = √3/2·xx -√3/2·yy
        #   S₂₂ = xy
        sqrt3 = math.sqrt(3.0)
        return np.array([
            [-0.5,     -0.5,      1.0,  0.0,   0.0,   0.0],
            [ 0.0,      0.0,      0.0,  0.0,   1.0,   0.0],
            [ 0.0,      0.0,      0.0,  0.0,   0.0,   1.0],
            [ 0.5*sqrt3, -0.5*sqrt3, 0.0, 0.0,   0.0,   0.0],
            [ 0.0,      0.0,      0.0,  1.0,   0.0,   0.0],
        ], dtype=np.float64)

    elif lbl == 'F':
        # Spherical order: C₃₀, C₃₁, S₃₁, C₃₂, S₃₂, C₃₃, S₃₃
        # Our Cartesian order: xxx(0), yyy(1), zzz(2), xxy(3), xxz(4),
        #                       xyy(5), yyz(6), xzz(7), yzz(8), xyz(9)
        # From HORTON (normalized), permuted to our column order.
        # HORTON col order: xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz
        # Permutation to our order:
        perm = [0, 6, 9, 1, 2, 3, 7, 5, 8, 4]
        s5 = math.sqrt(5.0)
        s6 = math.sqrt(6.0)
        s30 = math.sqrt(30.0)
        s3 = math.sqrt(3.0)
        s10 = math.sqrt(10.0)
        s2 = math.sqrt(2.0)
        Th = np.array([
            # C₃₀
            [ 0.0,  0.0, -0.3*s5,  0.0,  0.0,  0.0,  0.0, -0.3*s5,  0.0,  1.0],
            # C₃₁
            [-0.25*s6, 0.0, 0.0, -0.05*s30, 0.0, 0.2*s30, 0.0,  0.0,  0.0,  0.0],
            # S₃₁
            [ 0.0, -0.05*s30, 0.0, 0.0, 0.0,  0.0, -0.25*s6, 0.0, 0.2*s30, 0.0],
            # C₃₂
            [ 0.0,  0.0, 0.5*s3,  0.0,  0.0,  0.0,  0.0, -0.5*s3,  0.0,  0.0],
            # S₃₂
            [ 0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0],
            # C₃₃
            [ 0.25*s10, 0.0, 0.0, -0.75*s2, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0],
            # S₃₃
            [ 0.0, 0.75*s2, 0.0, 0.0, 0.0,  0.0, -0.25*s10, 0.0, 0.0,  0.0],
        ], dtype=np.float64)
        return Th[:, perm]

    else:
        raise ValueError(f"No transformation matrix for L={lbl}")


def _symbols_from_atomnos(atomnos):
    """Convert atomic numbers to element symbols."""
    _ATOMIC_SYMBOLS = {
        1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O',
        9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
        16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti',
        23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu',
        30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr',
        37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc',
        44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
        51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba',
        78: 'Pt', 79: 'Au', 80: 'Hg', 82: 'Pb',
    }
    return [_ATOMIC_SYMBOLS.get(int(z), f'Z{z}') for z in atomnos]


# ---------------------------------------------------------------------------
# NTO (Natural Transition Orbital) computation from fchk transition densities
# ---------------------------------------------------------------------------

def compute_ntos_from_fchk(filepath, state_idx, canon_wfn, homo_idx, nbasis,
                           mocoeffs_sph=None, T_sph_to_cart=None):
    """Compute NTO hole/particle wavefunctions for an excited state.

    Uses the "G to E trans densities" and "Orthonormal basis" from .fchk.
    The orthonormal basis X (nbasis × nindep) provides the metric to correctly
    transform the AO transition density to the MO basis without needing the
    explicit AO overlap matrix S.

    Parameters:
        filepath: path to .fchk file
        state_idx: 0-based excited state index
        canon_wfn: Wavefunction with canonical MO coefficients (Cartesian basis)
        homo_idx: 0-based HOMO index
        nbasis: number of Cartesian basis functions (for output verification)
        mocoeffs_sph: (nmo, nsph) array — original spherical MO coefficients
        T_sph_to_cart: (nsph, ncart) array — transformation matrix

    Returns:
        hole_wfn: Wavefunction (hole NTOs in Cartesian AO basis)
        part_wfn: Wavefunction (particle NTOs in Cartesian AO basis)
        eigenvalues: (nocc,) array of NTO amplitudes (Σ values)
    """
    with open(filepath, 'r') as f:
        fchk_text = f.read()

    # Read transition density matrices (full nbf×nbf, stored contiguously)
    g2e_all = _read_fchk_array(fchk_text, 'G to E trans densities', float)
    if g2e_all is None:
        raise ValueError("No 'G to E trans densities' found in .fchk file")

    # Read orthonormal basis X (nbasis × nindep)
    X_flat = _read_fchk_array(fchk_text, 'Orthonormal basis', float)
    if X_flat is None:
        raise ValueError("No 'Orthonormal basis' found in .fchk file")

    # Transition density is in the spherical basis (nsph × nsph)
    if mocoeffs_sph is None:
        raise ValueError("Spherical MO coefficients required for NTO computation")
    nsph = mocoeffs_sph.shape[1]  # number of spherical basis functions
    nmo = mocoeffs_sph.shape[0]   # number of MOs (independent functions)
    nindep = nmo  # for this file, nindep = nmo = 710

    full_size = nsph * nsph
    nmat = len(g2e_all) // full_size

    # Each state has 2 matrices (alpha, beta). Use alpha component (even index).
    mat_idx = state_idx * 2
    if mat_idx >= nmat:
        raise ValueError(f"State {state_idx} not found (only {nmat//2} states)")

    T_flat = g2e_all[mat_idx * full_size:(mat_idx + 1) * full_size]
    T_ao = T_flat.reshape(nsph, nsph).copy()  # (713, 713)

    # Orthonormal basis X: (nsph × nindep), stored as flat array of nindep*nsph
    # cclib/fchk store as row-major: first nsph values = column 0, etc.
    # Actually the fchk stores X flattened column-wise or row-wise?
    # Let's reshape to (nindep, nsph) then transpose if needed
    X_raw = X_flat.reshape(nindep, nsph)  # (710, 713) in C order
    X = X_raw.T.copy()  # (713, 710) — columns are orthonormal basis vectors

    # Verify X is orthonormal: X^T S X = I (but we can't check without S)
    # Instead check that X has full column rank
    XtX = X.T @ X  # (710, 710)
    XtX_inv = np.linalg.inv(XtX)

    # MO coefficients in the original AO basis: C_ao (nsph × nmo)
    C_ao = np.asarray(mocoeffs_sph, dtype=np.float64).T  # (713, 710)

    # C_ao = X @ C_ortho, so C_ortho = X^+ @ C_ao = (X^T X)^{-1} X^T @ C_ao
    C_ortho = XtX_inv @ X.T @ C_ao  # (710, 710)

    # Verify C_ortho is unitary (C_ortho^T C_ortho = I)
    err = np.max(np.abs(C_ortho.T @ C_ortho - np.eye(nmo)))
    if err > 1e-8:
        import warnings
        warnings.warn(f"C_ortho not perfectly unitary: max|C^TC - I| = {err:.2e}")

    # Transform transition density to orthonormal basis (with metric correction)
    # The correct transformation for a density-like matrix from AO to the
    # orthonormal basis (where X^T S X = I) includes the pseudoinverse:
    #   T_X = (X^T X)^{-1} X^T  T_ao  X (X^T X)^{-1}
    # Then transform to MO basis via the unitary C_ortho.
    XtX = X.T @ X  # (nindep, nindep)
    XtX_inv = np.linalg.inv(XtX)
    T_X = XtX_inv @ X.T @ T_ao @ X @ XtX_inv  # (710, 710) in orthonormal basis
    T_mo = C_ortho.T @ T_X @ C_ortho  # (710, 710) in MO basis

    nocc = homo_idx + 1
    nvirt = nmo - nocc

    # Extract occupied→virtual block
    T_ov = T_mo[:nocc, nocc:]  # (nocc, nvirt)

    # SVD
    U, sigma, Vt = np.linalg.svd(T_ov, full_matrices=False)
    # U: (nocc, nocc), sigma: (nocc,), Vt: (nocc, nvirt)

    # Build NTO coefficients in orthonormal basis, then back to AO
    C_occ_ortho = C_ortho[:, :nocc]  # (710, nocc)
    C_virt_ortho = C_ortho[:, nocc:]  # (710, nvirt)

    C_hole_ortho = C_occ_ortho @ U  # (710, nocc)
    C_part_ortho = C_virt_ortho @ Vt.T  # (710, nocc)

    # Transform back to spherical AO basis: C_ao = X @ C_ortho
    C_hole_sph = X @ C_hole_ortho  # (713, nocc)
    C_part_sph = X @ C_part_ortho  # (713, nocc)

    # Transform to Cartesian basis for evaluation
    if T_sph_to_cart is not None:
        # T_sph_to_cart: (nsph, ncart)
        C_hole_cart = T_sph_to_cart.T @ C_hole_sph  # (ncart, nocc)
        C_part_cart = T_sph_to_cart.T @ C_part_sph  # (ncart, nocc)
    else:
        C_hole_cart = C_hole_sph
        C_part_cart = C_part_sph

    # Create Wavefunction objects (nmo=nocc for both)
    hole_labels = [f'NTO-hole-{i+1}' for i in range(nocc)]
    part_labels = [f'NTO-particle-{i+1}' for i in range(nocc)]
    zero_energies = np.zeros(nocc, dtype=np.float64)

    hole_wfn = Wavefunction(C_hole_cart.T.copy(), zero_energies, hole_labels)
    part_wfn = Wavefunction(C_part_cart.T.copy(), zero_energies, part_labels)

    return hole_wfn, part_wfn, sigma


def get_nto_state_count(filepath):
    """Return the number of excited states with NTO data in the .fchk file."""
    with open(filepath, 'r') as f:
        text = f.read()
    nex = _read_fchk_array(text, 'Number of excited states', int)
    if nex is None:
        return 0
    return int(nex[0])


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
    
    Returns:
        vertices: (N, 3) array or None
        faces: (M, 3) array or None
    """
    from skimage import measure
    
    try:
        result = measure.marching_cubes(
            grid_values, level=isovalue, spacing=(spacing, spacing, spacing)
        )
        # skimage >=0.23 returns (verts, faces, normals, values)
        verts, faces = result[0], result[1]
        # Shift vertices to world coordinates
        verts = verts + origin
        return verts, faces
    except (ValueError, RuntimeError):
        # No surface at this isovalue
        return None, None


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
# Gaussian Cube file export (for Blender and other tools)
# ---------------------------------------------------------------------------

ATOMIC_NUMBERS = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
    'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,
    'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22,
    'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29,
    'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
    'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43,
    'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
    'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54,
}


def write_cube_file(path, atoms, grid_values, origin, spacing, mo_idx, energy, label):
    """Write a Gaussian cube file containing MO values on a grid.
    
    Parameters
    ----------
    path : Path or str — output .cube file path
    atoms : list of Atom namedtuples
    grid_values : (nx, ny, nz) numpy array — MO values on grid
    origin : (3,) array — grid lower-left-front corner (Å)
    spacing : float — grid step (Å)
    mo_idx : int — 0-based orbital index
    energy : float or None — orbital energy in Eh
    label : str — orbital type label ('canonical', 'localized', etc.)
    
    Format: Gaussian cube (Å, indexing: x outer, y middle, z inner)
    """
    nx, ny, nz = grid_values.shape
    
    with open(path, 'w') as f:
        # Comment lines
        energy_str = f"{energy:+.6f} Eh" if energy is not None else "N/A"
        f.write(f"  MO {mo_idx + 1}  {label}  {energy_str}  (units: Angstrom)\n")
        f.write(f"  Generated by Orbital Visualizer for Blender import\n")
        
        # Grid dimensions and origin
        f.write(f"{len(atoms):5d}{origin[0]:12.6f}{origin[1]:12.6f}{origin[2]:12.6f}\n")
        f.write(f"{nx:5d}{spacing:12.6f}    0.000000    0.000000\n")
        f.write(f"{ny:5d}    0.000000{spacing:12.6f}    0.000000\n")
        f.write(f"{nz:5d}    0.000000    0.000000{spacing:12.6f}\n")
        
        # Atom positions
        for atom in atoms:
            an = atom.atomic_number
            if an <= 0:
                an = ATOMIC_NUMBERS.get(atom.symbol, 6)
            f.write(f"{an:5d}    0.000000{atom.x:12.6f}{atom.y:12.6f}{atom.z:12.6f}\n")
        
        # Grid values — 6 per line, x outer loop, y middle, z inner
        count = 0
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    f.write(f"{grid_values[ix, iy, iz]:13.5E}")
                    count += 1
                    if count % 6 == 0:
                        f.write("\n")
        if count % 6 != 0:
            f.write("\n")


import json

def write_render_recipe(path, logpath, entries):
    """Write a render recipe JSON for Blender import.
    
    Parameters
    ----------
    path : Path or str — output .json file path
    logpath : Path or str — absolute path to the source GAMESS log
    entries : list of dict with keys:
        cube_file, mo_idx, wtype, energy, label, isovalue, grid_spacing
    """
    recipe = {
        "version": 1,
        "source_log": str(Path(logpath).resolve()),
        "orbitals": []
    }
    
    for entry in entries:
        recipe["orbitals"].append({
            "cube_file": entry["cube_file"],
            "mo_idx": entry["mo_idx"],
            "wtype": entry["wtype"],
            "energy": entry["energy"],
            "label": entry["label"],
            "isovalue": entry["isovalue"],
            "grid_spacing": entry.get("grid_spacing", 0.08),
        })
    
    with open(path, 'w') as f:
        json.dump(recipe, f, indent=2)


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


# ---------------------------------------------------------------------------
# Qt imports (needed by all GUI classes below)
# ---------------------------------------------------------------------------

from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QSplitter, QListWidget, QListWidgetItem, QTabWidget,
    QLabel, QSlider, QPushButton, QFileDialog, QStatusBar,
    QLineEdit,
    QMessageBox, QApplication,
    QDialog, QTreeWidget, QTreeWidgetItem, QAbstractItemView,
    QComboBox, QDoubleSpinBox
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal
from PyQt6.QtGui import QAction


class GridWorker(QThread):
    """Background thread for MO grid computation."""
    finished = pyqtSignal(object, object, float, int)  # grid_values, origin, spacing, mo_idx
    
    def __init__(self, session, mo_coeffs, grid_spacing, mo_idx, parent=None):
        super().__init__(parent)
        self.session = session
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
            self.session.atoms, self.session.basis_set, self.mo_coeffs,
            grid_spacing=self.grid_spacing
        )
        if not self._cancelled:
            self.finished.emit(grid_values, origin, spacing, self.mo_idx)


# ---------------------------------------------------------------------------
# Molecule session (one per loaded file)
# ---------------------------------------------------------------------------

class MoleculeSession:
    """Holds all data for one loaded molecule. Owns its basis cache."""
    __slots__ = ('atoms', 'basis_set', 'canon_wfn', 'local_wfn',
                 'nto_hole_wfns', 'nto_part_wfns', 'nto_sigmas',
                 'homo_idx', 'filepath', 'bonds', 'nto_state_count',
                 'mocoeffs_sph', 'T_sph_to_cart')
    
    def __init__(self, atoms, basis_set, canon_wfn, local_wfn, homo_idx, filepath,
                 mocoeffs_sph=None, T_sph_to_cart=None):
        self.atoms = atoms
        self.basis_set = basis_set
        self.canon_wfn = canon_wfn
        self.local_wfn = local_wfn
        self.homo_idx = homo_idx
        self.filepath = Path(filepath)
        self.bonds = detect_bonds(atoms)
        self.nto_hole_wfns = None   # list of Wavefunction per state (lazy)
        self.nto_part_wfns = None   # list of Wavefunction per state (lazy)
        self.nto_sigmas = None      # list of eigenvalue arrays per state
        self.nto_state_count = 0
        self.mocoeffs_sph = mocoeffs_sph    # original spherical MO coefficients (nmo, nsph)
        self.T_sph_to_cart = T_sph_to_cart  # (nsph, ncart) transformation matrix
    
    def ensure_ntos_loaded(self):
        """Lazy-load NTOs from .fchk file if available."""
        if self.nto_hole_wfns is not None:
            return
        ext = self.filepath.suffix.lower()
        if ext != '.fchk':
            self.nto_state_count = 0
            self.nto_hole_wfns = []
            self.nto_part_wfns = []
            self.nto_sigmas = []
            return
        try:
            nstates = get_nto_state_count(str(self.filepath))
        except Exception:
            nstates = 0
        self.nto_state_count = nstates
        self.nto_hole_wfns = []
        self.nto_part_wfns = []
        self.nto_sigmas = []
        if nstates > 0 and self.canon_wfn is not None:
            for s in range(nstates):
                try:
                    hwf, pwf, sig = compute_ntos_from_fchk(
                        str(self.filepath), s, self.canon_wfn,
                        self.homo_idx, self.basis_set.nbasis,
                        self.mocoeffs_sph, self.T_sph_to_cart)
                    self.nto_hole_wfns.append(hwf)
                    self.nto_part_wfns.append(pwf)
                    self.nto_sigmas.append(sig)
                except Exception:
                    # If one state fails, stop loading more
                    self.nto_state_count = len(self.nto_hole_wfns)
                    break


# ---------------------------------------------------------------------------
# Viewport widget (one OrbitalCanvas + its orbital state)
# ---------------------------------------------------------------------------

class ViewportWidget(QWidget):
    """One 3D viewport showing an orbital. Manages its own grid computation."""
    clicked = pyqtSignal(object)
    orbital_changed = pyqtSignal()
    
    def __init__(self, session, parent=None):
        super().__init__(parent)
        self.session = session
        self.mo_idx = -1
        self.wtype = 'canonical'
        self._grid_worker = None
        self._generation = 0
        self._current_grid_values = None
        self._current_origin = None
        self._current_spacing = None
        self._custom_wfn = None
        self._active = False
        self._build_ui()
    
    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.label = QLabel("No orbital")
        self.label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.label.setStyleSheet("background: #222; color: #aaa; padding: 2px;")
        self.label.mousePressEvent = lambda e: self.clicked.emit(self)
        layout.addWidget(self.label)
        self.canvas = OrbitalCanvas()
        layout.addWidget(self.canvas.native_widget, 1)
        if self.session is not None:
            self.canvas.add_atoms_and_bonds(self.session.atoms, self.session.bonds)
            self.canvas.set_camera_center(self.session.atoms)
    
    @property
    def active(self):
        return self._active
    
    @active.setter
    def active(self, val):
        self._active = val
        self.setStyleSheet("border: 2px solid #4a9eff;" if val else "border: 2px solid #333;")
    
    @property
    def current_wfn(self):
        if hasattr(self, '_custom_wfn') and self._custom_wfn is not None:
            return self._custom_wfn
        if self.wtype == 'canonical':
            return self.session.canon_wfn
        return self.session.local_wfn
    
    def set_orbital(self, mo_idx, wtype, isovalue, grid_spacing):
        # Use current_wfn property (handles canonical, localized, and NTO custom)
        wfn = self.current_wfn
        if wfn is None or mo_idx < 0 or mo_idx >= wfn.nmo:
            return
        self.mo_idx = mo_idx
        self.wtype = wtype
        self._cancel_computation()
        self._generation += 1
        gen = self._generation
        label = f"MO {mo_idx + 1} ({wfn.labels[mo_idx]})"
        if wtype == 'canonical' and mo_idx < len(wfn.energies):
            label += f"  —  {wfn.energies[mo_idx]:+.4f} Eh"
        elif 'nto' in wtype:
            label = f"{wfn.labels[mo_idx]}"
        self.label.setText(label)
        self.label.setStyleSheet("background: #222; color: #fff; padding: 2px;")
        self._start_compute(wfn.get_mo(mo_idx), grid_spacing, gen)
    
    def _start_compute(self, mo_coeffs, spacing, generation):
        self._grid_worker = GridWorker(self.session, mo_coeffs, spacing, self.mo_idx)
        self._grid_worker.finished.connect(
            lambda gv, o, s, mi: self._on_grid_done(gv, o, s, mi, generation))
        self._grid_worker.finished.connect(lambda *a: self._cleanup_worker(self._grid_worker))
        self._grid_worker.start()
    
    def _on_grid_done(self, grid_values, origin, spacing, mo_idx, generation):
        if generation != self._generation or mo_idx != self.mo_idx:
            return
        self._current_grid_values = grid_values
        self._current_origin = origin
        self._current_spacing = spacing
        self.orbital_changed.emit()
    
    def update_surface(self, isovalue):
        if self._current_grid_values is None:
            return
        vp, fp = extract_isosurface(self._current_grid_values, +isovalue,
                                        self._current_origin, self._current_spacing)
        vn, fn = extract_isosurface(self._current_grid_values, -isovalue,
                                        self._current_origin, self._current_spacing)
        self.canvas.set_orbital_surface(vp, fp, vn, fn)
    
    def _cancel_computation(self):
        if self._grid_worker is not None:
            if self._grid_worker.isRunning():
                self._grid_worker.cancel()
                self._grid_worker.wait(3000)
            self._grid_worker = None
    
    def _cleanup_worker(self, worker):
        if self._grid_worker is worker:
            worker.wait(500)
            if self._grid_worker is worker:
                self._grid_worker = None
    
    def shutdown(self):
        self._cancel_computation()
        self.canvas.canvas.close()
    
    def get_overlay_info(self):
        lines = []
        if self.mo_idx >= 0 and self.current_wfn is not None:
            wfn = self.current_wfn
            lines.append(f"Orbital {self.mo_idx + 1} ({wfn.labels[self.mo_idx]})")
            if self.wtype == 'canonical' and self.mo_idx < len(wfn.energies):
                e = wfn.energies[self.mo_idx]
                lines.append(f"Energy: {e:+.6f} Eh  ({e * 27.2114:+.2f} eV)")
        return lines


# ---------------------------------------------------------------------------
# Main GUI window
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Molecule tab (gallery + viewports + controls for one molecule)
# ---------------------------------------------------------------------------

class MoleculeTab(QWidget):
    """One tab per loaded molecule. Contains gallery, N viewports, control bar."""
    REFINEMENT_STEPS = [0.35, 0.20, 0.10]
    
    def __init__(self, session, parent=None):
        super().__init__(parent)
        self.session = session
        self._viewports = []
        self._active_viewport = None
        self._isovalue = 0.05
        self._grid_spacing = 0.35
        self._refining = {}
        self._refine_generations = {}
        
        self._build_ui()
        self._populate_gallery()
        self._set_viewport_count(1)
    
    def _build_ui(self):
        main_layout = QHBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        
        gallery_panel = QWidget()
        gallery_layout = QVBoxLayout(gallery_panel)
        gallery_layout.setContentsMargins(2, 2, 2, 2)
        
        self.gallery_tabs = QTabWidget()
        self.occupied_list = QListWidget()
        self.occupied_list.currentRowChanged.connect(self._on_occupied_selected)
        self.gallery_tabs.addTab(self.occupied_list, "Occupied")
        self.virtual_list = QListWidget()
        self.virtual_list.currentRowChanged.connect(self._on_virtual_selected)
        self.gallery_tabs.addTab(self.virtual_list, "Virtual")
        self.localized_list = QListWidget()
        self.localized_list.currentRowChanged.connect(self._on_localized_selected)
        gallery_layout.addWidget(self.gallery_tabs)
        
        # NTO gallery (shown only for .fchk files with excited states)
        self.nto_panel = QWidget()
        nto_layout = QVBoxLayout(self.nto_panel)
        nto_layout.setContentsMargins(0, 0, 0, 0)
        nto_state_layout = QHBoxLayout()
        nto_state_layout.addWidget(QLabel("State:"))
        self.nto_state_combo = QComboBox()
        self.nto_state_combo.currentIndexChanged.connect(self._on_nto_state_changed)
        nto_state_layout.addWidget(self.nto_state_combo)
        nto_layout.addLayout(nto_state_layout)
        self.nto_hole_list = QListWidget()
        self.nto_hole_list.currentRowChanged.connect(self._on_nto_hole_selected)
        nto_layout.addWidget(QLabel("Hole NTOs (occupied→virtual):"))
        nto_layout.addWidget(self.nto_hole_list)
        self.nto_part_list = QListWidget()
        self.nto_part_list.currentRowChanged.connect(self._on_nto_part_selected)
        nto_layout.addWidget(QLabel("Particle NTOs (virtual→occupied):"))
        nto_layout.addWidget(self.nto_part_list)
        self.nto_panel.setVisible(False)
        gallery_layout.addWidget(self.nto_panel)
        
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        right_layout.setContentsMargins(0, 0, 0, 0)
        
        self.viewport_container = QWidget()
        self.viewport_layout = QGridLayout(self.viewport_container)
        self.viewport_layout.setContentsMargins(0, 0, 0, 0)
        self.viewport_layout.setSpacing(2)
        right_layout.addWidget(self.viewport_container, 1)
        
        ctrl = QWidget()
        ctrl_layout = QHBoxLayout(ctrl)
        ctrl_layout.setContentsMargins(4, 2, 4, 2)
        
        ctrl_layout.addWidget(QLabel("Isovalue:"))
        self.isovalue_slider = QSlider(Qt.Orientation.Horizontal)
        self.isovalue_slider.setRange(10, 200)
        self.isovalue_slider.setValue(50)
        self.isovalue_slider.valueChanged.connect(self._on_isovalue_changed)
        ctrl_layout.addWidget(self.isovalue_slider)
        self.isovalue_edit = QLineEdit("0.050")
        self.isovalue_edit.setFixedWidth(50)
        self.isovalue_edit.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.isovalue_edit.setStyleSheet("QLineEdit { border: 1px solid #555; background: #1a1a1a; color: #fff; padding: 1px; }")
        self.isovalue_edit.editingFinished.connect(self._on_isovalue_edited)
        self.isovalue_edit.returnPressed.connect(self._on_isovalue_edited)
        ctrl_layout.addWidget(self.isovalue_edit)
        
        ctrl_layout.addSpacing(10)
        ctrl_layout.addWidget(QLabel("Grid:"))
        self.grid_slider = QSlider(Qt.Orientation.Horizontal)
        self.grid_slider.setRange(5, 40)
        self.grid_slider.setValue(35)
        self.grid_slider.valueChanged.connect(self._on_grid_changed)
        ctrl_layout.addWidget(self.grid_slider)
        self.grid_edit = QLineEdit("0.35")
        self.grid_edit.setFixedWidth(50)
        self.grid_edit.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.grid_edit.setStyleSheet("QLineEdit { border: 1px solid #555; background: #1a1a1a; color: #fff; padding: 1px; }")
        self.grid_edit.editingFinished.connect(self._on_grid_edited)
        self.grid_edit.returnPressed.connect(self._on_grid_edited)
        ctrl_layout.addWidget(self.grid_edit)
        
        for label, sp, sv in [("Coarse", 0.35, 35), ("Medium", 0.20, 20), ("Fine", 0.10, 10)]:
            btn = QPushButton(label)
            btn.clicked.connect(lambda checked, s=sp, v=sv: self._set_grid_preset(s, v))
            ctrl_layout.addWidget(btn)
        
        ctrl_layout.addStretch()
        
        self.split_btn = QPushButton("Split ▸")
        self.split_btn.clicked.connect(self._toggle_split)
        ctrl_layout.addWidget(self.split_btn)
        
        self.export_blender_btn = QPushButton("Export Blender")
        self.export_blender_btn.clicked.connect(self._export_to_blender)
        self.export_blender_btn.setToolTip("Export active viewport's orbital as Gaussian cube for Blender")
        ctrl_layout.addWidget(self.export_blender_btn)
        
        self.homo_btn = QPushButton("HOMO")
        self.homo_btn.clicked.connect(self._go_to_homo)
        ctrl_layout.addWidget(self.homo_btn)
        self.lumo_btn = QPushButton("LUMO")
        self.lumo_btn.clicked.connect(self._go_to_lumo)
        ctrl_layout.addWidget(self.lumo_btn)
        
        right_layout.addWidget(ctrl)
        
        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(gallery_panel)
        splitter.addWidget(right_panel)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([220, 980])
        main_layout.addWidget(splitter)
    
    def _populate_gallery(self):
        wfn = self.session.canon_wfn
        if wfn is None:
            return
        self.occupied_list.clear()
        for i in range(self.session.homo_idx + 1):
            e = wfn.energies[i]
            item = QListWidgetItem(f"MO {i+1}  ({e:+.3f} Eh)")
            item.setData(Qt.ItemDataRole.UserRole, i)
            self.occupied_list.addItem(item)
        self.virtual_list.clear()
        for i in range(self.session.homo_idx + 1, wfn.nmo):
            e = wfn.energies[i]
            item = QListWidgetItem(f"MO {i+1}  ({e:+.3f} Eh)")
            item.setData(Qt.ItemDataRole.UserRole, i)
            self.virtual_list.addItem(item)
        if self.session.local_wfn is not None:
            self.localized_list.clear()
            for i in range(self.session.local_wfn.nmo):
                item = QListWidgetItem(f"Loc {i+1}")
                item.setData(Qt.ItemDataRole.UserRole, i)
                self.localized_list.addItem(item)
            self.gallery_tabs.addTab(self.localized_list, "Localized")
        
        # NTO setup (lazy-loaded on first tab access or via state combo)
        if self.session.filepath.suffix.lower() == '.fchk':
            self._setup_nto_panel()
    
    def _set_viewport_count(self, n):
        for vp in self._viewports:
            vp.shutdown()
            vp.setParent(None)
        self._viewports = []
        self._refining = {}
        self._refine_generations = {}
        
        while self.viewport_layout.count():
            item = self.viewport_layout.takeAt(0)
            if item.widget():
                item.widget().setParent(None)
        
        for i in range(n):
            vp = ViewportWidget(self.session)
            vp.clicked.connect(lambda v=vp: self._activate_viewport(v))
            vp.orbital_changed.connect(lambda v=vp: self._on_viewport_updated(v))
            self._viewports.append(vp)
        
        if n == 1:
            self.viewport_layout.addWidget(self._viewports[0], 0, 0)
        elif n == 2:
            self.viewport_layout.addWidget(self._viewports[0], 0, 0)
            self.viewport_layout.addWidget(self._viewports[1], 0, 1)
        
        self._activate_viewport(self._viewports[-1])
        self.split_btn.setText("Split ▸" if n == 1 else "Unsplit ◂")
    
    def _toggle_split(self):
        n = 2 if len(self._viewports) == 1 else 1
        self._set_viewport_count(n)
    
    def _activate_viewport(self, vp):
        for v in self._viewports:
            v.active = (v is vp)
        self._active_viewport = vp
    
    def _on_viewport_updated(self, vp):
        vp.update_surface(self._isovalue)
        step = self._refining.get(id(vp), 0)
        if step < len(self.REFINEMENT_STEPS):
            next_sp = max(self.REFINEMENT_STEPS[step], self._grid_spacing)
            if vp._current_spacing is not None and next_sp < vp._current_spacing:
                self._refining[id(vp)] = step + 1
                gen = self._refine_generations.get(id(vp), 0) + 1
                self._refine_generations[id(vp)] = gen
                vp._start_compute(vp.current_wfn.get_mo(vp.mo_idx), next_sp, gen)
                return
        self._refining[id(vp)] = 999
    
    def _assign_orbital_to_active(self, mo_idx, wtype='canonical'):
        if self._active_viewport is None:
            return
        # Clear any NTO custom wavefunction when switching to canonical/localized
        self._active_viewport._custom_wfn = None
        self._refining[id(self._active_viewport)] = 0
        self._refine_generations[id(self._active_viewport)] = 0
        self._active_viewport.set_orbital(mo_idx, wtype, self._isovalue, self.REFINEMENT_STEPS[0])
    
    def _on_occupied_selected(self, row):
        if row < 0: return
        self._assign_orbital_to_active(self.occupied_list.item(row).data(Qt.ItemDataRole.UserRole), 'canonical')
    
    def _on_virtual_selected(self, row):
        if row < 0: return
        self._assign_orbital_to_active(self.virtual_list.item(row).data(Qt.ItemDataRole.UserRole), 'canonical')
    
    def _on_localized_selected(self, row):
        if row < 0: return
        self._assign_orbital_to_active(self.localized_list.item(row).data(Qt.ItemDataRole.UserRole), 'localized')
    
    # --- NTO handlers ---
    
    def _setup_nto_panel(self):
        """Set up NTO panel: load NTOs lazily and populate state combo."""
        self.session.ensure_ntos_loaded()
        nc = self.session.nto_state_count
        if nc > 0:
            self.nto_state_combo.blockSignals(True)
            self.nto_state_combo.clear()
            for s in range(nc):
                self.nto_state_combo.addItem(f"State {s+1}")
            self.nto_state_combo.blockSignals(False)
            self.nto_panel.setVisible(True)
            self._on_nto_state_changed(0)
    
    def _on_nto_state_changed(self, idx):
        if idx < 0 or idx >= self.session.nto_state_count:
            return
        self._populate_nto_lists(idx)
    
    def _populate_nto_lists(self, state_idx):
        hwfn = self.session.nto_hole_wfns[state_idx]
        pwfn = self.session.nto_part_wfns[state_idx]
        sigmas = self.session.nto_sigmas[state_idx]
        
        self.nto_hole_list.clear()
        for i in range(hwfn.nmo):
            lam = sigmas[i]
            pct = lam*lam / np.sum(sigmas**2) * 100
            item = QListWidgetItem(f"Hole {i+1}  λ={lam:.3f} ({pct:.1f}%)")
            item.setData(Qt.ItemDataRole.UserRole, i)
            self.nto_hole_list.addItem(item)
        
        self.nto_part_list.clear()
        for i in range(pwfn.nmo):
            lam = sigmas[i]
            pct = lam*lam / np.sum(sigmas**2) * 100
            item = QListWidgetItem(f"Part {i+1}  λ={lam:.3f} ({pct:.1f}%)")
            item.setData(Qt.ItemDataRole.UserRole, i)
            self.nto_part_list.addItem(item)
    
    def _on_nto_hole_selected(self, row):
        if row < 0: return
        state_idx = self.nto_state_combo.currentIndex()
        if state_idx < 0:
            return
        self._assign_nto_orbital_to_active(state_idx, 'hole', row)
    
    def _on_nto_part_selected(self, row):
        if row < 0: return
        state_idx = self.nto_state_combo.currentIndex()
        if state_idx < 0:
            return
        self._assign_nto_orbital_to_active(state_idx, 'particle', row)
    
    def _assign_nto_orbital_to_active(self, state_idx, nto_type, orb_idx):
        """Assign an NTO to the active viewport."""
        if self._active_viewport is None:
            return
        wfn = (self.session.nto_hole_wfns[state_idx] if nto_type == 'hole'
               else self.session.nto_part_wfns[state_idx])
        if orb_idx < 0 or orb_idx >= wfn.nmo:
            return
        self._refining[id(self._active_viewport)] = 0
        self._refine_generations[id(self._active_viewport)] = 0
        # Use a custom wtype to distinguish NTO from canonical
        wtype = f'nto_{nto_type}_{state_idx}'
        # Store the NTO wavefunction on the viewport so current_wfn finds it
        self._active_viewport._custom_wfn = wfn
        self._active_viewport.set_orbital(orb_idx, wtype, self._isovalue, self.REFINEMENT_STEPS[0])
    
    def _on_isovalue_changed(self, value):
        self._isovalue = value / 1000.0
        self.isovalue_edit.blockSignals(True)
        self.isovalue_edit.setText(f"{self._isovalue:.3f}")
        self.isovalue_edit.blockSignals(False)
        for vp in self._viewports:
            if vp._current_grid_values is not None:
                vp.update_surface(self._isovalue)
    
    def _on_isovalue_edited(self):
        """Handle manual isovalue input."""
        try:
            val = float(self.isovalue_edit.text())
            val = max(0.001, min(0.200, val))
            # Convert to slider value: slider = val * 1000
            slider_val = int(round(val * 1000))
            slider_val = max(10, min(200, slider_val))
            if slider_val != self.isovalue_slider.value():
                self.isovalue_slider.setValue(slider_val)
            else:
                # Force update display even if slider didn't move
                self._isovalue = slider_val / 1000.0
                self.isovalue_edit.setText(f"{self._isovalue:.3f}")
        except ValueError:
            self.isovalue_edit.setText(f"{self._isovalue:.3f}")
    
    def _on_grid_changed(self, value):
        self._grid_spacing = value / 100.0
        self.grid_edit.blockSignals(True)
        self.grid_edit.setText(f"{self._grid_spacing:.2f}")
        self.grid_edit.blockSignals(False)
        for vp in self._viewports:
            if vp.mo_idx >= 0:
                self._refining[id(vp)] = 0
                self._refine_generations[id(vp)] = 0
                vp._cancel_computation()
                vp._generation += 1
                vp._start_compute(vp.current_wfn.get_mo(vp.mo_idx), self._grid_spacing, vp._generation)
    
    def _on_grid_edited(self):
        """Handle manual grid spacing input."""
        try:
            val = float(self.grid_edit.text())
            val = max(0.05, min(0.40, val))
            slider_val = max(5, min(40, int(round(val * 100))))
            if slider_val != self.grid_slider.value():
                self.grid_slider.setValue(slider_val)
            else:
                self._grid_spacing = slider_val / 100.0
                self.grid_edit.setText(f"{self._grid_spacing:.2f}")
        except ValueError:
            self.grid_edit.setText(f"{self._grid_spacing:.2f}")
    
    def _set_grid_preset(self, spacing, slider_value):
        self.grid_slider.setValue(slider_value)
    
    def _go_to_homo(self):
        self.gallery_tabs.setCurrentIndex(0)
        self.occupied_list.setCurrentRow(self.session.homo_idx)
    
    def _go_to_lumo(self):
        self.gallery_tabs.setCurrentIndex(1)
        self.virtual_list.setCurrentRow(0)
    
    def _export_to_blender(self):
        """Export the active viewport's current orbital as a Gaussian cube + recipe JSON."""
        from PyQt6.QtWidgets import QFileDialog
        
        vp = self._active_viewport
        if vp is None or vp.mo_idx < 0 or vp._current_grid_values is None:
            QMessageBox.warning(self, "No Orbital", 
                                "Select an orbital first (click in gallery).")
            return
        
        dirpath = QFileDialog.getExistingDirectory(
            self, "Select export directory", "",
            QFileDialog.Option.ShowDirsOnly)
        if not dirpath:
            return
        
        dirpath = Path(dirpath)
        try:
            dirpath.mkdir(parents=True, exist_ok=True)
        except OSError as e:
            QMessageBox.critical(self, "Error", f"Cannot create directory:\n{e}")
            return
        
        self._write_blender_export(dirpath, vp)
    
    def _write_blender_export(self, dirpath, vp):
        """Write cube file and recipe for a single viewport's orbital."""
        wfn = vp.current_wfn
        energy = wfn.energies[vp.mo_idx] if vp.wtype == 'canonical' and vp.mo_idx < len(wfn.energies) else 0.0
        label = wfn.labels[vp.mo_idx] if vp.mo_idx < len(wfn.labels) else 'orbital'
        
        cube_name = f"mo_{vp.mo_idx + 1}_{vp.wtype}.cube"
        cube_path = dirpath / cube_name
        
        write_cube_file(cube_path, self.session.atoms, vp._current_grid_values,
                        vp._current_origin, vp._current_spacing,
                        vp.mo_idx, energy, label)
        
        recipe_path = dirpath / "render_recipe.json"
        write_render_recipe(recipe_path, self.session.filepath, [{
            "cube_file": cube_name,
            "mo_idx": vp.mo_idx,
            "wtype": vp.wtype,
            "energy": energy,
            "label": label,
            "isovalue": self._isovalue,
            "grid_spacing": vp._current_spacing,
        }])
        
        from PyQt6.QtWidgets import QApplication
        QApplication.instance().activeWindow().statusBar().showMessage(
            f"Exported {cube_name} to {dirpath}")
    
    def shutdown(self):
        for vp in self._viewports:
            vp.shutdown()

# ---------------------------------------------------------------------------
# Main GUI window (manages tabs)
# ---------------------------------------------------------------------------


class OrbitalViewer(QMainWindow):
    """Main application window with tabbed molecules and split viewports."""
    
    def __init__(self, logpath=None):
        super().__init__()
        self.setWindowTitle("Orbital Visualizer")
        self.resize(1200, 800)
        self._sessions = []
        
        self._build_ui()
        self._build_menu()
        
        if logpath is not None:
            self._open_file(Path(logpath))
    
    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)
        layout.setContentsMargins(4, 4, 4, 4)
        
        self.tab_widget = QTabWidget()
        self.tab_widget.setTabsClosable(True)
        self.tab_widget.tabCloseRequested.connect(self._close_tab)
        self.tab_widget.currentChanged.connect(self._on_tab_changed)
        layout.addWidget(self.tab_widget)
        
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready — Open a file (Ctrl+O): GAMESS .log / Gaussian .fchk")
    
    def _build_menu(self):
        mb = self.menuBar()
        fm = mb.addMenu("&File")
        a = QAction("&Open...", self); a.setShortcut("Ctrl+O"); a.triggered.connect(self._file_open); fm.addAction(a)
        a = QAction("&New Tab", self); a.setShortcut("Ctrl+T"); a.triggered.connect(lambda: self._file_open()); fm.addAction(a)
        fm.addSeparator()
        a = QAction("&Export Image...", self); a.setShortcut("Ctrl+E"); a.triggered.connect(self._file_export); fm.addAction(a)
        a = QAction("Export for &Blender...", self); a.setShortcut("Ctrl+B"); a.triggered.connect(self._export_blender_multi); fm.addAction(a)
        fm.addSeparator()
        a = QAction("&Close Tab", self); a.setShortcut("Ctrl+W"); a.triggered.connect(lambda: self._close_tab(self.tab_widget.currentIndex())); fm.addAction(a)
        a = QAction("&Quit", self); a.setShortcut("Ctrl+Q"); a.triggered.connect(self.close); fm.addAction(a)
    
    def _file_open(self):
        filters = (
            "Supported Files (*.log *.out *.fchk);;"
            "GAMESS Log (*.log *.out);;"
            "Gaussian FChk (*.fchk);;"
            "All Files (*)"
        )
        path, _ = QFileDialog.getOpenFileName(
            self, "Open Calculation File", "", filters)
        if path:
            self._open_file(Path(path))
    
    def _open_file(self, filepath):
        ext = filepath.suffix.lower()
        if ext == '.chk':
            QMessageBox.information(
                self, "Binary checkpoint",
                f"Binary .chk files cannot be parsed directly.\n\n"
                f"Convert to formatted checkpoint first:\n"
                f"  formchk {filepath.name} {filepath.stem}.fchk\n\n"
                f"Then open the .fchk file."
            )
            return
        
        self.status_bar.showMessage(f"Loading {filepath.name} ...")
        QApplication.processEvents()
        mocoeffs_sph = None
        T_s2c = None
        try:
            if ext in ('.log', '.out'):
                atoms, basis_set, canon_wfn, local_wfn, homo_idx = parse_gamess_log(filepath)
            elif ext == '.fchk':
                atoms, basis_set, canon_wfn, local_wfn, homo_idx, mocoeffs_sph, T_s2c = \
                    parse_gaussian_fchk(filepath)
            else:
                QMessageBox.warning(
                    self, "Unknown format",
                    f"Unrecognised file extension '{ext}'.\n"
                    f"Expected: .log, .out (GAMESS) or .fchk (Gaussian)."
                )
                self.status_bar.showMessage("Unrecognised file format")
                return
        except Exception as e:
            import traceback
            traceback.print_exc()
            QMessageBox.critical(self, "Error", f"Failed to parse file:\n{e}")
            self.status_bar.showMessage("Error loading file")
            return
        
        clear_basis_cache()
        session = MoleculeSession(atoms, basis_set, canon_wfn, local_wfn, homo_idx, filepath,
                                 mocoeffs_sph, T_s2c)
        self._sessions.append(session)
        
        tab = MoleculeTab(session)
        idx = self.tab_widget.addTab(tab, filepath.name)
        self.tab_widget.setCurrentIndex(idx)
        
        self.status_bar.showMessage(
            f"Loaded: {len(atoms)} atoms, {basis_set.nbasis} bf, "
            f"{canon_wfn.nmo} MOs, HOMO={homo_idx + 1}")
        
        # Pre-warm: auto-load HOMO so JIT compilation happens now,
        # not on first user click. GridWorker runs in background thread.
        tab._go_to_homo()
    
    def _close_tab(self, index):
        if index < 0 or index >= len(self._sessions):
            return
        tab = self.tab_widget.widget(index)
        if hasattr(tab, 'shutdown'):
            tab.shutdown()
        self.tab_widget.removeTab(index)
        del self._sessions[index]
    
    def _on_tab_changed(self, index):
        if index >= 0 and index < len(self._sessions):
            session = self._sessions[index]
            self.setWindowTitle(f"Orbital Visualizer — {session.filepath.name}")
    
    def _file_export(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Image", "orbital.png", "PNG Images (*.png);;All Files (*)")
        if path:
            self._export_with_overlay(path)
            self.status_bar.showMessage(f"Saved: {path}")
    
    def _export_with_overlay(self, path):
        from PyQt6.QtGui import QImage, QPainter, QColor, QFont
        
        tab = self.tab_widget.currentWidget()
        if tab is None or not hasattr(tab, '_active_viewport'):
            return
        vp = tab._active_viewport
        if vp is None:
            return
        
        img_array = vp.canvas.canvas.render()
        h, w = img_array.shape[:2]
        img_8bit = (np.clip(img_array, 0, 1) * 255).astype(np.uint8)
        qimg = QImage(img_8bit.data, w, h, w * 4, QImage.Format.Format_RGBA8888).copy()
        
        lines = vp.get_overlay_info()
        lines.append(f"Isovalue: ±{tab._isovalue:.3f}")
        if vp._current_spacing is not None:
            lines.append(f"Grid: {vp._current_spacing:.2f} Å")
        if self._sessions:
            s = self._sessions[self.tab_widget.currentIndex()]
            lines.insert(0, f"{s.filepath.name}")
        
        if not lines:
            qimg.save(path, 'PNG')
            return
        
        painter = QPainter(qimg)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        font = QFont("monospace", 12)
        font.setBold(True)
        painter.setFont(font)
        
        pad = 8
        lh = 20
        bh = len(lines) * lh + pad * 2
        metrics = painter.fontMetrics()
        bw = max(metrics.horizontalAdvance(l) for l in lines) + pad * 2
        
        painter.fillRect(0, 0, bw, bh, QColor(0, 0, 0, 180))
        painter.setPen(QColor(255, 255, 255, 255))
        y = pad + metrics.ascent()
        for line in lines:
            painter.drawText(pad, y, line)
            y += lh
        painter.end()
        qimg.save(path, 'PNG')
    
    def _export_blender_multi(self):
        """Open dialog to select multiple orbitals and export as cube files + recipe."""
        tab = self.tab_widget.currentWidget()
        if tab is None or not hasattr(tab, 'session'):
            QMessageBox.warning(self, "No Molecule", "Open a molecule first.")
            return
        
        session = tab.session
        wfn = session.canon_wfn
        if wfn is None:
            return
        
        dirpath = QFileDialog.getExistingDirectory(
            self, "Select export directory", "",
            QFileDialog.Option.ShowDirsOnly)
        if not dirpath:
            return
        dirpath = Path(dirpath)
        try:
            dirpath.mkdir(parents=True, exist_ok=True)
        except OSError as e:
            QMessageBox.critical(self, "Error", f"Cannot create directory:\n{e}")
            return
        
        # Build selection dialog
        dlg = QDialog(self)
        dlg.setWindowTitle("Select Orbitals for Blender Export")
        dlg.resize(500, 400)
        dlg_layout = QVBoxLayout(dlg)
        
        dlg_layout.addWidget(QLabel("Select orbitals to export:"))
        
        tree = QTreeWidget()
        tree.setHeaderLabels(["Orbital", "Energy (Eh)", "Type"])
        tree.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
        
        # Canonical orbitals
        canon_root = QTreeWidgetItem(tree, [f"Canonical ({wfn.nmo} MOs)", "", ""])
        canon_root.setFlags(canon_root.flags() & ~Qt.ItemFlag.ItemIsUserCheckable)
        for i in range(wfn.nmo):
            energy = wfn.energies[i] if i < len(wfn.energies) else 0.0
            item = QTreeWidgetItem(canon_root, [
                f"MO {i + 1} ({wfn.labels[i]})",
                f"{energy:+.4f}",
                "canonical"
            ])
            item.setData(0, Qt.ItemDataRole.UserRole, i)
            item.setData(2, Qt.ItemDataRole.UserRole + 1, "canonical")
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(0, Qt.CheckState.Unchecked)
        canon_root.setExpanded(True)
        
        # Localized orbitals
        if session.local_wfn is not None:
            local_wfn = session.local_wfn
            local_root = QTreeWidgetItem(tree, [f"Localized ({local_wfn.nmo} MOs)", "", ""])
            local_root.setFlags(local_root.flags() & ~Qt.ItemFlag.ItemIsUserCheckable)
            for i in range(local_wfn.nmo):
                item = QTreeWidgetItem(local_root, [
                    f"Loc {i + 1}",
                    "N/A",
                    "localized"
                ])
                item.setData(0, Qt.ItemDataRole.UserRole, i)
                item.setData(2, Qt.ItemDataRole.UserRole + 1, "localized")
                item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                item.setCheckState(0, Qt.CheckState.Unchecked)
            local_root.setExpanded(True)
        
        dlg_layout.addWidget(tree)
        
        # Grid spacing selector
        grid_layout = QHBoxLayout()
        grid_layout.addWidget(QLabel("Grid spacing:"))
        grid_combo = QComboBox()
        grid_combo.addItems(["0.10 (Fine)", "0.15 (Medium-fine)", "0.20 (Medium)", "0.35 (Coarse)"])
        grid_combo.setCurrentIndex(2)  # 0.20 default
        grid_layout.addWidget(grid_combo)
        grid_layout.addStretch()
        dlg_layout.addLayout(grid_layout)
        
        # Isovalue
        iso_layout = QHBoxLayout()
        iso_layout.addWidget(QLabel("Isovalue (for recipe only):"))
        iso_spin = QDoubleSpinBox()
        iso_spin.setRange(0.001, 0.5)
        iso_spin.setValue(0.05)
        iso_spin.setDecimals(3)
        iso_spin.setSingleStep(0.005)
        iso_layout.addWidget(iso_spin)
        iso_layout.addStretch()
        dlg_layout.addLayout(iso_layout)
        
        # Buttons
        btn_layout = QHBoxLayout()
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(dlg.reject)
        export_btn = QPushButton("Export Selected")
        export_btn.setDefault(True)
        btn_layout.addStretch()
        btn_layout.addWidget(cancel_btn)
        btn_layout.addWidget(export_btn)
        dlg_layout.addLayout(btn_layout)
        
        def do_export():
            spacing_str = grid_combo.currentText()
            spacing = float(spacing_str.split()[0])
            isovalue = iso_spin.value()
            
            entries = []
            # Collect checked items
            for i in range(canon_root.childCount()):
                child = canon_root.child(i)
                if child.checkState(0) == Qt.CheckState.Checked:
                    mo_idx = child.data(0, Qt.ItemDataRole.UserRole)
                    energy = wfn.energies[mo_idx] if mo_idx < len(wfn.energies) else 0.0
                    entries.append((mo_idx, "canonical", energy, wfn.labels[mo_idx]))
            
            if hasattr(session, 'local_wfn') and session.local_wfn is not None:
                for i in range(local_root.childCount()):
                    child = local_root.child(i)
                    if child.checkState(0) == Qt.CheckState.Checked:
                        mo_idx = child.data(0, Qt.ItemDataRole.UserRole)
                        entries.append((mo_idx, "localized", 0.0, "localized"))
            
            if not entries:
                QMessageBox.warning(dlg, "No Selection", "Check at least one orbital.")
                return
            
            dlg.accept()
            
            self.status_bar.showMessage(f"Exporting {len(entries)} orbitals...")
            QApplication.processEvents()
            
            recipe_entries = []
            for mo_idx, wtype, energy, label in entries:
                target_wfn = session.canon_wfn if wtype == "canonical" else session.local_wfn
                mo_coeffs = target_wfn.get_mo(mo_idx)
                
                grid_values, origin, grid_sp = eval_mo_on_grid(
                    session.atoms, session.basis_set, mo_coeffs,
                    grid_spacing=spacing, padding=5.0
                )
                
                cube_name = f"mo_{mo_idx + 1}_{wtype}.cube"
                cube_path = dirpath / cube_name
                write_cube_file(cube_path, session.atoms, grid_values, origin,
                               grid_sp, mo_idx, energy, label)
                
                recipe_entries.append({
                    "cube_file": cube_name,
                    "mo_idx": mo_idx,
                    "wtype": wtype,
                    "energy": energy,
                    "label": label,
                    "isovalue": isovalue,
                    "grid_spacing": grid_sp,
                })
            
            write_render_recipe(dirpath / "render_recipe.json", session.filepath, recipe_entries)
            self.status_bar.showMessage(
                f"Exported {len(entries)} orbitals to {dirpath}")
        
        export_btn.clicked.connect(do_export)
        dlg.exec()
    
    def closeEvent(self, event):
        for i in range(self.tab_widget.count()):
            tab = self.tab_widget.widget(i)
            if hasattr(tab, 'shutdown'):
                tab.shutdown()
        super().closeEvent(event)



# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def cli_main():
    """Command-line mode for headless rendering."""
    parser = argparse.ArgumentParser(description='Orbital Visualizer')
    parser.add_argument('logfile', nargs='?', default='spval.log',
                        help='GAMESS .log / .out or Gaussian .fchk file to visualize')
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
    parser.add_argument('--export-cube', type=str, default=None, metavar='DIR',
                        help='Export Gaussian cube file + recipe JSON to DIR (skips rendering)')
    parser.add_argument('--cli', action='store_true', default=True,
                        help='Run in CLI mode (default when --output given)')
    args = parser.parse_args()
    
    logpath = Path(args.logfile)
    if not logpath.exists():
        print(f"Error: file not found: {args.logfile}")
        sys.exit(1)
    
    print(f"Parsing {args.logfile} ...")
    ext = logpath.suffix.lower()
    if ext in ('.log', '.out'):
        atoms, basis_set, canon_wfn, local_wfn, homo_idx = parse_gamess_log(logpath)
    elif ext == '.fchk':
        atoms, basis_set, canon_wfn, local_wfn, homo_idx, mocoeffs_sph, T_s2c = parse_gaussian_fchk(logpath)
    elif ext == '.chk':
        print("Error: binary .chk files cannot be parsed directly.")
        print(f"  Convert first: formchk {logpath.name} {logpath.stem}.fchk")
        sys.exit(1)
    else:
        print(f"Error: unrecognised file extension '{ext}'")
        print("  Expected: .log, .out (GAMESS) or .fchk (Gaussian)")
        sys.exit(1)
    
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
    
    if args.export_cube is not None:
        # Export cube file + recipe JSON, skip rendering
        dirpath = Path(args.export_cube)
        dirpath.mkdir(parents=True, exist_ok=True)
        
        energy = wfn.energies[mo_idx] if mo_idx < len(wfn.energies) else 0.0
        label = wfn.labels[mo_idx] if mo_idx < len(wfn.labels) else 'orbital'
        
        cube_name = f"mo_{mo_idx + 1}_{args.type}.cube"
        cube_path = dirpath / cube_name
        write_cube_file(cube_path, atoms, grid_values, origin, spacing,
                       mo_idx, energy, label)
        
        recipe_path = dirpath / "render_recipe.json"
        write_render_recipe(recipe_path, logpath, [{
            "cube_file": cube_name,
            "mo_idx": mo_idx,
            "wtype": args.type,
            "energy": energy,
            "label": label,
            "isovalue": args.isovalue,
            "grid_spacing": spacing,
        }])
        
        print(f"Exported cube file: {cube_path}")
        print(f"Exported recipe:    {recipe_path}")
        return
    
    verts_pos, faces_pos = extract_isosurface(grid_values, +args.isovalue, origin, spacing)
    verts_neg, faces_neg = extract_isosurface(grid_values, -args.isovalue, origin, spacing)
    
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
        # Forward the positional logfile as part of remaining args for cli_main()
        if args.logfile:
            remaining = [args.logfile] + remaining
        sys.argv = [sys.argv[0]] + remaining
        cli_main()
        return
    
    # GUI mode: initialize vispy with PyQt6 backend before creating canvas
    from vispy.app import use_app
    use_app('pyqt6')
    
    # Check OpenGL availability before creating any widgets
    from PyQt6.QtWidgets import QApplication
    app = QApplication(sys.argv)
    from PyQt6.QtGui import QOpenGLContext, QSurfaceFormat
    
    # Request OpenGL 3.3 core profile (available on Mesa 10+, any distro from 2015+)
    fmt = QSurfaceFormat()
    fmt.setVersion(3, 3)
    fmt.setProfile(QSurfaceFormat.OpenGLContextProfile.CoreProfile)
    QSurfaceFormat.setDefaultFormat(fmt)
    
    # Quick smoke test: create temporary GL context
    temp_ctx = QOpenGLContext()
    if not temp_ctx.create():
        from PyQt6.QtWidgets import QMessageBox
        QMessageBox.critical(
            None, "OpenGL Error",
            "Cannot create OpenGL 3.3 context. Install system GPU drivers.\n\n"
            "Ubuntu/Debian:   sudo apt install libgl1-mesa-glx libegl1-mesa mesa-utils\n"
            "Fedora/RHEL:     sudo dnf install mesa-libGL mesa-libEGL glx-utils\n"
            "Arch:            sudo pacman -S mesa libglvnd\n"
            "openSUSE:        sudo zypper install Mesa-libGL1 Mesa-libEGL1\n\n"
            f"Qt platform: {app.platformName()}\n"
            "On Wayland, ensure qt6-wayland is installed (or try QT_QPA_PLATFORM=xcb)."
        )
        sys.exit(1)
    del temp_ctx
    
    app.setApplicationName("Orbital Visualizer")
    
    viewer = OrbitalViewer(logpath=args.logfile)
    viewer.show()
    
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
