#!/usr/bin/env python3

import os
import re
import numpy as np
import argparse
import sys

### General parameters ##############
try:
    import qcelemental as qcel
    amu2au = qcel.constants.unified_atomic_mass_unit/qcel.constants.electron_mass
    bohr2angs = qcel.constants.bohr2angstroms
    clight = qcel.constants.c
    hbar = qcel.constants.hbar
    bohr2m = qcel.constants.bohr2m
    electron_mass = qcel.constants.electron_mass
except ImportError:
    amu2au = 1822.8884853323707
    bohr2angs = 0.52917721067
    clight = 299792458.0
    hbar = 1.0545718e-34
    bohr2m = 5.2917721067e-11
    electron_mass = 9.10938356e-31
#####################################

ATOMIC_SYMBOLS = {
    1: 'H',  6: 'C',  7: 'N',  8: 'O',  9: 'F',
   11: 'Na', 12: 'Mg', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl',
   19: 'K',  20: 'Ca', 26: 'Fe', 28: 'Ni', 29: 'Cu', 30: 'Zn',
   35: 'Br', 53: 'I',
}

### Extractor functions ################################################
def extract_S(fname, target):
    """Read second-order transition moment tensor from DALTON output.
    Returns symmetric 3x3 matrix for the given target root."""
    Sind = {'XDIPLEN': 0, 'YDIPLEN': 1, 'ZDIPLEN': 2}
    S = np.zeros([3, 3])
    with open(fname) as f:
        for line in f:
            if 'Second order moment in a.u. for' in line:
                opA  = f.readline().split()[-3]
                opB  = f.readline().split()[-3]
                root = int(f.readline().split()[-3])
                if root != target:
                    continue
                f.readline()
                val = float(f.readline().split()[-1])
                S[Sind[opA], Sind[opB]] = val
                S[Sind[opB], Sind[opA]] = val  # fix: symmetrize
    return S


def extract_natoms(fname):
    natoms = -1
    with open(fname) as f:
        for line in f:
            if 'Total number of atoms' in line:
                natoms = int(line.split()[-1])
                break
    return natoms


def ltvector_to_matrix(Alt):
    """Convert lower-triangle vector to full symmetric matrix."""
    # fix: analytic triangular root n*(n+1)/2 = len => n = (-1 + sqrt(1+8*len)) / 2
    n = int((-1 + np.sqrt(1 + 8 * len(Alt))) / 2)
    A = np.zeros((n, n))
    idx = np.tril_indices(n)
    A[idx] = Alt
    d = A.diagonal()
    A += A.T - np.diag(d)
    return A


def read_array(fname, is_lt=False):
    with open(fname) as f:
        line1 = f.readline()
        dim = line1.split()[0]
        A = []
        for line in f:
            A += [float(x) for x in line.split()]
    if dim == '2D':
        i, j = [int(x) for x in line1.split()[1:]]
        # Revert Fortran order
        A = np.array(A).reshape(j, i).flatten('F').reshape(i, j)
    elif dim == '1D':
        i = int(line1.split()[-1])
        A = np.array(A)
        if is_lt:
            A = ltvector_to_matrix(A)
    return A


def read_fcc_geom(fname):
    """Read geometry from FCclasses3 state file (GEOM section, coords in Angstrom)."""
    atnames = None
    xyz = None
    with open(fname) as f:
        for line in f:
            if line.strip() == 'GEOM':   # fix: exact match, not prefix match
                natoms = int(f.readline())
                f.readline()  # title line
                atnames = []
                xyz = []
                for i in range(natoms):
                    a, x, y, z = f.readline().split()
                    atnames.append(a)
                    xyz.append(float(x))
                    xyz.append(float(y))
                    xyz.append(float(z))
                break
    # fix: raise meaningful error if GEOM section was not found
    if atnames is None:
        raise ValueError(f"No 'GEOM' section found in '{fname}'. "
                         "Check that the file is a valid FCclasses3 state file.")
    return atnames, np.array(xyz)


def _parse_gaussian_orientation_block(lines, start):
    """Parse one Standard/Input orientation table from a Gaussian log.
    Returns list of (symbol, x, y, z) in Angstrom, or None if malformed."""
    i = start + 5  # skip: header line, separator, 2 column headers, separator
    atoms = []
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('---'):
            break
        parts = line.split()
        if len(parts) < 6:
            break
        atomic_num = int(parts[1])
        x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
        atoms.append((ATOMIC_SYMBOLS.get(atomic_num, str(atomic_num)), x, y, z))
        i += 1
    return atoms if atoms else None


def read_gaussian_geom(fname):
    """Read the final optimized geometry from a Gaussian log file.
    Returns (atnames, xyz) where xyz is a flat Angstrom array, same format
    as read_fcc_geom."""
    with open(fname) as f:
        lines = f.readlines()
    last_atoms = None
    for i, line in enumerate(lines):
        if re.search(r'(Standard|Input) orientation:', line):
            atoms = _parse_gaussian_orientation_block(lines, i)
            if atoms:
                last_atoms = atoms
    if last_atoms is None:
        raise ValueError(f"No orientation block found in Gaussian log '{fname}'.")
    atnames = [a[0] for a in last_atoms]
    xyz = [c for a in last_atoms for c in a[1:]]
    return atnames, np.array(xyz)


def read_geom(fname):
    """Dispatch geometry reader based on file extension.
    Supports .fcc (FCclasses3 state file) and .log (Gaussian log file)."""
    ext = os.path.splitext(fname)[1].lower()
    if ext == '.log':
        return read_gaussian_geom(fname)
    elif ext == '.fcc':
        return read_fcc_geom(fname)
    else:
        raise ValueError(f"Unrecognised file extension '{ext}'. "
                         "Supported: .fcc (FCclasses3) or .log (Gaussian).")
#######################################################################

########### File writers ##############################################
def write_xyz(atnames, xyz, fname, title):
    natoms = int(len(xyz) / 3)
    with open(fname, 'w') as f:
        print(natoms, file=f)
        print(title, file=f)
        for i in range(natoms):
            j = 3 * i
            print(f'{atnames[i]:<5}   {xyz[j]:10.5f} {xyz[j+1]:10.5f} {xyz[j+2]:10.5f}', file=f)


def write_com(atnames, xyz, fname, title, header):
    natoms = int(len(xyz) / 3)
    # use only the basename so %chk= never contains a directory path
    chkname = os.path.splitext(os.path.basename(fname))[0] + '.chk'
    with open(fname, 'w') as f:
        new_header = header.strip().replace('CHKFILE', chkname).replace('TITLE', title) + '\n'
        print(new_header, file=f, end='')
        for i in range(natoms):
            j = 3 * i
            print(f'{atnames[i]:<5}   {xyz[j]:10.5f} {xyz[j+1]:10.5f} {xyz[j+2]:10.5f}', file=f)
        print('', file=f)
#######################################################################


def main(args):
    fname   = args.f
    disp    = args.disp
    dertype = args.der.upper()

    # fix: validate Q-specific arguments up front with clear messages
    if dertype == 'Q':
        for argname, val in [('lmat', args.lmat), ('mass', args.mass), ('freq', args.freq)]:
            if val is None:
                raise SystemExit(f"error: -{argname} is required when -der Q")
        Lfile  = args.lmat
        Ltype  = args.ltype
        mfile  = args.mass
        Frqfile = args.freq

    # fix: S-type is not implemented — fail explicitly
    if dertype == 'S':
        raise NotImplementedError("Displacement type 'S' (internal coordinates) is not yet implemented.")

    # fix: robust basename extraction (handles absolute paths and multi-dot filenames)
    basename = os.path.splitext(os.path.basename(fname))[0]

    # Output directory: create it if it does not exist
    outdir = args.outdir
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    outbase = os.path.join(outdir, basename) if outdir else basename

    # Load Gaussian header: if the argument is an existing file, read it;
    # otherwise treat the argument value itself as the literal header string.
    if args.gauhead is not None:
        if os.path.isfile(args.gauhead):
            with open(args.gauhead) as f:
                header = f.read()
        else:
            header = args.gauhead
    else:
        header = '''\
%chk=CHKFILE
%nproc=8
%mem=10gb

#p cam-b3lyp/def2TZVP td=(root=1,nstates=3) freq nosymm iop(7/33=1)

TITLE

0 1
'''

    ### Main ###

    atnames, xyz = read_geom(fname)
    natoms = len(atnames)

    if dertype == 'Q':
        # fix: removed dead 'Mass' variable; only Minv is needed
        mass = read_array(mfile) * amu2au
        mass3 = np.array([mass[int(i / 3)] for i in range(3 * natoms)])
        Minv = np.diag(1 / mass3)

        if Ltype == 'MWC':
            L = read_array(Lfile)
            L = np.sqrt(Minv).dot(L)
        else:
            raise NotImplementedError('Modes in internal coordinates not yet supported.')

        Freq = read_array(Frqfile)
        factor = np.sqrt(abs(Freq) * 1.e2 * clight * 2.0 * np.pi / hbar)
        factor *= bohr2m * np.sqrt(electron_mass)

        dQ = disp / factor
        for i, q in enumerate(dQ):
            dX = L[:, i] * q * bohr2angs
            xyz_fw = xyz.copy() + dX
            xyz_bw = xyz.copy() - dX
            write_com(atnames, xyz_fw, outbase + f'_Mode{i+1}_fw.com', f'Disp = +{q:<9.5f} +{disp:<9.5f}', header)
            write_com(atnames, xyz_bw, outbase + f'_Mode{i+1}_bw.com', f'Disp = -{q:<9.5f} -{disp:<9.5f}', header)

    else:  # type X: Cartesian displacements (disp in Angstrom)
        for i in range(natoms):
            for coord, label in enumerate(['xyz1', 'xyz2', 'xyz3']):
                xyz_fw = xyz.copy(); xyz_fw[3*i + coord] += disp
                xyz_bw = xyz.copy(); xyz_bw[3*i + coord] -= disp
                write_com(atnames, xyz_fw, outbase + f'_at{i+1}_{label}_fw.com', f'Disp = +{disp:<9.5f}', header)
                write_com(atnames, xyz_bw, outbase + f'_at{i+1}_{label}_bw.com', f'Disp = -{disp:<9.5f}', header)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generate displaced geometry input files for numerical TDM derivatives.')
    parser.add_argument('-f',       metavar='data.fcc',  help='Geometry source: FCclasses3 state file (.fcc) or Gaussian log file (.log)', required=True)
    parser.add_argument('-bmat',    metavar='B.dat',     help='B matrix file', required=False, default=None)
    parser.add_argument('-gmat',    metavar='G.dat',     help='G matrix file', required=False, default=None)
    parser.add_argument('-mass',    metavar='mass.dat',  help='Mass file (required for -der Q)', required=False)
    parser.add_argument('-freq',    metavar='Freq.dat',  help='Frequency file (required for -der Q)', required=False)
    parser.add_argument('-ltype',   help='Type of L matrix: MWC (default)', required=False, default='MWC')
    parser.add_argument('-lmat',    metavar='L.dat',     help='L matrix file (required for -der Q)', required=False)
    parser.add_argument('-der',     help='Displacement type: X (Cartesian, disp in Angstrom) or Q (normal mode, disp dimensionless)',
                        required=True)
    parser.add_argument('-disp',    help='Displacement magnitude: Angstrom for -der X, dimensionless for -der Q',
                        required=True, type=float)
    parser.add_argument('-outdir',  metavar='DIR',
                        help='Output directory for generated files (created if absent). '
                             'Default: current working directory.',
                        required=False, default=None)
    parser.add_argument('-gauhead', metavar='gau.txt|"string"',
                        help='Gaussian header block: path to a file, or a literal string in quotes. '
                             'Must contain the CHKFILE and TITLE placeholders. '
                             'If omitted, a default B3LYP/6-31G(d) header is used.',
                        required=False, default=None)
    args = parser.parse_args()
    main(args)
