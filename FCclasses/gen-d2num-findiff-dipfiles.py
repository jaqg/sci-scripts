#!/usr/bin/env python3

import numpy as np
import sys
import os
import re
import subprocess

def parse_fchk(fchkname,section):

    def totype(v,vtype):
        if vtype == 'R':
            return np.float64(v)
        elif vtype == 'I':
            return int(v)
        else:
            raise Exception('Unknown variable type: '+vtype)

    # Initialize number of items
    N=0

    # Read data from file
    with open(fchkname) as f:

        found = False
        for line in f:
            if section in line:
                found = True
                break

        if not found:
            raise Exception(f'Section {section} not found in FCHK file {fchkname}')

        line=line.split(section)[1]
        if 'N=' in line:
            # Vector data
            data=line.split()
            vartype=data[0]
            N=int(data[2])
            if vartype == 'R':
                items_per_line = 5
            elif vartype == 'I':
                items_per_line = 6
            nlines = int(N/items_per_line)
            if N%items_per_line != 0:
                nlines += 1
            aa=[]
            for i in range(nlines):
                line = f.readline()
                aa += [totype(x,vartype) for x in line.split()]

            if (len(aa) != N):
                raise Exception('Wrong number of data elements from fchk section: '+section)
            else:
                aa=np.array(aa)

        else:
            # Scalar data
            vartype=line.split()[-2]
            aa = totype(line.split()[-1],vartype)
            N=1


    if N==0:
        raise Exception('Section not found in fchk file')

    return aa

def Hltvector_to_H(Hlt):
    # Exact inverse of n*(n+1)/2 = len  =>  n = (-1 + sqrt(1+8*len)) / 2
    n = int((-1 + np.sqrt(1 + 8*len(Hlt))) / 2)
    A = np.zeros((n,n))
    # Lower
    idx = np.tril_indices(n)
    A[idx] = Hlt
    # Upper
    d = A.diagonal()
    A += A.T - np.diag(d)

    return A


def _sign_flip(vec, ref, imax):
    """Return vec (1D) flipped if its sign at imax disagrees with ref, or unchanged
    if ref[imax] is (near) zero (phase correction cannot be determined)."""
    if abs(ref[imax]) < 1e-12:
        return vec
    if vec[imax] / ref[imax] < 0.:
        return -vec
    return vec


def _sign_flip_grad(grad, dip0_disp, dip0_eq, imax):
    """Return gradient array (any shape) flipped if the displaced TDM (dip0_disp)
    has opposite sign to the equilibrium TDM (dip0_eq) at component imax."""
    if abs(dip0_eq[imax]) < 1e-12 or abs(dip0_disp[imax]) < 1e-12:
        return grad
    if dip0_disp[imax] * dip0_eq[imax] < 0.:
        return -grad
    return grad


def read_etran_from_fchk(fchk_fname, target_state=None):
    """Extract transition dipole moments and their Cartesian gradients from a
    Gaussian fchk file produced by a TD-DFT freq job with iop(7/33=1).

    Layout of 'ETran state values' array (Fortran 1-based, 16 values per state):
      offset+1        : total energy
      offset+2..+4    : electric TDM length gauge (x, y, z)
      offset+5..+7    : P operator / velocity TDM (x, y, z)
      offset+8..+10   : magnetic TDM * (-2)  (x, y, z)
      offset+11..+16  : unknowns
    After Nes*16 state values + 48 padding values, derivatives follow for the
    target state: 3*Nat blocks of 16, same sub-structure as above.

    Returns a dict with keys:
      nat, Nes, target,
      eldip_l (3,), eldip_p (3,), magdip (3,)  -- zeroth order
      grad_eldip_l (3*nat, 3) or None           -- d(eldip_l)/dxi
      grad_eldip_p  (3*nat, 3) or None
      grad_magdip   (3*nat, 3) or None
    """
    nat     = int(parse_fchk(fchk_fname, 'Number of atoms'))
    esc     = parse_fchk(fchk_fname, 'ETran scalars')
    Nes     = int(esc[0])
    Ntarget = int(esc[4]) if target_state is None else int(target_state)

    A = parse_fchk(fchk_fname, 'ETran state values')

    # Zeroth-order properties (convert Fortran 1-based to Python 0-based)
    # Fortran: j = (Ntarget-1)*16 + 2;  A(j:j+2)  →  Python A[j-1:j+2]
    off = (Ntarget - 1) * 16
    eldip_l = A[off + 1 : off + 4].copy()          # length-gauge TDM
    eldip_p = A[off + 4 : off + 7].copy()          # P operator TDM
    magdip  = A[off + 7 : off + 10].copy() / (-2.) # magnetic TDM (FCclasses conv.)

    # Check whether derivatives are present
    # Expected total size: Nes*16 + 48 + 3*Nat*16
    has_ders = (len(A) == Nes * 16 + 48 + 3 * nat * 16)
    if not has_ders:
        return dict(nat=nat, Nes=Nes, target=Ntarget,
                    eldip_l=eldip_l, eldip_p=eldip_p, magdip=magdip,
                    grad_eldip_l=None, grad_eldip_p=None, grad_magdip=None)

    grad_eldip_l = np.zeros((3 * nat, 3))
    grad_eldip_p = np.zeros((3 * nat, 3))
    grad_magdip  = np.zeros((3 * nat, 3))

    for j in range(1, 3 * nat + 1):              # j: 1-based Cartesian index
        kf = 16 * Nes + 48 + 16 * (j - 1)        # Fortran 1-based k
        # Fortran A(kf+2:kf+4) → Python A[kf+1:kf+4]
        grad_eldip_l[j - 1, :] = A[kf + 1 : kf + 4]
        grad_eldip_p[j - 1, :] = A[kf + 4 : kf + 7]
        grad_magdip [j - 1, :] = A[kf + 7 : kf + 10] / (-2.)

    return dict(nat=nat, Nes=Nes, target=Ntarget,
                eldip_l=eldip_l, eldip_p=eldip_p, magdip=magdip,
                grad_eldip_l=grad_eldip_l, grad_eldip_p=grad_eldip_p,
                grad_magdip=grad_magdip)


def _write_dipfile(fname, dip0, dip1):
    """Write an eldip_* / magdip_* file in the FCclasses format:
      line 1 : TDM value (3 components)
      line 2 : TDM value repeated
      lines 3+: TDM gradient, one Cartesian displacement per row (3 components)
    """
    with open(fname, 'w') as f:
        print(f' {dip0[0]:18.9e} {dip0[1]:18.9e} {dip0[2]:18.9e}', file=f)
        print(f' {dip0[0]:18.9e} {dip0[1]:18.9e} {dip0[2]:18.9e}', file=f)
        if dip1 is not None:
            for row in dip1:
                print(f' {row[0]:18.9e} {row[1]:18.9e} {row[2]:18.9e}', file=f)


def extract_eldip_files_from_fchk(ref_fchk, disp_pattern, out_basename=None, out_dir='', gauge='vel'):
    """Generate dipole files directly from Gaussian fchk files, replacing the
    need to run gen_fcc_dipfile for each displaced geometry.

    Args:
        ref_fchk      : path to the reference fchk (e.g. S1 freq at S0 equil.)
        disp_pattern  : format string for displaced fchk paths.
                        Must contain the placeholders {iat}, {ixyz}, {dir}:
                          {iat}  -- atom index (1-based integer)
                          {ixyz} -- Cartesian coordinate index (1, 2, or 3)
                          {dir}  -- displacement direction ('fw' or 'bw')
                        Example:
                          'numerical-derivatives/mol_S0_at{iat}_xyz{ixyz}_{dir}.fchk'
        out_basename  : prefix for the reference output file name
                        (default: stem of ref_fchk)
        out_dir       : directory where output files are written (default: current dir)
        gauge         : 'vel' (default) writes P_eldip_* files;
                        'len' writes eldip_* files (length gauge).
                        magdip_* is always written.

    Output files written to out_dir (gauge='vel'):
        eldip_{out_basename}_fchk                (reference, length-gauge TDM)
        P_eldip_{out_basename}_fchk              (reference, raw P operator)
        eldip_{disp_stem}_fchk                   (each displacement, length-gauge TDM)
        P_eldip_{disp_stem}_fchk                 (each displacement, raw P operator)
        magdip_{disp_stem}_fchk                  (each displacement)

    Output files written to out_dir (gauge='len'):
        eldip_{out_basename}_fchk                (reference, length-gauge TDM)
        eldip_{disp_stem}_fchk                   (each displacement, length-gauge TDM)
        magdip_{disp_stem}_fchk                  (each displacement)

    Note: eldip_* always contains length-gauge TDM (eldip_l from ETran), regardless
    of gauge. For vel gauge, P_eldip_* additionally stores the raw P operator used
    by the P-method compute functions.

    where disp_stem = stem of the displaced fchk file
    (e.g. mol_S0_at1_xyz1_bw for mol_S0_at1_xyz1_bw.fchk)
    """
    if out_basename is None:
        out_basename = os.path.splitext(os.path.basename(ref_fchk))[0]

    def _opath(fname):
        return os.path.join(out_dir, fname) if out_dir else fname

    # --- Reference geometry ---
    print(f'Reading reference: {ref_fchk}')
    ref = read_etran_from_fchk(ref_fchk)
    nat = ref['nat']
    # Always write eldip with length-gauge TDM (eldip_l).
    # For vel gauge also write P_eldip with the raw P operator.
    _write_dipfile(_opath(f'eldip_{out_basename}_fchk'), ref['eldip_l'], ref['grad_eldip_l'])
    if gauge == 'vel':
        _write_dipfile(_opath(f'P_eldip_{out_basename}_fchk'), ref['eldip_p'], ref['grad_eldip_p'])
    print(f'  nat={nat}  Nes={ref["Nes"]}  target state={ref["target"]}')

    # --- Displaced geometries ---
    n_written = 0
    for iat1 in range(1, nat + 1):
        for ii1 in range(1, 4):
            for direction in ('fw', 'bw'):
                try:
                    fpath = disp_pattern.format(iat=iat1, ixyz=ii1, dir=direction)
                except KeyError as e:
                    raise ValueError(
                        f'-disp_pattern is missing placeholder {e}. '
                        f'Required: {{iat}}, {{ixyz}}, {{dir}}') from e
                if not os.path.isfile(fpath):
                    raise FileNotFoundError(
                        f'Displaced fchk not found: {fpath}\n'
                        f'Check the -disp_pattern argument.')
                d = read_etran_from_fchk(fpath)
                # Name files after the stem of the displaced fchk
                tag = os.path.splitext(os.path.basename(fpath))[0]
                # Always write eldip with length-gauge TDM.
                # For vel gauge also write P_eldip with the raw P operator.
                _write_dipfile(_opath(f'eldip_{tag}_fchk'), d['eldip_l'], d['grad_eldip_l'])
                if gauge == 'vel':
                    _write_dipfile(_opath(f'P_eldip_{tag}_fchk'), d['eldip_p'], d['grad_eldip_p'])
                _write_dipfile(_opath(f'magdip_{tag}_fchk'), d['magdip'], d['grad_magdip'])
                n_written += 1

    extra = ', P_eldip' if gauge == 'vel' else ''
    print(f'Wrote reference + {n_written} displaced dipole files '
          f'(eldip{extra}, magdip) for basename "{out_basename}" in "{out_dir or "."}".')


def extract_eldip_files_using_gen(ref_fchk, disp_pattern, out_dir='', gauge='vel',
                                  gen_dipfile='gen_fcc_dipfile'):
    """Generate dipole files by calling gen_fcc_dipfile for each fchk file.

    Calls gen_fcc_dipfile with the requested gauge and writes eldip_* and
    magdip_* files for the reference and every displaced geometry.

    Args:
        ref_fchk     : path to the reference fchk (e.g. S1 freq at S0 equil.)
        disp_pattern : format string for displaced fchk paths, with {iat},
                       {ixyz}, {dir} placeholders.
        out_dir      : directory where output dipole files are written.
        gauge        : 'vel' (default) or 'len'.
        gen_dipfile  : path to the gen_fcc_dipfile executable.
    """
    gauge_str   = 'velocity' if gauge == 'vel' else 'lenght'
    gen_cmd     = os.path.abspath(gen_dipfile) if os.sep in gen_dipfile else gen_dipfile
    out_abs     = os.path.abspath(out_dir) if out_dir else os.path.abspath('.')
    os.makedirs(out_abs, exist_ok=True)

    def _run(fchk_path, eldip_fname, magdip_fname, dest_dir=None):
        # Run gen_fcc_dipfile from the fchk's own directory, passing only the
        # basename as -i so the binary writes output (and the auto-generated
        # P_eldip_* companion for velocity gauge) there without any path issues.
        # Afterwards, move all produced files to dest_dir if it differs from fchk_dir.
        # dest_dir defaults to out_abs; pass fchk_dir explicitly to keep files in place.
        fchk_abs  = os.path.abspath(fchk_path)
        fchk_dir  = os.path.dirname(fchk_abs)
        fchk_base = os.path.basename(fchk_abs)
        target    = dest_dir if dest_dir is not None else out_abs
        cmd = [gen_cmd, '-i', fchk_base, '-gauge', gauge_str,
               '-oe', eldip_fname, '-om', magdip_fname]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=fchk_dir)
        if result.returncode != 0:
            raise RuntimeError(
                f'{gen_dipfile} failed for {fchk_path}:\n'
                f'stdout: {result.stdout}\nstderr: {result.stderr}')
        if fchk_dir != target:
            for fname in [eldip_fname, magdip_fname, 'P_' + eldip_fname]:
                src = os.path.join(fchk_dir, fname)
                if os.path.isfile(src):
                    os.replace(src, os.path.join(target, fname))

    # --- Reference geometry (keep DIP files next to the reference fchk) ---
    print(f'Reading reference: {ref_fchk}')
    nat      = int(parse_fchk(ref_fchk, 'Number of atoms'))
    ref_stem = os.path.splitext(os.path.basename(ref_fchk))[0]
    ref_fchk_dir = os.path.dirname(os.path.abspath(ref_fchk))
    _run(ref_fchk, f'eldip_{ref_stem}_fchk', f'magdip_{ref_stem}_fchk',
         dest_dir=ref_fchk_dir)
    print(f'  nat={nat}')

    # --- Displaced geometries ---
    n_written = 0
    for iat1 in range(1, nat + 1):
        for ii1 in range(1, 4):
            for direction in ('fw', 'bw'):
                try:
                    fpath = disp_pattern.format(iat=iat1, ixyz=ii1, dir=direction)
                except KeyError as e:
                    raise ValueError(
                        f'-disp_pattern is missing placeholder {e}. '
                        f'Required: {{iat}}, {{ixyz}}, {{dir}}') from e
                if not os.path.isfile(fpath):
                    raise FileNotFoundError(
                        f'Displaced fchk not found: {fpath}\n'
                        f'Check the -disp_pattern argument.')
                tag = os.path.splitext(os.path.basename(fpath))[0]
                _run(fpath, f'eldip_{tag}_fchk', f'magdip_{tag}_fchk')
                n_written += 1

    print(f'Wrote reference + {n_written} displaced dipole files '
          f'(eldip, magdip) for ref "{ref_stem}" in "{out_abs}".')


def compute_magdip_d1num(nat, disp_stem_pattern, fname, delta, dip0, imax, dip_dir=''):
    fmagdip = open(fname,'w')

    def _dpath(f):
        return os.path.join(dip_dir, f) if dip_dir else f

    # Print constant term
    print(f' {dip0[0]:18.9e} {dip0[1]:18.9e} {dip0[2]:18.9e}', file=fmagdip)
    print(f' {dip0[0]:18.9e} {dip0[1]:18.9e} {dip0[2]:18.9e}', file=fmagdip)

    # Compute derivatives
    for iat1 in range(1,nat+1):
        for ii1 in range(1,4):
            stem_bw = disp_stem_pattern.format(iat=iat1, ixyz=ii1, dir='bw')
            stem_fw = disp_stem_pattern.format(iat=iat1, ixyz=ii1, dir='fw')
            dip0BW = _sign_flip(np.loadtxt(_dpath(f'magdip_{stem_bw}_fchk'))[0], dip0, imax)
            dip0FW = _sign_flip(np.loadtxt(_dpath(f'magdip_{stem_fw}_fchk'))[0], dip0, imax)
            dd = (dip0FW - dip0BW) / 2.0 / delta
            print(f' {dd[0]:18.9e} {dd[1]:18.9e} {dd[2]:18.9e}', file=fmagdip)
    fmagdip.close()


def compute_eldip_d1num_P(nat, disp_stem_pattern, refS0, fname, delta, dip_dir='', refS1=None, ref_basename=None):
    feldip = open(fname,'w')

    def _dpath(f):
        return os.path.join(dip_dir, f) if dip_dir else f

    # Read additional data (energies, gradients and Hessians from S0 and S1 states)
    if refS1 is None:
        refS1 = ref_basename + '.fchk' if ref_basename else disp_stem_pattern + '.fchk'
    Egs = parse_fchk(refS0,'Total Energy')
    Ees = parse_fchk(refS1,'Total Energy')
    Ggs = parse_fchk(refS0,'Cartesian Gradient')
    Ges = parse_fchk(refS1,'Cartesian Gradient')
    # Transition values
    Ev = Ees - Egs
    Gv = Ges - Ggs

    # Read P at reference geometry (ref_basename = S1 stem)
    P0 = np.loadtxt(_dpath(f'P_eldip_{ref_basename}_fchk'))[0]
    # And transform to eldip
    dip0 = - P0 / Ev

    # Select larger element of eldip to check sign
    dip0_abs = np.abs(dip0)
    imax = dip0_abs.argmax()

    # Print constant term
    print(f' {dip0[0]:18.9e} {dip0[1]:18.9e} {dip0[2]:18.9e}', file=feldip)
    print(f' {dip0[0]:18.9e} {dip0[1]:18.9e} {dip0[2]:18.9e}', file=feldip)

    # Compute derivatives of P
    P1 = np.zeros((nat*3,3))
    for iat1 in range(1,nat+1):
        for ii1 in range(1,4):
            i = (iat1-1)*3 + ii1 - 1
            stem_bw = disp_stem_pattern.format(iat=iat1, ixyz=ii1, dir='bw')
            stem_fw = disp_stem_pattern.format(iat=iat1, ixyz=ii1, dir='fw')
            P0BW = _sign_flip(np.loadtxt(_dpath(f'P_eldip_{stem_bw}_fchk'))[0], P0, imax)
            P0FW = _sign_flip(np.loadtxt(_dpath(f'P_eldip_{stem_fw}_fchk'))[0], P0, imax)
            dd = (P0FW - P0BW) / 2.0 / delta
            P1[i,:] = dd[:]

    # Compute derivatives
    dip1 = np.zeros((nat*3,3))
    dip1[:,0] = (-P1[:,0] - dip0[0] * Gv)/Ev
    dip1[:,1] = (-P1[:,1] - dip0[1] * Gv)/Ev
    dip1[:,2] = (-P1[:,2] - dip0[2] * Gv)/Ev

    # Write data
    for v in dip1:
        print(f' {v[0]:18.9e} {v[1]:18.9e} {v[2]:18.9e}', file=feldip)

    feldip.close()


def compute_magdip_d2num(nat, disp_stem_pattern, fname, delta, imax, symmetrize=True, dip_dir='', dip0_eq=None):
    fmagdip = open(fname,'w')

    def _dpath(f):
        return os.path.join(dip_dir, f) if dip_dir else f

    # Compute derivatives
    dip11 = np.zeros((nat*3,nat*3,3))
    for iat1 in range(1,nat+1):
        for ii1 in range(1,4):
            i = (iat1-1)*3 + ii1 - 1
            stem_bw = disp_stem_pattern.format(iat=iat1, ixyz=ii1, dir='bw')
            stem_fw = disp_stem_pattern.format(iat=iat1, ixyz=ii1, dir='fw')
            dip0BW = np.loadtxt(_dpath(f'magdip_{stem_bw}_fchk'))[0]
            dip1BW = np.loadtxt(_dpath(f'magdip_{stem_bw}_fchk'),skiprows=2)
            if dip0_eq is not None:
                dip1BW = _sign_flip_grad(dip1BW, dip0BW, dip0_eq, imax)
            dip0FW = np.loadtxt(_dpath(f'magdip_{stem_fw}_fchk'))[0]
            dip1FW = np.loadtxt(_dpath(f'magdip_{stem_fw}_fchk'),skiprows=2)
            if dip0_eq is not None:
                dip1FW = _sign_flip_grad(dip1FW, dip0FW, dip0_eq, imax)
            dd = (dip1FW - dip1BW) / 2.0 / delta
            for j in range(nat*3):
                dip11[i,j,:] = dd[j,:]

    # Symmetrize
    if symmetrize:
        for i in range(nat*3):
            for j in range(i,nat*3):
                dip11[i,j,:] = 0.5 * (dip11[i,j,:] + dip11[j,i,:])
                dip11[j,i,:] = dip11[i,j,:]

    # Write data
    for vv in dip11:
        for v in vv:
            print(f' {v[0]:18.9e} {v[1]:18.9e} {v[2]:18.9e}', file=fmagdip)

    fmagdip.close()

def compute_magdip_d2num_5p(nat, disp_stem_pattern, fname, delta, imax, symmetrize=True, dip_dir='', dip0_eq=None):
    fmagdip = open(fname,'w')

    def _dpath(f):
        return os.path.join(dip_dir, f) if dip_dir else f

    # Compute derivatives
    dip11 = np.zeros((nat*3,nat*3,3))
    for iat1 in range(1,nat+1):
        for ii1 in range(1,4):
            i = (iat1-1)*3 + ii1 - 1
            stem_bw = disp_stem_pattern.format(iat=iat1, ixyz=ii1, dir='bw')
            stem_fw = disp_stem_pattern.format(iat=iat1, ixyz=ii1, dir='fw')
            dip0BW_1 = np.loadtxt(_dpath(f'HalfDelta/magdip_{stem_bw}_fchk'))[0]
            dip1BW_1 = np.loadtxt(_dpath(f'HalfDelta/magdip_{stem_bw}_fchk'),skiprows=2)
            if dip0_eq is not None:
                dip1BW_1 = _sign_flip_grad(dip1BW_1, dip0BW_1, dip0_eq, imax)
            dip0FW_1 = np.loadtxt(_dpath(f'HalfDelta/magdip_{stem_fw}_fchk'))[0]
            dip1FW_1 = np.loadtxt(_dpath(f'HalfDelta/magdip_{stem_fw}_fchk'),skiprows=2)
            if dip0_eq is not None:
                dip1FW_1 = _sign_flip_grad(dip1FW_1, dip0FW_1, dip0_eq, imax)
            dip0BW_2 = np.loadtxt(_dpath(f'magdip_{stem_bw}_fchk'))[0]
            dip1BW_2 = np.loadtxt(_dpath(f'magdip_{stem_bw}_fchk'),skiprows=2)
            if dip0_eq is not None:
                dip1BW_2 = _sign_flip_grad(dip1BW_2, dip0BW_2, dip0_eq, imax)
            dip0FW_2 = np.loadtxt(_dpath(f'magdip_{stem_fw}_fchk'))[0]
            dip1FW_2 = np.loadtxt(_dpath(f'magdip_{stem_fw}_fchk'),skiprows=2)
            if dip0_eq is not None:
                dip1FW_2 = _sign_flip_grad(dip1FW_2, dip0FW_2, dip0_eq, imax)
            dd = (-dip1FW_2 + 8*dip1FW_1 - 8*dip1BW_1 + dip1BW_2) / 12.0 / delta
            for j in range(nat*3):
                dip11[i,j,:] = dd[j,:]

    # Symmetrize
    if symmetrize:
        for i in range(nat*3):
            for j in range(i,nat*3):
                dip11[i,j,:] = 0.5 * (dip11[i,j,:] + dip11[j,i,:])
                dip11[j,i,:] = dip11[i,j,:]

    # Write data
    for vv in dip11:
        for v in vv:
            print(f' {v[0]:18.9e} {v[1]:18.9e} {v[2]:18.9e}', file=fmagdip)

    fmagdip.close()


def compute_eldip_d2num_E(nat, disp_stem_pattern, fname, delta, imax, symmetrize=True, dip_dir='', dip0_eq=None):
    feldip = open(fname,'w')

    def _dpath(f):
        return os.path.join(dip_dir, f) if dip_dir else f

    # Compute derivatives
    dip11 = np.zeros((nat*3,nat*3,3))
    for iat1 in range(1,nat+1):
        for ii1 in range(1,4):
            i = (iat1-1)*3 + ii1 - 1
            stem_bw = disp_stem_pattern.format(iat=iat1, ixyz=ii1, dir='bw')
            stem_fw = disp_stem_pattern.format(iat=iat1, ixyz=ii1, dir='fw')
            dip0BW = np.loadtxt(_dpath(f'eldip_{stem_bw}_fchk'))[0]
            dip1BW = np.loadtxt(_dpath(f'eldip_{stem_bw}_fchk'),skiprows=2)
            if dip0_eq is not None:
                dip1BW = _sign_flip_grad(dip1BW, dip0BW, dip0_eq, imax)
            dip0FW = np.loadtxt(_dpath(f'eldip_{stem_fw}_fchk'))[0]
            dip1FW = np.loadtxt(_dpath(f'eldip_{stem_fw}_fchk'),skiprows=2)
            if dip0_eq is not None:
                dip1FW = _sign_flip_grad(dip1FW, dip0FW, dip0_eq, imax)
            dd = (dip1FW - dip1BW) / 2.0 / delta
            for j in range(nat*3):
                dip11[i,j,:] = dd[j,:]

    # Symmetrize
    if symmetrize:
        for i in range(nat*3):
            for j in range(i,nat*3):
                dip11[i,j,:] = 0.5 * (dip11[i,j,:] + dip11[j,i,:])
                dip11[j,i,:] = dip11[i,j,:]

    # Write data
    for vv in dip11:
        for v in vv:
            print(f' {v[0]:18.9e} {v[1]:18.9e} {v[2]:18.9e}', file=feldip)

    feldip.close()


def compute_eldip_d2num_E_5p(nat, disp_stem_pattern, fname, delta, imax, symmetrize=True, dip_dir='', dip0_eq=None):
    feldip = open(fname,'w')

    def _dpath(f):
        return os.path.join(dip_dir, f) if dip_dir else f

    # Compute derivatives
    dip11 = np.zeros((nat*3,nat*3,3))
    for iat1 in range(1,nat+1):
        for ii1 in range(1,4):
            i = (iat1-1)*3 + ii1 - 1
            stem_bw = disp_stem_pattern.format(iat=iat1, ixyz=ii1, dir='bw')
            stem_fw = disp_stem_pattern.format(iat=iat1, ixyz=ii1, dir='fw')
            dip0BW_1 = np.loadtxt(_dpath(f'HalfDelta/eldip_{stem_bw}_fchk'))[0]
            dip1BW_1 = np.loadtxt(_dpath(f'HalfDelta/eldip_{stem_bw}_fchk'),skiprows=2)
            if dip0_eq is not None:
                dip1BW_1 = _sign_flip_grad(dip1BW_1, dip0BW_1, dip0_eq, imax)
            dip0FW_1 = np.loadtxt(_dpath(f'HalfDelta/eldip_{stem_fw}_fchk'))[0]
            dip1FW_1 = np.loadtxt(_dpath(f'HalfDelta/eldip_{stem_fw}_fchk'),skiprows=2)
            if dip0_eq is not None:
                dip1FW_1 = _sign_flip_grad(dip1FW_1, dip0FW_1, dip0_eq, imax)
            dip0BW_2 = np.loadtxt(_dpath(f'eldip_{stem_bw}_fchk'))[0]
            dip1BW_2 = np.loadtxt(_dpath(f'eldip_{stem_bw}_fchk'),skiprows=2)
            if dip0_eq is not None:
                dip1BW_2 = _sign_flip_grad(dip1BW_2, dip0BW_2, dip0_eq, imax)
            dip0FW_2 = np.loadtxt(_dpath(f'eldip_{stem_fw}_fchk'))[0]
            dip1FW_2 = np.loadtxt(_dpath(f'eldip_{stem_fw}_fchk'),skiprows=2)
            if dip0_eq is not None:
                dip1FW_2 = _sign_flip_grad(dip1FW_2, dip0FW_2, dip0_eq, imax)
            dd = (-dip1FW_2 + 8*dip1FW_1 - 8*dip1BW_1 + dip1BW_2) / 12.0 / delta
            for j in range(nat*3):
                dip11[i,j,:] = dd[j,:]

    # Symmetrize
    if symmetrize:
        for i in range(nat*3):
            for j in range(i,nat*3):
                dip11[i,j,:] = 0.5 * (dip11[i,j,:] + dip11[j,i,:])
                dip11[j,i,:] = dip11[i,j,:]

    # Write data
    for vv in dip11:
        for v in vv:
            print(f' {v[0]:18.9e} {v[1]:18.9e} {v[2]:18.9e}', file=feldip)

    feldip.close()


def compute_eldip_d2num_P(nat, disp_stem_pattern, refS0, fname, delta, symmetrize=True, debug_fname=None, dip_dir='', refS1=None, ref_basename=None):
    feldip = open(fname,'w')

    def _dpath(f):
        return os.path.join(dip_dir, f) if dip_dir else f

    # Read additional data (energies, gradients and Hessians from S0 and S1 states)
    if refS1 is None:
        refS1 = ref_basename + '.fchk' if ref_basename else disp_stem_pattern + '.fchk'
    Egs = parse_fchk(refS0,'Total Energy')
    Ees = parse_fchk(refS1,'Total Energy')
    Ggs = parse_fchk(refS0,'Cartesian Gradient')
    Ges = parse_fchk(refS1,'Cartesian Gradient')
    Hgs = parse_fchk(refS0,'Cartesian Force Constants')
    Hes = parse_fchk(refS1,'Cartesian Force Constants')
    # Transition values
    Ev = Ees - Egs
    Gv = Ges - Ggs
    Hv = Hes - Hgs
    # Get complete Hessian matrix
    HHv = Hltvector_to_H(Hv)

    # Read P at reference geometry (ref_basename = S1 stem)
    P0 = np.loadtxt(_dpath(f'P_eldip_{ref_basename}_fchk'))[0]
    P1 = np.loadtxt(_dpath(f'P_eldip_{ref_basename}_fchk'),skiprows=2)
    # And transform to eldip
    dip0 = - P0 / Ev
    dip1 = - P1 / Ev
    for i in range(3*nat):
        dip1[i,:] += P0[:] * Gv[i] / Ev**2
    if debug_fname is not None:
        ftmp = open(debug_fname,'w')
        print(f'{dip0[0]:12.5e} {dip0[1]:12.5e} {dip0[2]:12.5e}', file=ftmp)
        print(f'{dip0[0]:12.5e} {dip0[1]:12.5e} {dip0[2]:12.5e}', file=ftmp)
        for v in dip1:
            print(f'{v[0]:12.5e} {v[1]:12.5e} {v[2]:12.5e}', file=ftmp)
        ftmp.close()

    # Select larger element of P to check sign
    P0_abs = np.abs(P0)
    imax = P0_abs.argmax()

    # Compute derivatives of P
    P11 = np.zeros((nat*3,nat*3,3))
    for iat1 in range(1,nat+1):
        for ii1 in range(1,4):
            i = (iat1-1)*3 + ii1 - 1
            stem_bw = disp_stem_pattern.format(iat=iat1, ixyz=ii1, dir='bw')
            stem_fw = disp_stem_pattern.format(iat=iat1, ixyz=ii1, dir='fw')
            P0BW = np.loadtxt(_dpath(f'P_eldip_{stem_bw}_fchk'))[0]
            P1BW = np.loadtxt(_dpath(f'P_eldip_{stem_bw}_fchk'),skiprows=2)
            P1BW = _sign_flip_grad(P1BW, P0BW, P0, imax)
            P0FW = np.loadtxt(_dpath(f'P_eldip_{stem_fw}_fchk'))[0]
            P1FW = np.loadtxt(_dpath(f'P_eldip_{stem_fw}_fchk'),skiprows=2)
            P1FW = _sign_flip_grad(P1FW, P0FW, P0, imax)
            dd = (P1FW - P1BW) / 2.0 / delta
            for j in range(nat*3):
                P11[i,j,:] = dd[j,:]

    # Compute intermediate elements
    dip1Gv = np.zeros((nat*3,nat*3,3))
    Gvdip1 = np.zeros((nat*3,nat*3,3))
    for i in range(nat*3):
        for j in range(nat*3):
            dip1Gv[i,j,0] = dip1[i,0] * Gv[j]
            dip1Gv[i,j,1] = dip1[i,1] * Gv[j]
            dip1Gv[i,j,2] = dip1[i,2] * Gv[j]
            Gvdip1[i,j,0] = Gv[i] * dip1[j,0]
            Gvdip1[i,j,1] = Gv[i] * dip1[j,1]
            Gvdip1[i,j,2] = Gv[i] * dip1[j,2]

    # Compute derivatives
    dip11 = np.zeros((nat*3,nat*3,3))
    dip11[:,:,0] = (-P11[:,:,0] - dip1Gv[:,:,0] - Gvdip1[:,:,0] - dip0[0] * HHv)/Ev
    dip11[:,:,1] = (-P11[:,:,1] - dip1Gv[:,:,1] - Gvdip1[:,:,1] - dip0[1] * HHv)/Ev
    dip11[:,:,2] = (-P11[:,:,2] - dip1Gv[:,:,2] - Gvdip1[:,:,2] - dip0[2] * HHv)/Ev

    # Symmetrize
    if symmetrize:
        for i in range(nat*3):
            for j in range(i,nat*3):
                dip11[i,j,:] = 0.5 * (dip11[i,j,:] + dip11[j,i,:])
                dip11[j,i,:] = dip11[i,j,:]

    # Write data
    for vv in dip11:
        for v in vv:
            print(f' {v[0]:18.9e} {v[1]:18.9e} {v[2]:18.9e}', file=feldip)

    feldip.close()


def compute_eldip_d2num_P_5p(nat, disp_stem_pattern, refS0, fname, delta, symmetrize=True, debug_fname=None, dip_dir='', refS1=None, ref_basename=None):
    feldip = open(fname,'w')

    def _dpath(f):
        return os.path.join(dip_dir, f) if dip_dir else f

    # Read additional data (energies, gradients and Hessians from S0 and S1 states)
    if refS1 is None:
        refS1 = ref_basename + '.fchk' if ref_basename else disp_stem_pattern + '.fchk'
    Egs = parse_fchk(refS0,'Total Energy')
    Ees = parse_fchk(refS1,'Total Energy')
    Ggs = parse_fchk(refS0,'Cartesian Gradient')
    Ges = parse_fchk(refS1,'Cartesian Gradient')
    Hgs = parse_fchk(refS0,'Cartesian Force Constants')
    Hes = parse_fchk(refS1,'Cartesian Force Constants')
    # Transition values
    Ev = Ees - Egs
    Gv = Ges - Ggs
    Hv = Hes - Hgs
    # Get complete Hessian matrix
    HHv = Hltvector_to_H(Hv)

    # Read P at reference geometry (ref_basename = S1 stem)
    P0 = np.loadtxt(_dpath(f'P_eldip_{ref_basename}_fchk'))[0]
    P1 = np.loadtxt(_dpath(f'P_eldip_{ref_basename}_fchk'),skiprows=2)
    # And transform to eldip
    dip0 = - P0 / Ev
    dip1 = - P1 / Ev
    for i in range(3*nat):
        dip1[i,:] += P0[:] * Gv[i] / Ev**2
    if debug_fname is not None:
        ftmp = open(debug_fname,'w')
        print(f'{dip0[0]:12.5e} {dip0[1]:12.5e} {dip0[2]:12.5e}', file=ftmp)
        print(f'{dip0[0]:12.5e} {dip0[1]:12.5e} {dip0[2]:12.5e}', file=ftmp)
        for v in dip1:
            print(f'{v[0]:12.5e} {v[1]:12.5e} {v[2]:12.5e}', file=ftmp)
        ftmp.close()

    # Select larger element of P to check sign
    P0_abs = np.abs(P0)
    imax = P0_abs.argmax()

    # Compute derivatives of P
    P11 = np.zeros((nat*3,nat*3,3))
    for iat1 in range(1,nat+1):
        for ii1 in range(1,4):
            i = (iat1-1)*3 + ii1 - 1
            stem_bw = disp_stem_pattern.format(iat=iat1, ixyz=ii1, dir='bw')
            stem_fw = disp_stem_pattern.format(iat=iat1, ixyz=ii1, dir='fw')
            P0BW_1 = np.loadtxt(_dpath(f'HalfDelta/P_eldip_{stem_bw}_fchk'))[0]
            P1BW_1 = np.loadtxt(_dpath(f'HalfDelta/P_eldip_{stem_bw}_fchk'),skiprows=2)
            P1BW_1 = _sign_flip_grad(P1BW_1, P0BW_1, P0, imax)
            P0FW_1 = np.loadtxt(_dpath(f'HalfDelta/P_eldip_{stem_fw}_fchk'))[0]
            P1FW_1 = np.loadtxt(_dpath(f'HalfDelta/P_eldip_{stem_fw}_fchk'),skiprows=2)
            P1FW_1 = _sign_flip_grad(P1FW_1, P0FW_1, P0, imax)
            P0BW_2 = np.loadtxt(_dpath(f'P_eldip_{stem_bw}_fchk'))[0]
            P1BW_2 = np.loadtxt(_dpath(f'P_eldip_{stem_bw}_fchk'),skiprows=2)
            P1BW_2 = _sign_flip_grad(P1BW_2, P0BW_2, P0, imax)
            P0FW_2 = np.loadtxt(_dpath(f'P_eldip_{stem_fw}_fchk'))[0]
            P1FW_2 = np.loadtxt(_dpath(f'P_eldip_{stem_fw}_fchk'),skiprows=2)
            P1FW_2 = _sign_flip_grad(P1FW_2, P0FW_2, P0, imax)
            dd = (-P1FW_2 + 8*P1FW_1 - 8*P1BW_1 + P1BW_2) / 12.0 / delta
            for j in range(nat*3):
                P11[i,j,:] = dd[j,:]

    # Compute intermediate elements
    dip1Gv = np.zeros((nat*3,nat*3,3))
    Gvdip1 = np.zeros((nat*3,nat*3,3))
    for i in range(nat*3):
        for j in range(nat*3):
            dip1Gv[i,j,0] = dip1[i,0] * Gv[j]
            dip1Gv[i,j,1] = dip1[i,1] * Gv[j]
            dip1Gv[i,j,2] = dip1[i,2] * Gv[j]
            Gvdip1[i,j,0] = Gv[i] * dip1[j,0]
            Gvdip1[i,j,1] = Gv[i] * dip1[j,1]
            Gvdip1[i,j,2] = Gv[i] * dip1[j,2]

    # Compute derivatives
    dip11 = np.zeros((nat*3,nat*3,3))
    dip11[:,:,0] = (-P11[:,:,0] - dip1Gv[:,:,0] - Gvdip1[:,:,0] - dip0[0] * HHv)/Ev
    dip11[:,:,1] = (-P11[:,:,1] - dip1Gv[:,:,1] - Gvdip1[:,:,1] - dip0[1] * HHv)/Ev
    dip11[:,:,2] = (-P11[:,:,2] - dip1Gv[:,:,2] - Gvdip1[:,:,2] - dip0[2] * HHv)/Ev

    # Symmetrize
    if symmetrize:
        for i in range(nat*3):
            for j in range(i,nat*3):
                dip11[i,j,:] = 0.5 * (dip11[i,j,:] + dip11[j,i,:])
                dip11[j,i,:] = dip11[i,j,:]

    # Write data
    for vv in dip11:
        for v in vv:
            print(f' {v[0]:18.9e} {v[1]:18.9e} {v[2]:18.9e}', file=feldip)

    feldip.close()


def _detect_mol_prefix(disp_dir):
    """Scan disp_dir for files matching <mol>_at<N>_xyz<M>_(fw|bw).fchk and return <mol>.

    Raises FileNotFoundError if the directory does not exist, ValueError if no
    matching files are found or if multiple different prefixes are detected
    (in which case -disp_pattern must be given explicitly).
    """
    pat = re.compile(r'^(.+)_at\d+_xyz[123]_(fw|bw)\.fchk$')
    try:
        entries = os.listdir(disp_dir)
    except FileNotFoundError:
        raise FileNotFoundError(f'-disp_dir not found: {disp_dir}')
    prefixes = {m.group(1) for name in entries if (m := pat.match(name))}
    if not prefixes:
        raise ValueError(
            f'No displaced fchk files found in "{disp_dir}".\n'
            f'Expected files matching <mol>_at<N>_xyz<M>_(fw|bw).fchk')
    if len(prefixes) > 1:
        raise ValueError(
            f'Multiple mol prefixes found in "{disp_dir}": {sorted(prefixes)}.\n'
            f'Use -disp_pattern to specify the pattern explicitly.')
    return prefixes.pop()


if __name__ == "__main__":
    import argparse
    import warnings

    # Input parser. Set flags
    parser = argparse.ArgumentParser(
        description="Compute second TDM derivatives from displaced fchk files."
    )
    parser.add_argument(
        "-f", metavar="",
        help="FCHK for the reference geometry S1 data",
        required=True)
    parser.add_argument(
        "-f0", metavar="",
        help="FCHK for the reference geometry S0 data",
    )
    parser.add_argument(
        "-gauge",
        help='Representation of the electric TDM: velocity (or vel, default) or length (or len)',
        default="vel",
    )
    parser.add_argument(
        "-method",
        help='Method to compute derivatives in the velocity gauge: P (from P_eldip files) or E (from eldip files)',
        default="P",
    )
    parser.add_argument(
        "-oe", metavar="eldip_d2num",
        help="Output file with eldip second ders",
        default=None,
    )
    parser.add_argument(
        "-om", metavar="magdip_d2num",
        help="Output file with magdip second ders",
        default=None,
    )
    parser.add_argument(
        "-delta",
        metavar="<dx>",
        help="Length of the bwd and fwd displacements in Angs",
        type=float,
        default=0.001,
    )
    parser.add_argument(
        "-nosymmetrize",
        action="store_true",
        help="Do not symmetrize second der matrices",
        default=False,
    )
    parser.add_argument(
        "-d5points",
        action="store_true",
        help="Use 5-points numerical 1st derivatives",
        default=False,
    )
    parser.add_argument(
        "-withPders",
        action="store_true",
        help="Compute derivatives from P elements (default when -withEders is not given)",
        default=False,
    )
    parser.add_argument(
        "-withEders",
        action="store_true",
        help="Compute derivatives from eldip elements",
        default=False,
    )
    parser.add_argument(
        "-debug",
        action="store_true",
        help="Write eldip_from_P debug file when using P-method",
        default=False,
    )
    # --- fchk extraction mode ---
    parser.add_argument(
        "-extract",
        action="store_true",
        help="Extract dipole files by calling gen_fcc_dipfile for the reference "
             "and each displaced fchk. "
             "Requires -f (reference fchk) and -disp_pattern.",
        default=False,
    )
    parser.add_argument(
        "-disp_dir",
        metavar="DIR",
        help="Directory containing the displaced fchk files. "
             "If -disp_pattern is not given, the mol prefix is auto-detected "
             "from files matching <mol>_at<N>_xyz<M>_(fw|bw).fchk in this directory "
             "and the pattern is constructed automatically.",
        default=None,
    )
    parser.add_argument(
        "-disp_pattern",
        metavar="PATTERN",
        help="Format string for displaced fchk file paths. "
             "Must contain {iat} (atom index, 1-based), {ixyz} (coordinate 1/2/3), "
             "and {dir} (fw or bw). Example: "
             "'numerical-derivatives/mol_S0_at{iat}_xyz{ixyz}_{dir}.fchk'. "
             "If omitted and -disp_dir is given, the pattern is built automatically.",
        default=None,
    )
    parser.add_argument(
        "-gen_dipfile",
        metavar="PATH",
        help="Path to the gen_fcc_dipfile executable used for dipole extraction "
             "(default: gen_fcc_dipfile, assumed on PATH). "
             "Example: ../../gen_fcc_dipfile",
        default="gen_fcc_dipfile",
    )
    # Parse input
    args = parser.parse_args()

    # Select gauge (needed before extract mode too)
    _gauge_map = {'vel': 'vel', 'velocity': 'vel', 'len': 'len', 'length': 'len'}
    gauge = _gauge_map.get(args.gauge.lower())
    if gauge is None:
        print('Gauge choice: velocity (vel) or length (len)')
        sys.exit()

    # Resolve disp_pattern from disp_dir if not given explicitly
    if args.disp_pattern is None and args.disp_dir is not None:
        try:
            mol = _detect_mol_prefix(args.disp_dir)
        except (FileNotFoundError, ValueError) as e:
            print(f'error: {e}')
            sys.exit(1)
        print(f'Auto-detected mol prefix: "{mol}"')
        args.disp_pattern = os.path.join(
            args.disp_dir, f'{mol}_at{{iat}}_xyz{{ixyz}}_{{dir}}.fchk')
        print(f'Using disp_pattern: {args.disp_pattern}')

    # ==================
    # Extract mode: generate dipole files from fchk, then exit
    # ==================
    if args.extract:
        if args.disp_pattern is None:
            print('error: -extract requires -disp_pattern or -disp_dir')
            sys.exit(1)
        extract_eldip_files_using_gen(
            ref_fchk=args.f,
            disp_pattern=args.disp_pattern,
            out_dir=os.path.dirname(args.disp_pattern),
            gauge=gauge,
            gen_dipfile=args.gen_dipfile,
        )
        sys.exit(0)

    # ==================
    # Initial settings
    # ==================
    # Take reference FCHK
    ref = args.f
    ext = ref.split('.')[-1].upper()
    if ext != 'FCHK':
        print('Only fchk files supported (-f)')
        sys.exit()
    ref0 = args.f0
    if ref0 is not None:
        ext = ref0.split('.')[-1].upper()
        if ext != 'FCHK':
            print('Only fchk files supported (-f0)')
            sys.exit()
    # S1 stem — used to name the reference dipole file and look it up
    basename = os.path.splitext(os.path.basename(ref))[0]
    # S0 stem — used to name the output d2num files (derivatives at S0 geometry)
    out_basename = os.path.splitext(os.path.basename(ref0))[0] if ref0 else basename
    # Displaced stem pattern — format string used to locate intermediate dipole files.
    # Derived from the directory-less, extension-less form of disp_pattern.
    # Example: "def2TZVP-H2O-OR-S0_at{iat}_xyz{ixyz}_{dir}"
    if args.disp_pattern:
        disp_stem_pattern = os.path.splitext(os.path.basename(args.disp_pattern))[0]
    else:
        disp_stem_pattern = None

    # Check the ref file is accessible
    try:
        f = open(ref)
    except FileNotFoundError:
        print('Make sure files '+ref+' and ders exit')
        sys.exit()

    # Derive dip_dir: directory where extracted dipole files live (and will be written).
    # Inferred from the directory component of disp_pattern; defaults to '' (current dir).
    dip_dir = os.path.dirname(args.disp_pattern) if args.disp_pattern else ''

    # Auto-extract: run extraction step if reference dipole file is missing.
    # The reference file name depends on gauge:
    #   vel → P_eldip_{basename}_fchk
    #   len → eldip_{basename}_fchk
    ref_dipfile_name = f'eldip_{basename}_fchk'
    ref_dipfile = os.path.join(dip_dir, ref_dipfile_name) if dip_dir else ref_dipfile_name
    if not os.path.isfile(ref_dipfile):
        if args.disp_pattern is None:
            print(f'error: {ref_dipfile} not found and -disp_pattern was not given.')
            print('Either run -extract first or provide -disp_pattern to auto-extract.')
            sys.exit(1)
        print(f'{ref_dipfile} not found — running extraction step automatically.')
        extract_eldip_files_using_gen(ref_fchk=ref, disp_pattern=args.disp_pattern,
                                      out_dir=dip_dir, gauge=gauge,
                                      gen_dipfile=args.gen_dipfile)
    else:
        print(f'Found {ref_dipfile} — skipping extraction step.')

    # Displacement over coordinates, in Bohr
    delta = args.delta * 1.88973

    # Output files
    eldip_fname  = args.oe
    magdip_fname = args.om
    if eldip_fname is None:
        eldip_fname = f'eldip_{out_basename}_{gauge}_d2num'
    if magdip_fname is None:
        magdip_fname = f'magdip_{out_basename}_d2num'

    # Symmetrize (default is yes)
    symmetrize = not args.nosymmetrize

    # Get system info from ref file
    nat = parse_fchk(ref,'Number of atoms')

    use_E = args.withEders or gauge == 'len'
    use_P = args.withPders or (not use_E)  # P is default when neither flag nor len gauge
    debug_fname = f'eldip_from_P_{out_basename}' if args.debug else None

    # Read imax for phase correction directly from the reference fchk —
    # no intermediate eldip_*/magdip_* reference files are needed.
    ref_data = read_etran_from_fchk(ref)
    imax_eldip  = int(np.abs(ref_data['eldip_l']).argmax())
    imax_magdip = int(np.abs(ref_data['magdip']).argmax())
    dip0_magdip = ref_data['magdip']

    eldip0 = ref_data['eldip_l']

    # Compute ders
    if args.d5points:
        compute_magdip_d2num_5p(nat, disp_stem_pattern, magdip_fname, delta, imax_magdip,
                                dip_dir=dip_dir, dip0_eq=dip0_magdip)
        if use_E:
            compute_eldip_d2num_E_5p(nat, disp_stem_pattern, eldip_fname, delta, imax_eldip,
                                     dip_dir=dip_dir, dip0_eq=eldip0)
        elif use_P:
            compute_eldip_d2num_P_5p(nat, disp_stem_pattern, ref0, eldip_fname, delta,
                                     debug_fname=debug_fname, dip_dir=dip_dir,
                                     refS1=ref, ref_basename=basename)
    else:
        compute_magdip_d2num(nat, disp_stem_pattern, magdip_fname, delta, imax_magdip,
                             dip_dir=dip_dir, dip0_eq=dip0_magdip)
        if use_E:
            compute_eldip_d2num_E(nat, disp_stem_pattern, eldip_fname, delta, imax_eldip,
                                  dip_dir=dip_dir, dip0_eq=eldip0)
        elif use_P:
            compute_eldip_d2num_P(nat, disp_stem_pattern, ref0, eldip_fname, delta,
                                  debug_fname=debug_fname, dip_dir=dip_dir,
                                  refS1=ref, ref_basename=basename)

    # Numerical first derivatives (for checking — commented out, not needed for production)
    # eldip1_fname  = eldip_fname.replace('_d2num', '') + '_d1num'
    # magdip1_fname = magdip_fname.replace('_d2num', '') + '_d1num'
    # compute_magdip_d1num(nat, disp_stem_pattern, magdip1_fname, delta, dip0_magdip,
    #                      imax_magdip, dip_dir=dip_dir)
    # if ref0 is not None and not use_E:
    #     compute_eldip_d1num_P(nat, disp_stem_pattern, ref0, eldip1_fname, delta,
    #                           dip_dir=dip_dir, refS1=ref, ref_basename=basename)
