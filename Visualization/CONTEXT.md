# CONTEXT.md — orbital-visualizer

## Glossary

**Molecular Orbital (MO):** A one-electron wavefunction ψ(r) expressed as linear combination of atomic basis functions χ_i: ψ(r) = Σ c_i · χ_i(r).

**Canonical Orbital:** MO from a converged SCF calculation. Energy-ordered, delocalized over the molecule. Has orbital energies.

**Localized Orbital:** MO transformed via unitary rotation (Boys, Pipek-Mezey, or Edmiston-Ruedenberg) to maximize spatial locality. No meaningful energy ordering. Produced only for occupied orbitals.

**Basis Function (χ_i):** Contracted Gaussian-type orbital (CGTO): a sum of primitive Gaussians with fixed contraction coefficients. Characterized by center (atom position), angular momentum (S, P, D, F, G), exponent α, and contraction coefficient.

**Primitive Gaussian:** Single Cartesian Gaussian: (x-X)ˡ · (y-Y)ᵐ · (z-Z)ⁿ · exp(-α|r-R|²). Unnormalized in GAMESS output — normalization must be applied during evaluation.

**Isosurface:** 3D mesh extracted by marching cubes from the grid of ψ values at a chosen isovalue. Positive lobe (+isovalue) and negative lobe (−isovalue) are separate meshes.

**Isovalue:** The value of ψ at which the isosurface is drawn. Default ±0.05. Larger values = tighter surface (closer to maxima).

## Decided Design

### Scope
Interactive orbital viewer — equivalent to VMD's Graphical Representations → Orbital panel as a standalone application. Primary use: browsing orbitals to identify specific bonding/antibonding features, not publication figure generation.

### Input
- CLI argument: `python orbital-visualizer.py calculation.log`
- GUI: File → Open dialog
- GAMESS `.log` files only initially; adapter pattern for future formats

### Parsing strategy
- cclib for geometry, basis set (`data.gbasis`), and canonical MO coefficients (`data.mocoeffs`)
- Post-cclib scan of raw `.log` text for localized orbital markers and coefficient blocks
- Markers: `THE BOYS LOCALIZED ORBITALS ARE`, `THE PIPEK-MEZEY POPULATION LOCALIZED ORBITALS ARE`, `EDMISTON-RUEDENBERG ENERGY LOCALIZED ORBITALS`
- Future formats use adapter pattern: populate shared internal data structures (Molecule, BasisSet, Wavefunction)

### Basis function evaluation
- numba JIT-compiled kernel for production
- gbasis library used during development for validation only
- Normalization applied in-kernel (primitive normalization constants computed on-the-fly)
- Supports S, P, D, F, G shells (sufficient for cc-pVTZ and aug-cc-pVTZ)

### Memory & computation strategy
On-the-fly streaming evaluation (same approach as VMD). No precomputed basis-function grid. For each grid point, evaluate all basis functions and accumulate MO linear combination. Memory = output grid only (~few MB). Per-orbital latency traded for zero memory pressure. CPU-parallelized via numba across available cores.

### Grid
- Atom-centered bounding box with fixed padding (±4 Å)
- User-controlled grid spacing slider (default 0.10 Å)
- Coarse browsing: 0.35–0.40 Å recommended for orbital identification
- Fine detail: 0.10–0.15 Å for inspection

### Progressive refinement
Coarse grid computed first → displayed immediately → medium grid refines → fine grid refines. User sees usable result at ~0.3 seconds. Clicking another orbital cancels pending refinement.

### Isosurface extraction
Two-pass marching cubes via scikit-image: one for +isovalue mesh, one for −isovalue mesh. Separate MeshVisual objects with independent colors.

### Thumbnail gallery
- Occupied/virtual toggle (tabs)
- Static cached snapshots rendered once at 0.4 Å grid, stored as images
- Thumbnails show |ψ| only (absolute value, single surface, no phase sign)
- Clicking a thumbnail triggers progressive refinement in main viewport
- Full-resolution view shows both lobes with distinct colors

### Rendering & GUI
- vispy for OpenGL 3D rendering (meshes, atom spheres, bond cylinders)
- PyQt6 for window, menus, sliders, gallery, file dialog
- Ball-and-stick atom representation (distance-based bond detection, covalent radii × 1.2 cutoff)
- Atom colors: CPK convention

### Orbital numbering
As printed in the `.log` file (GAMESS output order). No reordering.

### Defaults
- Isovalue: ±0.05 (fixed, user adjusts via slider)
- Grid spacing: 0.10 Å (slider default)
- Material: transparent for orbital surfaces
- Two ColorIDs for positive/negative lobes

### Threading
Async coarse-first. Background thread runs grid computation. GUI stays responsive. No progress bar — refinement is visible as the orbital sharpens. Cancel on next orbital click.

### Export
File → Export Image saves current viewport as PNG at screen resolution.

### Validation during development
gbasis used to verify numba kernel correctness. Removed from final script.

## Out of Scope
- Publication-quality figure generation (use PyMOL/VMD for that)
- CUDA/GPU acceleration (no discrete GPU on target hardware)
- Multiple molecules simultaneously
- Animation/trajectory playback
- Normal mode visualization
- Electron density / ESP / ELF maps
- Formats other than GAMESS `.log` (adapter ready, not implemented)
