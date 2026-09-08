# ADR 0001 — Grid evaluation boundary

## Status
Accepted (2026-06-12)

## Context
The Orbital Visualizer (`orbital-visualizer.py`) evaluates molecular orbitals on a 3D grid using numba-accelerated Cartesian Gaussian basis function kernels. This requires numba + LLVM, which are not part of Blender's bundled Python.

The Blender add-on could either:
1. Parse the `.log` file directly inside Blender and evaluate MOs on the grid (requires numba in Blender)
2. Receive pre-computed volumetric grids from the visualizer and only perform isosurface extraction and rendering

## Decision
**All grid evaluation happens in the Orbital Visualizer.** The Blender add-on receives `.cube` files containing pre-computed volumetric data.

The Blender add-on stores the grid array in memory and re-runs marching cubes when the user changes isovalue. It never evaluates basis functions.

## Consequences

### Pros
- **No numba dependency in Blender.** The add-on requires only numpy (ships with Blender) and pure Python marching cubes (included). No pip install, no LLVM, no version conflicts.
- **Clean separation of concerns.** Visualizer = computation, Blender = presentation.
- **Cube files are universal.** Users of other QM codes (Gaussian, ORCA) can generate cube files independently and still use the add-on.
- **Faster Blender startup.** No log parsing overhead.

### Cons
- **Grid spacing baked at export time.** Cannot be changed in Blender without re-exporting from visualizer. Mitigation: export at fine spacing (0.05–0.08 Å) for final work.
- **Orbital selection baked at export time.** Cannot browse all MOs in Blender — only those exported. Mitigation: this is by design (focused publication workflow).
- **File size.** 200³ cube file ≈ 32 MB. For 5 orbitals ≈ 160 MB. Acceptable for publication projects.

## Alternatives considered
- **Direct `.log` parsing in Blender with numba:** rejected due to numba installation complexity in Blender's bundled Python (LLVM dependency, varies across Blender versions and platforms).
- **Pure numpy MO evaluation in Blender:** feasible but 3–5× slower than numba. Would make isovalue changes sluggish (10+ seconds per update on fine grids). Rejected in favor of pre-computation.
