# ADR 0002 — Marching cubes implementation

## Status
Accepted (2026-06-12)

## Context
Isosurface extraction from the volumetric grid requires a marching cubes algorithm. Two options:

1. **scikit-image** (`skimage.measure.marching_cubes`) — well-tested, C-accelerated, but requires installing scikit-image into Blender's bundled Python.
2. **Pure Python implementation** — zero dependencies, slower for large grids, but self-contained.

## Decision
Use the **pure Python marching cubes implementation** from [Molecular Blender](https://github.com/smparker/molecular-blender) (originally by Paul Bourke, adapted by Tom Sapiens, Robert Forsman, Shane Parker).

## Rationale
- Blender ships its own Python. Installing scikit-image into it requires user action (`pip install scikit-image`), varies across platforms (Windows needs admin, Linux needs matching ABI), and can break across Blender version upgrades.
- The pure Python implementation is ~200 lines, well-understood (classic Bourke algorithm with `edgetable` + `tritable` lookup tables), and fast enough for our use case.
- Isosurface extraction runs on isovalue change only, not every frame. For a 200³ grid, the pure Python implementation takes ~1-3 seconds — acceptable for an interactive slider.
- Zero-dependency deployment: the add-on ships as a single `.zip` with no external requirements beyond Blender's built-in Python + numpy.

## Consequences

### Pros
- **Zero dependencies.** Users install the add-on and it works immediately.
- **Portable across Blender versions.** No ABI concerns, no pip install scripts.
- **Simple to audit and maintain.** Small code footprint, well-documented algorithm.

### Cons
- **Slower than scikit-image.** Acceptable: 1–3 seconds per isovalue change on typical grids vs ~0.1s with scikit-image.
- **Not optimized for very large grids.** >300³ grids may be slow (>10s). Mitigation: typical publication grids are 150³–250³.

## Alternatives considered
- **scikit-image:** rejected due to deployment complexity in Blender's Python environment.
- **PyMCubes:** same deployment issue as scikit-image.
- **Cython-compiled marching cubes** (from Molecular Blender): would require Cython + C compiler in user's environment. Rejected in favor of the pure Python fallback (same codebase, different execution path).
