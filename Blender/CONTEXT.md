# CONTEXT.md — Blender Orbital Renderer

> **MIGRATED** — canonical copy at `~/Documents/Universidad/theoria/CONTEXT.md`.
> This file kept for reference; the add-on now lives in `sci-blender-addons/` monorepo.

Domain glossary for the Blender add-on and its companion visualizer.  
No implementation details — only terms and their canonical meanings.

---

## Projects

### Orbital Visualizer
The PyQt6/vispy desktop application (`Visualization/orbital-visualizer.py`).  
Used for fast interactive browsing of molecular orbitals from GAMESS `.log` files, selection of orbitals for publication, and export of volumetric data.

### Blender Orbital Renderer (this add-on)
The Blender add-on that consumes exported data and produces publication-quality rendered figures.  
Handles isosurface extraction, mesh construction, materials, lighting, camera, and rendering.

---

## Core concepts

### Molecular Orbital (MO)
A linear combination of basis functions with associated coefficients. Identified by a 0-based index (`mo_idx`) and a type (`wtype`: canonical or localized). Has an energy (Hartree) for canonical orbitals.

### Orbital Lobe
The positive or negative phase region of a molecular orbital, visualized as an isosurface.  
Each orbital produces two lobes: **positive** (marching cubes at `+isovalue`) and **negative** (at `-isovalue`).

### Isovalue
The scalar threshold for isosurface extraction. A single value applied symmetrically: the positive lobe extracts at `+isovalue`, the negative lobe at `-isovalue`. Stored per-orbital.

### Volumetric Grid
A 3D numpy array of MO values evaluated on a uniform grid over the molecular bounding box.  
Defined by grid spacing (Å), origin, and shape (nx, ny, nz).  
**Evaluated in the Orbital Visualizer, consumed in Blender.**

### Grid Spacing
The distance between adjacent grid points in Ångström. Chosen at export time in the visualizer.  
Coarse (0.20–0.35 Å) for preview, fine (0.05–0.08 Å) for final export.  
Cannot be changed in Blender — baked into the `.cube` file.

### Marching Cubes
The algorithm that extracts a triangle mesh from a volumetric grid at a given isovalue.  
**Re-run in Blender whenever the isovalue changes** — operates on the stored grid, produces a new mesh.

---

## File formats

### Gaussian Cube File (`.cube`)
Text-based volumetric data format. Carries atom positions + atomic numbers, grid axes + spacing, and a 3D array of scalar values.  
One cube file per orbital, containing the full signed grid.  
All coordinates in Ångström (documented in header comment).

### Render Recipe (`.json`)
JSON file produced by the Orbital Visualizer, consumed by the Blender add-on.  
Contains: source log path, atom data (symbol, atomic number, position), bond pairs, and per-orbital entries (cube file path, MO index, type, energy, suggested isovalue, grid spacing).

### GAMESS Log File (`.log`)
Output file from a GAMESS quantum chemistry calculation.  
Parsed by the Orbital Visualizer via cclib. Contains atomic coordinates, basis set definition, MO coefficients, orbital energies.  
**Not read by the Blender add-on.**

---

## Blender scene model

### Atom Object
A UV sphere mesh at the atom's position, colored by CPK convention, sized by covalent radius (×0.3).

### Bond Object
A cylinder mesh between two bonded atoms, fixed radius (0.12 Å), neutral gray material.

### Orbital Lobe Object
A mesh produced by marching cubes, with material controls (color, roughness, metallic, alpha).  
Named `MO_{idx}_{positive|negative}`.

### Scene Singleton
Only one recipe is loaded at a time. Loading a new recipe clears all previous orbital and molecule objects.

---

## Rendering

### Default Camera
Orthographic. Auto-positioned: pointing at molecule centroid from +Z, distance = 2× longest bounding box axis.  
User-adjustable in the Render tab.

### Default Lighting
Three-point studio setup (key + fill + rim lights), auto-positioned relative to camera.

### Default Render Engine
Cycles (path tracer) with transparent film, 300+ samples, denoising enabled.  
User can switch to EEVEE in the Render tab.
