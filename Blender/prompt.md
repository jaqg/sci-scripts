# Session prompt — Blender Orbital Renderer (coding start)

This session builds the first working version of the Blender add-on.

## Your mission

Build a Blender add-on (`Blender/orbital_renderer/`) that:

1. Reads a `render_recipe.json` produced by the companion visualizer
2. Loads the referenced `.cube` files
3. Creates orbital isosurface meshes (positive + negative lobes) with live isovalue sliders
4. Creates atom spheres + bond cylinders from recipe data
5. Sets up camera + lights for publication-quality rendering
6. Provides a two-tab UI panel in Blender's sidebar

## Premise

Read these files for full context — they are **authoritative**:

- `Blender/CONTEXT.md` — glossary of every canonical term (workflow, file formats, scene model, rendering defaults)
- `Blender/docs/adr/0001-grid-evaluation-boundary.md` — why grid evaluation stays in the visualizer, not Blender
- `Blender/docs/adr/0002-marching-cubes-implementation.md` — why pure Python marching cubes over skimage
- `Blender/README.md` — attribution to Molecular Blender and batoms (GPLv3 sources)

**Critical rule:** `CONTEXT.md` is the glossary. Use its terms exactly. If you need a new term, add it to CONTEXT.md.

## Architecture

```
Blender/orbital_renderer/
├── __init__.py              # Add-on registration, panel classes
├── recipe_loader.py         # Reads render_recipe.json + .cube files
├── marching_cubes.py        # Pure Python marching cubes (from Molecular Blender)
├── isosurface.py            # Builds Blender mesh from grid + isovalue
├── molecule.py              # Atom spheres + bond cylinders
├── materials.py             # Principled BSDF with presets (from batoms)
├── render_setup.py          # Camera + lights auto-placement
├── operators.py             # Load recipe, adjust isovalue, render
└── preferences.py           # Add-on preferences (file paths, defaults)
```

## Reusable code (already cloned in Blender/)

| Need | Source | File |
|------|--------|------|
| Marching cubes | Molecular Blender | `Blender/molecular-blender/molecular_blender/marching_cube_py.py` |
| Cube file reader | Molecular Blender | `Blender/molecular-blender/molecular_blender/importers.py` (`molecule_from_cube`) |
| Periodic table | Molecular Blender | `Blender/molecular-blender/molecular_blender/periodictable.py` |
| Material creation | batoms | `Blender/beautiful-atoms/batoms/material/__init__.py` |
| UI patterns | batoms | `Blender/beautiful-atoms/batoms/plugins/isosurface/` |

**License:** GPLv3. Include the attribution header from README.md in `__init__.py`.

## Render recipe format

```json
{
  "version": 1,
  "source_log": "/absolute/path/to/calc.log",
  "atoms": [
    {"symbol": "C", "atomic_number": 6, "x": 0.0, "y": 0.0, "z": 0.0}
  ],
  "bonds": [[0, 1], [1, 2]],
  "orbitals": [
    {
      "cube_file": "mo_15.cube",
      "mo_idx": 15,
      "wtype": "canonical",
      "energy_eh": -0.3421,
      "label": "canonical",
      "isovalue": 0.05,
      "grid_spacing": 0.08
    }
  ]
}
```

`cube_file` paths are relative to the recipe JSON directory.

## Cube file format specifics

- Ångström units throughout (documented in header comment line 2)
- `importers.py` → `molecule_from_cube()` handles the format: comment lines, atom count, grid axes, atom positions, volumetric data
- Axes: positive `npoints` means Bohr (multiply by `bohr2ang`), negative `npoints` means Å (use absolute value). Handle both.
- Volumetric data: 6 floats per line, inner loop order: z → y → x
- One cube file = one orbital, full signed grid

## Blender scene model (from CONTEXT.md)

### On recipe load:
1. **Clear previous scene objects** (single recipe at a time)
2. **Create atoms** from recipe: UV spheres at positions, radius = 0.3 × covalent radius (from periodic table), CPK colors, ceramic material
3. **Create bonds** from recipe: cylinders between atom pairs, fixed radius 0.12 Å, neutral gray material
4. **For each orbital in recipe:**
   - Read `.cube` file → store 3D numpy grid array in memory (only current orbital keeps grid; others unload)
   - Run marching cubes at `+isovalue` → positive lobe mesh
   - Run marching cubes at `-isovalue` → negative lobe mesh
   - Create Blender mesh objects: `MO_{idx}_positive`, `MO_{idx}_negative`
   - Apply `shade_smooth`
   - Assign default materials (red for positive, blue for negative, translucent 0.6 alpha)
5. **Create camera** (orthographic, auto-framed: centroid +Z, distance = 2× longest axis)
6. **Create lights** (3-point studio: key + fill + rim, auto-positioned relative to camera)

### When user changes isovalue (slider in UI):
- Re-run marching cubes on the stored grid at new isovalue
- Update positive and negative lobe meshes
- Fast: 1-3s for 200³ grid

### When user switches orbital (via list):
- Unload current grid from memory
- Read new `.cube` file → store grid
- Run marching cubes at current isovalue
- Update meshes

## UI — two tabs (3D Viewport sidebar, "Orbitals" panel)

### Tab 1: Orbitals
- **"Load Recipe" button** — file browser, picks `render_recipe.json`
- **Orbital list** (UIList): each row shows MO index, type, energy, isovalue spinner
- **Per-orbital controls** (appear when an orbital is selected):
  - Isovalue slider (float, 0.001–0.500)
  - Positive lobe: color picker (default: red), alpha slider
  - Negative lobe: color picker (default: blue), alpha slider
  - Material style dropdown (default, metallic, ceramic, plastic, glass, emission — from batoms presets)
  - Subdivision checkbox (off by default)
  - Show/Hide checkbox per lobe

### Tab 2: Render
- **Camera**: ortho/perspective toggle, distance, lens/ortho_scale
- **Lighting**: intensity sliders for key/fill/rim, toggle visibility
- **Render engine**: Cycles/EEVEE dropdown
- **Samples** (Cycles only)
- **Resolution**: width × height
- **Transparent background** checkbox (default: on)
- **Show/Hide atoms** checkbox
- **Show/Hide bonds** checkbox
- **"Render" button** → opens render dialog or renders to file

## Materials (from batoms, adapted)

Default presets to expose via dropdown:
```python
material_styles = {
    "default":   {"Metallic": 0.10, "Roughness": 0.20, "IOR": 1.4},
    "metallic":  {"Metallic": 1.00, "Roughness": 0.20, "IOR": 1.4},
    "ceramic":   {"Metallic": 0.02, "Roughness": 0.00, "IOR": 1.4},
    "plastic":   {"Metallic": 0.00, "Roughness": 1.00, "IOR": 1.4},
    "mirror":    {"Metallic": 0.99, "Roughness": 0.01, "IOR": 2.0},
    "glass":     {"type": "Glass BSDF", "Roughness": 0.5, "IOR": 1.45},
}
```

For orbital lobes, wrap a Principled BSDF with alpha control. `blend_method = "BLEND"`, `show_transparent_back = False`.

## CPK colors

Use the same table as `Visualization/orbital-visualizer.py` (lines ~200-230: `CPK_COLORS` dict, keyed by atomic number). Duplicate into this add-on — keep Blender add-on self-contained.

## Periodic table

Extract from `molecular-blender/molecular_blender/periodictable.py`. Need: element symbol → covalent radius, van der Waals radius. Already keyed by lowercase symbol.

## Blender version

Target Blender 4.0+ (uses `blender_manifest.toml` for extensions, if 4.2+). Use `bpy` API compatible with 4.0+.

## Dependencies (inside Blender)

- `numpy` — ships with Blender
- `mathutils` — ships with Blender
- Nothing else. No pip install needed.

## First milestone

A minimal working add-on that:
1. Loads a recipe → atoms + bonds + one orbital appear with meshes
2. Isovalue slider updates the meshes
3. UI has both tabs (wired up but minimal)

Test with a manually-created `render_recipe.json` pointing to a sample `.cube` file. The `Blender/molecular-blender/examples/water.cube` file can serve as test data.

## File structure

Working directory: `/home/jose/Documents/GitHub/sci-scripts/Blender/`

```
Blender/
├── CONTEXT.md                           # Glossary (read first)
├── README.md                            # Attribution
├── prompt.md                            # This file
├── docs/adr/
│   ├── 0001-grid-evaluation-boundary.md
│   └── 0002-marching-cubes-implementation.md
├── orbital_renderer/                    # ← Build this
│   ├── __init__.py
│   ├── recipe_loader.py
│   ├── marching_cubes.py
│   ├── isosurface.py
│   ├── molecule.py
│   ├── materials.py
│   ├── render_setup.py
│   ├── operators.py
│   └── preferences.py
├── molecular-blender/                   # Reference only (GPLv3 source)
└── beautiful-atoms/                     # Reference only (GPLv3 source)
```

Do NOT modify files in `molecular-blender/` or `beautiful-atoms/`. Extract code into `orbital_renderer/` with attribution comments.

## Style

- Blender add-on conventions: register/unregister, PropertyGroup, operators, panels
- All strings user-visible in the UI
- Type hints encouraged
- Attribution header in every `.py` file per `README.md`
