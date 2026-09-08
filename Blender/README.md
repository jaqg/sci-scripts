# Blender Orbital Renderer — MIGRATED

> **Canonical repo:** `~/Documents/Universidad/theoria/` (add-on: `eidolon/`)
>
> This directory is a development workspace — the add-on, CLI script, and docs have moved.
> Remaining: `molecular-blender/` and `beautiful-atoms/` (third-party reference clones, GPLv3).

---

## Historical content below

# Blender Orbital Renderer

Blender add-on for high-quality molecular orbital rendering from GAMESS calculations, targeting publication/journal-frontmatter figures.

> **Status:** early development — first working milestone. Installation and usage subject to change.

---

## Installation

### Requirements

- **Blender 4.0 or later** (tested on 4.0+; may work on 3.6 with minor adjustments)
- **numpy** — ships with Blender, no extra install needed

### Install the add-on

**Option A — symlink (recommended for development):**

```bash
# Link the add-on into Blender's user scripts directory
ln -s /path/to/sci-scripts/Blender/orbital_renderer \
      ~/.config/blender/4.X/scripts/addons/orbital_renderer
```

**Option B — copy:**

```bash
cp -r /path/to/sci-scripts/Blender/orbital_renderer \
      ~/.config/blender/4.X/scripts/addons/
```

**Option C — install from ZIP (Blender 4.2+):**

1. Zip the `orbital_renderer/` folder
2. In Blender: Edit → Preferences → Add-ons → Install (top-right dropdown) → pick the `.zip`
3. Enable "Orbital Renderer" in the add-on list

### Enable the add-on

1. Open Blender → Edit → Preferences → Add-ons
2. Search "Orbital Renderer"
3. Check the box to enable
4. The "Orbitals" tab appears in 3D Viewport → Sidebar (press **N** to toggle sidebar)

---

## Quick start

### 1. Prepare input files

You need two things:

- **A `render_recipe.json`** — produced by the companion [orbital-visualizer.py](../Visualization/orbital-visualizer.py), or hand-written (see format below)
- **One or more `.cube` files** — volumetric grid data for each molecular orbital

A test recipe and cube file are included:

```
orbital_renderer/test_recipe.json          ← recipe (points to water.cube)
molecular-blender/examples/water.cube      ← example cube file (H₂O, MO 5, 32³ grid)
```

### 2. Load the recipe

1. Open Blender, go to **3D Viewport → Sidebar (N) → Orbitals tab**
2. Click **"Load Recipe"** — a file browser opens
3. Pick `orbital_renderer/test_recipe.json` (or your own recipe)

Blender creates:
- **Atom spheres** (CPK colors, ceramic material, sized by covalent radius × 0.3)
- **Bond cylinders** between bonded atom pairs (gray plastic)
- **Isosurface meshes** for the first orbital (positive lobe in red, negative in blue)
- **Camera** (orthographic, auto-framed from +Z direction)
- **3-point studio lights** (key + fill + rim)

### 3. Adjust the isosurface

- In the **Orbitals** tab, drag the **Isovalue** slider — the mesh updates in real time (1–3 s for typical 200³ grids)
- Toggle **Show Positive / Show Negative** to hide lobe meshes
- Change lobe colors with the color pickers, adjust **Alpha** for transparency
- Pick a **Material** style from the dropdown (default / metallic / ceramic / plastic / mirror / glass / emission)
- Enable **Subdivision** for smoother meshes at render time

### 4. Switch orbitals

If your recipe has multiple orbitals, click a row in the **orbital list** to switch. The previous grid is unloaded from memory and the new `.cube` file is read.

### 5. Render

1. Switch to the **Render** tab
2. Choose **render engine** (Cycles for publication quality, EEVEE for fast preview)
3. Set **resolution**, toggle **transparent background**, adjust **samples**
4. Control camera (ortho/perspective toggle, distance factor)
5. Adjust **light intensities** (key / fill / rim)
6. Toggle **Show Atoms / Show Bonds** to include/exclude the molecular model
7. Click **"Render"** → pick output file → Blender renders a still frame

### Render recipe format

```json
{
  "version": 1,
  "source_log": "/path/to/calc.log",
  "atoms": [
    {"symbol": "O", "atomic_number": 8, "x": 0.0, "y": 0.0, "z": 0.0},
    {"symbol": "H", "atomic_number": 1, "x": -0.757, "y": 0.586, "z": 0.0}
  ],
  "bonds": [[0, 1]],
  "orbitals": [
    {
      "cube_file": "mo_5.cube",
      "mo_idx": 5,
      "wtype": "canonical",
      "energy_eh": -0.3421,
      "label": "canonical",
      "isovalue": 0.05,
      "grid_spacing": 0.08
    }
  ]
}
```

- `cube_file` paths are **relative to the recipe JSON directory**
- Atom positions in **Ångström**
- `bonds` is a list of `[atom_index_a, atom_index_b]` pairs (0-based)
- One orbital entry per `.cube` file

---

## CLI — headless batch rendering

Use `render_orbital.py` to produce N figures from the command line without opening Blender's GUI. Combine with a style file for identical look across all figures.

```bash
# Render all orbitals from a recipe with a style file
blender --background --python render_orbital.py -- \
    --recipe recipe.json \
    --style pub_style.json \
    --output ./renders/ \
    --prefix figure

# Render specific orbitals
blender --background --python render_orbital.py -- \
    --recipe recipe.json \
    --style pub_style.json \
    --output ./renders/ \
    --orbital 0 2 5

# Override resolution
blender --background --python render_orbital.py -- \
    --recipe recipe.json \
    --output ./renders/ \
    --resolution 3840 2160
```

### CLI arguments

| Flag | Description |
|------|-------------|
| `--recipe`, `-r` | Path to `render_recipe.json` (required) |
| `--style`, `-s` | Path to orbital style JSON (falls back to defaults) |
| `--output`, `-o` | Output directory (required) |
| `--orbital`, `-m` | Space-separated orbital indices to render (0-based). Renders all if omitted. |
| `--resolution`, `-res` | Width Height (default: 1920 1080) |
| `--prefix`, `-p` | Output filename prefix (default: `orbital`) |

Output files are named `{prefix}_mo{N}_{type}.png`.

### Style files

A style file (`orbital_style.json`) captures all render decisions so you can reproduce the same look across figures. Like matplotlib's `mplstyle` sheets.

**From the GUI:** Render tab → Style box → **Save** (exports current settings) / **Load** (applies a style file to the scene).

**Style file format** (see `orbital_renderer/example_style.json`):

```json
{
  "version": 1,
  "camera":      { "ortho": true, "distance_factor": 2.0 },
  "lighting":    { "key_intensity": 100.0, "fill_intensity": 50.0, "rim_intensity": 80.0 },
  "render":      { "engine": "CYCLES", "samples": 300, "resolution_x": 1920, "resolution_y": 1080, "transparent": true },
  "orbitals":    { "default_isovalue": 0.05, "positive_color": [1.0,0.2,0.2], "negative_color": [0.2,0.4,1.0], "alpha": 0.6, "material_style": "default", "subdivision": false },
  "molecule":    { "show_atoms": true, "show_bonds": true }
}
```

Missing keys fall back to the built-in defaults.

---

## Architecture

```
orbital_renderer/
├── __init__.py          # Add-on registration, panels, PropertyGroups, UIList
├── recipe_loader.py     # Reads render_recipe.json + Gaussian cube files
├── marching_cubes.py    # Pure Python marching cubes (edgetable, tritable, polygonise)
├── isosurface.py        # Builds Blender mesh from grid + isovalue
├── molecule.py          # Atom spheres + bond cylinders
├── materials.py         # Principled BSDF presets + CPK color table
├── render_setup.py      # Camera + lights auto-placement
├── operators.py         # Load recipe, adjust isovalue, switch orbital, render
└── preferences.py       # Add-on preferences
```

---

## License

**GNU General Public License v3.0 (GPLv3).**  
This add-on incorporates code from the following GPLv3 projects.

---

## Reused modules & attribution

### 1. Marching cubes algorithm

**Source:** [Molecular Blender](https://github.com/smparker/molecular-blender) — `molecular_blender/marching_cube_py.py`

**Authors:** Paul Bourke (original algorithm), Tom Sapiens (Blender adaptation), Robert Forsman (optimizations), Shane Parker (molecular blender integration)

**License:** GPLv3

**What we take:** `edgetable`, `tritable`, `polygonise()`, `marching_cube_box()`, `marching_cube_outline()` — the full pure-Python marching cubes implementation with zero external dependencies.

**Why:** This removes the need for `scikit-image` inside Blender's bundled Python. The pure-Python fallback is fast enough for publication workflows (computed once per isovalue change, not real-time).

---

### 2. Cube file reader

**Source:** [Molecular Blender](https://github.com/smparker/molecular-blender) — `molecular_blender/importers.py` (function `molecule_from_cube`)

**Author:** Shane Parker, Joshua Szekely

**License:** GPLv3

**What we take:** The Gaussian cube file parser: comment lines, atom count, grid axes (Bohr↔Å conversion), atom positions, volumetric data block (z-y-x inner-loop ordering).

**Why:** Standard cube format reader tested across many quantum chemistry codes.

---

### 3. Periodic table data

**Source:** [Molecular Blender](https://github.com/smparker/molecular-blender) — `molecular_blender/periodictable.py`

**Author:** Shane Parker, Joshua Szekely

**License:** GPLv3

**What we take:** `Element` class and `elements` dictionary containing van der Waals radii, covalent radii, atomic masses, element names, and symbols for elements 1–109.

**Why:** Needed for bond detection (covalent radii) and atom sphere sizing (VDW radii). Complements the CPK color table already present in our companion `orbital-visualizer.py`.

---

### 4. Material creation system

**Source:** [Beautiful Atoms (batoms)](https://github.com/beautiful-atoms/beautiful-atoms) — `batoms/material/__init__.py`

**Author:** Xing Wang, Beautiful Atoms Team

**License:** GPLv3

**What we take:** `create_material()` function and `material_styles_dict` presets (default, metallic, ceramic, plastic, mirror, glass, emission). Sets up Blender's Principled BSDF node tree with configurable roughness, metallic, IOR, alpha, and optional vertex-color or attribute-driven coloring.

**Why:** Provides publication-ready material presets out of the box. We extend with a "translucent" style for semi-transparent orbital lobes.

---

### 5. UI patterns (inspiration, not copy-paste)

**Source:** [Beautiful Atoms (batoms)](https://github.com/beautiful-atoms/beautiful-atoms) — `batoms/plugins/isosurface/`

**Author:** Xing Wang, Beautiful Atoms Team

**License:** GPLv3

**What we reference:** The `PropertyGroup` + `UIList` + `Panel` pattern for managing multiple isosurfaces with per-surface isovalue, color, and material controls.

**Why:** This is the standard Blender add-on pattern for lists of configurable items. We adapt it for orbitals (positive/negative lobes per MO).

---

### 6. Render settings (inspiration)

**Source:** [Beautiful Atoms (batoms)](https://github.com/beautiful-atoms/beautiful-atoms) — `batoms/render/render.py`, `camera.py`, `light.py`

**Author:** Xing Wang, Beautiful Atoms Team

**License:** GPLv3

**What we reference:** Engine selection (Cycles/EEVEE), transparent film, resolution, sample count, camera setup (ortho/perspective), 3-point lighting auto-positioned relative to camera.

**Why:** These are the essential Blender-scene settings for producing publication-quality renders.

---

## Dependencies of this add-on

| Dependency | Needed for | Bundled? |
|-----------|-----------|----------|
| numpy | Array operations, mesh construction | Ships with Blender |
| (none else) | Marching cubes, cube reading, materials all use Blender's built-in `bpy` + `mathutils` + numpy | — |

**Note:** No `scikit-image`, no `numba`, no Cython compilation required inside Blender. The marching cubes algorithm is included directly. Grid evaluation happens externally in the companion `orbital-visualizer.py` script, which writes `.cube` files.

---

## Companion tool

**[orbital-visualizer.py](../Visualization/orbital-visualizer.py)** — PyQt6/vispy desktop application for fast orbital browsing, grid evaluation (numba-accelerated), and `.cube` file export. Not part of this Blender add-on; the add-on consumes its exported `.cube` files + `render_recipe.json`.

---

## Full attribution block (for add-on `__init__.py`)

```
This add-on incorporates code from:
  - Molecular Blender (c) 2014-2025 Shane Parker, Joshua Szekely — GPLv3
    marching cubes, cube file reader, periodic table data
  - Beautiful Atoms / batoms (c) Xing Wang, Beautiful Atoms Team — GPLv3
    material creation system, UI patterns, render settings patterns
```
