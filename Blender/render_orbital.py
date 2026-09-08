#!/usr/bin/env -S blender --background --python
# -*- coding: utf-8 -*-
#
#  Blender Orbital Renderer — render_orbital.py
#  Headless CLI entry point for batch orbital figure generation.
#
#  Usage:
#    blender --background --python render_orbital.py -- \
#        --recipe recipe.json \
#        --style style.json \
#        --output output_dir/ \
#        --orbital 0 1 2 \
#        --resolution 1920 1080
#
#  License: GPLv3
#

"""
Headless batch renderer for the Orbital Renderer add-on.

Run from the command line:
  blender --background --python render_orbital.py -- --recipe recipe.json --style pub.json --output ./renders/

All arguments after '--' are passed to this script.
"""

import sys
import os
import json
import argparse
import traceback


def _parse_args():
    """Parse command line arguments. Blender eats the first few, so we handle '--'."""
    # Blender passes everything after '--' as regular args to the Python script.
    # The script name may be at argv[0] or earlier. Find where our args start.
    try:
        sep = sys.argv.index("--")
        our_args = sys.argv[sep + 1:]
    except ValueError:
        # No '--' separator; try to find the script name
        our_args = []
        for i, arg in enumerate(sys.argv):
            if arg.endswith("render_orbital.py") or "render_orbital" == os.path.basename(arg):
                our_args = sys.argv[i + 1:]
                break

    parser = argparse.ArgumentParser(
        description="Headless orbital renderer for Blender"
    )
    parser.add_argument(
        "--recipe", "-r", required=True,
        help="Path to render_recipe.json",
    )
    parser.add_argument(
        "--style", "-s",
        help="Path to orbital style JSON (falls back to defaults)",
    )
    parser.add_argument(
        "--output", "-o", required=True,
        help="Output directory for rendered images",
    )
    parser.add_argument(
        "--orbital", "-m", type=int, nargs="*",
        help="Orbital indices to render (0-based). Renders all if omitted.",
    )
    parser.add_argument(
        "--resolution", "-res", type=int, nargs=2, default=[1920, 1080],
        metavar=("W", "H"),
        help="Override output resolution (default: 1920 1080)",
    )
    parser.add_argument(
        "--prefix", "-p", default="orbital",
        help="Output filename prefix (default: 'orbital')",
    )

    return parser.parse_args(our_args)


def _setup_addon():
    """Enable the orbital_renderer add-on and import its modules."""
    import bpy

    # Ensure the parent directory of this script is on sys.path
    # so that 'orbital_renderer' package is importable.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)

    # Also add the parent of script_dir (in case the add-on is a symlink)
    parent_dir = os.path.dirname(script_dir)
    if parent_dir not in sys.path:
        sys.path.insert(0, parent_dir)

    # Try to enable the add-on (works if installed in Blender's add-on dir)
    addon_module = "orbital_renderer"
    if addon_module not in bpy.context.preferences.addons:
        try:
            bpy.ops.preferences.addon_enable(module=addon_module)
        except Exception:
            pass

    # Import our modules
    from orbital_renderer import recipe_loader
    from orbital_renderer import style_loader
    from orbital_renderer import isosurface as iso
    from orbital_renderer import molecule as molmod
    from orbital_renderer import materials as mats
    from orbital_renderer import render_setup
    from orbital_renderer.operators import (
        _clear_objects,
        _get_or_create_collection,
        _load_orbital_grid,
        _update_camera_lights,
    )

    return {
        "recipe_loader": recipe_loader,
        "style_loader": style_loader,
        "iso": iso,
        "molmod": molmod,
        "mats": mats,
        "render_setup": render_setup,
        "_clear_objects": _clear_objects,
        "_get_or_create_collection": _get_or_create_collection,
        "_load_orbital_grid": _load_orbital_grid,
        "_update_camera_lights": _update_camera_lights,
    }


def render_orbital(
    mods,
    recipe: dict,
    style: dict,
    orbital_index: int,
    output_path: str,
    resolution: tuple,
):
    """
    Render a single orbital to file.

    Parameters
    ----------
    mods : dict
        Module references from _setup_addon().
    recipe : dict
        Loaded recipe.
    style : dict
        Loaded style.
    orbital_index : int
        Index into recipe["orbitals"].
    output_path : str
        Full output file path (e.g. '/tmp/renders/mo_5.png').
    resolution : tuple[int, int]
        (width, height).
    """
    import bpy
    scene = bpy.context.scene

    # Clear previous scene
    mods["_clear_objects"](["OR_Atoms", "OR_Bonds", "OR_Orbitals"])
    mods["iso"].clear_grid_cache()

    # Populate scene properties for this orbital
    props = scene.orbital_renderer

    # Apply style to scene properties
    mods["style_loader"].apply_style_to_scene(props, style)

    # Override resolution from CLI
    props.render_resolution_x = resolution[0]
    props.render_resolution_y = resolution[1]

    # Create atoms and bonds from recipe
    atoms_coll = mods["_get_or_create_collection"]("OR_Atoms")
    bonds_coll = mods["_get_or_create_collection"]("OR_Bonds")

    for i, atom in enumerate(recipe.get("atoms", [])):
        pos = (atom["x"], atom["y"], atom["z"])
        obj = mods["molmod"].create_atom(
            f"Atom_{i}_{atom['symbol']}",
            atom["atomic_number"],
            pos,
            atoms_coll,
        )
        mat = mods["mats"].create_atom_material(
            f"OR_AtomMat_{atom['symbol']}_{i}",
            atom["atomic_number"],
        )
        obj.data.materials.append(mat)

    for i, (a_idx, b_idx) in enumerate(recipe.get("bonds", [])):
        atoms = recipe["atoms"]
        if a_idx < len(atoms) and b_idx < len(atoms):
            pos_a = (atoms[a_idx]["x"], atoms[a_idx]["y"], atoms[a_idx]["z"])
            pos_b = (atoms[b_idx]["x"], atoms[b_idx]["y"], atoms[b_idx]["z"])
            obj = mods["molmod"].create_bond(f"Bond_{i}", pos_a, pos_b, bonds_coll)
            if obj:
                mat = mods["mats"].create_bond_material(f"OR_BondMat_{i}")
                obj.data.materials.append(mat)

    # Load the single orbital
    orbitals = recipe.get("orbitals", [])
    if orbital_index >= len(orbitals):
        print(f"  Orbital index {orbital_index} out of range (max {len(orbitals)-1})")
        return

    orb_data = orbitals[orbital_index]

    # Populate orbital properties
    props.orbitals.clear()
    item = props.orbitals.add()
    item.cube_file = orb_data.get("_cube_path", "")
    item.mo_idx = orb_data.get("mo_idx", orbital_index)
    item.wtype = orb_data.get("wtype", "canonical")
    item.energy_eh = orb_data.get("energy_eh", 0.0)
    item.label = orb_data.get("label", orb_data.get("wtype", "canonical"))
    item.grid_spacing = orb_data.get("grid_spacing", 0.08)

    # Apply orbital style defaults, then override with recipe values
    mods["style_loader"].apply_style_to_orbital(item, style)
    recipe_isoval = orb_data.get("isovalue")
    if recipe_isoval is not None:
        item.isovalue = recipe_isoval

    # Load the orbital mesh
    props.active_orbital_index = 0
    mods["_load_orbital_grid"](bpy.context, 0)

    # Set molecule visibility from style
    molecule_style = style.get("molecule", {})
    if atoms_coll:
        atoms_coll.hide_render = not molecule_style.get("show_atoms", True)
    if bonds_coll:
        bonds_coll.hide_render = not molecule_style.get("show_bonds", True)

    # Render
    mods["render_setup"].setup_render_settings(
        scene,
        engine=props.render_engine,
        samples=props.render_samples,
        resolution_x=resolution[0],
        resolution_y=resolution[1],
        transparent=props.render_transparent,
    )
    scene.render.filepath = output_path
    bpy.ops.render.render(write_still=True)

    print(f"  → {output_path}")


def main():
    """Entry point for headless rendering."""
    args = _parse_args()
    print(f"Orbital Renderer — batch mode")
    print(f"  Recipe:  {args.recipe}")
    print(f"  Style:   {args.style or 'default'}")
    print(f"  Output:  {args.output}")
    print(f"  Orbitals: {args.orbital or 'all'}")

    # Resolve absolute paths
    recipe_path = os.path.abspath(args.recipe)
    output_dir = os.path.abspath(args.output)
    os.makedirs(output_dir, exist_ok=True)

    # Load recipe
    mods = _setup_addon()
    recipe = mods["recipe_loader"].load_recipe(recipe_path)

    # Load style
    style = mods["style_loader"].DEFAULT_STYLE
    if args.style:
        style_path = os.path.abspath(args.style)
        if os.path.isfile(style_path):
            style = mods["style_loader"].load_style(style_path)
        else:
            print(f"  Warning: style file not found: {style_path}")

    # Determine which orbitals to render
    num_orbitals = len(recipe.get("orbitals", []))
    if args.orbital:
        orbital_indices = [i for i in args.orbital if 0 <= i < num_orbitals]
    else:
        orbital_indices = list(range(num_orbitals))

    if not orbital_indices:
        print("  No orbitals to render.")
        return

    print(f"  Rendering {len(orbital_indices)} orbital(s)...")

    resolution = tuple(args.resolution)
    prefix = args.prefix

    for idx in orbital_indices:
        orb_data = recipe["orbitals"][idx]
        mo_idx = orb_data.get("mo_idx", idx)
        label = orb_data.get("label", orb_data.get("wtype", "canonical"))
        filename = f"{prefix}_mo{mo_idx}_{label}.png"
        output_path = os.path.join(output_dir, filename)

        print(f"  [{idx+1}/{len(orbital_indices)}] MO {mo_idx} ({label})")
        try:
            render_orbital(mods, recipe, style, idx, output_path, resolution)
        except Exception as e:
            print(f"    ERROR: {e}")
            traceback.print_exc()

    print("Done.")


if __name__ == "__main__":
    main()
