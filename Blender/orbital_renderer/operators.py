# -*- coding: utf-8 -*-
#
#  Blender Orbital Renderer — operators.py
#  Blender operators for recipe loading, isovalue updates, and rendering.
#
#  License: GPLv3
#

"""Operators for the Orbital Renderer add-on."""

import os
import bpy
from bpy.props import StringProperty, FloatProperty, IntProperty
from bpy_extras.io_utils import ImportHelper, ExportHelper

from . import recipe_loader
from . import isosurface as iso
from . import molecule as molmod
from . import materials as mats
from . import render_setup


# ---------------------------------------------------------------------------
# Collection helper
# ---------------------------------------------------------------------------

def _get_or_create_collection(name: str, parent_collection=None):
    """Get or create a Blender collection under the given parent."""
    if parent_collection is None:
        parent_collection = bpy.context.scene.collection

    coll = bpy.data.collections.get(name)
    if coll is None:
        coll = bpy.data.collections.new(name)
        parent_collection.children.link(coll)
    return coll


def _clear_objects(collection_names):
    """Remove all objects from named collections, then remove collections.
    Also cleans up orphaned mesh data from removed objects."""
    for name in collection_names:
        coll = bpy.data.collections.get(name)
        if coll:
            # Remove objects (their mesh data will be orphaned)
            for obj in list(coll.objects):
                mesh = obj.data if obj.type == "MESH" else None
                bpy.data.objects.remove(obj, do_unlink=True)
                # Clean up orphaned mesh from this object
                if mesh and mesh.users == 0:
                    try:
                        bpy.data.meshes.remove(mesh)
                    except Exception:
                        pass
            bpy.data.collections.remove(coll)


# ---------------------------------------------------------------------------
# Operator: Load Recipe
# ---------------------------------------------------------------------------

class OR_OT_LoadRecipe(bpy.types.Operator, ImportHelper):
    """Load a render recipe JSON file"""
    bl_idname = "orbital_renderer.load_recipe"
    bl_label = "Load Recipe"
    bl_options = {"REGISTER", "UNDO"}

    filename_ext = ".json"
    filter_glob: StringProperty(default="*.json", options={"HIDDEN"})

    def execute(self, context):
        scene = context.scene
        props = scene.orbital_renderer

        recipe_path = self.filepath

        try:
            recipe = recipe_loader.load_recipe(recipe_path)
        except Exception as e:
            self.report({"ERROR"}, f"Failed to load recipe: {e}")
            return {"CANCELLED"}

        # Store recipe path
        props.recipe_path = recipe_path

        # Clear previous scene
        _clear_objects(["OR_Atoms", "OR_Bonds", "OR_Orbitals"])
        iso.clear_grid_cache()

        # Create collections
        atoms_coll = _get_or_create_collection("OR_Atoms")
        bonds_coll = _get_or_create_collection("OR_Bonds")
        orbitals_coll = _get_or_create_collection("OR_Orbitals")

        # Create atoms
        for i, atom in enumerate(recipe.get("atoms", [])):
            pos = (atom["x"], atom["y"], atom["z"])
            obj = molmod.create_atom(
                f"Atom_{i}_{atom['symbol']}",
                atom["atomic_number"],
                pos,
                atoms_coll,
            )
            # Assign CPK material
            mat = mats.create_atom_material(
                f"OR_AtomMat_{atom['symbol']}_{i}",
                atom["atomic_number"],
            )
            obj.data.materials.append(mat)

        # Create bonds
        for i, (a_idx, b_idx) in enumerate(recipe.get("bonds", [])):
            atoms = recipe["atoms"]
            if a_idx < len(atoms) and b_idx < len(atoms):
                pos_a = (atoms[a_idx]["x"], atoms[a_idx]["y"], atoms[a_idx]["z"])
                pos_b = (atoms[b_idx]["x"], atoms[b_idx]["y"], atoms[b_idx]["z"])
                obj = molmod.create_bond(
                    f"Bond_{i}",
                    pos_a,
                    pos_b,
                    bonds_coll,
                )
                if obj:
                    mat = mats.create_bond_material(f"OR_BondMat_{i}")
                    obj.data.materials.append(mat)

        # Store orbital data as Blender property group items
        props.orbitals.clear()
        for orb_data in recipe.get("orbitals", []):
            item = props.orbitals.add()
            item.cube_file = orb_data.get("_cube_path", "")
            item.mo_idx = orb_data.get("mo_idx", 0)
            item.wtype = orb_data.get("wtype", "canonical")
            item.energy_eh = orb_data.get("energy_eh", 0.0)
            item.label = orb_data.get("label", orb_data.get("wtype", "canonical"))
            item.isovalue = orb_data.get("isovalue", 0.05)
            item.grid_spacing = orb_data.get("grid_spacing", 0.08)

        # Load first orbital if any
        if len(props.orbitals) > 0:
            props.active_orbital_index = 0
            _load_orbital_grid(context, 0)

        self.report({"INFO"}, f"Loaded recipe with {len(props.orbitals)} orbitals")
        return {"FINISHED"}


# ---------------------------------------------------------------------------
# Grid loading helper
# ---------------------------------------------------------------------------

def _load_orbital_grid(context, orbital_index: int):
    """
    Load the .cube file for an orbital at the given index
    and create/update its isosurface meshes.

    Stores the grid data in the isosurface module cache under the key
    "current_grid".
    """
    props = context.scene.orbital_renderer
    if orbital_index < 0 or orbital_index >= len(props.orbitals):
        return

    orbital = props.orbitals[orbital_index]
    cube_path = orbital.cube_file

    if not cube_path or not os.path.isfile(cube_path):
        return

    # Load cube file
    try:
        cube_data = recipe_loader.load_cube_file(cube_path)
    except Exception as e:
        print(f"Failed to load cube file: {e}")
        return

    grid = cube_data["grid"]
    axes = cube_data["axes"]
    origin = cube_data["origin"]

    # Store grid in module cache
    iso.store_grid("current_grid", grid)

    # Also store metadata for mesh building
    iso.store_grid("current_axes", axes)
    iso.store_grid("current_origin", origin)

    # Create or update positive lobe
    orbitals_coll = _get_or_create_collection("OR_Orbitals")
    pos_name = f"MO_{orbital.mo_idx}_positive"
    neg_name = f"MO_{orbital.mo_idx}_negative"

    isoval = orbital.isovalue

    # Remove old lobe objects
    for old_name in [pos_name, neg_name]:
        old_obj = bpy.data.objects.get(old_name)
        if old_obj:
            bpy.data.objects.remove(old_obj, do_unlink=True)

    # Create positive lobe
    pos_obj = iso.create_isosurface(
        grid, isoval, pos_name, axes, origin, orbitals_coll
    )
    if pos_obj:
        pos_obj.hide_viewport = not orbital.show_positive
        pos_obj.hide_render = not orbital.show_positive
        bpy.ops.object.select_all(action="DESELECT")
        pos_obj.select_set(True)
        bpy.context.view_layer.objects.active = pos_obj
        bpy.ops.object.shade_smooth()
        # Assign material
        mat = mats.create_orbital_material(
            pos_name + "_mat",
            orbital.pos_color,
            alpha=orbital.pos_alpha,
            material_style=orbital.material_style,
        )
        pos_obj.data.materials.append(mat)

    # Create negative lobe
    neg_obj = iso.create_isosurface(
        grid, -isoval, neg_name, axes, origin, orbitals_coll
    )
    if neg_obj:
        neg_obj.hide_viewport = not orbital.show_negative
        neg_obj.hide_render = not orbital.show_negative
        bpy.ops.object.select_all(action="DESELECT")
        neg_obj.select_set(True)
        bpy.context.view_layer.objects.active = neg_obj
        bpy.ops.object.shade_smooth()
        mat = mats.create_orbital_material(
            neg_name + "_mat",
            orbital.neg_color,
            alpha=orbital.neg_alpha,
            material_style=orbital.material_style,
        )
        neg_obj.data.materials.append(mat)

    # Setup camera and lights based on all visible objects
    _update_camera_lights(context)


def _update_camera_lights(context):
    """Auto-position camera and lights based on current scene objects."""
    scene = context.scene
    props = scene.orbital_renderer

    # Collect only our scene objects (atoms, bonds, orbitals)
    our_collections = ["OR_Atoms", "OR_Bonds", "OR_Orbitals"]
    all_objects = []
    for coll_name in our_collections:
        coll = bpy.data.collections.get(coll_name)
        if coll:
            all_objects.extend([
                obj for obj in coll.objects
                if obj.type in {"MESH", "SURFACE", "META", "CURVE"}
            ])

    if not all_objects:
        return

    # Setup camera
    render_setup.setup_camera(
        scene,
        all_objects,
        ortho=props.camera_ortho,
        distance_factor=props.camera_distance,
    )

    # Setup lights
    cam_obj = bpy.data.objects.get("OR_Camera")
    if cam_obj:
        render_setup.setup_lights(
            scene,
            cam_obj,
            key_intensity=props.light_key_intensity,
            fill_intensity=props.light_fill_intensity,
            rim_intensity=props.light_rim_intensity,
        )

    # Apply render settings
    render_setup.setup_render_settings(
        scene,
        engine=props.render_engine,
        samples=props.render_samples,
        resolution_x=props.render_resolution_x,
        resolution_y=props.render_resolution_y,
        transparent=props.render_transparent,
    )


# ---------------------------------------------------------------------------
# Operator: Update Isovalue
# ---------------------------------------------------------------------------

class OR_OT_UpdateIsovalue(bpy.types.Operator):
    """Re-run marching cubes at the current isovalue"""
    bl_idname = "orbital_renderer.update_isovalue"
    bl_label = "Update Isovalue"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        props = context.scene.orbital_renderer
        idx = props.active_orbital_index

        if idx < 0 or idx >= len(props.orbitals):
            return {"CANCELLED"}

        orbital = props.orbitals[idx]
        grid = iso.get_grid("current_grid")
        axes = iso.get_grid("current_axes")
        origin = iso.get_grid("current_origin")

        if grid is None:
            self.report({"WARNING"}, "No grid data loaded. Reload recipe.")
            return {"CANCELLED"}

        isoval = orbital.isovalue
        pos_name = f"MO_{orbital.mo_idx}_positive"
        neg_name = f"MO_{orbital.mo_idx}_negative"

        # Update positive lobe
        pos_obj = bpy.data.objects.get(pos_name)
        if orbital.show_positive:
            if pos_obj:
                ok = iso.update_isosurface(pos_obj, grid, isoval, axes, origin)
                if not ok:
                    pos_obj.hide_viewport = True
                    pos_obj.hide_render = True
                else:
                    pos_obj.hide_viewport = False
                    pos_obj.hide_render = False
                    bpy.ops.object.select_all(action="DESELECT")
                    pos_obj.select_set(True)
                    bpy.context.view_layer.objects.active = pos_obj
                    bpy.ops.object.shade_smooth()
            else:
                # Create it fresh
                orbitals_coll = _get_or_create_collection("OR_Orbitals")
                pos_obj = iso.create_isosurface(
                    grid, isoval, pos_name, axes, origin, orbitals_coll
                )
                if pos_obj:
                    bpy.ops.object.shade_smooth()
                    mat = mats.create_orbital_material(
                        pos_name + "_mat",
                        orbital.pos_color,
                        alpha=orbital.pos_alpha,
                        material_style=orbital.material_style,
                    )
                    pos_obj.data.materials.append(mat)
        elif pos_obj:
            pos_obj.hide_viewport = True
            pos_obj.hide_render = True

        # Update negative lobe
        neg_obj = bpy.data.objects.get(neg_name)
        if orbital.show_negative:
            if neg_obj:
                ok = iso.update_isosurface(neg_obj, grid, -isoval, axes, origin)
                if not ok:
                    neg_obj.hide_viewport = True
                    neg_obj.hide_render = True
                else:
                    neg_obj.hide_viewport = False
                    neg_obj.hide_render = False
                    bpy.ops.object.select_all(action="DESELECT")
                    neg_obj.select_set(True)
                    bpy.context.view_layer.objects.active = neg_obj
                    bpy.ops.object.shade_smooth()
            else:
                orbitals_coll = _get_or_create_collection("OR_Orbitals")
                neg_obj = iso.create_isosurface(
                    grid, -isoval, neg_name, axes, origin, orbitals_coll
                )
                if neg_obj:
                    bpy.ops.object.shade_smooth()
                    mat = mats.create_orbital_material(
                        neg_name + "_mat",
                        orbital.neg_color,
                        alpha=orbital.neg_alpha,
                        material_style=orbital.material_style,
                    )
                    neg_obj.data.materials.append(mat)
        elif neg_obj:
            neg_obj.hide_viewport = True
            neg_obj.hide_render = True

        return {"FINISHED"}


# ---------------------------------------------------------------------------
# Operator: Select Orbital
# ---------------------------------------------------------------------------

class OR_OT_SelectOrbital(bpy.types.Operator):
    """Switch to the selected orbital and load its grid"""
    bl_idname = "orbital_renderer.select_orbital"
    bl_label = "Select Orbital"
    bl_options = {"REGISTER", "UNDO"}

    index: IntProperty(default=-1)

    def execute(self, context):
        props = context.scene.orbital_renderer
        idx = self.index

        if idx < 0 or idx >= len(props.orbitals):
            return {"CANCELLED"}

        props.active_orbital_index = idx
        _load_orbital_grid(context, idx)

        return {"FINISHED"}


# ---------------------------------------------------------------------------
# Operator: Render
# ---------------------------------------------------------------------------

class OR_OT_Render(bpy.types.Operator):
    """Render the current scene to file"""
    bl_idname = "orbital_renderer.render"
    bl_label = "Render"
    bl_options = {"REGISTER"}

    filepath: StringProperty(
        name="Output",
        description="Output file path",
        subtype="FILE_PATH",
    )

    def execute(self, context):
        scene = context.scene
        props = scene.orbital_renderer

        # Set up render settings
        render_setup.setup_render_settings(
            scene,
            engine=props.render_engine,
            samples=props.render_samples,
            resolution_x=props.render_resolution_x,
            resolution_y=props.render_resolution_y,
            transparent=props.render_transparent,
        )

        # Set output path
        if self.filepath:
            scene.render.filepath = self.filepath
        elif not scene.render.filepath:
            scene.render.filepath = "//orbital_render.png"

        # Show/hide atoms and bonds
        atoms_coll = bpy.data.collections.get("OR_Atoms")
        bonds_coll = bpy.data.collections.get("OR_Bonds")
        if atoms_coll:
            atoms_coll.hide_render = not props.show_atoms
        if bonds_coll:
            bonds_coll.hide_render = not props.show_bonds

        # Render
        bpy.ops.render.render(write_still=True)

        self.report({"INFO"}, f"Rendered to {scene.render.filepath}")
        return {"FINISHED"}


class OR_OT_RenderDialog(bpy.types.Operator, ExportHelper):
    """Render the current scene (with file dialog)"""
    bl_idname = "orbital_renderer.render_dialog"
    bl_label = "Render to File"

    filename_ext = ".png"
    filter_glob: StringProperty(default="*.png;*.jpg;*.exr;*.tiff", options={"HIDDEN"})

    def execute(self, context):
        # Delegate to OR_OT_Render with the chosen filepath
        bpy.ops.orbital_renderer.render(filepath=self.filepath)
        return {"FINISHED"}


# ---------------------------------------------------------------------------
# Operator: Load Style
# ---------------------------------------------------------------------------

class OR_OT_LoadStyle(bpy.types.Operator, ImportHelper):
    """Load an orbital style JSON file and apply it to the scene"""
    bl_idname = "orbital_renderer.load_style"
    bl_label = "Load Style"
    bl_options = {"REGISTER", "UNDO"}

    filename_ext = ".json"
    filter_glob: StringProperty(default="*.json", options={"HIDDEN"})

    def execute(self, context):
        from . import style_loader

        try:
            style = style_loader.load_style(self.filepath)
        except Exception as e:
            self.report({"ERROR"}, f"Failed to load style: {e}")
            return {"CANCELLED"}

        props = context.scene.orbital_renderer

        # Apply to scene
        style_loader.apply_style_to_scene(props, style)

        # Apply orbital defaults to all orbitals
        for orbital in props.orbitals:
            style_loader.apply_style_to_orbital(orbital, style)

        # Re-apply camera and lights if we have objects
        _update_camera_lights(context)

        self.report({"INFO"}, f"Loaded style: {self.filepath}")
        return {"FINISHED"}


# ---------------------------------------------------------------------------
# Operator: Save Style
# ---------------------------------------------------------------------------

class OR_OT_SaveStyle(bpy.types.Operator, ExportHelper):
    """Save current render settings as an orbital style JSON file"""
    bl_idname = "orbital_renderer.save_style"
    bl_label = "Save Style"
    bl_options = {"REGISTER"}

    filename_ext = ".json"
    filter_glob: StringProperty(default="*.json", options={"HIDDEN"})

    def execute(self, context):
        from . import style_loader

        props = context.scene.orbital_renderer
        style = style_loader.collect_style_from_scene(props)
        style_loader.save_style(style, self.filepath)

        self.report({"INFO"}, f"Saved style to {self.filepath}")
        return {"FINISHED"}
