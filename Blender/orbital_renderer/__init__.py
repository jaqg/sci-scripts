# -*- coding: utf-8 -*-
#
#  Blender Orbital Renderer — __init__.py
#  Add-on registration, property groups, UI panels, and orbital list.
#
#  This add-on incorporates code from:
#    - Molecular Blender (c) 2014-2025 Shane Parker, Joshua Szekely — GPLv3
#      marching cubes, cube file reader, periodic table data
#    - Beautiful Atoms / batoms (c) Xing Wang, Beautiful Atoms Team — GPLv3
#      material creation system, UI patterns, render settings patterns
#
#  License: GPLv3
#

"""
Blender Orbital Renderer — publication-quality molecular orbital visualization.

Two-tab UI panel in the 3D Viewport sidebar ("Orbitals"):
  - Orbitals tab: load recipe, browse orbitals, adjust isovalues and materials
  - Render tab: camera, lighting, render engine, and output settings
"""

import bpy
from bpy.props import (
    StringProperty,
    IntProperty,
    FloatProperty,
    BoolProperty,
    EnumProperty,
    FloatVectorProperty,
    CollectionProperty,
    PointerProperty,
)

# ---------------------------------------------------------------------------
# Blender metadata
# ---------------------------------------------------------------------------

bl_info = {
    "name": "Orbital Renderer",
    "author": "Orbital Renderer contributors",
    "version": (0, 1, 0),
    "blender": (4, 0, 0),
    "location": "View3D > Sidebar > Orbitals",
    "description": "Molecular orbital isosurface rendering from cube files and render recipes",
    "category": "3D View",
}


# ---------------------------------------------------------------------------
# Property Groups
# ---------------------------------------------------------------------------

class OR_OrbitalItem(bpy.types.PropertyGroup):
    """Per-orbital settings stored in the scene."""
    cube_file: StringProperty(
        name="Cube File",
        description="Absolute path to the cube file",
        default="",
    )
    mo_idx: IntProperty(name="MO Index", default=0)
    wtype: StringProperty(name="Type", default="canonical")
    energy_eh: FloatProperty(name="Energy (Eh)", default=0.0)
    label: StringProperty(name="Label", default="")

    isovalue: FloatProperty(
        name="Isovalue",
        description="Isosurface threshold",
        default=0.05,
        min=0.001,
        max=0.500,
        precision=4,
        update=lambda self, ctx: _on_isovalue_change(self, ctx),
    )
    grid_spacing: FloatProperty(name="Grid Spacing (Å)", default=0.08)

    pos_color: FloatVectorProperty(
        name="Positive Lobe Color",
        subtype="COLOR",
        default=(1.0, 0.2, 0.2, 1.0),
        size=4,
        update=lambda self, ctx: _on_lobe_color_change(self, ctx),
    )
    neg_color: FloatVectorProperty(
        name="Negative Lobe Color",
        subtype="COLOR",
        default=(0.2, 0.4, 1.0, 1.0),
        size=4,
        update=lambda self, ctx: _on_lobe_color_change(self, ctx),
    )
    pos_alpha: FloatProperty(
        name="Positive Alpha",
        default=0.6,
        min=0.0,
        max=1.0,
        update=lambda self, ctx: _on_lobe_color_change(self, ctx),
    )
    neg_alpha: FloatProperty(
        name="Negative Alpha",
        default=0.6,
        min=0.0,
        max=1.0,
        update=lambda self, ctx: _on_lobe_color_change(self, ctx),
    )

    material_style: EnumProperty(
        name="Material",
        items=[
            ("default",  "Default",  "Soft diffuse with slight metallic"),
            ("metallic", "Metallic", "Full metallic reflection"),
            ("ceramic",  "Ceramic",  "Smooth glossy"),
            ("plastic",  "Plastic",  "Fully rough diffuse"),
            ("mirror",   "Mirror",   "Near-perfect mirror"),
            ("glass",    "Glass",    "Transparent glass BSDF"),
            ("emission", "Emission", "Self-illuminating"),
        ],
        default="default",
        update=lambda self, ctx: _on_lobe_color_change(self, ctx),
    )
    subdivision: BoolProperty(
        name="Subdivision",
        description="Apply subdivision surface modifier",
        default=False,
        update=lambda self, ctx: _on_subdivision_change(self, ctx),
    )
    show_positive: BoolProperty(
        name="Show Positive",
        default=True,
        update=lambda self, ctx: _on_lobe_visibility_change(self, ctx),
    )
    show_negative: BoolProperty(
        name="Show Negative",
        default=True,
        update=lambda self, ctx: _on_lobe_visibility_change(self, ctx),
    )


class OR_SceneProperties(bpy.types.PropertyGroup):
    """Scene-level state for the Orbital Renderer."""
    recipe_path: StringProperty(
        name="Recipe Path",
        description="Path to the loaded render recipe",
        default="",
    )
    active_orbital_index: IntProperty(
        name="Active Orbital",
        default=0,
        update=lambda self, ctx: _on_orbital_select(self, ctx),
    )

    # Orbitals collection
    orbitals: CollectionProperty(type=OR_OrbitalItem)

    # Camera settings
    camera_ortho: BoolProperty(name="Orthographic", default=True)
    camera_distance: FloatProperty(
        name="Distance Factor",
        default=2.0,
        min=0.5,
        max=5.0,
    )

    # Light settings
    light_key_intensity: FloatProperty(
        name="Key Intensity", default=100.0, min=0.0, max=1000.0,
    )
    light_fill_intensity: FloatProperty(
        name="Fill Intensity", default=50.0, min=0.0, max=1000.0,
    )
    light_rim_intensity: FloatProperty(
        name="Rim Intensity", default=80.0, min=0.0, max=1000.0,
    )

    # Render settings
    render_engine: EnumProperty(
        name="Engine",
        items=[
            ("CYCLES", "Cycles", "Path-traced renderer"),
            ("BLENDER_EEVEE", "EEVEE", "Real-time renderer"),
        ],
        default="CYCLES",
    )
    render_samples: IntProperty(
        name="Samples", default=300, min=1, max=4096,
    )
    render_resolution_x: IntProperty(
        name="Width", default=1920, min=1, max=8192,
    )
    render_resolution_y: IntProperty(
        name="Height", default=1080, min=1, max=8192,
    )
    render_transparent: BoolProperty(name="Transparent BG", default=True)
    show_atoms: BoolProperty(name="Show Atoms", default=True)
    show_bonds: BoolProperty(name="Show Bonds", default=True)


# ---------------------------------------------------------------------------
# Property update callbacks
# ---------------------------------------------------------------------------

def _on_isovalue_change(self, context):
    """Called when isovalue slider changes — re-run marching cubes."""
    # Only update if this is the active orbital
    props = context.scene.orbital_renderer
    idx = props.active_orbital_index
    if 0 <= idx < len(props.orbitals) and props.orbitals[idx] == self:
        bpy.ops.orbital_renderer.update_isovalue()


def _on_lobe_color_change(self, context):
    """Called when lobe color/alpha/material changes."""
    from . import materials as mats

    props = context.scene.orbital_renderer
    idx = props.active_orbital_index
    if not (0 <= idx < len(props.orbitals) and props.orbitals[idx] == self):
        return

    mo_idx = self.mo_idx

    # Update by replacing material on objects
    for lobe_suffix, color, alpha in [
        ("positive", self.pos_color, self.pos_alpha),
        ("negative", self.neg_color, self.neg_alpha),
    ]:
        obj_name = f"MO_{mo_idx}_{lobe_suffix}"
        obj = bpy.data.objects.get(obj_name)
        if obj is None:
            continue

        mat_name = f"{obj_name}_mat"
        mat = mats.create_orbital_material(
            mat_name, color, alpha=alpha, material_style=self.material_style,
        )
        # Replace old material on the object
        if obj.data.materials:
            old_mat = obj.data.materials[0]
            obj.data.materials.clear()
            if old_mat and old_mat.users == 0:
                bpy.data.materials.remove(old_mat)
        obj.data.materials.append(mat)


def _on_lobe_visibility_change(self, context):
    """Called when show_positive/show_negative toggles."""
    props = context.scene.orbital_renderer
    idx = props.active_orbital_index
    if not (0 <= idx < len(props.orbitals) and props.orbitals[idx] == self):
        return

    mo_idx = self.mo_idx
    pos_obj = bpy.data.objects.get(f"MO_{mo_idx}_positive")
    neg_obj = bpy.data.objects.get(f"MO_{mo_idx}_negative")

    if pos_obj:
        pos_obj.hide_viewport = not self.show_positive
        pos_obj.hide_render = not self.show_positive
    if neg_obj:
        neg_obj.hide_viewport = not self.show_negative
        neg_obj.hide_render = not self.show_negative


def _on_orbital_select(self, context):
    """Called when active orbital index changes via UI list selection."""
    bpy.ops.orbital_renderer.select_orbital(index=self.active_orbital_index)


def _on_subdivision_change(self, context):
    """Called when subdivision checkbox toggles."""
    props = context.scene.orbital_renderer
    idx = props.active_orbital_index
    if not (0 <= idx < len(props.orbitals) and props.orbitals[idx] == self):
        return

    mo_idx = self.mo_idx
    for lobe_suffix in ("positive", "negative"):
        obj_name = f"MO_{mo_idx}_{lobe_suffix}"
        obj = bpy.data.objects.get(obj_name)
        if obj is None:
            continue

        # Find existing subdivision modifier
        mod = obj.modifiers.get("OR_Subdivision")
        if self.subdivision:
            if mod is None:
                mod = obj.modifiers.new("OR_Subdivision", "SUBSURF")
                mod.levels = 1
                mod.render_levels = 3
            mod.show_viewport = True
            mod.show_render = True
        else:
            if mod:
                obj.modifiers.remove(mod)


# ---------------------------------------------------------------------------
# UI List — Orbital List
# ---------------------------------------------------------------------------

class OR_UL_OrbitalList(bpy.types.UIList):
    """Draws the list of orbitals in the recipe."""
    bl_idname = "OR_UL_OrbitalList"

    def draw_item(self, context, layout, data, item, icon, active_data, active_propname):
        if self.layout_type in {"DEFAULT", "COMPACT"}:
            row = layout.row(align=True)
            row.label(text=f"MO {item.mo_idx}", icon="MESH_DATA")
            row.label(text=item.label)
            row.label(text=f"{item.energy_eh:+.4f} Eh")
            row.prop(item, "isovalue", text="")
        elif self.layout_type == "GRID":
            layout.alignment = "CENTER"
            layout.label(text=f"MO {item.mo_idx}")


# ---------------------------------------------------------------------------
# Panel: Orbitals Tab
# ---------------------------------------------------------------------------

class OR_PT_Orbitals(bpy.types.Panel):
    """Orbitals panel — load recipe, browse orbitals, adjust isovalues."""
    bl_label = "Orbitals"
    bl_idname = "OR_PT_Orbitals"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Orbitals"
    bl_order = 1

    def draw(self, context):
        layout = self.layout
        props = context.scene.orbital_renderer

        # Load Recipe button
        row = layout.row()
        row.operator("orbital_renderer.load_recipe", text="Load Recipe", icon="FILE_FOLDER")
        if props.recipe_path:
            row = layout.row()
            row.label(text=props.recipe_path, icon="FILE_TEXT")

        layout.separator()

        # Orbital list
        if len(props.orbitals) > 0:
            row = layout.row()
            row.template_list(
                "OR_UL_OrbitalList",
                "orbital_list",
                props,
                "orbitals",
                props,
                "active_orbital_index",
                rows=4,
            )

            layout.separator()

            # Active orbital controls
            idx = props.active_orbital_index
            if 0 <= idx < len(props.orbitals):
                orbital = props.orbitals[idx]

                box = layout.box()
                box.label(text=f"MO {orbital.mo_idx} — {orbital.label}", icon="ORBIT_DATA")
                box.label(text=f"Energy: {orbital.energy_eh:+.4f} Eh")
                box.label(text=f"Grid spacing: {orbital.grid_spacing:.3f} Å")

                # Isovalue slider
                box.prop(orbital, "isovalue", slider=True)

                # Positive lobe
                sub = box.box()
                sub.label(text="Positive Lobe (+)", icon="ADD")
                sub.prop(orbital, "show_positive", text="Visible")
                row = sub.row()
                row.prop(orbital, "pos_color", text="")
                sub.prop(orbital, "pos_alpha", slider=True)

                # Negative lobe
                sub = box.box()
                sub.label(text="Negative Lobe (−)", icon="REMOVE")
                sub.prop(orbital, "show_negative", text="Visible")
                row = sub.row()
                row.prop(orbital, "neg_color", text="")
                sub.prop(orbital, "neg_alpha", slider=True)

                # Material
                box.prop(orbital, "material_style")

                # Subdivision
                box.prop(orbital, "subdivision")

        else:
            layout.label(text="No recipe loaded", icon="INFO")


# ---------------------------------------------------------------------------
# Panel: Render Tab
# ---------------------------------------------------------------------------

class OR_PT_Render(bpy.types.Panel):
    """Render panel — camera, lights, render settings."""
    bl_label = "Render"
    bl_idname = "OR_PT_Render"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Orbitals"
    bl_order = 2

    def draw(self, context):
        layout = self.layout
        props = context.scene.orbital_renderer

        # Camera
        box = layout.box()
        box.label(text="Camera", icon="CAMERA_DATA")
        box.prop(props, "camera_ortho")
        box.prop(props, "camera_distance")

        # Lighting
        box = layout.box()
        box.label(text="Lighting", icon="LIGHT_AREA")
        box.prop(props, "light_key_intensity")
        box.prop(props, "light_fill_intensity")
        box.prop(props, "light_rim_intensity")

        # Render engine
        box = layout.box()
        box.label(text="Render Engine", icon="SCENE")
        box.prop(props, "render_engine")

        if props.render_engine == "CYCLES":
            box.prop(props, "render_samples")

        # Resolution
        box = layout.box()
        box.label(text="Resolution", icon="RENDER_RESULT")
        row = box.row(align=True)
        row.prop(props, "render_resolution_x", text="W")
        row.prop(props, "render_resolution_y", text="H")

        # Options
        box = layout.box()
        box.label(text="Options", icon="SETTINGS")
        box.prop(props, "render_transparent")
        box.prop(props, "show_atoms")
        box.prop(props, "show_bonds")

        # Render button
        layout.separator()
        row = layout.row(align=True)
        row.scale_y = 1.5
        row.operator("orbital_renderer.render_dialog", text="Render", icon="RENDER_STILL")

        # Style management
        layout.separator()
        box = layout.box()
        box.label(text="Style", icon="FILE_SCRIPT")
        row = box.row(align=True)
        row.operator("orbital_renderer.load_style", text="Load", icon="IMPORT")
        row.operator("orbital_renderer.save_style", text="Save", icon="EXPORT")


# ---------------------------------------------------------------------------
# Registration
# ---------------------------------------------------------------------------

classes = [
    OR_OrbitalItem,
    OR_SceneProperties,
    OR_UL_OrbitalList,
    OR_PT_Orbitals,
    OR_PT_Render,
]


def _import_operators():
    """Lazy import operator classes for registration."""
    from . import operators
    return [
        operators.OR_OT_LoadRecipe,
        operators.OR_OT_UpdateIsovalue,
        operators.OR_OT_SelectOrbital,
        operators.OR_OT_Render,
        operators.OR_OT_RenderDialog,
    ]


def _import_preferences():
    from . import preferences
    return [preferences.OR_AddonPreferences]


def register():
    """Register all classes and properties."""
    # Import and collect all classes
    from . import operators as _op
    from . import preferences as _pref

    all_classes = classes + [
        _op.OR_OT_LoadRecipe,
        _op.OR_OT_UpdateIsovalue,
        _op.OR_OT_SelectOrbital,
        _op.OR_OT_Render,
        _op.OR_OT_RenderDialog,
        _op.OR_OT_LoadStyle,
        _op.OR_OT_SaveStyle,
        _pref.OR_AddonPreferences,
    ]

    for cls in all_classes:
        bpy.utils.register_class(cls)

    # Register scene property
    bpy.types.Scene.orbital_renderer = PointerProperty(type=OR_SceneProperties)


def unregister():
    """Unregister all classes and properties."""
    # Unregister scene property
    del bpy.types.Scene.orbital_renderer

    from . import operators as _op
    from . import preferences as _pref

    all_classes = classes + [
        _op.OR_OT_LoadRecipe,
        _op.OR_OT_UpdateIsovalue,
        _op.OR_OT_SelectOrbital,
        _op.OR_OT_Render,
        _op.OR_OT_RenderDialog,
        _op.OR_OT_LoadStyle,
        _op.OR_OT_SaveStyle,
        _pref.OR_AddonPreferences,
    ]

    for cls in reversed(all_classes):
        bpy.utils.unregister_class(cls)


if __name__ == "__main__":
    register()
