# -*- coding: utf-8 -*-
#
#  Blender Orbital Renderer — preferences.py
#  Add-on preferences panel.
#
#  License: GPLv3
#

"""Add-on preferences for the Orbital Renderer."""

import bpy
from bpy.props import StringProperty, FloatProperty, IntProperty


class OR_AddonPreferences(bpy.types.AddonPreferences):
    """Orbital Renderer preferences"""
    bl_idname = __package__

    default_recipe_dir: StringProperty(
        name="Default Recipe Directory",
        description="Default path for loading render recipe JSON files",
        subtype="DIR_PATH",
        default="",
    )

    default_render_dir: StringProperty(
        name="Default Render Output",
        description="Default path for rendered images",
        subtype="DIR_PATH",
        default="",
    )

    def draw(self, context):
        layout = self.layout
        layout.prop(self, "default_recipe_dir")
        layout.prop(self, "default_render_dir")
