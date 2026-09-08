# -*- coding: utf-8 -*-
#
#  Blender Orbital Renderer — render_setup.py
#  Camera, lighting, and render settings for publication-quality output.
#
#  Patterns adapted from Beautiful Atoms / batoms:
#    (c) Xing Wang, Beautiful Atoms Team — GPLv3
#
#  License: GPLv3
#

"""Setup camera, lights, and render settings."""

import math
import bpy
import mathutils
import numpy as np


# ---------------------------------------------------------------------------
# Camera
# ---------------------------------------------------------------------------

def compute_centroid_and_bbox(objects):
    """
    Compute centroid and bounding box of a list of Blender objects.

    Parameters
    ----------
    objects : list of bpy.types.Object

    Returns
    -------
    centroid : mathutils.Vector
    size : mathutils.Vector (max extents)
    """
    if not objects:
        return mathutils.Vector((0, 0, 0)), mathutils.Vector((5, 5, 5))

    min_corner = mathutils.Vector((float("inf"),) * 3)
    max_corner = mathutils.Vector((float("-inf"),) * 3)

    for obj in objects:
        for corner in obj.bound_box:
            world_corner = obj.matrix_world @ mathutils.Vector(corner)
            min_corner.x = min(min_corner.x, world_corner.x)
            min_corner.y = min(min_corner.y, world_corner.y)
            min_corner.z = min(min_corner.z, world_corner.z)
            max_corner.x = max(max_corner.x, world_corner.x)
            max_corner.y = max(max_corner.y, world_corner.y)
            max_corner.z = max(max_corner.z, world_corner.z)

    if min_corner.x == float("inf"):
        return mathutils.Vector((0, 0, 0)), mathutils.Vector((5, 5, 5))

    centroid = (min_corner + max_corner) / 2.0
    size = max_corner - min_corner

    return centroid, size


def setup_camera(
    scene,
    objects,
    camera_name: str = "OR_Camera",
    ortho: bool = True,
    distance_factor: float = 2.0,
):
    """
    Create or update camera, auto-framing the given objects.

    Parameters
    ----------
    scene : bpy.types.Scene
    objects : list of bpy.types.Object
        Objects to frame.
    camera_name : str
    ortho : bool
        True for orthographic, False for perspective.
    distance_factor : float
        Camera distance = distance_factor × longest bounding box axis.

    Returns
    -------
    bpy.types.Object
        The camera object.
    """
    centroid, size = compute_centroid_and_bbox(objects)
    longest_axis = max(size.x, size.y, size.z)
    if longest_axis < 0.01:
        longest_axis = 5.0

    distance = longest_axis * distance_factor

    # Position camera at centroid + Z
    camera_pos = centroid + mathutils.Vector((0, 0, distance))

    # Look-at direction (camera points toward -Z in local space)
    look_dir = (centroid - camera_pos).normalized()

    # Create or reuse camera
    cam_obj = bpy.data.objects.get(camera_name)
    if cam_obj is None:
        cam_data = bpy.data.cameras.new(camera_name)
        cam_obj = bpy.data.objects.new(camera_name, cam_data)
        scene.collection.objects.link(cam_obj)

    cam_obj.location = camera_pos
    cam_obj.data.type = "ORTHO" if ortho else "PERSP"

    # Point camera at centroid
    direction = centroid - cam_obj.location
    rot_quat = direction.to_track_quat("-Z", "Y")
    cam_obj.rotation_euler = rot_quat.to_euler()

    # Set ortho scale to frame the object
    if ortho:
        cam_obj.data.ortho_scale = longest_axis * 1.5
    else:
        # Adjust lens or FOV
        cam_obj.data.lens = 50.0

    # Make active
    scene.camera = cam_obj

    return cam_obj


# ---------------------------------------------------------------------------
# Lighting (3-point studio)
# ---------------------------------------------------------------------------

def setup_lights(
    scene,
    camera_obj,
    key_intensity: float = 100.0,
    fill_intensity: float = 50.0,
    rim_intensity: float = 80.0,
):
    """
    Create or update a 3-point studio light setup.

    Parameters
    ----------
    scene : bpy.types.Scene
    camera_obj : bpy.types.Object
        Camera to position lights relative to.
    key_intensity : float
    fill_intensity : float
    rim_intensity : float

    Returns
    -------
    tuple of (key, fill, rim) light objects
    """
    camera_dir = mathutils.Vector((0, 0, -1))
    camera_dir.rotate(camera_obj.rotation_euler)

    camera_right = mathutils.Vector((1, 0, 0))
    camera_right.rotate(camera_obj.rotation_euler)

    camera_up = mathutils.Vector((0, 1, 0))
    camera_up.rotate(camera_obj.rotation_euler)

    cam_pos = camera_obj.location.copy()

    # Key light: 45° to the right, 30° up
    key_dir = (camera_dir + camera_right * 0.7 + camera_up * 0.5).normalized()
    key_pos = cam_pos + key_dir * 10.0

    # Fill light: 45° to the left, 15° up, half intensity
    fill_dir = (camera_dir - camera_right * 0.7 + camera_up * 0.3).normalized()
    fill_pos = cam_pos + fill_dir * 10.0

    # Rim light: behind the subject, above
    rim_dir = (-camera_dir + camera_up * 0.8).normalized()
    rim_pos = cam_pos + rim_dir * 10.0

    key_light = _ensure_light("OR_KeyLight", key_pos, "AREA", key_intensity, scene)
    fill_light = _ensure_light("OR_FillLight", fill_pos, "AREA", fill_intensity, scene)
    rim_light = _ensure_light("OR_RimLight", rim_pos, "AREA", rim_intensity, scene)

    # Point lights at origin (approximate subject position)
    for light in [key_light, fill_light, rim_light]:
        direction = mathutils.Vector((0, 0, 0)) - light.location
        if direction.length > 0.001:
            rot_quat = direction.to_track_quat("-Z", "Y")
            light.rotation_euler = rot_quat.to_euler()

    return key_light, fill_light, rim_light


def _ensure_light(name, position, light_type, intensity, scene):
    """Get or create a light object."""
    obj = bpy.data.objects.get(name)
    if obj is None:
        light_data = bpy.data.lights.new(name, light_type)
        obj = bpy.data.objects.new(name, light_data)
        scene.collection.objects.link(obj)

    obj.location = position
    obj.data.energy = intensity
    obj.data.type = light_type
    if light_type == "AREA":
        obj.data.size = 5.0

    return obj


# ---------------------------------------------------------------------------
# Render settings
# ---------------------------------------------------------------------------

def setup_render_settings(
    scene,
    engine: str = "CYCLES",
    samples: int = 300,
    resolution_x: int = 1920,
    resolution_y: int = 1080,
    transparent: bool = True,
    use_denoising: bool = True,
):
    """
    Configure render settings for publication output.

    Parameters
    ----------
    scene : bpy.types.Scene
    engine : str
        "CYCLES", "BLENDER_EEVEE", or "BLENDER_WORKBENCH".
    samples : int
        Render samples (Cycles only).
    resolution_x : int
    resolution_y : int
    transparent : bool
        Render with transparent film.
    use_denoising : bool
        Enable denoising (Cycles only).
    """
    scene.render.engine = engine
    scene.render.resolution_x = resolution_x
    scene.render.resolution_y = resolution_y
    scene.render.resolution_percentage = 100
    scene.render.film_transparent = transparent
    scene.render.image_settings.file_format = "PNG"
    scene.render.image_settings.color_mode = "RGBA" if transparent else "RGB"
    scene.render.image_settings.color_depth = "16"

    if engine == "CYCLES":
        scene.cycles.samples = samples
        scene.cycles.use_denoising = use_denoising
        scene.cycles.device = "GPU"
    elif engine == "BLENDER_EEVEE":
        scene.eevee.taa_render_samples = max(16, samples // 4)
        scene.eevee.use_taa_reprojection = False
