# -*- coding: utf-8 -*-
#
#  Blender Orbital Renderer — style_loader.py
#  Style file format, default style, load/save/apply for render presets.
#
#  Like matplotlib's mplstyle, an orbital_style.json captures all render
#  decisions (camera, lights, engine, materials, colors) so multiple figures
#  can be produced with identical look from the CLI.
#
#  License: GPLv3
#

"""Orbital render style — load, save, apply, defaults."""

import json
import os


# ---------------------------------------------------------------------------
# Default style
# ---------------------------------------------------------------------------

DEFAULT_STYLE = {
    "version": 1,
    "camera": {
        "ortho": True,
        "distance_factor": 2.0,
    },
    "lighting": {
        "key_intensity": 100.0,
        "fill_intensity": 50.0,
        "rim_intensity": 80.0,
    },
    "render": {
        "engine": "CYCLES",
        "samples": 300,
        "resolution_x": 1920,
        "resolution_y": 1080,
        "transparent": True,
    },
    "orbitals": {
        "default_isovalue": 0.05,
        "positive_color": [1.0, 0.2, 0.2],
        "negative_color": [0.2, 0.4, 1.0],
        "alpha": 0.6,
        "material_style": "default",
        "subdivision": False,
    },
    "molecule": {
        "show_atoms": True,
        "show_bonds": True,
    },
}


# ---------------------------------------------------------------------------
# Load / save
# ---------------------------------------------------------------------------

def load_style(style_path: str) -> dict:
    """
    Load a style JSON file. Missing keys fall back to defaults.

    Parameters
    ----------
    style_path : str
        Path to the style JSON file.

    Returns
    -------
    dict
        Merged style (loaded values over DEFAULT_STYLE).

    Raises
    ------
    FileNotFoundError
    ValueError
    """
    if not os.path.isfile(style_path):
        raise FileNotFoundError(f"Style file not found: {style_path}")

    with open(style_path, "r") as f:
        user_style = json.load(f)

    return _deep_merge(DEFAULT_STYLE, user_style)


def save_style(style: dict, style_path: str):
    """Save a style dict to a JSON file."""
    with open(style_path, "w") as f:
        json.dump(style, f, indent=2)


def get_default_style() -> dict:
    """Return a deep copy of the default style."""
    return json.loads(json.dumps(DEFAULT_STYLE))


# ---------------------------------------------------------------------------
# Apply to Blender scene
# ---------------------------------------------------------------------------

def apply_style_to_scene(scene_props, style: dict):
    """
    Apply a style dict to the add-on's scene property group.

    Parameters
    ----------
    scene_props : OR_SceneProperties
        The bpy scene property group (context.scene.orbital_renderer).
    style : dict
        Style dict (from load_style or DEFAULT_STYLE).
    """
    cam = style.get("camera", {})
    scene_props.camera_ortho = cam.get("ortho", True)
    scene_props.camera_distance = cam.get("distance_factor", 2.0)

    light = style.get("lighting", {})
    scene_props.light_key_intensity = light.get("key_intensity", 100.0)
    scene_props.light_fill_intensity = light.get("fill_intensity", 50.0)
    scene_props.light_rim_intensity = light.get("rim_intensity", 80.0)

    render = style.get("render", {})
    scene_props.render_engine = render.get("engine", "CYCLES")
    scene_props.render_samples = render.get("samples", 300)
    scene_props.render_resolution_x = render.get("resolution_x", 1920)
    scene_props.render_resolution_y = render.get("resolution_y", 1080)
    scene_props.render_transparent = render.get("transparent", True)

    molecule = style.get("molecule", {})
    scene_props.show_atoms = molecule.get("show_atoms", True)
    scene_props.show_bonds = molecule.get("show_bonds", True)


def apply_style_to_orbital(orbital_item, style: dict):
    """
    Apply orbital defaults from a style to an orbital property group item.

    Parameters
    ----------
    orbital_item : OR_OrbitalItem
        The orbital item in the scene properties.
    style : dict
        Style dict.
    """
    orb_style = style.get("orbitals", {})
    orbital_item.isovalue = orb_style.get("default_isovalue", 0.05)
    pos_col = orb_style.get("positive_color", [1.0, 0.2, 0.2])
    neg_col = orb_style.get("negative_color", [0.2, 0.4, 1.0])
    alpha = orb_style.get("alpha", 0.6)
    orbital_item.pos_color = (*pos_col, alpha)
    orbital_item.neg_color = (*neg_col, alpha)
    orbital_item.pos_alpha = alpha
    orbital_item.neg_alpha = alpha
    orbital_item.material_style = orb_style.get("material_style", "default")
    orbital_item.subdivision = orb_style.get("subdivision", False)


def collect_style_from_scene(scene_props) -> dict:
    """
    Build a style dict from the current scene property values.

    Parameters
    ----------
    scene_props : OR_SceneProperties

    Returns
    -------
    dict
    """
    return {
        "version": 1,
        "camera": {
            "ortho": scene_props.camera_ortho,
            "distance_factor": scene_props.camera_distance,
        },
        "lighting": {
            "key_intensity": scene_props.light_key_intensity,
            "fill_intensity": scene_props.light_fill_intensity,
            "rim_intensity": scene_props.light_rim_intensity,
        },
        "render": {
            "engine": scene_props.render_engine,
            "samples": scene_props.render_samples,
            "resolution_x": scene_props.render_resolution_x,
            "resolution_y": scene_props.render_resolution_y,
            "transparent": scene_props.render_transparent,
        },
        "orbitals": {
            "default_isovalue": 0.05,
            "positive_color": [1.0, 0.2, 0.2],
            "negative_color": [0.2, 0.4, 1.0],
            "alpha": 0.6,
            "material_style": "default",
            "subdivision": False,
        },
        "molecule": {
            "show_atoms": scene_props.show_atoms,
            "show_bonds": scene_props.show_bonds,
        },
    }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge override into base, returning a new dict."""
    result = {}
    for key in base:
        if key in override:
            if isinstance(base[key], dict) and isinstance(override[key], dict):
                result[key] = _deep_merge(base[key], override[key])
            else:
                result[key] = override[key]
        else:
            result[key] = base[key]
    # Include keys only in override
    for key in override:
        if key not in base:
            result[key] = override[key]
    return result
