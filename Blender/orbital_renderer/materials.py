# -*- coding: utf-8 -*-
#
#  Blender Orbital Renderer — materials.py
#  Material presets and CPK color table.
#
#  CPK colors from companion orbital-visualizer.py.
#  Material presets adapted from Beautiful Atoms / batoms:
#    (c) Xing Wang, Beautiful Atoms Team — GPLv3
#
#  License: GPLv3
#

"""Material creation system with Principled BSDF presets and CPK atom colors."""

import bpy


# ---------------------------------------------------------------------------
# CPK color table (by atomic number)
# ---------------------------------------------------------------------------

CPK_COLORS = {
    1:  (1.00, 1.00, 1.00),  # H - white
    2:  (0.85, 1.00, 1.00),  # He
    3:  (0.80, 0.50, 1.00),  # Li
    4:  (0.76, 1.00, 0.00),  # Be
    5:  (1.00, 0.71, 0.71),  # B
    6:  (0.35, 0.35, 0.35),  # C - dark gray
    7:  (0.14, 0.14, 0.82),  # N - blue
    8:  (1.00, 0.05, 0.05),  # O - red
    9:  (0.50, 0.70, 0.30),  # F
    10: (0.85, 1.00, 1.00),  # Ne
    11: (0.67, 0.36, 0.95),  # Na
    12: (0.54, 1.00, 0.00),  # Mg
    13: (0.75, 0.65, 0.65),  # Al
    14: (0.50, 0.60, 0.60),  # Si
    15: (1.00, 0.50, 0.00),  # P
    16: (1.00, 1.00, 0.19),  # S - yellow
    17: (0.12, 0.94, 0.12),  # Cl - green
    35: (0.65, 0.16, 0.16),  # Br
    53: (0.58, 0.00, 0.58),  # I
}

# Fallback color for elements not in the table
DEFAULT_CPK_COLOR = (0.70, 0.70, 0.70)


def get_cpk_color(atomic_number: int):
    """Return CPK RGB tuple (0-1) for an atomic number."""
    return CPK_COLORS.get(atomic_number, DEFAULT_CPK_COLOR)


# ---------------------------------------------------------------------------
# Material style presets (from batoms, adapted)
# ---------------------------------------------------------------------------

MATERIAL_STYLES = {
    "default": {
        "type": "Principled BSDF",
        "inputs": {"Metallic": 0.10, "Roughness": 0.20, "IOR": 1.4},
    },
    "metallic": {
        "type": "Principled BSDF",
        "inputs": {"Metallic": 1.00, "Roughness": 0.20, "IOR": 1.4},
    },
    "ceramic": {
        "type": "Principled BSDF",
        "inputs": {"Metallic": 0.02, "Roughness": 0.00, "IOR": 1.4},
    },
    "plastic": {
        "type": "Principled BSDF",
        "inputs": {"Metallic": 0.00, "Roughness": 1.00, "IOR": 1.4},
    },
    "mirror": {
        "type": "Principled BSDF",
        "inputs": {"Metallic": 0.99, "Roughness": 0.01, "IOR": 2.0},
    },
    "glass": {
        "type": "Glass BSDF",
        "inputs": {"Roughness": 0.5, "IOR": 1.45},
    },
    "emission": {
        "type": "Emission",
        "inputs": {"Strength": 5.0},
    },
}

MATERIAL_STYLE_ITEMS = [
    ("default",  "Default",  "Soft diffuse with slight metallic"),
    ("metallic", "Metallic", "Full metallic reflection"),
    ("ceramic",  "Ceramic",  "Smooth glossy with no roughness"),
    ("plastic",  "Plastic",  "Fully rough diffuse"),
    ("mirror",   "Mirror",   "Near-perfect mirror reflection"),
    ("glass",    "Glass",    "Transparent glass BSDF"),
    ("emission", "Emission", "Self-illuminating emission shader"),
]


# ---------------------------------------------------------------------------
# Material creation
# ---------------------------------------------------------------------------

def create_material(
    name: str,
    color,
    material_style: str = "default",
    alpha = None,
):
    """
    Create a Blender material with a Principled BSDF (or specialty shader).

    Parameters
    ----------
    name : str
        Material name in bpy.data.materials.
    color : tuple[float, float, float, float]
        RGBA base color. Alpha channel used if alpha is None.
    material_style : str
        Key from MATERIAL_STYLES.
    alpha : float or None
        Explicit alpha override. If None, uses color[3] when available.

    Returns
    -------
    bpy.types.Material
    """
    if material_style not in MATERIAL_STYLES:
        material_style = "default"

    style = MATERIAL_STYLES[material_style]
    node_type = style["type"]
    node_inputs = style["inputs"].copy()

    # Create material
    material = bpy.data.materials.new(name)
    material.diffuse_color = color[:3] if len(color) >= 3 else color
    material.blend_method = "BLEND"
    material.use_nodes = True
    material.show_transparent_back = False

    # Set up nodes
    nodes = material.node_tree.nodes
    links = material.node_tree.links

    # Clear default nodes
    nodes.clear()

    # Create output node
    output_node = nodes.new("ShaderNodeOutputMaterial")
    output_node.location = (300, 0)

    # Create shader node
    if node_type == "Glass BSDF":
        shader_node = nodes.new("ShaderNodeBsdfGlass")
    elif node_type == "Emission":
        shader_node = nodes.new("ShaderNodeEmission")
    else:
        shader_node = nodes.new("ShaderNodeBsdfPrincipled")

    shader_node.location = (0, 0)

    # Set base color (handle both Principled BSDF and other shaders)
    color_input_name = "Base Color" if "Base Color" in shader_node.inputs else "Color"
    base_color = color[:3] if len(color) >= 3 else color[:3]
    alpha_val = alpha if alpha is not None else (color[3] if len(color) >= 4 else 1.0)

    shader_node.inputs[color_input_name].default_value = (*base_color, 1.0)

    # Apply material style inputs
    for key, value in node_inputs.items():
        if key in shader_node.inputs:
            shader_node.inputs[key].default_value = value

    # Handle alpha for Principled BSDF
    if "Alpha" in shader_node.inputs:
        shader_node.inputs["Alpha"].default_value = alpha_val

    # Link shader to output
    if node_type == "Glass BSDF" or node_type == "Emission":
        links.new(shader_node.outputs.get("BSDF") or shader_node.outputs[0], output_node.inputs["Surface"])
    else:
        links.new(shader_node.outputs["BSDF"], output_node.inputs["Surface"])

    return material


def create_orbital_material(
    name: str,
    color,
    alpha: float = 0.6,
    material_style: str = "default",
):
    """
    Create a material for an orbital isosurface lobe.

    Orbital materials use Principled BSDF with alpha transparency.

    Parameters
    ----------
    name : str
        Material name.
    color : tuple[float, float, float]
        RGB color for the lobe.
    alpha : float
        Transparency (0 = fully transparent, 1 = opaque).
    material_style : str
        Material style preset.

    Returns
    -------
    bpy.types.Material
    """
    rgba = (*color, alpha)
    return create_material(name, rgba, material_style=material_style, alpha=alpha)


def create_atom_material(name: str, atomic_number: int):
    """
    Create a ceramic CPK-colored material for an atom sphere.

    Parameters
    ----------
    name : str
        Material name.
    atomic_number : int
        Atomic number for CPK color lookup.

    Returns
    -------
    bpy.types.Material
    """
    color = get_cpk_color(atomic_number)
    rgba = (*color, 1.0)
    return create_material(name, rgba, material_style="ceramic")


def create_bond_material(name: str = "OR_Bond"):
    """
    Create a neutral gray material for bond cylinders.

    Parameters
    ----------
    name : str
        Material name.

    Returns
    -------
    bpy.types.Material
    """
    gray = (1.0, 1.0, 1.0, 1.0)
    return create_material(name, gray, material_style="plastic")
