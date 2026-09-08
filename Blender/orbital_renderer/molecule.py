# -*- coding: utf-8 -*-
#
#  Blender Orbital Renderer — molecule.py
#  Creates atom sphere and bond cylinder meshes.
#
#  Periodic table data adapted from Molecular Blender:
#    (c) Shane Parker, Joshua Szekely — GPLv3
#
#  License: GPLv3
#

"""Build atom (UV sphere) and bond (cylinder) Blender mesh objects."""

import bpy
import mathutils


# ---------------------------------------------------------------------------
# Periodic table — covalent radii (Å)
# ---------------------------------------------------------------------------

COVALENT_RADII = {
    1:  0.31,  2:  0.28,  3:  1.28,  4:  0.96,  5:  0.84,
    6:  0.76,  7:  0.71,  8:  0.66,  9:  0.57,  10: 0.58,
    11: 1.66, 12: 1.41, 13: 1.21, 14: 1.11, 15: 1.07,
    16: 1.05, 17: 1.02, 18: 1.06, 19: 2.03, 20: 1.76,
    21: 1.70, 22: 1.36, 23: 1.53, 24: 1.39, 25: 1.39,
    26: 1.32, 27: 1.26, 28: 1.24, 29: 1.32, 30: 1.22,
    31: 1.22, 32: 1.20, 33: 1.19, 34: 1.20, 35: 1.20,
    36: 1.16, 37: 2.20, 38: 1.95, 39: 1.90, 40: 1.75,
    41: 1.64, 42: 1.54, 43: 1.47, 44: 1.46, 45: 1.42,
    46: 1.39, 47: 1.45, 48: 1.44, 49: 1.42, 50: 1.39,
    51: 1.39, 52: 1.38, 53: 1.39, 54: 1.40, 55: 2.44,
    56: 2.15, 57: 2.07, 58: 2.04, 59: 2.03, 60: 2.01,
    61: 1.99, 62: 1.98, 63: 1.98, 64: 1.96, 65: 1.94,
    66: 1.92, 67: 1.92, 68: 1.89, 69: 1.90, 70: 1.87,
    71: 1.87, 72: 1.75, 73: 1.70, 74: 1.62, 75: 1.51,
    76: 1.44, 77: 1.41, 78: 1.36, 79: 1.36, 80: 1.32,
    81: 1.45, 82: 1.46, 83: 1.48, 84: 1.40, 85: 1.50,
    86: 1.50, 87: 2.60, 88: 2.21, 89: 2.15, 90: 2.06,
    91: 2.00, 92: 1.96, 93: 1.90, 94: 1.87, 95: 1.80,
    96: 1.69,
}

ATOM_RADIUS_SCALE = 0.3
ATOM_SEGMENTS = 32
ATOM_RINGS = 16

BOND_RADIUS = 0.12
BOND_SEGMENTS = 16


# ---------------------------------------------------------------------------
# Atom creation
# ---------------------------------------------------------------------------

def get_covalent_radius(atomic_number: int) -> float:
    """Return covalent radius in Å, or a reasonable default."""
    return COVALENT_RADII.get(atomic_number, 1.50)


def create_atom(
    name: str,
    atomic_number: int,
    position,
    collection: bpy.types.Collection,
):
    """
    Create a UV sphere at atom position, sized by covalent radius.

    Parameters
    ----------
    name : str
        Object name in Blender.
    atomic_number : int
        Atomic number for radius lookup.
    position : tuple[float, float, float]
        Atom position in Å.
    collection : bpy.types.Collection
        Blender collection to link the object into.

    Returns
    -------
    bpy.types.Object
        The created UV sphere object.
    """
    radius = get_covalent_radius(atomic_number) * ATOM_RADIUS_SCALE

    # Create mesh
    bpy.ops.mesh.primitive_uv_sphere_add(
        radius=radius,
        location=position,
        segments=ATOM_SEGMENTS,
        ring_count=ATOM_RINGS,
    )
    obj = bpy.context.active_object
    obj.name = name

    # Link to collection
    if obj.name not in collection.objects:
        collection.objects.link(obj)

    # Unlink from default collection if needed
    for coll in obj.users_collection:
        if coll != collection:
            coll.objects.unlink(obj)

    # Smooth shading
    obj.data.use_auto_smooth = False
    bpy.ops.object.shade_smooth()

    return obj


# ---------------------------------------------------------------------------
# Bond creation
# ---------------------------------------------------------------------------

def create_bond(
    name: str,
    pos1,
    pos2,
    collection: bpy.types.Collection,
    radius: float = BOND_RADIUS,
):
    """
    Create a cylinder between two atom positions.

    Parameters
    ----------
    name : str
        Object name.
    pos1 : tuple[float, float, float]
        Start position in Å.
    pos2 : tuple[float, float, float]
        End position in Å.
    collection : bpy.types.Collection
        Blender collection to link into.
    radius : float
        Cylinder radius in Å.

    Returns
    -------
    bpy.types.Object
        The created cylinder object.
    """
    p1 = mathutils.Vector(pos1)
    p2 = mathutils.Vector(pos2)
    mid = (p1 + p2) / 2.0
    direction = p2 - p1
    length = direction.length

    if length < 1e-6:
        return None

    # Create cylinder
    bpy.ops.mesh.primitive_cylinder_add(
        radius=radius,
        depth=length,
        location=mid,
        vertices=BOND_SEGMENTS,
    )
    obj = bpy.context.active_object
    obj.name = name

    # Orient cylinder along the bond direction
    # Default cylinder is along Z; rotate to align with direction
    z_axis = mathutils.Vector((0, 0, 1))
    rot_quat = z_axis.rotation_difference(direction)
    obj.rotation_mode = "QUATERNION"
    obj.rotation_quaternion = rot_quat

    # Link to collection
    if obj.name not in collection.objects:
        collection.objects.link(obj)

    for coll in obj.users_collection:
        if coll != collection:
            coll.objects.unlink(obj)

    # Smooth shading
    bpy.ops.object.shade_smooth()

    return obj
