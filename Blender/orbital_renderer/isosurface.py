# -*- coding: utf-8 -*-
#
#  Blender Orbital Renderer — isosurface.py
#  Builds Blender meshes from volumetric grids via marching cubes.
#
#  License: GPLv3
#

"""Convert volumetric grid data to Blender mesh objects via marching cubes."""

import numpy as np
import bpy

from . import marching_cubes


# Module-level grid storage: holds the currently loaded orbital's grid data.
# Keyed by a unique id string; cleared on recipe load.
_grid_cache = {}


def store_grid(key: str, grid_data: np.ndarray):
    """Store a volumetric grid in the module cache."""
    _grid_cache[key] = grid_data


def get_grid(key: str):
    """Retrieve a volumetric grid from cache, or None."""
    return _grid_cache.get(key)


def clear_grid_cache():
    """Clear all cached grids."""
    _grid_cache.clear()


# ---------------------------------------------------------------------------
# Mesh building
# ---------------------------------------------------------------------------

def _triangles_to_mesh(triangles, axes, origin):
    """
    Convert marching-cubes triangles to a Blender mesh.

    Parameters
    ----------
    triangles : list of list of (x, y, z)
        Triangles from marching_cubes, vertices in grid-index space.
    axes : np.ndarray, shape (3, 3)
        Grid axis step vectors in Å.
    origin : np.ndarray, shape (3,)
        Grid origin in Å.

    Returns
    -------
    bpy.types.Mesh
    """
    n_tris = len(triangles)
    if n_tris == 0:
        return None

    n_verts = n_tris * 3
    vertices = np.zeros((n_verts, 3), dtype=np.float32)
    faces = np.zeros((n_tris, 3), dtype=np.int32)

    for i, tri in enumerate(triangles):
        for j in range(3):
            # Grid-space coordinates
            gi = tri[j]
            # Transform to physical Å coordinates:
            # physical = origin + grid_x * axis[0] + grid_y * axis[1] + grid_z * axis[2]
            phys = (
                origin
                + gi[0] * axes[0, :]
                + gi[1] * axes[1, :]
                + gi[2] * axes[2, :]
            )
            vertices[i * 3 + j] = phys
        faces[i] = [i * 3, i * 3 + 1, i * 3 + 2]

    # Build Blender mesh
    mesh = bpy.data.meshes.new("isosurface_mesh")
    mesh.from_pydata(vertices.tolist(), [], faces.tolist())
    mesh.update(calc_edges=True)
    mesh.validate()

    return mesh


def create_isosurface(
    grid: np.ndarray,
    isovalue: float,
    name: str,
    axes: np.ndarray,
    origin: np.ndarray,
    collection: bpy.types.Collection,
):
    """
    Create a Blender mesh object from a volumetric grid and an isovalue.

    Parameters
    ----------
    grid : np.ndarray, shape (nx, ny, nz)
        Volumetric data.
    isovalue : float
        Isosurface threshold.
    name : str
        Object name.
    axes : np.ndarray, shape (3, 3)
        Grid axis step vectors in Å.
    origin : np.ndarray, shape (3,)
        Grid origin in Å.
    collection : bpy.types.Collection
        Collection to link the object into.

    Returns
    -------
    bpy.types.Object or None
        The created Blender object, or None if no surface was found.
    """
    triangles = marching_cubes.marching_cubes(grid, isovalue)

    if not triangles:
        return None

    mesh = _triangles_to_mesh(triangles, axes, origin)
    if mesh is None:
        return None

    obj = bpy.data.objects.new(name, mesh)

    # Link to collection
    collection.objects.link(obj)

    return obj


def update_isosurface(
    obj: bpy.types.Object,
    grid: np.ndarray,
    isovalue: float,
    axes: np.ndarray,
    origin: np.ndarray,
):
    """
    Update an existing Blender mesh object with new isovalue.

    Parameters
    ----------
    obj : bpy.types.Object
        Existing Blender object to update.
    grid : np.ndarray
        Volumetric data.
    isovalue : float
        New isosurface threshold.
    axes : np.ndarray, shape (3, 3)
        Grid axis step vectors in Å.
    origin : np.ndarray, shape (3,)
        Grid origin in Å.

    Returns
    -------
    bool
        True if mesh was updated, False if no surface at this isovalue.
    """
    triangles = marching_cubes.marching_cubes(grid, isovalue)

    if not triangles:
        return False

    new_mesh = _triangles_to_mesh(triangles, axes, origin)
    if new_mesh is None:
        return False

    old_mesh = obj.data
    obj.data = new_mesh
    if old_mesh and old_mesh.users == 0:
        bpy.data.meshes.remove(old_mesh)

    return True
