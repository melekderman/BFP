# mesh.py
import numpy as np
import mfem.ser as mfem
import os

__all__ = ["create_2D_mesh", "create_3D_mesh"]


def create_2D_mesh(nx, ny, x_start, x_end, y_start, y_end):
    """Return a structured quadrilateral MFEM mesh with nx-by-ny elements.

    ``x_start``/``x_end`` and ``y_start``/``y_end`` give coordinate bounds.
    For transport, use increasing position bounds in cm and physical kinetic
    energy bounds in MeV; ``ny`` counts energy elements. The solver uses this
    single coupled (x, E) mesh, not separate spatial and energy meshes.

    Creates ``mesh/usr/`` in the working directory and writes
    ``{nx}x{ny}_2D.mesh`` if absent, independently of solver output controls.
    """

    x_coords = np.linspace(x_start, x_end, nx + 1)
    y_coords = np.linspace(y_start, y_end, ny + 1)

    mesh = mfem.Mesh(nx, ny, "QUADRILATERAL", True, 0.0, 0.0)

    verts = mesh.GetVertexArray()
    expected_num = (nx + 1) * (ny + 1)
    num_verts = mesh.GetNV()

    if num_verts != expected_num:
        print(
            "Warning: Unexpected number of vertices! ({} != {})".format(
                num_verts, expected_num
            )
        )

    k = 0
    for j in range(ny + 1):
        for i in range(nx + 1):
            verts[k][0] = x_coords[i]
            verts[k][1] = y_coords[j]
            k += 1

    target_dir = os.path.join(os.getcwd(), "mesh", "usr")
    file_name = f"{nx}x{ny}_2D.mesh"
    file_path = os.path.join(target_dir, file_name)

    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
        print(f"Directory '{target_dir}' was created.")

    if not os.path.exists(file_path):
        mesh.Print(file_path)
        print(f"File '{file_path}' was successfully created.")
    else:
        print(f"File '{file_path}' already exists.")

    return mesh


"""
# Example usage:
nx = 20
ny = 50
x_min = 0.0
x_max = 0.3
y_start = 1.0
y_end = 0.01

# Create the mesh
mesh = create_custom_mesh(nx, ny, x_min, x_max, y_start, y_end)

# Plot the mesh using glvis
glvis(mesh)

# Print updated vertices for verification
verts = mesh.GetVertexArray()
for k, v in enumerate(verts):
    print("Vertex {}: {}".format(k, v))
"""


def create_3D_mesh(nx, ny, nz, x_start, x_end, y_start, y_end, z_start, z_end):
    """Return a structured hexahedral MFEM mesh with nx-by-ny-by-nz elements.

    The start/end arguments bound each coordinate. This is a geometry helper;
    BFP's transport solvers currently support only the two-coordinate (x, E)
    slab formulation. Creates ``mesh/usr/{nx}x{ny}_3D.mesh`` if absent in the
    working directory; the legacy filename does not include nz.
    """

    x_coords = np.linspace(x_start, x_end, nx + 1)
    y_coords = np.linspace(y_start, y_end, ny + 1)
    z_coords = np.linspace(z_start, z_end, nz + 1)

    mesh = mfem.Mesh(nx, ny, nz, "HEXAHEDRON", True, 0.0, 0.0, 0.0)

    verts = mesh.GetVertexArray()

    expected_num = (nx + 1) * (ny + 1) * (nz + 1)
    num_verts = mesh.GetNV()
    if num_verts != expected_num:
        print(
            "Warning: Unexpected number of vertices! ({} != {})".format(
                num_verts, expected_num
            )
        )

    k = 0
    for k_z in range(nz + 1):
        for k_y in range(ny + 1):
            for k_x in range(nx + 1):
                verts[k][0] = x_coords[k_x]
                verts[k][1] = y_coords[k_y]
                verts[k][2] = z_coords[k_z]
                k += 1

    target_dir = os.path.join(os.getcwd(), "mesh", "usr")
    file_name = f"{nx}x{ny}_3D.mesh"
    file_path = os.path.join(target_dir, file_name)

    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
        print(f"Directory '{target_dir}' was created.")

    if not os.path.exists(file_path):

        mesh.Print(file_path)
        print(f"File '{file_path}' was successfully created.")
    else:
        print(f"File '{file_path}' already exists.")

    return mesh


"""
# Example usage:
nx = 20
ny = 50
nz = 30
x_start = 0.0
x_end = 0.3
y_start = 1.0
y_end = 0.01
z_start = 0.0
z_end = 0.3

# Create the mesh
mesh2 = create_3D_mesh(nx, ny, nz, x_start, x_end, y_start, y_end, z_start, z_end)

# Plot the mesh using glvis
glvis(mesh)

# Print updated vertices for verification:
verts2 = mesh2.GetVertexArray()
for k, v in enumerate(verts2):
    print("Vertex {}: {}".format(k, v))

"""
