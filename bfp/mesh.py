import mfem.par as mfem

def E_1d(L, NE):
    """
    1D Mesh (in the E direction)
    Domain: E ∈ [0, L].
    NE: Number of elements.
    """
    mesh = mfem.Mesh.MakeCartesian1D(NE, L)
    return mesh

def x_1d(Lx, Nx):
    """
    1D Mesh (in the x direction)
    Domain: x ∈ [0, Lx].
    Nx: Number of elements.
    """
    mesh = mfem.Mesh.MakeCartesian1D(Nx, Lx)
    return mesh

def t_1d(Lt, Nt):
    """
    1D Mesh (in the t direction)
    Domain: t ∈ [0, Lt].
    Nt: Number of elements.
    """
    mesh = mfem.Mesh.MakeCartesian1D(Nt, Lt)
    return mesh

def xy_2d(Lx, Ly, Nx, Ny):
    """
    2D Mesh (x-y plane)
    Domain: x ∈ [0, Lx], y ∈ [0, Ly].
    Nx, Ny: Number of elements in the x and y directions respectively.
    """
    mesh = mfem.Mesh.MakeCartesian2D([Nx, Ny], mfem.Element.QUADRILATERAL, Lx, Ly)
    return mesh

def xyE_3d(Lx, Ly, E_max, Nx, Ny, NE):
    """
    3D Mesh (x, y, E)
    Domain: x ∈ [0, Lx], y ∈ [0, Ly], E ∈ [0, E_max].
    Nx, Ny, NE: Number of elements in the x, y, and E directions respectively.
    """
    mesh = mfem.Mesh.MakeCartesian3D([Nx, Ny, NE], mfem.Element.HEXAHEDRON, Lx, Ly, E_max)
    return mesh

def xyz_3d(Lx, Ly, Lz, Nx, Ny, Nz):
    """
    3D Mesh (x, y, z)
    Domain: x ∈ [0, Lx], y ∈ [0, Ly], z ∈ [0, Lz].
    Nx, Ny, Nz: Number of elements in the x, y, and z directions respectively.
    """
    mesh = mfem.Mesh.MakeCartesian3D([Nx, Ny, Nz], mfem.Element.HEXAHEDRON, Lx, Ly, Lz)
    return mesh