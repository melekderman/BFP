import mfem.par as mfem

def E_1d(L: float, NE: int) -> mfem.Mesh:
    """Generates a 1D mesh in the E direction.
    
    Args:
        L (float): Length of the domain in the E direction.
        NE (int): Number of elements.

    Returns:
        mfem.Mesh: A 1D Cartesian mesh.
    """
    return mfem.Mesh.MakeCartesian1D(NE, L)

def x_1d(Lx: float, Nx: int) -> mfem.Mesh:
    """Generates a 1D mesh in the x direction.
    
    Args:
        Lx (float): Length of the domain in the x direction.
        Nx (int): Number of elements.

    Returns:
        mfem.Mesh: A 1D Cartesian mesh.
    """
    return mfem.Mesh.MakeCartesian1D(Nx, Lx)

def t_1d(Lt: float, Nt: int) -> mfem.Mesh:
    """Generates a 1D mesh in the t direction.
    
    Args:
        Lt (float): Length of the domain in the t direction.
        Nt (int): Number of elements.

    Returns:
        mfem.Mesh: A 1D Cartesian mesh.
    """
    return mfem.Mesh.MakeCartesian1D(Nt, Lt)

def xy_2d(Lx: float, Ly: float, Nx: int, Ny: int) -> mfem.Mesh:
    """Generates a 2D mesh in the x-y plane.
    
    Args:
        Lx (float): Length of the domain in the x direction.
        Ly (float): Length of the domain in the y direction.
        Nx (int): Number of elements in the x direction.
        Ny (int): Number of elements in the y direction.

    Returns:
        mfem.Mesh: A 2D Cartesian mesh.
    """
    return mfem.Mesh.MakeCartesian2D([Nx, Ny], mfem.Element.QUADRILATERAL, Lx, Ly)

def xyE_3d(Lx: float, Ly: float, E_max: float, Nx: int, Ny: int, NE: int) -> mfem.Mesh:
    """Generates a 3D mesh in the x, y, and E directions.
    
    Args:
        Lx (float): Length of the domain in the x direction.
        Ly (float): Length of the domain in the y direction.
        E_max (float): Length of the domain in the E direction.
        Nx (int): Number of elements in the x direction.
        Ny (int): Number of elements in the y direction.
        NE (int): Number of elements in the E direction.

    Returns:
        mfem.Mesh: A 3D Cartesian mesh.
    """
    return mfem.Mesh.MakeCartesian3D([Nx, Ny, NE], mfem.Element.HEXAHEDRON, Lx, Ly, E_max)

def xyz_3d(Lx: float, Ly: float, Lz: float, Nx: int, Ny: int, Nz: int) -> mfem.Mesh:
    """Generates a 3D mesh in the x, y, and z directions.
    
    Args:
        Lx (float): Length of the domain in the x direction.
        Ly (float): Length of the domain in the y direction.
        Lz (float): Length of the domain in the z direction.
        Nx (int): Number of elements in the x direction.
        Ny (int): Number of elements in the y direction.
        Nz (int): Number of elements in the z direction.

    Returns:
        mfem.Mesh: A 3D Cartesian mesh.
    """
    return mfem.Mesh.MakeCartesian3D([Nx, Ny, Nz], mfem.Element.HEXAHEDRON, Lx, Ly, Lz)