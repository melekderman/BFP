# Parametreler: eleman sayıları (her eksendeki eleman sayısı)
# Örnek: 2×2×2 küp için (verilen örneğe eşdeğer)
nelx = 2
nely = 2
nelz = 2

# Türevler:
num_elems = nelx * nely * nelz
num_vertices = (nelx + 1) * (nely + 1) * (nelz + 1)
nnx = 2 * nelx + 1
nny = 2 * nely + 1
nnz = 2 * nelz + 1

def vertex_index(i, j, k, nx, ny):
    return i + (nx) * j + (nx) * (ny) * k

with open("cube.mesh", "w") as f:
    f.write("MFEM mesh v1.0\n\n")
    f.write("dimension\n3\n\n")
    # elements section
    f.write("elements\n")
    f.write(f"{num_elems}\n")
    for k in range(nelz):
        for j in range(nely):
            for i in range(nelx):
                nx_global = nelx + 1
                ny_global = nely + 1
                base = vertex_index(i, j, k, nx_global, ny_global)
                v0 = base
                v1 = base + 1
                v3 = base + nx_global
                v2 = v1 + nx_global  # = base + 1 + nx_global
                v4 = base + nx_global * ny_global
                v5 = v4 + 1
                v7 = v4 + nx_global
                v6 = v5 + nx_global
                # Yaz: "1 5 v0 v1 v2 v3 v4 v5 v6 v7"
                f.write(f"1 5 {v0} {v1} {v2} {v3} {v4} {v5} {v6} {v7}\n")
    f.write("\n")
    # boundary section
    # Toplam boundary face sayısı = 2*(nelx*nely + nelx*nelz + nely*nelz)
    boundary_faces = []
    # (a) Z-boundaries (marker 3)
    # z = 0 face: for i=0..nelx-1, j=0..nely-1 on k=0
    for j in range(nely):
        for i in range(nelx):
            # Use ordering: [v0, v3, v2, v1]
            v0 = vertex_index(i, j, 0, nelx+1, nely+1)
            v1 = vertex_index(i+1, j, 0, nelx+1, nely+1)
            v2 = vertex_index(i+1, j+1, 0, nelx+1, nely+1)
            v3 = vertex_index(i, j+1, 0, nelx+1, nely+1)
            boundary_faces.append(f"3 3 {v0} {v3} {v2} {v1}")
    # z = max face: k = nelz
    for j in range(nely):
        for i in range(nelx):
            # Use ordering: [v0, v1, v2, v3]
            v0 = vertex_index(i, j, nelz, nelx+1, nely+1)
            v1 = vertex_index(i+1, j, nelz, nelx+1, nely+1)
            v2 = vertex_index(i+1, j+1, nelz, nelx+1, nely+1)
            v3 = vertex_index(i, j+1, nelz, nelx+1, nely+1)
            boundary_faces.append(f"3 3 {v0} {v1} {v2} {v3}")
    # (b) X-boundaries (marker 1)
    # x = 0 face: i = 0; loop j=0..nely-1, k=0..nelz-1
    for k in range(nelz):
        for j in range(nely):
            # Use ordering: [v0, v4, v7, v3]
            v0 = vertex_index(0, j, k, nelx+1, nely+1)
            v3 = vertex_index(0, j+1, k, nelx+1, nely+1)
            v4 = vertex_index(0, j, k+1, nelx+1, nely+1)
            v7 = vertex_index(0, j+1, k+1, nelx+1, nely+1)
            boundary_faces.append(f"1 3 {v0} {v4} {v7} {v3}")
    # x = max face: i = nelx; loop j and k; use ordering: [v0, v3, v7, v4] with indices computed at i=nelx
    for k in range(nelz):
        for j in range(nely):
            v0 = vertex_index(nelx, j, k, nelx+1, nely+1)
            v3 = vertex_index(nelx, j+1, k, nelx+1, nely+1)
            v7 = vertex_index(nelx, j+1, k+1, nelx+1, nely+1)
            v4 = vertex_index(nelx, j, k+1, nelx+1, nely+1)
            # To match sample, reverse order to: [v0, v3, v7, v4] → output as: v0, v3, v7, v4? 
            # Sample x=max first face: "1 3 2 5 14 11" corresponds to v0=2, v3=5, v7=14, v4=11.
            boundary_faces.append(f"1 3 {v0} {v3} {v7} {v4}")
    # (c) Y-boundaries (marker 2)
    # y = 0 face: j = 0; loop i=0..nelx-1, k=0..nelz-1; use ordering: [v0, v1, v5, v4]
    for k in range(nelz):
        for i in range(nelx):
            v0 = vertex_index(i, 0, k, nelx+1, nely+1)
            v1 = vertex_index(i+1, 0, k, nelx+1, nely+1)
            v5 = vertex_index(i+1, 0, k+1, nelx+1, nely+1)
            v4 = vertex_index(i, 0, k+1, nelx+1, nely+1)
            boundary_faces.append(f"2 3 {v0} {v1} {v5} {v4}")
    # y = max face: j = nely; loop i and k; use ordering: [v0, v4, v5, v1]
    for k in range(nelz):
        for i in range(nelx):
            v0 = vertex_index(i, nely, k, nelx+1, nely+1)
            v4 = vertex_index(i, nely, k+1, nelx+1, nely+1)
            v5 = vertex_index(i+1, nely, k+1, nelx+1, nely+1)
            v1 = vertex_index(i+1, nely, k, nelx+1, nely+1)
            # To match sample, order as: [v0, v4, v5, v1] becomes [v0, v4, v5, v1] ? 
            # Sample y=max first face: "2 3 6 15 16 7" where 6,15,16,7 are obtained by [6,15,16,7].
            boundary_faces.append(f"2 3 {v0} {v4} {v5} {v1}")
    f.write("boundary\n")
    f.write(f"{len(boundary_faces)}\n")
    for face in boundary_faces:
        f.write(face + "\n")
    f.write("\n")
    # vertices section
    f.write("vertices\n")
    f.write(f"{num_vertices}\n")
    for k in range(nelz+1):
        for j in range(nely+1):
            for i in range(nelx+1):
                x = i / nelx
                y = j / nely
                z = k / nelz
                f.write(f"{x}\n{y}\n{z}\n")
    f.write("\n")
    # nodes section (quadratic FE nodes)
    f.write("nodes\n")
    f.write("FiniteElementSpace\n")
    f.write("FiniteElementCollection: H1_3D_P2\n")
    f.write("VDim: 3\n")
    f.write("Ordering: 0\n\n")
    total_nodes = nnx * nny * nnz
    for k in range(nnz):
        for j in range(nny):
            for i in range(nnx):
                x = i / (nnx - 1)
                y = j / (nny - 1)
                z = k / (nnz - 1)
                f.write(f"{x}\n{y}\n{z}\n")
