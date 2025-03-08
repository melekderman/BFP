import unittest
import os
import shutil
import numpy as np
from mesh import create_2D_mesh, create_3D_mesh


class TestMeshCreation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.target_dir = os.path.join(os.getcwd(), 'mesh', 'usr')
        if os.path.exists(cls.target_dir):
            shutil.rmtree(cls.target_dir)

    def test_create_2D_mesh(self):
        nx, ny = 10, 20
        x_start, x_end = 0.0, 1.0
        y_start, y_end = 0.0, 2.0

        mesh = create_2D_mesh(nx, ny, x_start, x_end, y_start, y_end)
        
        self.assertEqual(mesh.GetNV(), (nx + 1) * (ny + 1))

        verts = mesh.GetVertexArray()
        x_coords = np.linspace(x_start, x_end, nx + 1)
        y_coords = np.linspace(y_start, y_end, ny + 1)

        k = 0
        for j in range(ny + 1):
            for i in range(nx + 1):
                self.assertAlmostEqual(verts[k][0], x_coords[i])
                self.assertAlmostEqual(verts[k][1], y_coords[j])
                k += 1

        expected_file = os.path.join(self.target_dir, f'{nx}x{ny}_2D.mesh')
        self.assertTrue(os.path.exists(expected_file))

    def test_create_3D_mesh(self):
        nx, ny, nz = 5, 5, 5
        x_start, x_end = 0.0, 1.0
        y_start, y_end = 0.0, 1.0
        z_start, z_end = 0.0, 1.0

        mesh = create_3D_mesh(nx, ny, nz, x_start, x_end, y_start, y_end, z_start, z_end)

        self.assertEqual(mesh.GetNV(), (nx + 1) * (ny + 1) * (nz + 1))

        verts = mesh.GetVertexArray()
        x_coords = np.linspace(x_start, x_end, nx + 1)
        y_coords = np.linspace(y_start, y_end, ny + 1)
        z_coords = np.linspace(z_start, z_end, nz + 1)

        k = 0
        for k_z in range(nz + 1):
            for k_y in range(ny + 1):
                for k_x in range(nx + 1):
                    self.assertAlmostEqual(verts[k][0], x_coords[k_x])
                    self.assertAlmostEqual(verts[k][1], y_coords[k_y])
                    self.assertAlmostEqual(verts[k][2], z_coords[k_z])
                    k += 1

        expected_file = os.path.join(self.target_dir, f'{nx}x{ny}_3D.mesh')
        self.assertTrue(os.path.exists(expected_file))


if __name__ == '__main__':
    unittest.main()
