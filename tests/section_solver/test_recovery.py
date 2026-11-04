import unittest
import numpy as np
from section_solver import Profile, adapter
from section_solver.recovery import point


class PointTest(unittest.TestCase):
    def test_units_and_mixed(self):
        mesh = adapter.generate(Profile("circle",2,2),128,.008)
        k = adapter.solve_fields(mesh,.3)
        y,z = .31,-.18
        loads = np.eye(6)
        def components(f):
            r = point(k,f,y,z)
            return np.array([r[key] for key in ("sigmaX","tauXY","tauXZ")])
        s = np.array([components(f) for f in loads])
        np.testing.assert_allclose(s[:,0], [1/np.pi,0,0,0,4*z/np.pi,-4*y/np.pi], atol=1e-12)
        np.testing.assert_allclose(s[3,1:],np.array([-z,y])*2/np.pi, atol=.001)
        # Circular profile: rotating both point and unit shear by 90 degrees.
        rotated = point(k,loads[2],-z,y)
        np.testing.assert_allclose([rotated["tauXY"],rotated["tauXZ"]],[-s[1,2],s[1,1]],atol=.001)
        mixed = [.1,1,-.31,.17,.09,-.06]
        np.testing.assert_allclose(components(mixed), mixed @ s, atol=1e-12)
        np.testing.assert_allclose(components(np.multiply(-2,mixed)), -2*components(mixed))
        self.assertEqual(point(k,mixed,2,0)["status"],"unavailable")

    def test_hole_and_edge(self):
        mesh = adapter.generate(Profile("rectangular_tube",3,1,.15),32,.02)
        k = adapter.solve_fields(mesh,.3)
        f = [0,1,0,0,0,0]
        self.assertEqual(point(k,f,0,0)["status"],"unavailable")
        self.assertTrue(point(k,f,1.35,.35)["singularPeak"])
        edges = {}
        for i, row in enumerate(mesh.triangles):
            for index in row[3:]:
                if index in edges:
                    a,b = edges[index],i
                    xy = mesh.nodes[index]
                    self.assertEqual(point(k,f,*xy)["element"],a)
                    self.assertEqual(point(k,f,*xy,b)["element"],b)
                    return
                edges[index] = i
        self.fail("Expected shared edge")

