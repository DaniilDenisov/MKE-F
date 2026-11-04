import unittest
from section_solver import Profile, adapter
from section_solver.peak import maximum


class PeakTest(unittest.TestCase):
    def test_sharp_corner(self):
        k = adapter.solve_fields(adapter.generate(Profile("square_tube",1,1,.1),32,.01),.3)
        shear = maximum(k,[0,0,1,0,0,0])
        self.assertEqual(shear["status"],"unavailable")
        self.assertGreater(shear["mesh_peak"]["value"],0)
        axial = maximum(k,[-1,0,0,0,0,0])
        self.assertEqual(axial["status"],"ok")
        self.assertAlmostEqual(axial["value"],1/.36)

