import unittest

import numpy as np
from shapely.geometry import Point

from section_solver import Profile, validate_nu


class GeometryTest(unittest.TestCase):
    def test_profiles(self):
        for p in (Profile("circle",1,1), Profile("rectangle",3,1), Profile("circular_tube",1,1,.15),
                  Profile("square_tube",1,1,.15), Profile("rectangular_tube",3,1,.15),
                  Profile("rectangular_tube",3,1,.15,.15), Profile("square_tube",1,1,.05,.05)):
            h = p.h
            shape = p.polygon(512)
            np.testing.assert_allclose(shape.bounds, [-h/2, -.5, h/2, .5])
            np.testing.assert_allclose(shape.centroid.coords[0], [0, 0], atol=1e-13)
            self.assertEqual(shape.contains(Point(0, 0)), not p.t)
            # Polygon boundary integrals, with signed inner rings.
            moments = np.zeros(3)
            for ring in [shape.exterior, *shape.interiors]:
                a, b = np.asarray(ring.coords)[:-1], np.asarray(ring.coords)[1:]
                cross = a[:,0]*b[:,1]-a[:,1]*b[:,0]
                iz, iy = (cross @ (a*a+a*b+b*b))/12
                moments += [cross.sum()/2, iy, iz]
            np.testing.assert_allclose(abs(moments), p.properties, rtol=1e-4)

    def test_invalid(self):
        for args in [("square_tube", 2, 1, .1), ("rectangle", -1, 1),
                     ("circle", float("nan"), 1), ("circular_tube", 1, 1, .5),
                     ("square_tube",1,1,.2,.3), ("rectangle",2,1,0,.1),
                     ("square_tube",.1,.1,.005,.045), ("square_tube",.1,.1,.005,-.001)]:
            with self.assertRaises(ValueError):
                Profile(*args)
        for nu in [-1, .5, float("nan"), True]:
            with self.assertRaises(ValueError):
                validate_nu(nu)

