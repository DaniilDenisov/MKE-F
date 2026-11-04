import unittest
import numpy as np
from section_solver import Profile
from section_solver.engine import Engine, Options
from section_solver.peak import maximum
from section_solver.member import _length_peak, coefficients, member, planar, truss


class MemberTest(unittest.TestCase):
    def test_interval_certificate(self):
        c = coefficients([[.8,1,-1],[.1,2,-1.8],[0],[0],[0],[0]])
        def at(x):
            value = np.linalg.norm(np.polynomial.polynomial.polyval(x,c.T))
            return dict(value=value,upper=value)
        result = _length_peak(c,at,np.ones(6),lambda:None,257,1e-6)
        dense = np.linalg.norm(np.polynomial.polynomial.polyval(np.linspace(0,1,10001),c.T),axis=0)
        self.assertGreater(len(result["samples"]),42)
        self.assertGreaterEqual(result["upper"],dense.max())
        self.assertLessEqual(result["relative_gap"],1e-6)

    def test_polynomial_extrema(self):
        self.assertEqual(truss([-1,2],2)["maximum"]["value"],.5)
        self.assertAlmostEqual(truss([0,4,-4],2)["maximum"]["xi"],.5)
        p = Profile("rectangular_tube",3,1,.15)
        e = Engine()
        rows = planar([0],[0],[0,4,-4],90)
        self.assertTrue(np.array_equal(rows[4],[0,4,-4,0]))
        self.assertFalse(rows[5].any())
        r = member(e,p,.3,rows)
        self.assertEqual(r["maximumStatus"],"ok")
        best = r["maximum"]
        self.assertAlmostEqual(best["xi"],.5)
        self.assertAlmostEqual(best["value"],p.b/2/p.properties[1])
        self.assertGreaterEqual(len(best["samples"]),33)
        self.assertEqual(member(e,p,.3,rows,max_samples=10)["status"],"resource_limit")
        rows[1,0]=1
        self.assertEqual(member(e,p,.3,rows)["status"],"unavailable")

    def test_spatial_end(self):
        r = member(Engine(),Profile("rectangle",2,1),.3,[[0,1],[.2],[0],[.1],[0],[0]])
        self.assertEqual(r["status"],"ok",r.get("reason"))
        self.assertEqual(r["maximum"]["xi"],1.)

    def test_singular_member_reports_only_mesh_peak(self):
        engine = Engine()
        p = Profile("square_tube",.1,.1,.005)
        for torque,xi in (([0,1],1.),([0,4,-4],.5)):
            result = member(engine,p,.25,[[0],[0],[0],torque,[0],[0]],
                            Options(levels=1),include_mesh_peak=True)
            self.assertEqual(result["maximumStatus"],"unavailable")
            self.assertNotIn("value",result["maximum"])
            diagnostic = result["maximum"]["mesh_peak"]
            self.assertEqual(diagnostic["xi"],xi)
            direct = maximum(result["kernel"],[0,0,0,1,0,0],tolerance=.000625)
            self.assertAlmostEqual(diagnostic["value"],direct["mesh_peak"]["value"])

    def test_rounded_trial(self):
        p = Profile.from_section({"type":"square_tube","dimensions":{"b":.1,"t":.005,"inner_radius":.005}})
        r = member(Engine(),p,.25,[[0],[0],[-200],[-400],[400,-400],[0]])
        self.assertEqual(r["status"],"ok",r.get("reason"))
        self.assertEqual(r["maximum"]["xi"],0.)
        self.assertTrue(all(h["peak_change"]<=.005 for h in r["history"][-2:]))
        self.assertGreater(r["maximum"]["value"],400*.05/p.properties[1])
        self.assertAlmostEqual(r["maximum"]["value"],maximum(r["kernel"],[0,0,-200,-400,400,0])["value"])

