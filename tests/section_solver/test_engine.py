import unittest
from dataclasses import replace
from unittest.mock import patch

import numpy as np
from section_solver import adapter, Profile
from section_solver.engine import Engine, Options


class EngineTest(unittest.TestCase):
    def test_cache_and_limits(self):
        engine = Engine()
        p = Profile("rectangle",2,1)
        f = np.array([1.,.2,0,.1,0,0])
        o = Options(levels=1)
        with patch.object(adapter,"generate",wraps=adapter.generate) as mesh, patch.object(
                adapter,"solve_fields",wraps=adapter.solve_fields) as solve:
            a = engine.solve(p,.3,f,o)
            self.assertEqual(a["status"],"nonconverged")
            b = engine.solve(p,.3,2*f,o)
            self.assertAlmostEqual(b["maximum"]["value"],2*a["maximum"]["value"])
            self.assertEqual((mesh.call_count,solve.call_count),(1,1))
            engine.solve(p,.2,f,o)
            self.assertEqual((mesh.call_count,solve.call_count),(1,2))
            engine.solve(p,.3,f,replace(o,levels=2))
            self.assertEqual((mesh.call_count,solve.call_count),(2,3))
            engine.solve(Profile("rectangle",3,1),.3,f,o)
            self.assertEqual((mesh.call_count,solve.call_count),(3,4))
            for limit in [dict(max_bytes=engine.resident_bytes-1),dict(max_elements=1),
                          dict(seconds=1e-12),dict(max_cells=1)]:
                self.assertEqual(engine.solve(p,.3,f,replace(o,**limit))["status"],"resource_limit")
            self.assertEqual(engine.solve(p,.3,f,o,cancelled=lambda: True)["status"],"cancelled")
        fresh = Engine().solve(p,.3,f,o)
        self.assertEqual(fresh["maximum"],a["maximum"])

    def test_convergence(self):
        r = Engine().solve(Profile("rectangle",2,1),.3,[1,.2,0,.1,0,0])
        self.assertEqual(r["status"],"ok",r["history"])
        self.assertGreaterEqual(len(r["history"]),3)

    def test_dimensional_scaling(self):
        # Preserve stresses when lengths, forces and moments are scaled together.
        ref = None
        for s in (1.,1e-60,1e60):
            f = [s**2,.2*s**2,0,.1*s**3,0,0]
            result = Engine().solve(Profile("rectangle",2*s,s),.3,f,Options(levels=2))
            self.assertEqual(len(result["history"]),2,result.get("reason"))
            vals = [result["maximum"]["value"],result["J_FE"]/s**4,
                    *result["history"][-1]["field_changes"]]
            if ref is None:
                ref = vals
            else:
                np.testing.assert_allclose(vals,ref,rtol=.2,atol=.002)

