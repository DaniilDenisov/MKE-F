import json
from pathlib import Path
import unittest

import numpy as np
from section_solver import Engine, Options, Profile, member, point
from section_solver.peak import search

BASE = json.loads(Path(__file__).with_name("migration.json").read_text())


class MigrationTest(unittest.TestCase):
    def test_both_element_sides_and_mesh_peaks(self):
        for case in BASE["profiles"]:
            with self.subTest(section=case["section"]):
                r = Engine().solve(Profile.from_section(case["section"]), BASE["nu"], BASE["load"],
                                   Options(levels=case["levels"]))
                self.assertAlmostEqual(r["J_FE"], case["J"], places=9)
                for ref in case["sides"]:
                    p = point(r["kernel"], BASE["load"], ref["y"], ref["z"], ref["element"])
                    for field in ("sigmaX","tauXY","tauXZ","value"):
                        self.assertAlmostEqual(p[field],ref[field],delta=BASE["point_atol"])
                peak, ref = r["maximum"], case["maximum"]
                self.assertEqual(peak["status"],ref["status"])
                if ref.get("singularPeak"):
                    self.assertNotIn("value",peak)
                    peak,ref = peak["mesh_peak"],ref["mesh_peak"]
                self.assertLess(abs(peak["value"]/ref["value"]-1),BASE["peak_rtol"])

    def test_old_js_search(self):
        vertices,stress = BASE["adaptive"]["element"]
        result = search(np.array([vertices]),np.array([stress]),BASE["peak_rtol"])
        self.assertGreater(result["cells"],1)
        self.assertLess(abs(result["value"]/BASE["adaptive"]["peak"]-1),BASE["peak_rtol"])

    def test_old_js_member(self):
        ref = BASE["member"]
        r = member(Engine(),Profile.from_section(ref["section"]),BASE["nu"],ref["coefficients"])
        self.assertEqual(r["maximum"]["xi"],ref["xi"])
        self.assertAlmostEqual(r["maximum"]["value"],ref["value"])
