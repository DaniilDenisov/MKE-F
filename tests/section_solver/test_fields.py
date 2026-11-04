import unittest
from unittest.mock import patch
import numpy as np
from section_solver import adapter, Profile


class FieldTest(unittest.TestCase):
    def test_resultants_and_local_recovery(self):
        mesh = adapter.generate(Profile("rectangle", 3, 1), 32, .012)
        extract = adapter._extract
        def inspect(section, *args):
            k = extract(section, *args)
            el = section.elements[0]
            p = [.17, .29, .54]
            a,b,c = p
            weights = [a*(2*a-1),b*(2*b-1),c*(2*c-1),4*a*b,4*b*c,4*c*a]
            xy = p @ el.coords.T[:3]
            for i, name in enumerate(("vx", "vy", "mzz")):
                loads = dict(adapter._arguments(section), vx=0, vy=0, mzz=0)
                loads[name] = 1
                direct = el.local_element_stress(tuple(xy), **loads, **adapter._local_args(section, el))
                np.testing.assert_allclose(weights @ k.units[0,:,i], direct[1:], atol=1e-10)
            return k
        with patch.object(adapter, "_extract", inspect):
            k = adapter.solve_fields(mesh, .3)
        v = mesh.vertices
        area = abs(np.linalg.det(v[:,1:]-v[:,:1]))/2
        resultant = np.zeros((3,3))
        for bary, w in [([1/3]*3,-27), ([.6,.2,.2],25), ([.2,.6,.2],25), ([.2,.2,.6],25)]:
            a,b,c = bary
            shape = [a*(2*a-1),b*(2*b-1),c*(2*c-1),4*a*b,4*b*c,4*c*a]
            tau = np.einsum("n,enfc->efc",shape,k.units)
            y,z = np.einsum("n,enc->ce",bary,v)
            moment = y[:,None]*tau[:,:,1]-z[:,None]*tau[:,:,0]
            resultant += w/48*np.einsum("e,efc->fc",area,np.dstack((tau,moment)))
        np.testing.assert_allclose(resultant, np.eye(3), atol=1e-3)
        other = adapter.solve_fields(mesh, .1)
        self.assertGreater(abs(k.units[:,:,:2]-other.units[:,:,:2]).max(), .001)
        self.assertAlmostEqual(k.torsion_report(k.J_FE*2)["relative_difference"], -.5)
        nodes = mesh.triangles[:,[5,3,4]].ravel()
        values = k.units[:,3:].reshape(-1,3,2)
        order = np.argsort(nodes)
        jumps = np.diff(values[order],axis=0)[np.diff(nodes[order])==0]
        self.assertGreater(abs(jumps).max(), 1e-6)

