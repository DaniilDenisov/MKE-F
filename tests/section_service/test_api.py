"""Exercise the real HTTP server and independently call the public package."""
import json
import os
import subprocess
import sys
import time
import unittest
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import psutil
import numpy as np
from section_solver import Profile, Engine, Options, member, point


def request(path, body=None, method=None):
    req = Request("http://127.0.0.1:8301/api/v1/sections/"+path,
                  None if body is None else json.dumps(body).encode(),
                  {"Content-Type":"application/json"},method=method)
    with urlopen(req,timeout=3) as response:
        return json.load(response)


def wait(key):
    for _ in range(3000):
        row = request("jobs/"+key)
        if row["state"] not in ("running","queued"):
            return row
        time.sleep(.02)
    raise AssertionError("Job did not finish")


class ServiceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.server = subprocess.Popen([sys.executable,"-m","uvicorn","section_service.app:app",
                                      "--port","8301","--no-access-log"],
                                      env=dict(os.environ,SECTION_QUEUE_CAPACITY="2"))
        cls.addClassCleanup(cls.shutdown)
        for _ in range(200):
            if cls.server.poll() is not None:
                raise RuntimeError("Server exited")
            try:
                request("health")
                return
            except URLError:
                time.sleep(.05)
        raise RuntimeError("No health response")

    @classmethod
    def shutdown(cls):
        cls.server.terminate()
        cls.server.wait(timeout=10)

    def test_operations(self):
        p = Profile("rectangle",2,1)
        e = Engine()
        opts = Options(levels=2)
        args = dict(section={"type":"rectangle","dimensions":{"b":1,"h":2}},
                    nu=.3,options={"levels":2},J_model=.5)
        f = [.1,1,-.31,.17,.09,-.06]
        rows = [[v] for v in f]
        for op,body in [("solve",args),("maximum",dict(args,forces=f)),
                        ("point",dict(args,forces=f,y=.1,z=.2)),
                        ("member",dict(args,coefficients=rows))]:
            job = wait(request(op,body)["id"])
            self.assertEqual(job["state"],"completed",job)
            r = job["result"]
            if op=="member":
                direct = member(e,p,.3,rows,opts)
            else:
                direct = e.solve(p,.3,[0.]*6 if op=="solve" else f,opts)
            self.assertEqual(r["fieldStatus"],direct["fieldStatus"])
            self.assertAlmostEqual(r["J_FE"],direct["J_FE"])
            if op=="point":
                direct_point = point(direct["kernel"],f,.1,.2)
                self.assertAlmostEqual(r["point"]["value"],direct_point["value"])
                self.assertEqual(r["point"]["status"],direct["fieldStatus"])
            if op in ("maximum","member"):
                self.assertAlmostEqual(r["maximum"]["value"],direct["maximum"]["value"])
            self.assertNotIn("kernel",r)
        with self.assertRaises(HTTPError) as ex:
            request("solve",dict(args,nu=.5))
        self.assertEqual(ex.exception.code,422)

    def test_diagnostic_field(self):
        f = [.1,1,-.31,.17,.09,-.06]
        body = dict(section={"type":"rectangle","dimensions":{"b":1,"h":2}},
                    nu=.3,forces=f,options={"levels":1},J_model=.4)
        job = wait(request("diagnostic",body)["id"])
        r = job["result"]
        direct = Engine().solve(Profile("rectangle",2,1),.3,f,Options(levels=1))
        self.assertEqual(r["maximum"],direct["maximum"])
        self.assertEqual(r["level"],r["history"][-1]["level"])
        self.assertEqual(r["J_model"],.4)
        path = "jobs/"+job["id"]
        with urlopen("http://127.0.0.1:8301/api/v1/sections/"+path+"/field") as resp:
            data = resp.read()
        self.assertEqual(len(data),r["display"]["bytes"])
        version,nn,ne,div = map(int,np.frombuffer(data,"<u4",4))
        self.assertEqual((version,ne,div),(1,r["elements"],4))
        nodes = np.frombuffer(data,"<f8",2*nn,16).reshape(-1,2)
        tris = np.frombuffer(data,"<u4",6*ne,16+16*nn).reshape(-1,6)
        values = np.frombuffer(data,"<f8",offset=16+16*nn+24*ne).reshape(ne,15,4)
        np.testing.assert_equal(nodes,direct["kernel"].mesh.nodes)
        for i in (0,ne-1):
            v = nodes[tris[i,:3]]
            sample = 0
            for a in range(5):
                for b in range(5-a):
                    p = point(direct["kernel"],f,*((4-a-b)*v[0]/4+a*v[1]/4+b*v[2]/4),i)
                    np.testing.assert_allclose(values[i,sample],[p[k] for k in ("sigmaX","tauXY","tauXZ","value")],atol=1e-11)
                    sample += 1
        self.assertEqual(request(path,method="DELETE")["result"]["display"]["status"],"released")
        with self.assertRaises(HTTPError) as err:
            request(path+"/field")
        self.assertEqual(err.exception.code,404)

    def test_queue_cancel_and_deadline(self):
        slow = dict(section={"type":"circle","dimensions":{"d":2}},nu=.3,
                    options={"area_fraction":.0002})
        first = request("solve",slow)["id"]
        second = request("member",{"axial":[1],"area":2})["id"]
        with self.assertRaises(HTTPError) as ex:
            request("solve",slow)
        self.assertEqual(ex.exception.code,429)
        self.assertEqual(request("jobs/"+second,method="DELETE")["state"],"cancelled")
        self.assertEqual(request("jobs/"+first)["state"],"running")
        workers = [p for p in psutil.Process(self.server.pid).children() if
                   "spawn_main" in " ".join(p.cmdline())]
        self.assertEqual(len(workers),1)
        self.assertEqual(request("jobs/"+first,method="DELETE")["state"],"cancelled")
        self.assertFalse(workers[0].is_running())
        tiny = dict(slow,options={"seconds":.001})
        self.assertEqual(wait(request("solve",tiny)["id"])["result"]["status"],"resource_limit")
        next_job = wait(request("member",{"axial":[0,4,-4],"area":2})["id"])
        self.assertEqual(next_job["result"]["maximum"]["value"],.5)

