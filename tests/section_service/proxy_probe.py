"""Run inside the rod solver container against the Nginx proxy."""
import json
import sys
from time import monotonic,sleep
from urllib.error import HTTPError
from urllib.request import Request,urlopen


def request(path,body=None):
    req = Request("http://web:8080/api/v1/"+path,
                  None if body is None else json.dumps(body).encode(),
                  headers={"Host":"localhost","Content-Type":"application/json"})
    with urlopen(req,timeout=5) as r:
        return json.load(r)


request("health")
if sys.argv[1]=="down":
    try:
        request("sections/health")
        raise AssertionError("Stopped service is still available")
    except HTTPError as err:
        assert err.code in (502,504),err.code
    request("health")
else:
    deadline=monotonic()+30
    while True:
        try:
            assert request("sections/health")["status"]=="ok"
            break
        except HTTPError:
            if monotonic()>=deadline:
                raise
            sleep(.2)
    # Exercise the section FEM through Nginx, including job polling.
    job=request("sections/maximum",{"section":{"type":"rectangle","dimensions":{"h":2,"b":1}},
                                    "nu":.3,"forces":[1,0,0,0,0,0]})
    while True:
        job=request("sections/jobs/"+job["id"])
        if job["state"] not in ("queued","running"):
            break
        assert monotonic()<deadline
        sleep(.05)
    assert job["state"]=="completed" and job["result"]["status"]=="ok",job
    assert abs(job["result"]["maximum"]["value"]-.5)<1e-12
print("Proxy probe passed:",sys.argv[1])

