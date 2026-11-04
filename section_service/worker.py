from section_solver import Engine, Options, Profile, member, point, truss
from . import display


def compute(engine, operation, data):
    if operation == "member" and "axial" in data:
        return truss(data["axial"],data["area"])
    profile = Profile.from_section(data["section"])
    options = Options(**data["options"])
    args = dict(options=options,J_model=data["J_model"])
    if operation == "member":
        result = member(engine,profile,data["nu"],data["coefficients"],
                        max_samples=data["max_samples"],include_mesh_peak=data["include_mesh_peak"],**args)
    else:
        load = data.get("forces",[0.]*6)
        result = engine.solve(profile,data["nu"],load,**args)
    kernel = result.pop("kernel",None)
    if operation == "diagnostic":
        result["display"] = {"status":"unavailable"}
        if kernel is not None:
            result["level"] = result["history"][-1]["level"]
            result["elements"] = len(kernel.mesh.triangles)
            result["display"], result["_field"] = display.export(kernel,load)
    if operation == "solve":
        result.pop("maximum",None)
    if operation == "point" and kernel is not None and result["status"] not in ("cancelled","resource_limit"):
        p = point(kernel,load,data["y"],data["z"],data.get("element"))
        if p["status"]=="ok" and any(load[1:4]):
            p["status"] = result["fieldStatus"]
        result.update(point=p,status=p["status"])
        result.pop("maximum",None)
    return result


def serve(connection):
    import json
    engine = Engine()
    while True:
        try:
            operation,data = connection.recv()
        except EOFError:
            break
        try:
            result = compute(engine,operation,data)
            field = result.pop("_field",None)
            # Validate finite, HTTP-serializable output before reporting completion.
            json.dumps(result,allow_nan=False)
            connection.send(("completed",result,field))
            del field
            if result.get("status") == "resource_limit":
                engine = Engine()
        except Exception as err:
            connection.send(("failed",{"reason":str(err)},None))

