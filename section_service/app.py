import asyncio
from contextlib import asynccontextmanager, suppress
import os

from fastapi import FastAPI, Response

from .jobs import Jobs
from .models import Member, Maximum, Point, Solve, Truss


@asynccontextmanager
async def lifespan(app):
    jobs = Jobs(max(1,min(32,int(os.getenv("SECTION_QUEUE_CAPACITY","4")))))
    app.state.jobs = jobs
    task = asyncio.create_task(jobs.run())
    try:
        yield
    finally:
        task.cancel()
        with suppress(asyncio.CancelledError):
            await task


prefix = "/api/v1/sections"
app = FastAPI(title="MKE-F sections",lifespan=lifespan,docs_url=None,redoc_url=None,
              openapi_url=prefix+"/openapi.json")


@app.get(prefix+"/health")
async def health():
    return {"status":"ok"}


@app.post(prefix+"/solve",status_code=202)
async def solve(body: Solve):
    return app.state.jobs.submit("solve",body.model_dump())


@app.post(prefix+"/point",status_code=202)
async def point(body: Point):
    return app.state.jobs.submit("point",body.model_dump())


@app.post(prefix+"/maximum",status_code=202)
async def maximum(body: Maximum):
    return app.state.jobs.submit("maximum",body.model_dump())


@app.post(prefix+"/diagnostic",status_code=202)
async def diagnostic(body: Maximum):
    return app.state.jobs.submit("diagnostic",body.model_dump())


@app.get(prefix+"/jobs/{key}/field",response_class=Response)
async def field(key: str):
    return Response(app.state.jobs.field(key),media_type="application/octet-stream",
                    headers={"Cache-Control":"no-store"})


@app.post(prefix+"/member",status_code=202)
async def member(body: Member | Truss):
    return app.state.jobs.submit("member",body.model_dump())


@app.get(prefix+"/jobs/{key}")
async def status(key: str):
    return app.state.jobs.get(key)


@app.delete(prefix+"/jobs/{key}")
async def cancel(key: str):
    return app.state.jobs.cancel(key)

