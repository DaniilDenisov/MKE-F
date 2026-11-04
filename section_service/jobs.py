"""One bounded FIFO and one reusable numerical process per ASGI server."""
import asyncio
from collections import deque
import multiprocessing
from time import monotonic
from uuid import uuid4

from fastapi import HTTPException

from .worker import serve
from .display import RETAINED_LIMIT


class Jobs:
    def __init__(self, capacity):
        self.capacity = capacity
        self.records = {}
        self.fields = {}
        self.pending = deque()
        self.active = self.process = self.connection = None
        self.wake = asyncio.Event()

    def prune(self):
        done = sorted((k for k,r in self.records.items() if "finished" in r),
                      key=lambda k:self.records[k]["finished"])
        for i,k in enumerate(done):
            if i < len(done)-32 or monotonic()-self.records[k]["finished"]>900:
                del self.records[k]
                self.fields.pop(k,None)

    def field(self,key):
        self.get(key)
        if key not in self.fields:
            raise HTTPException(404,"Field unavailable or expired")
        return self.fields[key]

    def get(self,key):
        self.prune()
        if key not in self.records:
            raise HTTPException(404,"Unknown or expired job")
        return {k:v for k,v in self.records[key].items() if k != "finished"}

    def submit(self,operation,data):
        self.prune()
        if len(self.pending)+int(self.active is not None)>=self.capacity:
            raise HTTPException(429,"Section queue is full",headers={"Retry-After":"1"})
        key = uuid4().hex
        self.records[key] = dict(id=key,operation=operation,state="queued")
        self.pending.append((key,operation,data))
        self.wake.set()
        return self.get(key)

    def stop_worker(self):
        if self.process is not None:
            if self.process.is_alive():
                self.process.kill()
            self.process.join()
            self.process.close()
            self.connection.close()
            self.process = self.connection = None

    def finish(self,key,state,result,field=None):
        self.prune()
        if field is not None:
            while self.fields and sum(map(len,self.fields.values()))+len(field)>RETAINED_LIMIT:
                oldest = next(iter(self.fields))
                self.fields.pop(oldest)
                self.records[oldest]["result"]["display"]["status"] = "evicted"
            self.fields[key] = field
        self.records[key].update(state=state,result=result,finished=monotonic())

    def cancel(self,key):
        state = self.get(key)["state"]
        if state == "running":
            self.stop_worker()
            self.active = None
        elif state == "queued":
            self.pending = deque(row for row in self.pending if row[0]!=key)
        else:
            if self.fields.pop(key,None) is not None:
                self.records[key]["result"]["display"]["status"] = "released"
            return self.get(key)
        self.finish(key,"cancelled",dict(status="cancelled"))
        self.wake.set()
        return self.get(key)

    def dispatch(self):
        key,operation,data = self.pending.popleft()
        self.active = key
        self.deadline = monotonic()+data.get("options",{}).get("seconds",120)
        if self.process is None:
            ctx = multiprocessing.get_context("spawn")
            self.connection,child = ctx.Pipe()
            self.process = ctx.Process(target=serve,args=(child,),daemon=True)
            self.process.start()
            child.close()
        self.connection.send((operation,data))
        self.records[key]["state"] = "running"

    async def run(self):
        try:
            while True:
                try:
                    if self.active is None and self.pending:
                        self.dispatch()
                    if self.active is None:
                        self.wake.clear()
                        try:
                            await asyncio.wait_for(self.wake.wait(),timeout=5)
                        except asyncio.TimeoutError:
                            self.prune()
                        continue
                    if monotonic() >= self.deadline:
                        self.stop_worker()
                        self.finish(self.active,"completed",dict(status="resource_limit",reason="Job deadline"))
                    elif self.connection.poll():
                        state,result,field = self.connection.recv()
                        self.finish(self.active,state,result,field)
                        del field
                    elif not self.process.is_alive():
                        raise EOFError("Section process exited")
                    else:
                        await asyncio.sleep(.025)
                        continue
                except (OSError,EOFError) as err:
                    self.stop_worker()
                    self.finish(self.active,"failed",dict(reason=str(err) or "Section process exited"))
                self.active = None
        finally:
            self.stop_worker()
            self.fields.clear()

