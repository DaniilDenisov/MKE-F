from __future__ import annotations

import asyncio
import json
import os
import re
import shutil
import signal
import uuid
from contextlib import asynccontextmanager
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from fastapi import FastAPI, HTTPException, Request, status
from fastapi.middleware.trustedhost import TrustedHostMiddleware
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel, ConfigDict, Field, field_validator


TERMINAL_STATES = {"succeeded", "failed", "canceled", "timed_out"}
ERROR_PATTERN = re.compile(r"^MKEF_ERROR\t([^\t]+)\t(.*)$", re.MULTILINE)


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def env_int(name: str, default: int, minimum: int = 1) -> int:
    raw = os.getenv(name)
    if raw is None:
        return default
    value = int(raw)
    if value < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return value


@dataclass(frozen=True)
class Settings:
    jobs_root: Path = Path("/jobs")
    repository_root: Path = Path("/app")
    octave_executable: str = "octave-cli"
    runner_script: Path = Path("/app/scripts/run_solver_job.m")
    max_case_bytes: int = 10 * 1024 * 1024
    max_result_bytes: int = 100 * 1024 * 1024
    max_log_bytes: int = 256 * 1024
    max_active_jobs: int = 20
    job_timeout_seconds: int = 30 * 60
    job_ttl_seconds: int = 24 * 60 * 60
    cancel_grace_seconds: int = 5
    allowed_origins: tuple[str, ...] = (
        "http://localhost:8080",
        "http://127.0.0.1:8080",
    )
    check_octave: bool = True

    @classmethod
    def from_env(cls) -> "Settings":
        origins = tuple(
            value.strip()
            for value in os.getenv(
                "MKEF_ALLOWED_ORIGINS",
                "http://localhost:8080,http://127.0.0.1:8080",
            ).split(",")
            if value.strip()
        )
        return cls(
            jobs_root=Path(os.getenv("MKEF_JOBS_ROOT", "/jobs")),
            repository_root=Path(os.getenv("MKEF_REPOSITORY_ROOT", "/app")),
            octave_executable=os.getenv("MKEF_OCTAVE_EXECUTABLE", "octave-cli"),
            runner_script=Path(
                os.getenv("MKEF_RUNNER_SCRIPT", "/app/scripts/run_solver_job.m")
            ),
            max_case_bytes=env_int("MKEF_MAX_CASE_BYTES", 10 * 1024 * 1024),
            max_result_bytes=env_int(
                "MKEF_MAX_RESULT_BYTES", 100 * 1024 * 1024
            ),
            max_log_bytes=env_int("MKEF_MAX_LOG_BYTES", 256 * 1024),
            max_active_jobs=env_int("MKEF_MAX_ACTIVE_JOBS", 20),
            job_timeout_seconds=env_int("MKEF_JOB_TIMEOUT_SECONDS", 30 * 60),
            job_ttl_seconds=env_int("MKEF_JOB_TTL_SECONDS", 24 * 60 * 60),
            cancel_grace_seconds=env_int("MKEF_CANCEL_GRACE_SECONDS", 5),
            allowed_origins=origins,
        )


@dataclass
class Job:
    id: uuid.UUID
    name: str
    directory: Path
    input_path: Path
    result_path: Path
    status: str = "queued"
    created_at: datetime = field(default_factory=utc_now)
    started_at: datetime | None = None
    finished_at: datetime | None = None
    error: dict[str, str] | None = None
    process: asyncio.subprocess.Process | None = None
    log: bytearray = field(default_factory=bytearray)

    def append_log(self, chunk: bytes, limit: int) -> None:
        self.log.extend(chunk)
        overflow = len(self.log) - limit
        if overflow > 0:
            del self.log[:overflow]

    def log_text(self) -> str:
        return bytes(self.log).decode("utf-8", errors="replace")


class JobCreate(BaseModel):
    model_config = ConfigDict(extra="forbid")

    name: str = Field(min_length=1, max_length=128)
    caseText: str = Field(min_length=1)

    @field_validator("name")
    @classmethod
    def validate_name(cls, value: str) -> str:
        value = value.strip()
        if not value:
            raise ValueError("name must contain a visible character")
        return value


class QueueFullError(RuntimeError):
    pass


class JobManager:
    def __init__(self, settings: Settings):
        self.settings = settings
        self.jobs: dict[uuid.UUID, Job] = {}
        self.queue: asyncio.Queue[uuid.UUID] = asyncio.Queue(
            maxsize=settings.max_active_jobs
        )
        self.worker_task: asyncio.Task[None] | None = None
        self.cleanup_task: asyncio.Task[None] | None = None
        self.octave_version: str | None = None
        self.octave_error: str | None = None

    async def start(self) -> None:
        self.settings.jobs_root.mkdir(parents=True, exist_ok=True)
        if self.settings.check_octave:
            await self._detect_octave()
        else:
            self.octave_version = "test-runtime"
        self.worker_task = asyncio.create_task(self._worker(), name="mkef-worker")
        self.cleanup_task = asyncio.create_task(self._cleanup_loop(), name="mkef-cleanup")

    async def close(self) -> None:
        for job in list(self.jobs.values()):
            if job.status in {"queued", "running"}:
                await self.cancel(job.id)
        for task in (self.worker_task, self.cleanup_task):
            if task:
                task.cancel()
        await asyncio.gather(
            *(task for task in (self.worker_task, self.cleanup_task) if task),
            return_exceptions=True,
        )

    async def _detect_octave(self) -> None:
        try:
            process = await asyncio.create_subprocess_exec(
                self.settings.octave_executable,
                "--version",
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.STDOUT,
            )
            output, _ = await asyncio.wait_for(process.communicate(), timeout=10)
            if process.returncode != 0:
                raise RuntimeError(f"version command exited with {process.returncode}")
            first_line = output.decode("utf-8", errors="replace").splitlines()[0]
            self.octave_version = first_line.strip()
        except Exception as exception:  # health endpoint reports the exact startup issue
            self.octave_error = str(exception)

    async def submit(self, name: str, case_text: str) -> Job:
        active = sum(
            job.status not in TERMINAL_STATES for job in self.jobs.values()
        )
        if active >= self.settings.max_active_jobs or self.queue.full():
            raise QueueFullError("The local solver queue is full.")

        job_id = uuid.uuid4()
        directory = self.settings.jobs_root / str(job_id)
        directory.mkdir(mode=0o700)
        job = Job(
            id=job_id,
            name=name.strip(),
            directory=directory,
            input_path=directory / "case.txt",
            result_path=directory / "result.json",
        )
        try:
            await asyncio.to_thread(job.input_path.write_text, case_text, encoding="utf-8")
        except Exception:
            shutil.rmtree(directory, ignore_errors=True)
            raise
        self.jobs[job_id] = job
        self.queue.put_nowait(job_id)
        return job

    def get(self, job_id: uuid.UUID) -> Job:
        job = self.jobs.get(job_id)
        if job is None:
            raise KeyError(job_id)
        return job

    async def cancel(self, job_id: uuid.UUID) -> Job:
        job = self.get(job_id)
        if job.status in TERMINAL_STATES:
            return job
        job.status = "canceled"
        job.error = {
            "code": "MKEF:Canceled",
            "message": "The calculation was canceled.",
        }
        job.finished_at = utc_now()
        if job.process and job.process.returncode is None:
            await self._terminate_process(job.process)
        return job

    async def _worker(self) -> None:
        while True:
            job_id = await self.queue.get()
            try:
                job = self.jobs.get(job_id)
                if job and job.status == "queued":
                    await self._run(job)
            finally:
                self.queue.task_done()

    async def _run(self, job: Job) -> None:
        job.status = "running"
        job.started_at = utc_now()
        environment = os.environ.copy()
        environment.update(
            {
                "HOME": "/tmp/mkef-home",
                "OCTAVE_HISTFILE": "/dev/null",
                "GNUTERM": "dumb",
            }
        )
        command = (
            self.settings.octave_executable,
            "--no-gui",
            "--quiet",
            "--no-history",
            "--norc",
            str(self.settings.runner_script),
            str(job.input_path),
            str(job.result_path),
            job.name,
        )
        capture_task: asyncio.Task[None] | None = None
        try:
            job.process = await asyncio.create_subprocess_exec(
                *command,
                cwd=self.settings.repository_root,
                env=environment,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.STDOUT,
                start_new_session=True,
            )
            capture_task = asyncio.create_task(self._capture_output(job))
            try:
                await asyncio.wait_for(
                    job.process.wait(), timeout=self.settings.job_timeout_seconds
                )
            except asyncio.TimeoutError:
                job.status = "timed_out"
                job.error = {
                    "code": "MKEF:SolverTimeout",
                    "message": (
                        "The calculation exceeded the configured "
                        f"{self.settings.job_timeout_seconds} second timeout."
                    ),
                }
                await self._terminate_process(job.process)

            if capture_task:
                await capture_task
            if job.status == "canceled" or job.status == "timed_out":
                return
            if job.process.returncode != 0:
                job.status = "failed"
                job.error = self._process_error(job)
                return
            await asyncio.to_thread(self._validate_result, job)
            job.status = "succeeded"
        except Exception as exception:
            if job.status not in {"canceled", "timed_out"}:
                job.status = "failed"
                job.error = {
                    "code": "MKEF:RunnerFailure",
                    "message": str(exception),
                }
        finally:
            if capture_task and not capture_task.done():
                capture_task.cancel()
                await asyncio.gather(capture_task, return_exceptions=True)
            job.process = None
            if job.finished_at is None:
                job.finished_at = utc_now()

    async def _capture_output(self, job: Job) -> None:
        assert job.process and job.process.stdout
        while True:
            chunk = await job.process.stdout.read(8192)
            if not chunk:
                return
            job.append_log(chunk, self.settings.max_log_bytes)

    def _process_error(self, job: Job) -> dict[str, str]:
        match = ERROR_PATTERN.search(job.log_text())
        if match:
            return {"code": match.group(1), "message": match.group(2).strip()}
        return {
            "code": "MKEF:SolverProcessFailed",
            "message": f"Octave exited with status {job.process.returncode}.",
        }

    async def _terminate_process(self, process: asyncio.subprocess.Process) -> None:
        if process.returncode is not None:
            return
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            return
        try:
            await asyncio.wait_for(
                process.wait(), timeout=self.settings.cancel_grace_seconds
            )
        except asyncio.TimeoutError:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                return
            await process.wait()

    def _validate_result(self, job: Job) -> None:
        if not job.result_path.is_file():
            raise RuntimeError("Octave completed without creating result.json.")
        size = job.result_path.stat().st_size
        if size > self.settings.max_result_bytes:
            job.result_path.unlink(missing_ok=True)
            raise RuntimeError(
                "The generated result exceeds the configured "
                f"{self.settings.max_result_bytes} byte limit."
            )
        with job.result_path.open("r", encoding="utf-8") as result_file:
            result = json.load(result_file)
        if not isinstance(result, dict):
            raise RuntimeError("The generated result is not a JSON object.")
        if result.get("format") != "mkef-postprocessor" or result.get("version") not in (1, 2, 3):
            raise RuntimeError("The generated result has an unsupported schema.")

    async def _cleanup_loop(self) -> None:
        interval = min(60, max(1, self.settings.job_ttl_seconds // 4))
        while True:
            await asyncio.sleep(interval)
            threshold = utc_now().timestamp() - self.settings.job_ttl_seconds
            expired = [
                job_id
                for job_id, job in self.jobs.items()
                if job.status in TERMINAL_STATES
                and job.finished_at
                and job.finished_at.timestamp() < threshold
            ]
            for job_id in expired:
                job = self.jobs.pop(job_id, None)
                if job:
                    await asyncio.to_thread(shutil.rmtree, job.directory, True)

    def public_status(self, job: Job) -> dict[str, Any]:
        end = job.finished_at or utc_now()
        elapsed = 0.0
        if job.started_at:
            elapsed = max(0.0, (end - job.started_at).total_seconds())
        return {
            "id": str(job.id),
            "name": job.name,
            "status": job.status,
            "createdAt": job.created_at.isoformat(),
            "startedAt": job.started_at.isoformat() if job.started_at else None,
            "finishedAt": job.finished_at.isoformat() if job.finished_at else None,
            "elapsedSeconds": round(elapsed, 3),
            "error": job.error,
            "logTail": job.log_text(),
        }


def safe_download_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._ -]+", "_", name).strip(" .")
    return (cleaned or "MKE-F-result") + ".json"


def create_app(settings: Settings | None = None) -> FastAPI:
    active_settings = settings or Settings.from_env()
    manager = JobManager(active_settings)

    @asynccontextmanager
    async def lifespan(application: FastAPI):
        await manager.start()
        application.state.job_manager = manager
        try:
            yield
        finally:
            await manager.close()

    application = FastAPI(
        title="MKE-F Solver API",
        version="1.0.0",
        docs_url=None,
        redoc_url=None,
        openapi_url=None,
        lifespan=lifespan,
    )
    application.state.job_manager = manager
    application.add_middleware(
        TrustedHostMiddleware,
        allowed_hosts=["solver", "localhost", "127.0.0.1", "testserver"],
    )

    @application.middleware("http")
    async def reject_cross_origin(request: Request, call_next):
        origin = request.headers.get("origin")
        if (
            request.method in {"POST", "PUT", "PATCH", "DELETE"}
            and origin
            and origin not in active_settings.allowed_origins
        ):
            return JSONResponse(
                status_code=status.HTTP_403_FORBIDDEN,
                content={"detail": "Cross-origin requests are not allowed."},
            )
        return await call_next(request)

    def resolve_job(job_id: uuid.UUID) -> Job:
        try:
            return manager.get(job_id)
        except KeyError:
            raise HTTPException(status_code=404, detail="Job not found or expired.")

    @application.get("/api/v1/health")
    async def health():
        if manager.octave_version:
            return {"status": "ready", "octaveVersion": manager.octave_version}
        return JSONResponse(
            status_code=503,
            content={
                "status": "unavailable",
                "octaveVersion": None,
                "error": manager.octave_error or "Octave is unavailable.",
            },
        )

    @application.post("/api/v1/jobs", status_code=202)
    async def create_job(payload: JobCreate):
        encoded_size = len(payload.caseText.encode("utf-8"))
        if encoded_size > active_settings.max_case_bytes:
            raise HTTPException(
                status_code=413,
                detail=(
                    "Case text exceeds the configured "
                    f"{active_settings.max_case_bytes} byte limit."
                ),
            )
        if "\x00" in payload.caseText:
            raise HTTPException(status_code=422, detail="Case text contains a null byte.")
        try:
            job = await manager.submit(payload.name, payload.caseText)
        except QueueFullError as exception:
            raise HTTPException(status_code=429, detail=str(exception))
        return manager.public_status(job)

    @application.get("/api/v1/jobs/{job_id}")
    async def get_job(job_id: uuid.UUID):
        return manager.public_status(resolve_job(job_id))

    @application.post("/api/v1/jobs/{job_id}/cancel")
    async def cancel_job(job_id: uuid.UUID):
        try:
            job = await manager.cancel(job_id)
        except KeyError:
            raise HTTPException(status_code=404, detail="Job not found or expired.")
        return manager.public_status(job)

    @application.get("/api/v1/jobs/{job_id}/result")
    async def get_result(job_id: uuid.UUID):
        job = resolve_job(job_id)
        if job.status != "succeeded":
            raise HTTPException(
                status_code=409,
                detail={"status": job.status, "error": job.error},
            )
        return FileResponse(
            job.result_path,
            media_type="application/json",
            filename=safe_download_name(job.name),
        )

    return application


app = create_app()
