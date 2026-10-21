import asyncio
import stat
import time
import uuid
from datetime import timedelta
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from backend.app import Job, JobManager, LeasedFileResponse, Settings, create_app, utc_now


FAKE_RUNNER = r'''#!/usr/bin/env python3
import json
import pathlib
import sys
import time

input_path = pathlib.Path(sys.argv[-3])
output_path = pathlib.Path(sys.argv[-2])
text = input_path.read_text(encoding="utf-8")
if "SLOW" in text:
    time.sleep(5)
if "FAIL" in text:
    print("MKEF_ERROR\tMKEF:FakeFailure\tRequested fake failure.", flush=True)
    raise SystemExit(7)
if "NOISY" in text:
    print("x" * 4096, flush=True)
payload = {"format": "mkef-postprocessor", "version": 1, "metadata": {"title": sys.argv[-1]}}
if "LARGE" in text:
    payload["padding"] = "y" * 4096
output_path.write_text(json.dumps(payload), encoding="utf-8")
'''


@pytest.fixture
def client_factory(tmp_path):
    clients = []

    def factory(**overrides):
        runner = tmp_path / "fake-octave"
        runner.write_text(FAKE_RUNNER, encoding="utf-8")
        runner.chmod(runner.stat().st_mode | stat.S_IXUSR)
        jobs = tmp_path / ("jobs-" + str(len(clients)))
        settings = Settings(
            jobs_root=jobs,
            repository_root=tmp_path,
            octave_executable=str(runner),
            runner_script=tmp_path / "ignored-runner.m",
            max_case_bytes=overrides.get("max_case_bytes", 1024),
            max_result_bytes=overrides.get("max_result_bytes", 1024 * 1024),
            max_log_bytes=overrides.get("max_log_bytes", 1024),
            max_active_jobs=overrides.get("max_active_jobs", 20),
            job_timeout_seconds=overrides.get("job_timeout_seconds", 10),
            job_ttl_seconds=overrides.get("job_ttl_seconds", 3600),
            cancel_grace_seconds=1,
            allowed_origins=("http://localhost:8080",),
            check_octave=False,
        )
        context = TestClient(create_app(settings))
        client = context.__enter__()
        clients.append(context)
        return client

    yield factory
    for context in reversed(clients):
        context.__exit__(None, None, None)


def submit(client, case_text="OK", name="Case"):
    response = client.post("/api/v1/jobs", json={"name": name, "caseText": case_text})
    assert response.status_code == 202, response.text
    return response.json()["id"]


def wait_for_terminal(client, job_id, timeout=4):
    deadline = time.time() + timeout
    while time.time() < deadline:
        response = client.get(f"/api/v1/jobs/{job_id}")
        assert response.status_code == 200
        job = response.json()
        if job["status"] in {"succeeded", "failed", "canceled", "timed_out"}:
            return job
        time.sleep(0.03)
    raise AssertionError("job did not reach a terminal state")


def test_health_and_successful_result(client_factory):
    client = client_factory()
    assert client.get("/api/v1/health").json()["status"] == "ready"
    job_id = submit(client, name="Demo frame")
    job = wait_for_terminal(client, job_id)
    assert job["status"] == "succeeded"
    result = client.get(f"/api/v1/jobs/{job_id}/result")
    assert result.status_code == 200
    assert result.json()["format"] == "mkef-postprocessor"
    assert "Demo%20frame.json" in result.headers["content-disposition"]


def test_request_rejects_extra_fields_nulls_size_and_cross_origin(client_factory):
    client = client_factory(max_case_bytes=4)
    assert client.post("/api/v1/jobs", json={"name": "x", "caseText": "ok", "script": "evil.m"}).status_code == 422
    assert client.post("/api/v1/jobs", json={"name": "x", "caseText": "12345"}).status_code == 413
    assert client.post("/api/v1/jobs", json={"name": "x", "caseText": "a\u0000"}).status_code == 422
    assert client.post(
        "/api/v1/jobs",
        json={"name": "x", "caseText": "ok"},
        headers={"Origin": "https://attacker.invalid"},
    ).status_code == 403


def test_process_failure_is_structured_and_result_is_not_ready(client_factory):
    client = client_factory()
    job_id = submit(client, "FAIL")
    job = wait_for_terminal(client, job_id)
    assert job["status"] == "failed"
    assert job["error"] == {"code": "MKEF:FakeFailure", "message": "Requested fake failure."}
    assert client.get(f"/api/v1/jobs/{job_id}/result").status_code == 409


def test_timeout_and_cancel_are_terminal(client_factory):
    timeout_client = client_factory(job_timeout_seconds=1)
    timeout_id = submit(timeout_client, "SLOW")
    assert wait_for_terminal(timeout_client, timeout_id)["status"] == "timed_out"

    cancel_client = client_factory()
    cancel_id = submit(cancel_client, "SLOW")
    deadline = time.time() + 2
    while time.time() < deadline:
        if cancel_client.get(f"/api/v1/jobs/{cancel_id}").json()["status"] == "running":
            break
        time.sleep(0.02)
    response = cancel_client.post(f"/api/v1/jobs/{cancel_id}/cancel")
    assert response.status_code == 200
    assert response.json()["status"] == "canceled"
    assert cancel_client.post(f"/api/v1/jobs/{cancel_id}/cancel").json()["status"] == "canceled"


def test_fifo_queue_limit_and_queued_cancel(client_factory):
    client = client_factory(max_active_jobs=2)
    first = submit(client, "SLOW")
    second = submit(client, "OK")
    assert client.post("/api/v1/jobs", json={"name": "third", "caseText": "OK"}).status_code == 429
    canceled = client.post(f"/api/v1/jobs/{second}/cancel").json()
    assert canceled["status"] == "canceled"
    client.post(f"/api/v1/jobs/{first}/cancel")


def test_log_is_bounded_and_large_result_fails(client_factory):
    client = client_factory(max_log_bytes=64, max_result_bytes=512)
    noisy_id = submit(client, "NOISY")
    noisy = wait_for_terminal(client, noisy_id)
    assert noisy["status"] == "succeeded"
    assert len(noisy["logTail"].encode("utf-8")) <= 64

    large_id = submit(client, "LARGE")
    large = wait_for_terminal(client, large_id)
    assert large["status"] == "failed"
    assert "exceeds" in large["error"]["message"]


def test_terminal_jobs_expire(client_factory):
    client = client_factory(job_ttl_seconds=1)
    job_id = submit(client)
    assert wait_for_terminal(client, job_id)["status"] == "succeeded"
    deadline = time.time() + 3
    while time.time() < deadline:
        if client.get(f"/api/v1/jobs/{job_id}").status_code == 404:
            return
        time.sleep(0.1)
    raise AssertionError("expired job was not removed")


def test_repeated_queued_cancellation_releases_slots(client_factory):
    client = client_factory(max_active_jobs=2)
    first = submit(client, "SLOW")
    for _ in range(25):
        waiting = submit(client)
        assert client.post(f"/api/v1/jobs/{waiting}/cancel").json()["status"] == "canceled"
    last = submit(client)
    client.post(f"/api/v1/jobs/{first}/cancel")
    assert wait_for_terminal(client, last)["status"] == "succeeded"


@pytest.mark.parametrize("shutdown", [False, True])
def test_cancel_during_spawn_terminates_child(client_factory, monkeypatch, shutdown):
    client = client_factory()
    manager = client.app.state.job_manager
    original = asyncio.create_subprocess_exec

    async def scenario():
        spawned, release = asyncio.Event(), asyncio.Event()
        children = []

        async def delayed_spawn(*args, **kwargs):
            process = await original(*args, **kwargs)
            children.append(process)
            spawned.set()
            await release.wait()
            return process

        monkeypatch.setattr(asyncio, "create_subprocess_exec", delayed_spawn)
        job = await manager.submit("race", "SLOW")
        closing = None
        try:
            await asyncio.wait_for(spawned.wait(), 3)
            if shutdown:
                closing = asyncio.create_task(manager.close())
                await asyncio.sleep(0)
                assert not closing.done()
            else:
                await manager.cancel(job.id)
                assert job.executing and job.finished_at is None
                await manager.cleanup()
                assert job.directory.is_dir()
        finally:
            release.set()
        await asyncio.wait_for(closing if closing else manager.queue.join(), 3)
        assert job.status == "canceled"
        assert not job.executing and job.finished_at is not None
        assert children[0].returncode is not None
        if not shutdown:
            assert manager.healthy()

    client.portal.call(scenario)


def test_cancel_during_result_validation_stays_canceled(client_factory, monkeypatch):
    client = client_factory()
    manager = client.app.state.job_manager
    original = asyncio.to_thread

    async def scenario():
        validating, release = asyncio.Event(), asyncio.Event()

        async def delayed_validation(function, *args, **kwargs):
            if function == manager._validate_result:
                validating.set()
                await release.wait()
            return await original(function, *args, **kwargs)

        monkeypatch.setattr(asyncio, "to_thread", delayed_validation)
        job = await manager.submit("validation race", "OK")
        try:
            await asyncio.wait_for(validating.wait(), 3)
            await manager.cancel(job.id)
        finally:
            release.set()
        await asyncio.wait_for(manager.queue.join(), 3)
        assert job.status == "canceled"
        assert job.error["code"] == "MKEF:Canceled"

    client.portal.call(scenario)


def test_capture_failure_terminates_child_and_worker_recovers(client_factory, monkeypatch):
    client = client_factory()
    manager = client.app.state.job_manager
    original = manager._capture_output

    async def scenario():
        children = []

        async def fail_capture(job):
            children.append(job.process)
            raise RuntimeError("capture failed")

        monkeypatch.setattr(manager, "_capture_output", fail_capture)
        job = await manager.submit("capture failure", "SLOW")
        await asyncio.wait_for(manager.queue.join(), 3)
        assert job.status == "failed"
        assert children[0].returncode is not None
        assert manager.healthy()
        monkeypatch.setattr(manager, "_capture_output", original)
        next_job = await manager.submit("recovery", "OK")
        await asyncio.wait_for(manager.queue.join(), 3)
        assert next_job.status == "succeeded"

    client.portal.call(scenario)


def test_timeout_also_covers_output_capture(client_factory, monkeypatch):
    client = client_factory(job_timeout_seconds=1)
    manager = client.app.state.job_manager

    async def stalled_capture(job):
        await asyncio.Event().wait()

    monkeypatch.setattr(manager, "_capture_output", stalled_capture)
    job_id = submit(client)
    assert wait_for_terminal(client, job_id)["status"] == "timed_out"


@pytest.mark.parametrize("task_name", ["worker_task", "cleanup_task"])
def test_health_detects_stopped_background_task(client_factory, task_name):
    client = client_factory()
    manager = client.app.state.job_manager

    async def stop_task():
        task = getattr(manager, task_name)
        task.cancel()
        await asyncio.gather(task, return_exceptions=True)

    client.portal.call(stop_task)
    assert client.get("/api/v1/health").status_code == 503
    assert client.post("/api/v1/jobs", json={"name": "x", "caseText": "OK"}).status_code == 429


def test_concurrent_submission_preserves_existing_capacity(client_factory):
    client = client_factory(max_active_jobs=2)
    manager = client.app.state.job_manager

    async def scenario():
        from backend.app import QueueFullError
        results = await asyncio.gather(
            *(manager.submit(str(i), "SLOW") for i in range(12)), return_exceptions=True
        )
        assert sum(isinstance(item, Job) for item in results) == 2
        assert sum(isinstance(item, QueueFullError) for item in results) == 10

    client.portal.call(scenario)


@pytest.mark.parametrize("abort", [False, True])
def test_download_lease_survives_cleanup_and_releases_on_abort(tmp_path, abort):
    async def scenario():
        manager = JobManager(Settings(jobs_root=tmp_path, job_ttl_seconds=1))
        directory = tmp_path / "result"
        directory.mkdir()
        result = directory / "result.json"
        result.write_text('{"ok": true}')
        job = Job(uuid.uuid4(), "download", directory, directory / "case.txt", result,
                  status="succeeded", finished_at=utc_now() - timedelta(seconds=10))
        manager.jobs[job.id] = job
        response = LeasedFileResponse(job, media_type="application/json")
        chunks = []

        async def send(message):
            await manager.cleanup()
            assert result.is_file() and manager.get(job.id) is job
            if message["type"] == "http.response.body":
                if abort:
                    raise ConnectionError("client disconnected")
                chunks.append(message.get("body", b""))

        async def receive():
            return {"type": "http.disconnect"}

        if abort:
            with pytest.raises(ConnectionError):
                await response({"type": "http", "method": "GET", "headers": []}, receive, send)
        else:
            await response({"type": "http", "method": "GET", "headers": []}, receive, send)
            assert b''.join(chunks) == b'{"ok": true}'
        assert job.readers == 0
        await manager.cleanup()
        assert not directory.exists() and job.id not in manager.jobs

    asyncio.run(scenario())


def test_failed_cleanup_retries_without_losing_job(tmp_path, monkeypatch):
    async def scenario():
        manager = JobManager(Settings(jobs_root=tmp_path, job_ttl_seconds=1))
        directory = tmp_path / "expired"
        directory.mkdir()
        job = Job(uuid.uuid4(), "retry", directory, directory / "case.txt", directory / "result.json",
                  status="failed", finished_at=utc_now() - timedelta(seconds=10))
        manager.jobs[job.id] = job
        import shutil
        original = shutil.rmtree

        def fail_once(*args):
            raise PermissionError("busy")

        monkeypatch.setattr(shutil, "rmtree", fail_once)
        await manager.cleanup()
        assert manager.get(job.id) is job and not job.deleting
        monkeypatch.setattr(shutil, "rmtree", original)
        await manager.cleanup()
        assert job.id not in manager.jobs and not directory.exists()

    asyncio.run(scenario())
