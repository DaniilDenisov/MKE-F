import stat
import time
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from backend.app import Settings, create_app


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
