import json
import subprocess
import time
from pathlib import Path

import pytest
from fastapi.testclient import TestClient
from backend.app import Settings, create_app


REPOSITORY_ROOT = Path("/app")
RUNNER = REPOSITORY_ROOT / "scripts" / "run_solver_job.m"


@pytest.mark.parametrize(
    ("case_name", "analysis_type"),
    [
        ("CasePreprocessorStatic.txt", "static"),
        ("CasePreprocessorModal.txt", "modal"),
        ("CasePreprocessorTransient.txt", "transient"),
    ],
)
def test_fixed_octave_runner_exports_each_analysis(tmp_path, case_name, analysis_type):
    if not RUNNER.is_file():
        pytest.skip("container Octave runner is not available")
    result_path = tmp_path / f"{analysis_type}.json"
    completed = subprocess.run(
        [
            "octave-cli",
            "--no-gui",
            "--quiet",
            "--no-history",
            "--norc",
            str(RUNNER),
            str(REPOSITORY_ROOT / "examples" / "cases" / case_name),
            str(result_path),
            f"Integration {analysis_type}",
        ],
        cwd=REPOSITORY_ROOT,
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    result = json.loads(result_path.read_text(encoding="utf-8"))
    assert result["format"] == "mkef-postprocessor"
    assert result["version"] == 1
    assert result["analysis"]["type"] == analysis_type
    assert result["metadata"]["title"] == f"Integration {analysis_type}"


def test_uniform_load_through_real_api(tmp_path):
    if not RUNNER.is_file():
        pytest.skip("container Octave runner is not available")
    settings = Settings(jobs_root=tmp_path / "jobs", repository_root=REPOSITORY_ROOT,
                        runner_script=RUNNER, job_timeout_seconds=60)
    text = (REPOSITORY_ROOT / "examples/cases/CaseUniformFrame.txt").read_text()
    with TestClient(create_app(settings)) as client:
        response = client.post("/api/v1/jobs", json={"name": "Uniform beam", "caseText": text})
        assert response.status_code == 202, response.text
        job_id = response.json()["id"]
        deadline = time.monotonic() + 65
        while time.monotonic() < deadline:
            job = client.get(f"/api/v1/jobs/{job_id}").json()
            if job["status"] in {"succeeded", "failed", "timed_out", "canceled"}:
                break
            time.sleep(0.05)
        assert job["status"] == "succeeded", job
        result = client.get(f"/api/v1/jobs/{job_id}/result").json()
        assert result["version"] == 2
        assert result["model"]["nodalLoads"] == []
        assert len(result["model"]["elementLoads"]) == 1
        assert result["analysis"]["reactions"][:3] == pytest.approx([0, 2000, 2000])
        assert result["analysis"]["elementResults"][0]["localEndForces"][3:] == pytest.approx([0, 0, 0], abs=1e-8)
