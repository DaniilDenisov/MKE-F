import json
import math
import re
import subprocess
import time
from pathlib import Path

import pytest
from fastapi.testclient import TestClient
from backend.app import Settings, create_app


REPOSITORY_ROOT = Path("/app")
RUNNER = REPOSITORY_ROOT / "scripts" / "run_solver_job.m"


@pytest.mark.parametrize("section_type,invalid_wall", [
    ("rectangle", False),
    *[(kind, invalid) for kind in ("square_tube", "rectangular_tube", "circular_tube")
      for invalid in (False, True)],
])
def test_sections_through_real_api(tmp_path, section_type, invalid_wall):
    if not RUNNER.is_file():
        pytest.skip("container Octave runner is not available")
    settings = Settings(jobs_root=tmp_path / "jobs", repository_root=REPOSITORY_ROOT,
                        runner_script=RUNNER)
    text = (REPOSITORY_ROOT / "examples/cases/CaseSectionBeam.txt").read_text()
    dimensions = {"b": .1, "h": .2}
    area, inertia, extent = .02, .1 * .2**3 / 12, .1
    if section_type == "circular_tube":
        dimensions = {"d": .1, "t": .005}
        area = math.pi * (.1**2 - .09**2) / 4
        inertia, extent = math.pi * (.1**4 - .09**4) / 64, .05
    elif section_type in {"square_tube", "rectangular_tube"}:
        height = .1 if section_type == "square_tube" else .2
        dimensions = {"b": .1, "t": .005}
        if section_type == "rectangular_tube":
            dimensions["h"] = height
        area = .1 * height - .09 * (height - .01)
        inertia, extent = (.1 * height**3 - .09 * (height - .01)**3) / 12, height / 2
    if invalid_wall:
        dimensions["t"] = .05
    section = {"id": 1, "name": section_type, "type": section_type, "dimensions": dimensions}
    text = re.sub(r'\{"id":1[^\r\n]+', lambda _: json.dumps(section), text)
    with TestClient(create_app(settings)) as client:
        response = client.post("/api/v1/jobs", json={"name": "Sections", "caseText": text})
        assert response.status_code == 202, response.text
        path = "/api/v1/jobs/" + response.json()["id"]
        deadline = time.monotonic() + 65
        while time.monotonic() < deadline:
            job = client.get(path).json()
            if job["status"] in {"succeeded", "failed", "timed_out", "canceled"}:
                break
            time.sleep(0.05)
        if invalid_wall:
            assert job["status"] == "failed", job
            assert job["error"]["code"] == "MKEF:InvalidSection"
            assert client.get(path + "/result").status_code == 409
        else:
            assert job["status"] == "succeeded", job
            result = client.get(path + "/result").json()
            assert result["version"] == 6
            assert result["model"]["sections"] == [section]
            properties = result["model"]["elements"][0]["properties"]
            assert properties["area"] == pytest.approx(area)
            assert properties["momentOfInertia"] == pytest.approx(inertia)
            stress = result["analysis"]["elementResults"][0]["sectionStress"]
            assert stress["maximum"]["value"] == pytest.approx(500 * extent / inertia)
            assert stress["minimum"]["value"] == pytest.approx(-500 * extent / inertia)
            assert stress["maximum"]["xi"] == pytest.approx(.5)


@pytest.mark.parametrize("case_name", ["CaseReleaseStatic", "CaseReleaseModal", "CaseReleaseTransient"])
def test_releases_through_real_api(tmp_path, case_name):
    if not RUNNER.is_file():
        pytest.skip("container Octave runner is not available")
    settings = Settings(jobs_root=tmp_path / "jobs", repository_root=REPOSITORY_ROOT,
                        runner_script=RUNNER, job_timeout_seconds=60)
    text = (REPOSITORY_ROOT / "examples/cases" / f"{case_name}.txt").read_text()
    with TestClient(create_app(settings)) as client:
        response = client.post("/api/v1/jobs", json={"name": case_name, "caseText": text})
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
        assert result["version"] == 4
        model = result["model"]
        assert len(model["dofRegistry"]) == model["numberOfDOFs"]
        assert len(model["releases"]) == 2
        assert model["dofMap"][0][2] == 0
        if case_name == "CaseReleaseStatic":
            assert result["analysis"]["elementResults"][0]["localEndForces"][2::3] == pytest.approx([0, 0], abs=1e-8)


@pytest.mark.parametrize("case_name", ["CaseMPCFrame", "CaseMPCModal", "CaseMPCTransient", "CaseMPCUniform"])
def test_mpc_through_real_api(tmp_path, case_name):
    if not RUNNER.is_file():
        pytest.skip("container Octave runner is not available")
    settings = Settings(jobs_root=tmp_path / "jobs", repository_root=REPOSITORY_ROOT,
                        runner_script=RUNNER, job_timeout_seconds=60)
    case_text = (REPOSITORY_ROOT / "examples" / "cases" / f"{case_name}.txt").read_text()
    with TestClient(create_app(settings)) as client:
        submitted = client.post("/api/v1/jobs", json={"name": case_name, "caseText": case_text})
        assert submitted.status_code == 202, submitted.text
        job_id = submitted.json()["id"]
        deadline = time.monotonic() + 65
        while time.monotonic() < deadline:
            job = client.get(f"/api/v1/jobs/{job_id}").json()
            if job["status"] in {"succeeded", "failed", "timed_out", "canceled"}:
                break
            time.sleep(0.05)
        assert job["status"] == "succeeded", job
        result = client.get(f"/api/v1/jobs/{job_id}/result").json()
        assert result["version"] == 3
        assert len(result["model"]["mpcs"]) == 1
        assert isinstance(result["analysis"]["mpcMultipliers"], list)
        assert "supportReactions" in result["analysis"]
        assert "mpcForces" in result["analysis"]
        if case_name == "CaseMPCUniform":
            assert len(result["model"]["elementLoads"]) == 1


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


@pytest.mark.parametrize(("case_name", "version", "reactions"), [
    ("CaseUniformFrame", 2, [0, 2000, 2000]),
    ("CaseLinearFrame", 5, [0, 1200, 1400]),
    ("CaseTriangleFrame", 5, [0, 1000, 4000/3]),
    ("CaseLinearRelease", 5, None),
    ("CaseLinearMPC", 5, None),
    ("CaseSectionBeam", 6, None),
])
def test_distributed_load_through_real_api(tmp_path, case_name, version, reactions):
    if not RUNNER.is_file():
        pytest.skip("container Octave runner is not available")
    settings = Settings(jobs_root=tmp_path / "jobs", repository_root=REPOSITORY_ROOT,
                        runner_script=RUNNER, job_timeout_seconds=60)
    text = (REPOSITORY_ROOT / "examples/cases" / f"{case_name}.txt").read_text()
    with TestClient(create_app(settings)) as client:
        response = client.post("/api/v1/jobs", json={"name": case_name, "caseText": text})
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
        assert result["version"] == version
        if version == 6:
            assert result["model"]["sections"][0]["type"] == "rectangle"
            stress = result["analysis"]["elementResults"][0]["sectionStress"]
            assert stress["maximum"]["value"] == pytest.approx(750000)
            assert stress["minimum"]["value"] == pytest.approx(-750000)
            assert stress["maximum"]["xi"] == pytest.approx(0.5)
        if case_name != "CaseLinearMPC":
            assert result["model"]["nodalLoads"] == []
        assert len(result["model"]["elementLoads"]) == 1
        assert result["analysis"]["equilibriumResidual"] == pytest.approx([0, 0, 0], abs=1e-8)
        if reactions is not None:
            assert result["analysis"]["reactions"][:3] == pytest.approx(reactions)
            assert result["analysis"]["elementResults"][0]["localEndForces"][3:] == pytest.approx([0, 0, 0], abs=1e-8)
        if version == 5:
            assert "qx" not in result["model"]["elementLoads"][0]
            assert "qy2" in result["model"]["elementLoads"][0]
        if case_name == "CaseLinearRelease":
            assert len(result["model"]["dofRegistry"]) == result["model"]["numberOfDOFs"]
            assert result["analysis"]["elementResults"][0]["localEndForces"][2::3] == pytest.approx([0, 0], abs=1e-8)
        if case_name == "CaseLinearMPC":
            assert len(result["model"]["mpcs"]) == 1
            assert len(result["analysis"]["mpcMultipliers"]) == 1
