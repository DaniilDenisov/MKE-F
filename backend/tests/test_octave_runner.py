import json
import subprocess
from pathlib import Path

import pytest


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
