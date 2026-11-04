"""The downloadable 3D guide must reproduce its stated values through the API."""
import time

import pytest
from fastapi.testclient import TestClient

from backend.app import Settings, create_app
from backend.tests.test_spatial_runner import ROOT, RUNNER


def test_reference_quick_start_api(tmp_path):
    if not RUNNER.is_file():
        pytest.skip('requires container Octave runtime')
    settings = Settings(jobs_root=tmp_path / 'jobs', repository_root=ROOT,
                        runner_script=RUNNER, job_timeout_seconds=60)
    text = (ROOT / 'examples/cases/CaseSpatialQuickStart.txt').read_text()
    with TestClient(create_app(settings)) as client:
        response = client.post('/api/v1/jobs', json={'name': '3D quick start', 'caseText': text})
        assert response.status_code == 202
        url = '/api/v1/jobs/' + response.json()['id']
        deadline = time.monotonic() + 65
        while time.monotonic() < deadline:
            job = client.get(url).json()
            if job['status'] in {'succeeded', 'failed', 'timed_out', 'canceled'}:
                break
            time.sleep(.05)
        assert job['status'] == 'succeeded', job
        data = client.get(url + '/result').json()
        assert data['version'] == 8
        assert data['model']['elements'][0]['sectionId'] == 1
        result = data['analysis']
        assert result['displacements'][7] == pytest.approx(-.0002)
        assert result['nodalSupportReactions'][1] == pytest.approx(1000)
        stress = result['elementResults'][0]['sectionStress']
        assert stress['minimum']['value'] == pytest.approx(-3e6)
        assert stress['maximum']['value'] == pytest.approx(3e6)
