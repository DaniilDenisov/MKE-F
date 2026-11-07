"""Real Octave/API coverage, including singleton JSON dimensions."""
import json
import time
from pathlib import Path

import pytest
from fastapi.testclient import TestClient
from backend.app import Settings, create_app

ROOT = Path('/app')
RUNNER = ROOT / 'scripts/runtime/run_solver_job.m'
CASES = ['CaseSpatialTruss', 'CaseSpatialFrame', 'CaseSpatialTrussModal',
         'CaseSpatialFrameModal', 'CaseSpatialTrussTransient', 'CaseSpatialFrameTransient']


@pytest.mark.parametrize('name', CASES)
@pytest.mark.parametrize('with_section', [False, True])
def test_spatial_api(tmp_path, name, with_section):
    if not RUNNER.is_file():
        pytest.skip('requires container Octave runtime')
    settings = Settings(jobs_root=tmp_path / 'jobs', repository_root=ROOT,
                        runner_script=RUNNER, job_timeout_seconds=60)
    text = (ROOT / f'examples/cases/{name}.txt').read_text()
    section = {'id': 1, 'name': 'API section', 'type': 'channel' if 'Truss' in name else 'rectangle',
               'dimensions': {'b': .1, 'h': .2, 'tw': .008, 'tf': .012},
               'torsion': {'mode': 'manual', 'value': 1e-6}}
    if with_section:
        text += '\nsections\n1\n' + json.dumps(section) + '\nelement_sections\n1\n1,1,-37.5\n'
    with TestClient(create_app(settings)) as client:
        response = client.post('/api/v1/jobs', json={
            'name': name, 'caseText': text})
        assert response.status_code == 202, response.text
        path = '/api/v1/jobs/' + response.json()['id']
        deadline = time.monotonic() + 65
        while time.monotonic() < deadline:
            job = client.get(path).json()
            if job['status'] in {'succeeded', 'failed', 'timed_out', 'canceled'}:
                break
            time.sleep(.05)
        assert job['status'] == 'succeeded', job
        response = client.get(path + '/result')
        assert response.status_code == 200
        data = response.json()
        assert data['version'] == 8
        m, a = data['model'], data['analysis']
        assert m['dimension'] == 3
        assert len(m['nodes']) == 2 and len(m['elements']) == 1
        assert m['elements'][0]['type'] == (1133 if 'Truss' in name else 1436)
        assert len(m['dofRegistry']) == m['numberOfDOFs']
        assert isinstance(m['mpcs'], list) and isinstance(m['supports'], list)
        assert isinstance(m['nodeKinematics']['entries'][0], list)
        assert 'stiffness' not in m and 'mass' not in m
        if with_section:
            assert m['sections'] == [section]
            assert m['elements'][0]['sectionId'] == 1
            assert m['elements'][0]['sectionAngle'] == -37.5
            assert m['elements'][0]['properties']['area'] == pytest.approx(.003808 if 'Truss' in name else .02)
            if 'Frame' in name:
                assert m['elements'][0]['properties']['torsionConstant'] == 1e-6
        if a['type'] == 'modal':
            assert isinstance(a['frequenciesHz'], list)
            assert isinstance(a['modeShapes'][0], list)
            assert all(f > 0 for f in a['frequenciesHz'])
        elif a['type'] == 'static':
            assert len(a['equilibriumResidual']) == 6
            assert max(map(abs, a['equilibriumResidual'])) < 1e-8
            if with_section:
                assert isinstance(a['elementResults'][0]['sectionStress']['minimum'], dict)
            if 'Truss' in name and not with_section:
                assert a['displacements'][5] == pytest.approx(1e-6)
        else:
            assert isinstance(a['displacements'][0], list)
            assert len(a['displacements'][0]) == len(a['time'])


@pytest.mark.parametrize('kind', ['channel', 'tee', 'ibeam'])
@pytest.mark.parametrize('suffix', ['', 'Modal', 'Transient'])
def test_open_frame_section_rejected_through_api(tmp_path, kind, suffix):
    if not RUNNER.is_file():
        pytest.skip('requires container Octave runtime')
    settings = Settings(jobs_root=tmp_path / 'jobs', repository_root=ROOT,
                        runner_script=RUNNER, job_timeout_seconds=60)
    text = (ROOT / f'examples/cases/CaseSpatialFrame{suffix}.txt').read_text()
    section = {'id': 1, 'name': kind, 'type': kind,
               'dimensions': {'b': .1, 'h': .2, 'tw': .008, 'tf': .012}}
    text += '\nsections\n1\n' + json.dumps(section) + '\nelement_sections\n1\n1,1,0\n'
    with TestClient(create_app(settings)) as client:
        response = client.post('/api/v1/jobs', json={'name': kind, 'caseText': text})
        assert response.status_code == 202
        path = '/api/v1/jobs/' + response.json()['id']
        deadline = time.monotonic() + 65
        while time.monotonic() < deadline:
            job = client.get(path).json()
            if job['status'] in {'succeeded', 'failed', 'timed_out', 'canceled'}:
                break
            time.sleep(.05)
        assert job['status'] == 'failed', job
        assert job['error']['code'] == 'MKEF:UnsupportedFrameSection'
        assert f'1 ({kind})' in job['error']['message']
        assert client.get(path + '/result').status_code == 409


def test_spatial_singleton_partial_arrays():
    path = ROOT / 'output/CaseSpatialTrussTransientPartial.json'
    if not path.is_file():
        pytest.skip('run Octave spatial export fixtures first')
    data = json.loads(path.read_text())
    a = data['analysis']
    assert isinstance(a['globalDOFIds'], list) and len(a['globalDOFIds']) == 1
    assert len(a['displacements']) == 1 and isinstance(a['displacements'][0], list)
    assert a['mpcMultipliers'] == []
