(function (M) {
  'use strict';
  var tests = [];
  function test(name, operation) { tests.push({ name: name, operation: operation }); }
  function assert(condition, message) { if (!condition) throw new Error(message || 'Assertion failed'); }
  function close(actual, expected, tolerance) { assert(Math.abs(actual - expected) <= (tolerance || 1e-10), actual + ' != ' + expected); }
  function fixture() {
    return {
      format: 'mkef-postprocessor', version: 1,
      metadata: { title: 'Browser test', units: { length: 'm', force: 'N', moment: 'N*m', stress: 'Pa', time: 's' } },
      model: {
        dimension: 2, dofPerNode: 3, dofLabels: ['ux', 'uy', 'rz'],
        nodes: [{ id: 1, x: 1, y: 2 }, { id: 2, x: 4, y: 6 }],
        elements: [{ id: 1, type: 113, nodeIds: [1, 2], properties: { area: 1, youngsModulus: 2, density: 3, momentOfInertia: 4 } }],
        dofMap: [[1, 2, 3], [4, 5, 6]], supports: [{ type: 1, nodeId: 1 }]
      },
      analysis: {
        type: 'static', displacements: [0, 0, 0, 0, 0, 0], loadVector: [0, 0, 0, 0, 1, 0],
        reactions: [0, -1, -4, 0, 0, 0], equilibriumResidual: [0, 0, 0],
        elementResults: [{ elementId: 1, type: 113, localEndForces: [0, -1, -4, 0, 1, 0], axialStrain: 0, axialStress: 0, axialForce: 0 }]
      }
    };
  }
  function trussFixture() {
    var data = fixture();
    data.model.dofPerNode = 2; data.model.dofLabels = ['ux', 'uy'];
    data.model.dofMap = [[1, 2], [3, 4]];
    data.model.elements[0] = { id: 1, type: 112, nodeIds: [1, 2], properties: { area: 1, youngsModulus: 2, density: 3 } };
    data.analysis.displacements = [0, 0, .01, 0]; data.analysis.loadVector = [0, 0, 10, 0];
    data.analysis.reactions = [-10, 0, 0, 0];
    data.analysis.elementResults[0] = { elementId: 1, type: 112, localEndForces: [-10, 0, 10, 0], axialStrain: .002, axialStress: 20, axialForce: 10 };
    return data;
  }

  test('schema accepts a valid frame', function () { assert(M.validateDataset(fixture()).raw.analysis.type === 'static'); });
  test('schema reports a JSON path', function () {
    var data = fixture(); data.model.elements[0].nodeIds[1] = 99;
    try { M.validateDataset(data); } catch (error) { assert(error.message.indexOf('$.model.elements[0].nodeIds[1]') === 0); return; }
    throw new Error('Invalid node reference was accepted');
  });
  test('112 interpolation preserves endpoints', function () {
    var points = M.geometry.sample112({ x: 0, y: 0 }, { x: 2, y: 0 }, { x: 1, y: 2 }, { x: 3, y: 4 }, 2);
    close(points[0].x, 2); close(points[0].y, 4); close(points[1].x, 8); close(points[1].y, 8);
  });
  test('113 Hermite interpolation preserves endpoints and rotations', function () {
    var L = 5, count = 1001, theta1 = 0.2, theta2 = -0.1;
    var points = M.geometry.sample113({ x: 0, y: 0 }, { x: L, y: 0 }, { x: 0, y: 0, r: theta1 }, { x: 0, y: 0, r: theta2 }, 1, count);
    close(points[0].x, 0); close(points[0].y, 0); close(points[count - 1].x, L); close(points[count - 1].y, 0);
    var dx = L / (count - 1);
    close((points[1].y - points[0].y) / dx, theta1, 0.001);
    close((points[count - 1].y - points[count - 2].y) / dx, theta2, 0.001);
  });
  test('inclined frame transforms axial displacement', function () {
    var first = { x: 1, y: 2 }, second = { x: 4, y: 6 };
    var points = M.geometry.sample113(first, second, { x: 0, y: 0, r: 0 }, { x: 0.003, y: 0.004, r: 0 }, 1, 5);
    close(points[4].x, 4.003); close(points[4].y, 6.004);
  });
  test('mathematical Y is inverted exactly once', function () { var p = M.geometry.svgPoint({ x: 2, y: 3 }); close(p.x, 2); close(p.y, -3); });
  test('automatic scale detects zero deformation', function () {
    var dataset = M.validateDataset(fixture()); var info = M.geometry.automaticScale(dataset, [0, 0, 0, 0, 0, 0], 31);
    assert(info.zero && info.scale === 1);
  });
  test('support restraint mapping matches the solver', function () {
    assert(M.geometry.restrainedLocalDOFs(1, 3).join(',') === '0,1,2');
    assert(M.geometry.restrainedLocalDOFs(2, 3).join(',') === '1,2');
    assert(M.geometry.restrainedLocalDOFs(3, 3).join(',') === '0,2');
    assert(M.geometry.restrainedLocalDOFs(4, 3).join(',') === '0,1');
  });
  test('renderer creates stable finite SVG layers', function () {
    var dataset = M.validateDataset(fixture());
    var renderer = new M.Renderer(document.getElementById('test-svg'), document.getElementById('test-selection'));
    renderer.render(dataset, dataset.raw.analysis.displacements, {
      showOriginal: true, showDeformed: true, showNodes: true,
      showNodeLabels: true, showElementLabels: true,
      scaleMode: 'auto', manualScale: 1, samples: 31
    });
    assert(document.querySelectorAll('#test-svg [data-layer]').length === 9, 'Layer count changed');
    document.querySelectorAll('#test-svg path').forEach(function (path) {
      assert(!/NaN|Infinity/.test(path.getAttribute('d')), 'Non-finite SVG path');
    });
  });
  test('cantilever diagram convention is explicit', function () {
    var q = [0, -100, -50, 0, 100, 0];
    close(M.staticResults.frameDiagram(q, 'N', 0.4), 0);
    close(M.staticResults.frameDiagram(q, 'V', 0.4), 100);
    close(M.staticResults.frameDiagram(q, 'M', 0), 50);
    close(M.staticResults.frameDiagram(q, 'M', 1), 0);
  });
  test('static renderer creates supports, loads, reactions, and diagrams', function () {
    var dataset = M.validateDataset(fixture()), svg = document.getElementById('test-svg');
    var renderer = new M.Renderer(svg, document.getElementById('test-selection'));
    var geometryInfo = renderer.render(dataset, dataset.raw.analysis.displacements, {
      showOriginal: true, showDeformed: true, showNodes: true,
      showNodeLabels: false, showElementLabels: false,
      scaleMode: 'auto', manualScale: 1, samples: 31
    });
    M.staticResults.render(renderer, dataset, {
      showSupports: true, showLoads: true, showReactions: true,
      trussResult: 'none', diagram: 'M'
    }, geometryInfo);
    assert(svg.querySelectorAll('.support-symbol').length === 1, 'Support was not rendered');
    assert(svg.querySelectorAll('[data-load-node-id]').length > 0, 'Load was not rendered');
    assert(svg.querySelectorAll('[data-reaction-node-id]').length > 0, 'Reaction was not rendered');
    assert(svg.querySelectorAll('[data-diagram-element-id]').length === 1, 'Diagram was not rendered');
    assert(renderer.layers.labels.style.display !== 'none', 'Static legend is hidden');
  });
  test('standalone SVG export embeds styles and escapes metadata', function () {
    var svg = document.getElementById('test-svg');
    var text = M.exporting.serialize(svg, { title: '<unsafe & title>', analysisType: 'static' });
    assert(text.indexOf('xmlns="http://www.w3.org/2000/svg"') >= 0, 'SVG namespace missing');
    assert(text.indexOf('<style') >= 0, 'Embedded styles missing');
    assert(text.indexOf('&lt;unsafe &amp; title&gt;') >= 0, 'Metadata was not safely escaped');
    assert(text.indexOf('http://') === text.indexOf('http://www.w3.org/2000/svg'), 'Unexpected external URL');
  });
  test('truss axial result uses signed coloring and labels', function () {
    var dataset = M.validateDataset(trussFixture()), svg = document.getElementById('test-svg');
    var renderer = new M.Renderer(svg, document.getElementById('test-selection'));
    var geometryInfo = renderer.render(dataset, dataset.raw.analysis.displacements, {
      showOriginal: true, showDeformed: true, showNodes: true,
      showNodeLabels: false, showElementLabels: false,
      scaleMode: 'auto', manualScale: 1, samples: 31
    });
    M.staticResults.render(renderer, dataset, {
      showSupports: true, showLoads: true, showReactions: true,
      trussResult: 'axialForce', diagram: 'none'
    }, geometryInfo);
    assert(svg.querySelector('.deformed-element').style.stroke !== '', 'Tension color was not applied');
    assert(Array.from(svg.querySelectorAll('.result-label')).some(function (label) { return label.textContent.indexOf('10') >= 0; }), 'Axial result label missing');
  });
  test('PNG rasterization produces a nonempty blob', function () {
    return M.exporting.pngBlob(document.getElementById('test-svg'), { title: 'PNG test', analysisType: 'static' }, 1).then(function (blob) {
      assert(blob.type === 'image/png' && blob.size > 0, 'PNG blob is empty');
    });
  });

  document.addEventListener('DOMContentLoaded', async function () {
    var lines = [], failed = 0;
    for (var i = 0; i < tests.length; i += 1) {
      var entry = tests[i];
      try { await entry.operation(); lines.push('PASS ' + entry.name); }
      catch (error) { failed += 1; lines.push('FAIL ' + entry.name + ': ' + error.message); }
    }
    document.getElementById('test-output').textContent = lines.join('\n') + '\n\n' + (failed ? failed + ' failed' : tests.length + ' passed');
    document.body.setAttribute('data-test-status', failed ? 'failed' : 'passed');
  });
}(window.MKEFPost));
