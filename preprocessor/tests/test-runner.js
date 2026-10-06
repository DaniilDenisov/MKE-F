(function () {
  'use strict';
  var tests = [], passed = 0;
  function test(name, operation) { tests.push({ name: name, operation: operation }); }
  function assert(condition, message) { if (!condition) throw new Error(message || 'Assertion failed.'); }
  function throws(fragment, operation) { try { operation(); } catch (error) { assert(error.message.indexOf(fragment) >= 0, 'Unexpected error: ' + error.message); return; } throw new Error('Expected an error containing ' + fragment); }
  function fixture(record) { return ['# browser fixture', 'analysis', record || 'static', 'nodes', '2', '0,0,0', '1,0,0', 'elems_113', '1', '113,1,2,0.1,200,10,0.01', 'bcfix', '1', '1,1,0,0,0', 'bcforce_stat', '1', '10,2,5,-2,1', ''].join('\n'); }

  test('complete support catalog and family preservation', function () {
    var S = MKEFSupports;
    for (var type=1; type<=7; type++) assert(S.type(S.fromType(type,1,3),3) === type);
    assert(S.type(S.fromType(4,1,2),2) === 1);
    assert(S.type(S.fromType(5,1,2),2) === 3);
    throws('Unsupported', function () { S.fromType(7,1,2); });
    var model = {elementType:112,supports:[S.fromType(1,1,2)]};
    S.changeFamily(model,113); assert(S.type(model.supports[0],3) === 4);
    model.supports.push(S.fromType(7,2,3));
    throws('thetaZ',function () { S.changeFamily(model,112); });
    assert(model.elementType === 113);
  });
  test('MPC parser round trip and validation parity', function () {
    var text=fixture()+'mpc\n1\n2,2,0,1,1,2,1\n', model=MKEFPre.caseFormat.parse(text);
    assert(model.mpcs.length===1 && model.mpcs[0].masters[0].coefficient===1);
    assert(MKEFPre.caseFormat.serialize(MKEFPre.caseFormat.parse(MKEFPre.caseFormat.serialize(model)))===MKEFPre.caseFormat.serialize(model));
    ['2,2,1,1,1,2,1','2,2,0,2,1,2,1','2,2,0,1,1,2,NaN','2,2,0,1,1,2,0','2,2,0,1,2,2,1','2,4,0,1,1,2,1','2,2,0,2,1,2,1,1,2,2','1,2,0,1,2,2,1'].forEach(function (record) {
      throws('MPC',function () { MKEFPre.caseFormat.parse(fixture()+'mpc\n1\n'+record); });
    });
    throws('Cycle',function () { MKEFPre.caseFormat.parse(fixture()+'mpc\n2\n2,2,0,1,2,1,1\n2,1,0,1,2,2,1'); });
    throws('more than one',function () { MKEFPre.caseFormat.parse(text+'mpc\n1\n2,2,0,1,1,2,1'); });
  });
  test('axis helper respects rigid motion, orientation and dependent availability', function () {
    [[2,0],[0,2],[2,2],[3,4]].forEach(function (end) {
      var model={elementType:112,nodes:[{x:0,y:0},{x:end[0],y:end[1]},{x:end[0]*.25,y:end[1]*.25}],supports:[],mpcs:[]};
      var m=MKEFMPC.onAxis(model,1,2,3);
      function motion(node,dof) { var p=model.nodes[node-1]; return dof===1?2-.03*p.y:-4+.03*p.x; }
      var right=m.masters.reduce(function (v,a) { return v+a.coefficient*motion(a.node,a.dof); },0);
      assert(Math.abs(motion(m.depNode,m.depDOF)-right)<1e-12);
      if (end[0]===end[1]) { assert(m.depDOF===1); model.supports=[MKEFSupports.fromType(3,3,2)]; assert(MKEFMPC.onAxis(model,1,2,3).depDOF===2); }
      model.nodes[2].x+=100; throws('inside',function () { MKEFMPC.onAxis(model,1,2,3); });
    });
  });
  test('MPC node deletion, renumbering and geometry edits', function () {
    var model=MKEFPre.caseFormat.parse(fixture()); model.nodes.push({x:2,y:0});
    model.mpcs=[{depNode:3,depDOF:2,rhs:0,masters:[{node:2,dof:2,coefficient:.5}]}];
    var saved=JSON.stringify(model.mpcs); MKEFPre.modelEdit.updateNode(model,2,'x',3); assert(JSON.stringify(model.mpcs)===saved);
    MKEFPre.modelEdit.deleteNode(model,0); assert(model.mpcs[0].depNode===2 && model.mpcs[0].masters[0].node===1);
    MKEFPre.modelEdit.deleteNode(model,0); assert(model.mpcs.length===0);
  });
  test('parses and serializes a configured static case', function () { var model = MKEFPre.caseFormat.parse(fixture()), text = MKEFPre.caseFormat.serialize(model); assert(model.analysis.type === 'static'); assert(text.indexOf('analysis\nstatic\n') === 0); assert(text.indexOf('bcforce_stat\n1\n10,2,5,-2,1') > 0); });
  test('canonical round trip preserves semantic data', function () { var first = MKEFPre.caseFormat.parse(fixture()), second = MKEFPre.caseFormat.parse(MKEFPre.caseFormat.serialize(first)); assert(JSON.stringify(first) === JSON.stringify(second)); });
  test('legacy cases require a task only at export', function () { var model = MKEFPre.caseFormat.parse(fixture().replace('analysis\nstatic\n', '')); assert(model.analysis.type === ''); throws('Select static', function () { MKEFPre.caseFormat.serialize(model); }); });
  test('transient settings round trip exactly', function () { var text = fixture('transient,0.001,0.01,2,2').replace('bcforce_stat\n1\n10,2,5,-2,1\n', 'bcforce_harm\n1\n11,2,5,-2,1,3\n'), model = MKEFPre.caseFormat.parse(text), output = MKEFPre.caseFormat.serialize(model); assert(model.analysis.timeStep === .001 && model.analysis.monitorDOF === 2); assert(output.indexOf('transient,0.001,0.01,2,2') > 0); });
  test('modal cases reject loads', function () { throws('not valid for modal', function () { MKEFPre.caseFormat.parse(fixture('modal')); }); });
  test('static cases reject time-dependent loads', function () { throws('not valid for static', function () { MKEFPre.caseFormat.parse(fixture().replace('bcforce_stat\n1\n10,', 'bcforce_step\n1\n13,')); }); });
  test('unknown markers report their source line', function () { throws('Line 1', function () { MKEFPre.caseFormat.parse('mystery\n'); }); });
  test('duplicate constraints are detected by DOF', function () { var model = MKEFPre.caseFormat.parse(fixture()); model.supports.push({ type: 2, node: 1 }); assert(MKEFPre.caseFormat.validate(model).some(function (message) { return message.indexOf('duplicates') >= 0; })); });
  test('truss cases reject nonzero nodal moments', function () { var model = MKEFPre.caseFormat.newModel('static', 112); model.nodes = [{ x: 0, y: 0 }, { x: 1, y: 0 }]; model.elements = [{ node1: 1, node2: 2, area: 1, youngsModulus: 1, density: 1, momentOfInertia: null }]; model.loads = [{ type: 10, node: 2, fx: 0, fy: 0, mz: 1, frequency: null }]; assert(MKEFPre.caseFormat.validate(model).some(function (message) { return message.indexOf('truss node') >= 0; })); });
  test('zero-length elements are rejected', function () { var model = MKEFPre.caseFormat.parse(fixture()); model.nodes[1] = { x: 0, y: 0 }; assert(MKEFPre.caseFormat.validate(model).some(function (message) { return message.indexOf('zero length') >= 0; })); });
  test('canonical export groups repeated load types', function () { var model = MKEFPre.caseFormat.parse(fixture()); model.loads.push({ type: 10, node: 2, fx: 1, fy: 0, mz: 0, frequency: null }); var text = MKEFPre.caseFormat.serialize(model); assert(text.indexOf('bcforce_stat\n2\n') > 0); });
  test('exact node coordinate edits are serialized', function () { var model = MKEFPre.caseFormat.parse(fixture()); model.nodes[1].x = 2.25; model.nodes[1].y = -3.5; var text = MKEFPre.caseFormat.serialize(model); assert(text.indexOf('2.25,-3.5,0') > 0); });
  test('renderer exposes model-coordinate grid and axes', function () { var svg = document.getElementById('test-svg'), renderer = new MKEFPre.Renderer(svg); renderer.setGrid(true, .5); renderer.draw(MKEFPre.caseFormat.parse(fixture())); assert(svg.querySelectorAll('.coordinate-grid line').length > 0, 'Grid lines are missing.'); var labels = Array.prototype.map.call(svg.querySelectorAll('.coordinate-axes text'), function (item) { return item.textContent; }).join(' '); assert(labels.indexOf('X') >= 0 && labels.indexOf('Y') >= 0, 'Axis labels are missing.'); });
  test('renderer viewBox follows the full canvas aspect ratio', function () { var svg = document.getElementById('test-svg'), renderer = new MKEFPre.Renderer(svg), rect = svg.getBoundingClientRect(); renderer.view = { x: -1, y: -1, width: 2, height: 2 }; renderer.synchronizeAspectRatio(); assert(Math.abs(renderer.view.width / renderer.view.height - rect.width / rect.height) < 1e-6, 'Canvas and viewBox aspect ratios differ.'); });
  test('node tool shows a dashed placement preview at the normalized point', function () { var svg = document.getElementById('test-svg'), renderer = new MKEFPre.Renderer(svg, { placementPoint: function () { return { x: 2, y: 3 }; } }); renderer.draw(MKEFPre.caseFormat.parse(fixture())); renderer.setTool('node'); renderer.updateNodePreview({ x: 1.8, y: 3.2 }); var preview = svg.querySelector('.node-placement-preview'); assert(preview && preview.getAttribute('visibility') === 'visible', 'Node preview is not visible.'); assert(preview.getAttribute('cx') === '2' && preview.getAttribute('cy') === '-3', 'Node preview does not use normalized coordinates.'); });
  test('member tool can start with the selected node as its first endpoint', function () { var svg = document.getElementById('test-svg'), endpoints = null, renderer = new MKEFPre.Renderer(svg, { addElement: function (first, second) { endpoints = [first, second]; } }); renderer.draw(MKEFPre.caseFormat.parse(fixture())); renderer.setTool('member', 1); assert(svg.querySelector('[data-node-id="2"]').classList.contains('pending'), 'Selected first endpoint is not pending.'); svg.querySelector('[data-node-id="1"]').dispatchEvent(new MouseEvent('click', { bubbles: true })); assert(endpoints && endpoints[0] === 1 && endpoints[1] === 0, 'Member endpoints do not start with the selected node.'); });
  test('member tool retains the first endpoint when creation is rejected', function () { var svg = document.getElementById('test-svg'), renderer = new MKEFPre.Renderer(svg, { addElement: function () { return false; } }); renderer.draw(MKEFPre.caseFormat.parse(fixture())); renderer.setTool('member', 0); svg.querySelector('[data-node-id="2"]').dispatchEvent(new MouseEvent('click', { bubbles: true })); assert(svg.querySelector('[data-node-id="1"]').classList.contains('pending'), 'Rejected creation discarded the first endpoint.'); });
  test('support types use distinct marker geometry', function () { var svg = document.getElementById('test-svg'), model = MKEFPre.caseFormat.parse(fixture()), renderer = new MKEFPre.Renderer(svg); model.nodes = [{ x: 0, y: 0 }, { x: 1, y: 0 }, { x: 2, y: 0 }, { x: 3, y: 0 }]; model.supports = [{ type: 1, node: 1 }, { type: 2, node: 2 }, { type: 3, node: 3 }, { type: 4, node: 4 }]; renderer.draw(model); assert(svg.querySelector('.support-type-1 .support-ground') && !svg.querySelector('.support-type-1 path'), 'Type 1 fixed marker is missing.'); assert(svg.querySelectorAll('.support-type-2 circle').length === 1 && svg.querySelector('.support-type-2 path'), 'Type 2 horizontal roller marker is missing.'); assert(svg.querySelectorAll('.support-type-3 circle').length === 1 && svg.querySelector('.support-type-3 path'), 'Type 3 vertical roller marker is missing.'); assert(svg.querySelector('.support-type-4 path') && !svg.querySelector('.support-type-4 circle'), 'Type 4 pinned marker is missing.'); });
  test('loads render above every model layer', function () { var svg = document.getElementById('test-svg'), renderer = new MKEFPre.Renderer(svg); renderer.draw(MKEFPre.caseFormat.parse(fixture())); assert(svg.lastElementChild.classList.contains('loads-layer'), 'Loads layer is not topmost.'); });

  function uniformFixture() { return MKEFPre.caseFormat.parse(fixture() + 'eload_uniform\n1\n20,1,2,0,-7\n'); }
  test('uniform load sections round trip and accumulate', function () {
    var model = MKEFPre.caseFormat.parse(fixture() + 'eload_uniform\n1\n20,1,1,2,-7\neload_uniform\n1\n20,1,2,0,3\n');
    assert(model.elementLoads.length === 2);
    assert(JSON.stringify(model) === JSON.stringify(MKEFPre.caseFormat.parse(MKEFPre.caseFormat.serialize(model))));
  });
  test('uniform parser rejects missing, nonfinite, zero, invalid and unsupported records', function () {
    ['20,1,2,,7', '20,1,2,0,Infinity', '20,1,2,0,0', '20,2,2,0,7', '20,1,3,0,7', '21,1,2,0,7', '20,1,2,7'].forEach(function (record) {
      throws('Line', function () { MKEFPre.caseFormat.parse(fixture() + 'eload_uniform\n1\n' + record); });
    });
    throws('only for static', function () { MKEFPre.caseFormat.parse(fixture('modal').replace('bcforce_stat\n1\n10,2,5,-2,1\n', '') + 'eload_uniform\n1\n20,1,1,0,1'); });
    var model = uniformFixture(); model.elementType = 112; assert(MKEFPre.caseFormat.validate(model).some(function (s) { return s.indexOf('frame element 113') >= 0; }));
  });
  test('coordinate switch preserves physical force on inclined member', function () {
    var model = uniformFixture(); model.nodes[1] = {x:3,y:4};
    MKEFPre.modelEdit.changeCoordinates(model, 0, 1);
    assert(Math.abs(model.elementLoads[0].qx + 5.6) < 1e-12);
    assert(Math.abs(model.elementLoads[0].qy + 4.2) < 1e-12);
    MKEFPre.modelEdit.changeCoordinates(model, 0, 2);
    assert(Math.abs(model.elementLoads[0].qx) < 1e-12 && Math.abs(model.elementLoads[0].qy + 7) < 1e-12);
    model.nodes[1] = model.nodes[0]; throws('Valid element geometry', function () { MKEFPre.modelEdit.changeCoordinates(model, 0, 1); });
  });
  test('all actual element edits remove loads, unchanged values do not', function () {
    ['node1','node2','area','youngsModulus','density','momentOfInertia'].forEach(function (field) {
      var model = uniformFixture(), old = model.elements[0][field];
      MKEFPre.modelEdit.updateElement(model, 0, field, old); assert(model.elementLoads.length === 1);
      MKEFPre.modelEdit.updateElement(model, 0, field, old+1); assert(model.elementLoads.length === 0);
    });
  });
  test('moving a node preserves nodal conditions and removes incident element loads', function () {
    var model = uniformFixture(); MKEFPre.modelEdit.updateNode(model, 1, 'x', 1); assert(model.elementLoads.length === 1);
    MKEFPre.modelEdit.updateNode(model, 1, 'x', 2); assert(model.elementLoads.length === 0 && model.loads.length === 1 && model.supports.length === 1);
  });
  test('node deletion cascades and renumbers only surviving references', function () {
    var model = uniformFixture(); model.nodes.push({x:2,y:0}, {x:3,y:0});
    model.elements.push(Object.assign({}, model.elements[0], {node1:3,node2:4}));
    model.elementLoads.push({type:20,elementId:2,coordinateSystem:1,qx:2,qy:0});
    model.supports.push({type:4,node:3}); model.analysis.monitorNode=2;
    MKEFPre.modelEdit.deleteNode(model, 1);
    assert(model.nodes.length===3 && model.elements.length===1 && model.loads.length===0);
    assert(model.elements[0].node1===2 && model.elements[0].node2===3);
    assert(model.elementLoads.length===1 && model.elementLoads[0].elementId===1 && model.elementLoads[0].qx===2);
    assert(model.supports[0].node===1 && model.supports[1].node===2 && model.analysis.monitorNode===0);
    MKEFPre.modelEdit.deleteNode(model, 0); assert(model.supports.length===1 && model.supports[0].node===1);
  });
  test('deleting a member preserves nodes and nodal conditions', function () {
    var model=uniformFixture(); MKEFPre.modelEdit.deleteElement(model,0);
    assert(model.elements.length===0 && model.elementLoads.length===0 && model.nodes.length===2 && model.loads.length===1 && model.supports.length===1);
  });
  test('uniform arrows and local axes are visible for selected element', function () {
    var renderer=new MKEFPre.Renderer(document.getElementById('test-svg')); renderer.draw(uniformFixture()); renderer.setSelection({kind:'element',index:0});
    assert(document.querySelectorAll('#test-svg .element-load-symbol line').length===9);
    assert(document.querySelectorAll('#test-svg .local-axis').length===2);
  });

  function linearFixture() { return MKEFPre.caseFormat.parse(fixture()+'eload_linear\n1\n21,1,2,0,0,0,-8\n'); }
  test('mixed linear and uniform loads round trip in order', function () {
    var F=MKEFPre.caseFormat, model=F.parse(F.serialize(linearFixture())+'eload_uniform\n1\n20,1,1,2,0\neload_linear\n1\n21,1,1,3,-4,-3,4\n');
    assert(JSON.stringify(F.parse(F.serialize(model)).elementLoads)===JSON.stringify(model.elementLoads));
  });
  test('linear parser and validation reject invalid records and tasks', function () {
    ['21,1,2,0,0,0,0','20,1,2,0,1,0,2','21,1,2,0,NaN,0,1','21,1,2,0,Inf,0,1','21,1,2,0,,0,1','21,3,2,0,1,0,2','21,1,3,0,1,0,2','21,1,2,0,1,0'].forEach(function (record) {
      throws('Line',function () { MKEFPre.caseFormat.parse(fixture()+'eload_linear\n1\n'+record); });
    });
    var model=linearFixture(); model.analysis.type='modal'; assert(MKEFPre.caseFormat.validate(model).some(function (s) { return s.includes('only for static'); }));
    model.analysis.type='transient'; assert(MKEFPre.caseFormat.validate(model).some(function (s) { return s.includes('only for static'); }));
    model.analysis.type='static'; model.elementType=112; assert(MKEFPre.caseFormat.validate(model).some(function (s) { return s.includes('frame element 113'); }));
  });
  test('load type conversion and both endpoint coordinate transforms preserve intensity', function () {
    var E=MKEFPre.modelEdit, model=uniformFixture(); E.changeLoadType(model,0,21);
    assert(model.elementLoads[0].qy1===-7 && model.elementLoads[0].qy2===-7);
    model.elementLoads[0].qy1=0; model.nodes[1]={x:3,y:4}; E.changeCoordinates(model,0,1);
    var q=model.elementLoads[0]; assert(q.qx1===0 && q.qy1===0 && Math.abs(q.qx2+5.6)<1e-12 && Math.abs(q.qy2+4.2)<1e-12);
    E.changeCoordinates(model,0,2); E.changeLoadType(model,0,20);
    assert(Math.abs(model.elementLoads[0].qx)<1e-12 && Math.abs(model.elementLoads[0].qy+3.5)<1e-12);
  });
  test('linear arrows grow to the endpoint and include an envelope', function () {
    var renderer=new MKEFPre.Renderer(document.getElementById('test-svg')); renderer.draw(linearFixture());
    var lines=Array.from(document.querySelectorAll('#test-svg .element-load-symbol line'));
    assert(lines.length===8); var lengths=lines.map(function (line) { return Math.hypot(line.x2.baseVal.value-line.x1.baseVal.value,line.y2.baseVal.value-line.y1.baseVal.value); });
    assert(lengths[7]>lengths[0]*7.9); assert(document.querySelector('#test-svg .element-load-envelope'));
    assert(document.querySelector('#test-svg .element-load-label').textContent.includes('q1=(0, 0)'));
  });
  test('linear loads follow existing deletion and geometry invalidation rules', function () {
    var model=linearFixture(); MKEFPre.modelEdit.setRelease(model,0,1,true); assert(model.elementLoads.length===1);
    MKEFPre.modelEdit.updateNode(model,1,'x',2); assert(model.elementLoads.length===0);
    model=linearFixture(); MKEFPre.modelEdit.deleteElement(model,0); assert(model.elementLoads.length===0);
  });

  document.addEventListener('DOMContentLoaded', function () {
    var output = document.getElementById('test-output'), status = document.getElementById('test-status'), lines = [];
    tests.forEach(function (item) { try { item.operation(); passed += 1; lines.push('PASS  ' + item.name); } catch (error) { lines.push('FAIL  ' + item.name + '\n      ' + error.message); } });
    output.textContent = lines.join('\n'); status.textContent = passed + ' / ' + tests.length + ' tests passed'; status.setAttribute('data-status', passed === tests.length ? 'passed' : 'failed');
  });
}());
