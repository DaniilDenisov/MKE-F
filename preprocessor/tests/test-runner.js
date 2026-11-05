(function () {
  'use strict';
  var tests = [], passed = 0;
  function test(name, operation) { tests.push({ name: name, operation: operation }); }
  function assert(condition, message) { if (!condition) throw new Error(message || 'Assertion failed.'); }
  function throws(fragment, operation) { try { operation(); } catch (error) { assert(error.message.indexOf(fragment) >= 0, 'Unexpected error: ' + error.message); return; } throw new Error('Expected an error containing ' + fragment); }
  function fixture(record) { return ['# browser fixture', 'analysis', record || 'static', 'nodes', '2', '0,0,0', '1,0,0', 'elems_113', '1', '113,1,2,0.1,200,10,0.01', 'bcfix', '1', '1,1,0,0,0', 'bcforce_stat', '1', '10,2,5,-2,1', ''].join('\n'); }
  function editor(model) { return new MKEFPre.spatialEdit.History(MKEFPre.editorAdapter.fromNative(model)); }

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
    var history=editor(model); model=history.model;
    var saved=JSON.stringify(model.mpcs); history.moveNode(2,{x:3,y:0,z:0}); assert(JSON.stringify(model.mpcs)===saved);
    history.remove({kind:'node',index:0}); assert(model.mpcs[0].depNode===2 && model.mpcs[0].masters[0].node===1);
    history.remove({kind:'node',index:0}); assert(model.mpcs.length===0);
    history.travel(true); assert(model.mpcs.length===1);
  });
  test('parses and serializes a configured static case', function () { var model = MKEFPre.caseFormat.parse(fixture()), text = MKEFPre.caseFormat.serialize(model); assert(model.analysis.type === 'static'); assert(text.indexOf('analysis\nstatic\n') === 0); assert(text.indexOf('bcforce_stat\n1\n10,2,5,-2,1') > 0); });
  test('canonical round trip preserves semantic data', function () { var first = MKEFPre.caseFormat.parse(fixture()), second = MKEFPre.caseFormat.parse(MKEFPre.caseFormat.serialize(first)); assert(JSON.stringify(first) === JSON.stringify(second)); });
  test('TXT requires an explicit analysis', function () { throws('Select static', function () { MKEFPre.caseFormat.parse(fixture().replace('analysis\nstatic\n', '')); }); });
  test('transient settings round trip exactly', function () { var text = fixture('transient,0.001,0.01,2,2').replace('bcforce_stat\n1\n10,2,5,-2,1\n', 'bcforce_harm\n1\n11,2,5,-2,1,3\n'), model = MKEFPre.caseFormat.parse(text), output = MKEFPre.caseFormat.serialize(model); assert(model.analysis.timeStep === .001 && model.analysis.monitorDOF === 2); assert(output.indexOf('transient,0.001,0.01,2,2') > 0); });
  test('modal cases reject loads', function () { throws('not valid for modal', function () { MKEFPre.caseFormat.parse(fixture('modal')); }); });
  test('static cases reject time-dependent loads', function () { throws('not valid for static', function () { MKEFPre.caseFormat.parse(fixture().replace('bcforce_stat\n1\n10,', 'bcforce_step\n1\n13,')); }); });
  test('unknown markers report their source line', function () { throws('Line 1', function () { MKEFPre.caseFormat.parse('mystery\n'); }); });
  test('duplicate constraints are detected by DOF', function () { var model = MKEFPre.caseFormat.parse(fixture()); model.supports.push({ type: 2, node: 1 }); assert(MKEFPre.caseFormat.validate(model).some(function (message) { return message.indexOf('duplicates') >= 0; })); });
  test('truss cases reject nonzero nodal moments', function () { var model = MKEFPre.caseFormat.newModel('static', 112); model.nodes = [{ x: 0, y: 0 }, { x: 1, y: 0 }]; model.elements = [{ node1: 1, node2: 2, area: 1, youngsModulus: 1, density: 1, momentOfInertia: null }]; model.loads = [{ type: 10, node: 2, fx: 0, fy: 0, mz: 1, frequency: null }]; assert(MKEFPre.caseFormat.validate(model).some(function (message) { return message.indexOf('truss node') >= 0; })); });
  test('zero-length elements are rejected', function () { var model = MKEFPre.caseFormat.parse(fixture()); model.nodes[1] = { x: 0, y: 0 }; assert(MKEFPre.caseFormat.validate(model).some(function (message) { return message.indexOf('zero length') >= 0; })); });
  test('canonical export groups repeated load types', function () { var model = MKEFPre.caseFormat.parse(fixture()); model.loads.push({ type: 10, node: 2, fx: 1, fy: 0, mz: 0, frequency: null }); var text = MKEFPre.caseFormat.serialize(model); assert(text.indexOf('bcforce_stat\n2\n') > 0); });
  test('exact node coordinate edits are serialized', function () { var model = MKEFPre.caseFormat.parse(fixture()); model.nodes[1].x = 2.25; model.nodes[1].y = -3.5; var text = MKEFPre.caseFormat.serialize(model); assert(text.indexOf('2.25,-3.5,0') > 0); });










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
    var history=editor(model), C=MKEFPre.editorConditions; model=history.model;
    history.condition('elementLoads',0,C.convertLoadAxes(model,model.elementLoads[0],1));
    assert(Math.abs(model.elementLoads[0].components[0] + 5.6) < 1e-12);
    assert(Math.abs(model.elementLoads[0].components[1] + 4.2) < 1e-12);
    history.condition('elementLoads',0,C.convertLoadAxes(model,model.elementLoads[0],2));
    assert(Math.abs(model.elementLoads[0].components[0]) < 1e-12 && Math.abs(model.elementLoads[0].components[1] + 7) < 1e-12);
    var saved=JSON.stringify(model); throws('Another node', function () { history.moveNode(1,model.nodes[0]); });
    assert(JSON.stringify(model)===saved);
  });
  test('all actual element edits remove loads, unchanged values do not', function () {
    ['node1','node2','area','youngsModulus','density','momentOfInertia'].forEach(function (field) {
      var history=editor(uniformFixture()),model=history.model;
      history.addNode({x:2,y:1,z:0});
      var element=Object.assign({},model.elements[0]);
      history.member(0,element); assert(model.elementLoads.length === 1);
      element[field]=field==='node1'||field==='node2'?3:element[field]+1;
      history.member(0,element); assert(model.elementLoads.length === 0);
      history.travel(true); assert(model.elementLoads.length === 1);
    });
  });
  test('moving a node preserves nodal conditions and removes incident element loads', function () {
    var history=editor(uniformFixture()),model=history.model;
    history.moveNode(1,{x:1,y:0,z:0}); assert(model.elementLoads.length === 1);
    history.moveNode(1,{x:2,y:0,z:0}); assert(model.elementLoads.length === 0 && model.loads.length === 1 && model.supports.length === 1);
    history.travel(true); assert(model.elementLoads.length===1);
  });
  test('node deletion cascades and renumbers only surviving references', function () {
    var model = uniformFixture(); model.nodes.push({x:2,y:0}, {x:3,y:0});
    model.elements.push(Object.assign({}, model.elements[0], {node1:3,node2:4}));
    model.elementLoads.push({type:20,elementId:2,coordinateSystem:1,qx:2,qy:0});
    model.supports.push({type:4,node:3});
    var history=editor(model); model=history.model;
    history.remove({kind:'node',index:1});
    assert(model.nodes.length===3 && model.elements.length===1 && model.loads.length===0);
    assert(model.elements[0].node1===2 && model.elements[0].node2===3);
    assert(model.elementLoads.length===1 && model.elementLoads[0].elementId===1 && model.elementLoads[0].components[0]===2);
    assert(model.supports[0].node===1 && model.supports[1].node===2);
    history.remove({kind:'node',index:0}); assert(model.supports.length===1 && model.supports[0].node===1);
    history.travel(true); history.travel(true); assert(model.nodes.length===4 && model.elementLoads.length===2);
  });
  test('deleting a member preserves nodes and nodal conditions', function () {
    var history=editor(uniformFixture()),model=history.model; history.remove({kind:'element',index:0});
    assert(model.elements.length===0 && model.elementLoads.length===0 && model.nodes.length===2 && model.loads.length===1 && model.supports.length===1);
  });
  test('deleting a transient monitor requires a replacement and restores it on undo', function () {
    var native=MKEFPre.caseFormat.parse(fixture('transient,0.001,0.01,2,2').replace('bcforce_stat\n1\n10,2,5,-2,1\n','bcforce_harm\n1\n11,2,5,-2,1,3\n'));
    native.nodes.push({x:2,y:0});
    var history=editor(native),model=history.model,saved=JSON.stringify(model);
    throws('remaining monitor',function(){history.remove({kind:'node',index:1});});
    assert(JSON.stringify(model)===saved);
    history.remove({kind:'node',index:1},3); assert(model.analysis.monitorNode===2);
    history.travel(true); assert(JSON.stringify(model)===saved);
    history.travel(false); assert(model.nodes.length===2 && model.analysis.monitorNode===2);
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
    var native=uniformFixture(); native.nodes[1]={x:3,y:4};
    var history=editor(native),model=history.model,C=MKEFPre.editorConditions;
    history.condition('elementLoads',0,C.convertLoadType(model.elementLoads[0],21));
    assert(model.elementLoads[0].components[1]===-7 && model.elementLoads[0].components[4]===-7);
    var load=Object.assign({},model.elementLoads[0],{components:[0,0,0,0,-7,0]});
    history.condition('elementLoads',0,C.convertLoadAxes(model,load,1));
    var q=model.elementLoads[0].components; assert(q[0]===0 && q[1]===0 && Math.abs(q[3]+5.6)<1e-12 && Math.abs(q[4]+4.2)<1e-12);
    history.condition('elementLoads',0,C.convertLoadAxes(model,model.elementLoads[0],2));
    throws('Confirm',function(){C.convertLoadType(model.elementLoads[0],20,false);});
    history.condition('elementLoads',0,C.convertLoadType(model.elementLoads[0],20,true));
    assert(Math.abs(model.elementLoads[0].components[0])<1e-12 && Math.abs(model.elementLoads[0].components[1]+3.5)<1e-12);
  });

  test('linear loads follow existing deletion and geometry invalidation rules', function () {
    var history=editor(linearFixture()),model=history.model;
    history.setReleaseMask(1,1,['Mz']); assert(model.elementLoads.length===1);
    history.moveNode(1,{x:2,y:0,z:0}); assert(model.elementLoads.length===0);
    history.travel(true); assert(model.elementLoads.length===1 && model.releases.length===1);
    history.remove({kind:'element',index:0}); assert(model.elementLoads.length===0 && model.releases.length===0);
  });

  document.addEventListener('DOMContentLoaded', function () {
    var output = document.getElementById('test-output'), status = document.getElementById('test-status'), lines = [];
    tests.forEach(function (item) { try { item.operation(); passed += 1; lines.push('PASS  ' + item.name); } catch (error) { lines.push('FAIL  ' + item.name + '\n      ' + error.message); } });
    output.textContent = lines.join('\n'); status.textContent = passed + ' / ' + tests.length + ' tests passed'; status.setAttribute('data-status', passed === tests.length ? 'passed' : 'failed');
  });
}());
