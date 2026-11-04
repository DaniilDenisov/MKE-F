(function (M) {
  'use strict';
  var tests = [];
  function test(name, operation) { tests.push({ name: name, operation: operation }); }
  function assert(condition, message) { if (!condition) throw new Error(message || 'Assertion failed'); }
  function close(actual, expected, tolerance) { assert(Math.abs(actual - expected) <= (tolerance || 1e-10), actual + ' != ' + expected); }
  // Synthetic v8 display fixtures. Identity operators are deliberate UI data;
  // numerical K/M accuracy is verified against Octave by the recovery suites.
  function completeFixture(data){
    var m=data.model,d=m.dofPerNode,n=2*d,rows=new Map(m.nodes.map((node,i)=>[node.id,i]));
    data.metadata.exporter={name:'mkef-octave',contractVersion:1};
    m.numberOfDOFs=m.nodes.length*d;
    m.recovery={version:1,coordinates:'local-physical-end',stiffnessFormulation:d===2?'linear-bar':'euler-bernoulli',massFormulation:'consistent-v1'};
    m.elements.forEach(function(e,i){
      var a=m.nodes[rows.get(e.nodeIds[0])],b=m.nodes[rows.get(e.nodeIds[1])],L=Math.hypot(b.x-a.x,b.y-a.y),c=(b.x-a.x)/L,s=(b.y-a.y)/L;
      e.length=L;e.globalDOFs=e.nodeIds.flatMap(id=>m.dofMap[rows.get(id)]);
      var block=d===2?[[c,s],[-s,c]]:[[c,s,0],[-s,c,0],[0,0,1]];
      e.localFromSystem=Array.from({length:n},(_,r)=>Array.from({length:n},(_,k)=>Math.floor(r/d)===Math.floor(k/d)?block[r%d][k%d]:0));
      if(data.analysis.type==='static'){
        var r=data.analysis.elementResults[i],f=r.localEndForces;
        r.localDisplacements=Array(n).fill(0);r.equivalentLocalLoadVector=r.equivalentLocalLoadVector||Array(n).fill(0);
        r.sectionForceCoefficients=[[-f[0],0,0,0],[0,0,0,0],[0,0,0,0]];
        if(d===2){var point={value:-f[0]/e.properties.area,xi:0};r.sectionStress={minimum:point,maximum:point,positions:[0,1]};}
        if(d===3){r.sectionForceCoefficients[1]=[-f[1],0,0,0];r.sectionForceCoefficients[2]=[-f[2],f[2]+f[5],0,0];}
      }else{
        var identity=Array.from({length:n},(_,r)=>Array.from({length:n},(_,k)=>r===k?1:0));
        e.recoveryOperators={K:identity,M:identity};
      }
    });
    if(data.analysis.type==='modal')data.analysis.eigenvalues=data.analysis.angularFrequenciesRadPerSec.map(v=>v*v);
    return data;
  }
  function fixture() {
    return completeFixture({
      format: 'mkef-postprocessor', version: 8,
      metadata: { title: 'Browser test', units: { length: 'm', force: 'N', moment: 'N*m', stress: 'Pa', time: 's' } },
      model: { dimension: 2, dofPerNode: 3, dofLabels: ['ux', 'uy', 'rz'],
        nodes: [{ id: 1, x: 0, y: 0 }, { id: 2, x: .5, y: .5 }],
        elements: [{ id: 1, type: 113, nodeIds: [1, 2], properties: { area: 1, youngsModulus: 2, density: 3, momentOfInertia: 4 } }],
        dofMap: [[1, 2, 3], [4, 5, 6]], supports: [{ type: 1, nodeId: 1 }] },
      analysis: { type: 'static', displacements: [0, 0, 0, 0, 0, 0], loadVector: [0, 0, 0, 0, 1, 0], reactions: [0, -1, -4, 0, 0, 0], equilibriumResidual: [0, 0, 0],
        elementResults: [{ elementId: 1, type: 113, localEndForces: [0, -1, -4, 0, 1, 0], axialStrain: 0, axialStress: 0, axialForce: 0 }] }
    });
  }
  function trussFixture() {
    var data = fixture(); data.model.dofPerNode = 2; data.model.dofLabels = ['ux', 'uy']; data.model.dofMap = [[1, 2], [3, 4]];
    data.model.elements[0] = { id: 1, type: 112, nodeIds: [1, 2], properties: { area: 1, youngsModulus: 2, density: 3 } };
    data.analysis.displacements = [0, 0, .01, 0]; data.analysis.loadVector = [0, 0, 10, 0]; data.analysis.reactions = [-10, 0, 0, 0];
    data.analysis.elementResults[0] = { elementId: 1, type: 112, localEndForces: [-10, 0, 10, 0], axialStrain: .002, axialStress: 20, axialForce: 10 };
    return completeFixture(data);
  }
  function modalFixture() {
    var data = fixture(); data.analysis = { type: 'modal', frequenciesHz: [2, 5], angularFrequenciesRadPerSec: [4 * Math.PI, 10 * Math.PI],
      modeShapes: [[0, 0], [0, 0], [0, 0], [.1, -.1], [.2, .3], [.01, -.02]] };
    return completeFixture(data);
  }
  function transientFixture(reduced) {
    var data = fixture(), ids = reduced ? [2] : [1, 2, 3, 4, 5, 6];
    var rows = ids.map(function (id) { return [0, id * .001, id * .002]; });
    data.analysis = { type: 'transient', time: [0, .1, .2], globalDOFIds: ids, displacements: rows, velocities: rows,
      spectrumFrequencyHz: [0, 5, 10], displacementAmplitudeSpectrum: rows,
      sampling: { originalSampleCount: 3, exportedSampleCount: 3, timeStride: 1 } };
    return completeFixture(data);
  }
  function largeFixture(count) {
    var nodes = [], elements = [], dofMap = [], elementResults = [], dofCount = (count + 1) * 2;
    for (var i = 0; i <= count; i += 1) { nodes.push({ id: i + 1, x: i, y: 0 }); dofMap.push([2 * i + 1, 2 * i + 2]); }
    for (var e = 0; e < count; e += 1) { elements.push({ id: e + 1, type: 112, nodeIds: [e + 1, e + 2], properties: { area: 1, youngsModulus: 2, density: 3 } }); elementResults.push({ elementId: e + 1, type: 112, localEndForces: [0, 0, 0, 0], axialStrain: 0, axialStress: 0, axialForce: 0 }); }
    return completeFixture({ format: 'mkef-postprocessor', version: 8, metadata: { title: 'Large', units: { length: '', force: '', moment: '', stress: '', time: '' } },
      model: { dimension: 2, dofPerNode: 2, dofLabels: ['ux', 'uy'], nodes: nodes, elements: elements, dofMap: dofMap, supports: [] },
      analysis: { type: 'static', displacements: new Array(dofCount).fill(0), loadVector: new Array(dofCount).fill(0), reactions: new Array(dofCount).fill(0), equilibriumResidual: [0, 0, 0], elementResults: elementResults } });
  }

  test('schema accepts static, modal, and transient v8 data', function () { assert(M.validateDataset(fixture())); assert(M.validateDataset(modalFixture())); assert(M.validateDataset(transientFixture(false))); });
  test('schema reports the invalid JSON path', function () { var data = fixture(); data.model.elements[0].nodeIds[1] = 99; try { M.validateDataset(data); } catch (error) { assert(error.message.indexOf('$.model.elements[0].nodeIds[1]') === 0); return; } throw new Error('Invalid node reference was accepted'); });
  test('112 interpolation preserves endpoints', function () { var points = M.geometry.sample112({ x: 0, y: 0 }, { x: 2, y: 0 }, { x: 1, y: 2 }, { x: 3, y: 4 }, 2); close(points[0].x, 2); close(points[0].y, 4); close(points[1].x, 8); close(points[1].y, 8); });
  test('113 Hermite interpolation preserves endpoints and rotations', function () { var L = 5, count = 1001, points = M.geometry.sample113({ x: 0, y: 0 }, { x: L, y: 0 }, { x: 0, y: 0, r: .2 }, { x: 0, y: 0, r: -.1 }, 1, count), dx = L / (count - 1); close(points[0].x, 0); close(points[count - 1].x, L); close((points[1].y - points[0].y) / dx, .2, .001); close((points[count - 1].y - points[count - 2].y) / dx, -.1, .001); });
  test('displacement field sampling finds an interior frame maximum', function () { var data = fixture(); data.analysis.displacements = [0, 0, 1, 0, 0, -1]; var dataset = M.validateDataset(data), field = M.staticResults.displacementMagnitudeField(dataset, 31), values = field.byElementId.get(1); close(values[0], 0); close(values[values.length - 1], 0); assert(field.maximum > 0, 'Interior frame displacement was missed'); });
  test('mathematical Y is inverted exactly once', function () { var point = M.geometry.svgPoint({ x: 2, y: 3 }); close(point.x, 2); close(point.y, -3); });







  test('selected static node details focus on displacement and nonzero actions', function () { var dataset = M.validateDataset(trussFixture()), lines = M.staticResults.nodeDetails(dataset, 2).join('\n'); assert(lines.indexOf('|u| = 0.01 m') >= 0); assert(lines.indexOf('ux=0.01 m') >= 0 && lines.indexOf('uy=0 m') >= 0); assert(lines.indexOf('Load: Fx=10 N') >= 0); assert(lines.indexOf('Reaction:') < 0); assert(lines.indexOf('Global DOFs') < 0 && lines.indexOf('Coordinates') < 0); });
  test('selected node reports the displayed element-end result', function () { var frame = M.validateDataset(fixture()), moment = M.staticResults.nodeResultDetails(frame, 1, 'M').join('\n'), truss = M.staticResults.nodeResultDetails(M.validateDataset(trussFixture()), 2, 'axialForce').join('\n'); assert(moment.indexOf('Displayed bending moment M') >= 0 && moment.indexOf('E1 · end 1 = 4 N*m') >= 0, 'Frame diagram value is missing at the selected node.'); assert(truss.indexOf('Displayed axial force') >= 0 && truss.indexOf('E1 · end 2 = 10 N') >= 0, 'Truss result value is missing at the selected node.'); });


  test('modal canonicalization is deterministic and does not mutate JSON', function () { var data = modalFixture(), original = JSON.stringify(data.analysis.modeShapes), dataset = M.validateDataset(data), canonical = M.modalView.canonicalMode(dataset, 1); assert(JSON.stringify(data.analysis.modeShapes) === original, 'Source mode matrix was mutated'); assert(canonical.pivotDOF === 5 && canonical.sign === 1); close(canonical.vector[3], -1 / 3); close(canonical.vector[4], 1); close(canonical.vector[5], -1 / 15); assert(canonical.vector[0] === 0 && canonical.vector[1] === 0 && canonical.vector[2] === 0, 'Fixed DOFs changed'); });
  test('modal sign uses the earliest tied translational component', function () { var data = modalFixture(); data.analysis.modeShapes = [[0], [0], [0], [-.5], [.5], [.25]]; data.analysis.frequenciesHz = [2]; data.analysis.angularFrequenciesRadPerSec = [4 * Math.PI]; data.analysis.eigenvalues=[(4*Math.PI)**2]; var canonical = M.modalView.canonicalMode(M.validateDataset(data), 0); assert(canonical.pivotDOF === 4 && canonical.sign === -1); close(canonical.vector[3], 1); close(canonical.vector[4], -1); close(canonical.vector[5], -.5); });
  test('modal display has a stable sampled peak scale', function () { var dataset = M.validateDataset(modalFixture()), canonical = M.modalView.canonicalMode(dataset, 1), automatic = M.geometry.automaticScale(dataset, canonical.vector, 31), result = M.modalView.display(dataset, 1, 2, 31); close(result.peakScale, 2 * automatic.scale); close(result.scale, result.peakScale); close(result.frequencyHz, 5); assert(!result.zero); });
  test('rotation-only modal shapes render and all-zero modes are disabled', function () { var data = modalFixture(); data.analysis.modeShapes = [[0, 0], [0, 0], [.1, 0], [0, 0], [0, 0], [-.2, 0]]; var dataset = M.validateDataset(data), rotational = M.modalView.display(dataset, 0, 1, 31), zero = M.modalView.display(dataset, 1, 1, 31); assert(!rotational.zero && rotational.pivotDOF === 6 && Number.isFinite(rotational.peakScale) && rotational.peakScale > 0, 'Rotation-only mode was suppressed'); close(rotational.vector[2], -.5); close(rotational.vector[5], 1); assert(zero.zero); close(zero.peakScale, 0); });
  test('modal phase math covers peaks, zero, speed, and navigation bounds', function () { var vector = [1, -2, 0], negative = M.modalView.phaseVector(vector, -1), zero = M.modalView.phaseVector(vector, 0); assert(negative[0] === -1 && negative[1] === 2 && zero.every(function (value) { return value === 0; })); close(M.modalView.phaseAt(0, 0, 1).factor, 1); close(M.modalView.phaseAt(0, 500, 1).factor, 0); close(M.modalView.phaseAt(0, 1000, 1).factor, -1); close(M.modalView.phaseAt(0, 500, 2).factor, -1); close(M.modalView.phaseAngleForFactor(-1), Math.PI); assert(M.modalView.clampModeIndex(-1, 2) === 0 && M.modalView.clampModeIndex(5, 2) === 1); });
  test('transient full and reduced exports are distinguished', function () { var full = M.validateDataset(transientFixture(false)), reduced = M.validateDataset(transientFixture(true)); assert(M.transientView.hasFullDisplacements(full)); assert(!M.transientView.hasFullDisplacements(reduced)); close(M.transientView.vectorAt(full, 'displacements', 2)[5], .012); assert(M.transientView.vectorAt(reduced, 'displacements', 1)[0] === undefined, 'Missing DOF was represented as zero'); assert(M.transientView.series(reduced, 2, 'displacements').y.length === 3); });
  test('transient histories retain exact samples and report physical units', function () { var dataset = M.validateDataset(transientFixture(false)), translation = M.transientView.series(dataset, 4, 'velocities'), rotation = M.transientView.series(dataset, 6, 'displacements'), spectrum = M.transientView.series(dataset, 4, 'spectrum'); assert(translation.y === dataset.raw.analysis.velocities[3], 'History samples were copied or changed'); assert(translation.yLabel.indexOf('m/s') >= 0 && translation.yLabel.indexOf('ux') >= 0); assert(rotation.yLabel.indexOf('rad') >= 0 && spectrum.yLabel.indexOf('m') >= 0); close(translation.y[2], .008); });
  test('transient field and decimation metadata expose omissions', function () { var data = transientFixture(false); delete data.analysis.velocities; data.analysis.sampling = { originalSampleCount: 9, exportedSampleCount: 3, timeStride: 4 }; var dataset = M.validateDataset(data), fields = M.transientView.fieldStatus(dataset), sampling = M.transientView.samplingSummary(dataset); assert(fields.present.indexOf('displacements') >= 0 && fields.omitted.indexOf('velocities') >= 0 && fields.omitted.indexOf('reactions') >= 0); assert(sampling.decimated && sampling.original === 9 && sampling.exported === 3 && sampling.stride === 4); });
  test('transient playback uses exported time values, stride, and exact boundaries', function () { var time = [0, .1, 1, 2]; assert(M.transientView.stepFrameIndex(0, -1, 4, 2) === 0); assert(M.transientView.stepFrameIndex(0, 1, 4, 2) === 2); assert(M.transientView.stepFrameIndex(2, 1, 4, 2) === 3); assert(M.transientView.playbackFrames(4, 0, 2).join(',') === '0,2,3'); assert(M.transientView.playbackFrame(time, 0, 2500, 1, 2, 5000).index === 2, 'Nonuniform exported time was ignored'); var end = M.transientView.playbackFrame(time, 1, 5000, 1, 2, 5000); assert(end.index === 3 && end.complete); assert(M.transientView.playbackFrame(time, 2, 0, 1, 1, 5000).index === 2, 'Resume did not start at the retained frame'); });
  test('transient schema rejects incomplete spectrum pairs and impossible sampling metadata', function () { var missing = transientFixture(false); delete missing.analysis.displacementAmplitudeSpectrum; try { M.validateDataset(missing); throw new Error('Incomplete spectrum was accepted'); } catch (error) { assert(error.message.indexOf('must be provided together') >= 0); } var sampling = transientFixture(false); sampling.analysis.sampling.originalSampleCount = 2; try { M.validateDataset(sampling); throw new Error('Impossible sample counts were accepted'); } catch (error) { assert(error.message.indexOf('cannot be smaller') >= 0); } });
  test('history downsampling preserves extrema and point budget', function () { var x = [], y = []; for (var i = 0; i < 10000; i += 1) { x.push(i); y.push(i === 5001 ? 100 : Math.sin(i)); } var result = M.transientView.downsample({ x: x, y: y, xLabel: 'x', yLabel: 'y' }, 2000); assert(result.x.length <= 2000); assert(Math.max.apply(null, result.y) === 100); });
  test('chart renders cursor and retains original tooltip samples', function () { var chart = document.getElementById('test-chart'), x = [], y = []; for (var i = 0; i < 3000; i += 1) { x.push(i); y.push(i * 2); } M.charts.render(chart, { x: x, y: y, xLabel: 'Time', yLabel: 'u' }, 1); assert(chart.querySelector('.chart-series')); assert(chart.querySelector('.chart-cursor')); assert(chart.querySelector('.chart-hit-area')); assert(chart.__chartState.series.x.length === 3000, 'Original chart samples were discarded'); });
  test('chart keyboard selection reports an original sample', function () { var chart = document.getElementById('test-chart'), selected = null; M.charts.render(chart, { x: [0, .1, .2], y: [0, 2, 4], xLabel: 'Time', yLabel: 'u' }, .1, function (index, x, y) { selected = [index, x, y]; }); chart.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowRight', bubbles: true })); assert(selected && selected[0] === 0 && selected[1] === 0 && selected[2] === 0); });








  function uniformFixture() {
    var data=fixture(); data.model.nodes[1]={id:2,x:2,y:0};
    data.model.nodalLoads=[]; data.model.elementLoads=[{type:20,elementId:1,coordinateSystem:2,qx:0,qy:-1000}];
    data.analysis.loadVector=[0,-1000,-1000/3,0,-1000,1000/3];
    data.analysis.elementResults[0].localEndForces=[0,1000,0,0,1000,0];
    data.analysis.elementResults[0].equivalentLocalLoadVector=data.analysis.loadVector.slice();
    completeFixture(data);data.analysis.elementResults[0].sectionForceCoefficients=[[0,0,0,0],[-1000,2000,0,0],[0,2000,-2000,0]];
    return data;
  }
  function linearFixture() {
    var data=uniformFixture();
    data.model.elementLoads=[{type:21,elementId:1,coordinateSystem:2,qx1:0,qy1:0,qx2:0,qy2:-1000}];
    data.analysis.elementResults[0].localEndForces=[0,1000/3,0,0,2000/3,0];
    data.analysis.elementResults[0].sectionForceCoefficients=[[0,0,0,0],[-1000/3,0,1000,0],[0,2000/3,0,-2000/3]];
    return data;
  }
  test('v8 validates linear loads and rejects unsupported or incomplete data', function () {
    assert(M.validateDataset(linearFixture()));
    [function(d){delete d.model.elementLoads[0].qy2;},function(d){d.model.elementLoads[0].qy2=0;},function(d){d.model.elementLoads[0].qx1=Infinity;},function(d){d.model.releases=[{elementId:1,end:1,component:'Mz'}];},function(d){delete d.analysis.elementResults[0].equivalentLocalLoadVector;}].forEach(function (change) {
      var data=linearFixture(); change(data); var rejected=false; try { M.validateDataset(data); } catch (_) { rejected=true; } assert(rejected);
    });
  });
  test('triangle diagrams recover exact shear endpoints and cubic moment maximum', function () {
    var dataset=M.validateDataset(linearFixture()), e=dataset.raw.model.elements[0], S=M.staticResults;
    close(S.elementDiagram(dataset,e,'V',0),-1000/3); close(S.elementDiagram(dataset,e,'V',1),2000/3);
    close(S.elementDiagram(dataset,e,'M',0),0); close(S.elementDiagram(dataset,e,'M',1),0);
    var xi=1/Math.sqrt(3), maximum=4000/(9*Math.sqrt(3));
    close(S.elementDiagram(dataset,e,'M',xi),maximum);
    close(Math.max.apply(null,S.diagramSamples(dataset,e,'M',6).map(function (v) { return v.value; })),maximum);
  });
  test('sign-changing loads include both moment extrema and N/V stationary points', function () {
    var data=linearFixture(), q=data.model.elementLoads[0]; q.qx1=q.qy1=1000; q.qx2=q.qy2=-1000;
    data.analysis.elementResults[0].localEndForces=[0,-320,0,0,0,0];
    data.analysis.elementResults[0].sectionForceCoefficients=[[0,-2000,2000,0],[320,-2000,2000,0],[0,-640,2000,-4000/3]];
    var dataset=M.validateDataset(data), e=dataset.raw.model.elements[0], S=M.staticResults;
    var moments=S.diagramSamples(dataset,e,'M',6);
    assert(moments.some(function(v){return Math.abs(v.xi-.2)<1e-12;}) && moments.some(function(v){return Math.abs(v.xi-.8)<1e-12;}));
    ['N','V'].forEach(function (quantity) { assert(S.diagramSamples(dataset,e,quantity,6).some(function(v){return Math.abs(v.xi-.5)<1e-12;})); });
    close(S.elementDiagram(dataset,e,'N',.5),-500);
  });

  test('v8 validation rejects malformed load data and missing equivalent vectors', function () {
    assert(M.validateDataset(uniformFixture()));
    [function(d){d.model.elementLoads[0].elementId=7;},function(d){d.model.elementLoads[0].qy=0;},function(d){d.model.elementLoads[0].qx=Infinity;},function(d){delete d.analysis.elementResults[0].equivalentLocalLoadVector;},function(d){d.model.nodalLoads={};}].forEach(function (change) {
      var data=uniformFixture(); change(data); var rejected=false; try {M.validateDataset(data);} catch (_) {rejected=true;} assert(rejected);
    });
  });
  test('uniform beam diagrams recover the interior moment maximum and both end shears', function () {
    var dataset=M.validateDataset(uniformFixture()), e=dataset.raw.model.elements[0];
    close(M.staticResults.elementDiagram(dataset,e,'M',0),0); close(M.staticResults.elementDiagram(dataset,e,'M',1),0);
    close(M.staticResults.elementDiagram(dataset,e,'M',.5),500);
    close(M.staticResults.elementDiagram(dataset,e,'V',0),-1000); close(M.staticResults.elementDiagram(dataset,e,'V',1),1000);
    var values=M.staticResults.diagramSamples(dataset,e,'M',6); close(Math.max.apply(null,values.map(function(v){return v.value;})),500);
    dataset.raw.model.elementLoads[0].qx=20; dataset.resultsByElementId.get(1).localEndForces[0]=-40;
    dataset.resultsByElementId.get(1).sectionForceCoefficients[0]=[40,-40,0,0];
    close(M.staticResults.elementDiagram(dataset,e,'N',0),40); close(M.staticResults.elementDiagram(dataset,e,'N',1),0);
  });



  test('2000-element planar frames retain bounded objects during repeated updates',function(){var r=sceneRenderer(largeFixture(2000));r.updateDeformation(r.dataset.raw.analysis.displacements,1,31);r.flush();r.scene.capture();var objects=r.sceneNodeCount();for(var i=0;i<20;i++){r.updateDeformation(r.dataset.raw.analysis.displacements,i/10,31);r.flush();}assert(r.sceneNodeCount()===objects);assert(r.scene.labelCount<=200);assert(r.frame.original.owners.length===2000);});
  test('legacy nonconsecutive IDs survive shared geometry, selection and export',function(){var data=fixture();data.model.nodes[0].id=10;data.model.nodes[1].id=20;data.model.elements[0].id=7;data.model.elements[0].nodeIds=[10,20];data.model.supports[0].nodeId=10;data.analysis.elementResults[0].elementId=7;var r=sceneRenderer(data);r.setSelection({kind:'node',id:20});assert(r.scene.selection.index===1);r.setSelection({kind:'element',id:7});assert(r.scene.selection.index===0);r.flush();r.scene.capture();assert(r.scene.conditions.data.supports[0].node===1);assert(r.scene.labelFrame.some(l=>l.text==='E7'));assert(r.frame.original.positions.every(Number.isFinite));});
  // Retained WebGL layers replace the retired SVG DOM renderer. Numerical tests above are unchanged.
  var activeRenderer;
  function sceneRenderer(data) {if(activeRenderer)activeRenderer.dispose();var r=activeRenderer=new M.WebGLRenderer(document.getElementById('test-host'));r.mount(M.validateDataset(data));assert(r.scene,'WebGL unavailable');return r;}
  test('shared renderer retains buffers, camera and selection through deformation',function(){var r=sceneRenderer(fixture());r.updateDeformation(r.dataset.raw.analysis.displacements,1,31);r.flush();r.setSelection({kind:'element',id:1});var object=r.scene.results.objects.deformed,view=JSON.stringify(r.scene.getViewState()),count=r.sceneNodeCount();r.updateDeformation([0,0,0,.01,.02,0],2,31);r.setSelection({kind:'element',id:1});r.flush();assert(r.scene.results.objects.deformed===object);assert(r.sceneNodeCount()===count);var a=JSON.parse(view),b=r.scene.getViewState();close(a.height,b.height);a.target.forEach((v,i)=>close(v,b.target[i]));assert(r.scene.selection.index===0);});
  test('deformed nodes are pickable at displaced physical coordinates',function(){var r=sceneRenderer(fixture());r.updateDeformation([0,0,0,.01,.02,0],10,31);r.flush();var p=r.frame.resultNodes.positions;close(p[3],.6);close(p[4],.7);r.scene.capture();var point=r.scene.local({x:.6,y:.7,z:0}).project(r.scene.camera),rect=r.scene.canvas.getBoundingClientRect();r.scene.pick(rect.left+(point.x+1)*rect.width/2,rect.top+(1-point.y)*rect.height/2);assert(r.selection&&r.selection.kind==='node'&&r.selection.id===2);});
  test('screen labels and nodes remain bounded at different zoom levels',function(){var r=sceneRenderer(fixture());r.scene.capture();assert(r.scene.labelFrame.length<=200);assert(r.scene.nodes.material.size===11);r.scene.camera.zoom*=2;r.scene.capture();assert(r.scene.nodes.material.size===11);assert(r.scene.labelFrame.every(l=>Number.isFinite(l.x)&&Number.isFinite(l.y)));});
  test('planar supports map physical rotations to Z and preserve all restraint types',function(){for(var type=1;type<=7;type++){var d=fixture();d.model.supports=[{type:type,nodeId:1}];var ds=M.validateDataset(d),c=M.planarPresentation.conditions(ds),mask=MKEFSupports.mask(d.model.supports[0],3);assert(c.supports[0].mask.join(',')===[+mask[0],+mask[1],0,0,0,+mask[2]].join(','));}});
  test('forces and moments produce finite reaction geometry without equivalent load duplication',function(){var r=sceneRenderer(fixture());M.staticResults.mount(r,r.dataset);r.flush();assert(r.frame.loads.positions.length>0&&r.frame.reactions.positions.length>0);assert(r.frame.reactions.positions.every(Number.isFinite));r=sceneRenderer(uniformFixture());M.staticResults.mount(r,r.dataset);r.flush();assert(!r.frame.loads.positions.length);r.setDistributedLoads(1,true,true);r.flush();var count=r.frame.distributed.positions.length;assert(count>0);r.setDistributedLoads(1,false,true);r.flush();assert(r.frame.distributed.positions.length===count/2);r.setDistributedLoads(1,false,false);r.flush();assert(!r.frame.distributed.positions.length);});
  test('displacement gradients retain exact extrema and zero neutral color',function(){var r=sceneRenderer(trussFixture()),legend=r.applyResult(r.dataset,'displacementMagnitude',31);r.flush();assert(legend.kind==='color'&&legend.unit==='m'&&legend.ticks.length===3);close(legend.max,.01);assert(r.colors.get(1).length<=9);var values=JSON.stringify(r.frame.original.colors);r.updateDeformation(r.dataset.raw.analysis.displacements,25,31);r.flush();assert(JSON.stringify(r.frame.original.colors)===values);r=sceneRenderer(fixture());legend=r.applyResult(r.dataset,'displacementMagnitude',31);assert(legend.ticks.length===1&&legend.ticks[0]===0);assert(r.colors.get(1)==='#66788a');});
  test('diagrams preserve end-force convention, fill and signed element colors',function(){var data=uniformFixture();data.analysis.elementResults[0].sectionForceCoefficients[2]=[50,-50,0,0];var ds=M.validateDataset(data);close(M.staticResults.elementDiagram(ds,ds.raw.model.elements[0],'M',0),50);close(M.staticResults.elementDiagram(ds,ds.raw.model.elements[0],'M',1),0);var r=sceneRenderer(uniformFixture()),legend=r.applyResult(r.dataset,'M',6);r.flush();assert(legend.text.startsWith('M diagram'));assert(r.frame.diagrams.positions.length>12&&r.frame.diagramFill.positions.length>0);r=sceneRenderer(trussFixture());r.applyResult(r.dataset,'axialForce',31);r.flush();assert(r.frame.original.colors.length>0);});
  test('standalone vector export escapes text and includes labels, legend, chart and diagram fills',function(){var r=sceneRenderer(trussFixture()),legend=r.applyResult(r.dataset,'displacementMagnitude',31);r.updateDeformation(r.dataset.raw.analysis.displacements,1,31);r.flush();var text=M.exporting.serialize(M.exporting.sceneSource(r.scene,false),{title:'<unsafe & title>',analysisType:'static'},document.getElementById('test-chart'),legend);assert(text.includes('&lt;unsafe &amp; title&gt;'));assert(text.includes('data-export-legend')&&text.includes('Displacement magnitude |u| [m]')&&text.includes('chart-series'));assert(text.includes('data-layer="deformed"')&&!text.includes('<image'));r=sceneRenderer(uniformFixture());r.setDistributedLoads(1,true,false);r.applyResult(r.dataset,'M',6);r.flush();text=new XMLSerializer().serializeToString(M.exporting.sceneSource(r.scene,false));assert(text.includes('diagramFill')&&text.includes('distributed'));});
  for(const type of ['modal','transient'])test(type+' export freezes current geometry and includes snapshot metadata',function(){var r=sceneRenderer(type==='modal'?modalFixture():transientFixture(false)),vector=type==='modal'?M.modalView.display(r.dataset,0,1,31).vector:M.transientView.vectorAt(r.dataset,'displacements',1);r.updateDeformation(vector,-.5,31);r.flush();var source=M.exporting.sceneSource(r.scene,false),before=new XMLSerializer().serializeToString(source),context={title:type,analysisType:type,modalPhaseFactor:-.5,transientSampling:{original:3,exported:3,stride:1},captionLines:['Current sample / phase']};r.updateDeformation(vector,10,31);r.flush();assert(new XMLSerializer().serializeToString(source)===before);var text=M.exporting.serialize(source,context,null,null);assert(text.includes('modalPhaseFactor')&&text.includes('transientSampling')&&text.includes('data-export-caption'));});
  for(const type of ['static','modal','transient'])test(type+' WebGL PNG includes a frozen scene and export caption',async function(){var r=sceneRenderer(type==='static'?uniformFixture():type==='modal'?modalFixture():transientFixture(false));r.flush();var source=M.exporting.sceneSource(r.scene,true);assert(source.querySelector('image').getAttribute('href').startsWith('data:image/png'));var blob=await M.exporting.pngBlob(source,{title:type,analysisType:type,captionLines:['Frame snapshot']},1,document.getElementById('test-chart'),null);assert(blob.type==='image/png'&&blob.size>100);});

  document.addEventListener('DOMContentLoaded', async function () {
    var lines = [], failed = 0;
    for (var i = 0; i < tests.length; i += 1) { var entry = tests[i]; try { await entry.operation(); lines.push('PASS ' + entry.name); } catch (error) { failed += 1; lines.push('FAIL ' + entry.name + ': ' + error.message); } }
    document.getElementById('test-output').textContent = lines.join('\n') + '\n\n' + (failed ? failed + ' failed' : tests.length + ' passed');
    document.body.setAttribute('data-test-status', failed ? 'failed' : 'passed');
  });
}(window.MKEFPost));
