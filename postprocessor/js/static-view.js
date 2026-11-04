(function (M) {
  'use strict';
  var S = M.staticResults = {};
  var viridis = [
    { offset: 0, color: '#440154' },
    { offset: .25, color: '#3b528b' },
    { offset: .5, color: '#21918c' },
    { offset: .75, color: '#5ec962' },
    { offset: 1, color: '#fde725' }
  ];
  function unit(dataset, name) { return dataset.raw.metadata.units[name] || ''; }
  function format(value) { if (value === 0) return '0'; var a = Math.abs(value); return a >= 10000 || a < .001 ? value.toExponential(3) : Number(value.toPrecision(5)).toString(); }
  function withUnit(value, suffix) { return format(value) + (suffix ? ' ' + suffix : ''); }
  function tolerance(values) { return 1e-10 * Math.max(1, Math.max.apply(null, values.map(Math.abs).concat([0]))); }
  function rgb(hex) { return [parseInt(hex.slice(1, 3), 16), parseInt(hex.slice(3, 5), 16), parseInt(hex.slice(5, 7), 16)]; }
  function colorAt(amount) {
    var value = Math.max(0, Math.min(1, amount));
    var right = viridis.findIndex(function (stop) { return stop.offset >= value; });
    if (right <= 0) return viridis[0].color;
    var a = viridis[right - 1], b = viridis[right], span = b.offset - a.offset, fraction = span ? (value - a.offset) / span : 0;
    var first = rgb(a.color), second = rgb(b.color);
    return '#' + first.map(function (component, index) { return Math.round(component + (second[index] - component) * fraction).toString(16).padStart(2, '0'); }).join('');
  }
  function limitedStops(values, maximum) {
    var count = Math.min(9, values.length), stops = [], used = new Set();
    for (var i = 0; i < count; i += 1) {
      var index = Math.round(i * (values.length - 1) / Math.max(count - 1, 1));
      if (!used.has(index)) {
        used.add(index);
        stops.push({ offset: index / (values.length - 1), color: colorAt(maximum ? values[index] / maximum : 0) });
      }
    }
    return stops;
  }

  S.nodalLoadVector = function (dataset) {
    if (!dataset.raw.model.elementLoads) return dataset.raw.analysis.loadVector;
    var vector = new Array(dataset.raw.analysis.loadVector.length).fill(0);
    dataset.raw.model.nodalLoads.forEach(function (load) {
      var dofs = dataset.raw.model.dofMap[dataset.nodeIndexById.get(load.nodeId)];
      [load.fx, load.fy, load.mz].forEach(function (value, index) { if (dofs[index]) vector[dofs[index]-1] += value; });
    });
    return vector;
  };

  S.mount = function(renderer,dataset) {
    renderer.setForces('loads',S.nodalLoadVector(dataset),'#b63a3f',false);
    renderer.setForces('reactions',dataset.raw.analysis.supportReactions||dataset.raw.analysis.reactions,'#6d3ca5',true);
  };

  S.displacementMagnitudeField = function (dataset, sampleCount) {
    var byElementId = new Map(), maximum = 0, displacement = dataset.raw.analysis.displacements;
    dataset.raw.model.elements.forEach(function (element) {
      var values = M.geometry.sampleElementDisplacements(dataset, element, displacement, sampleCount).map(function (component) {
        var value = Math.hypot(component.x, component.y);
        maximum = Math.max(maximum, value);
        return value;
      });
      byElementId.set(element.id, values);
    });
    return { maximum: maximum, byElementId: byElementId };
  };

  S.applyDisplacementMagnitude = function (renderer, dataset, sampleCount) {
    var field = S.displacementMagnitudeField(dataset, sampleCount || M.config.defaultFrameSamples);
    dataset.raw.model.elements.forEach(function (element) {
      var values = field.byElementId.get(element.id);
      if (!field.maximum) renderer.setElementColor(element.id, '#66788a');
      else renderer.setElementGradient(element.id, dataset.nodesById.get(element.nodeIds[0]), dataset.nodesById.get(element.nodeIds[1]), limitedStops(values, field.maximum));
    });
    return {
      kind: 'color',
      label: 'Displacement magnitude |u|',
      unit: unit(dataset, 'length'),
      min: 0,
      max: field.maximum,
      ticks: field.maximum ? [0, field.maximum / 2, field.maximum] : [0],
      stops: viridis.map(function (stop) { return { offset: stop.offset, color: stop.color }; }),
      summary: field.maximum ? 'Maximum ' + withUnit(field.maximum, unit(dataset, 'length')) : 'Zero displacement'
    };
  };

  function compactVector(dataset, vector, nodeId, prefix) {
    var row = dataset.nodeIndexById.get(nodeId), dofs = dataset.raw.model.dofMap[row], threshold = tolerance(vector), values = [];
    dofs.forEach(function (id, index) {
      if (!id) return;
      var value = vector[id - 1];
      if (Math.abs(value) <= threshold) return;
      var label = index === 0 ? 'Fx' : index === 1 ? 'Fy' : 'Mz';
      values.push(label + '=' + withUnit(value, unit(dataset, index < 2 ? 'force' : 'moment')));
    });
    return values.length ? prefix + ': ' + values.join(', ') : '';
  }

  S.nodeDetails = function (dataset, nodeId) {
    var model = dataset.raw.model, row = dataset.nodeIndexById.get(nodeId), dofs = model.dofMap[row], displacement = dataset.raw.analysis.displacements;
    var values = dofs.map(function (id) { return displacement[id - 1]; });
    var lengthUnit = unit(dataset, 'length'), components = [];
    components.push((model.dofLabels[0] || 'ux') + '=' + withUnit(values[0], lengthUnit));
    components.push((model.dofLabels[1] || 'uy') + '=' + withUnit(values[1], lengthUnit));
    if (values.length > 2) components.push((model.dofLabels[2] || 'rz') + '=' + (dofs[2] ? withUnit(values[2], 'rad') : 'not defined'));
    var lines = ['Node ' + nodeId, '|u| = ' + withUnit(Math.hypot(values[0], values[1]), lengthUnit), 'Components: ' + components.join(', ')];
    var load = compactVector(dataset, S.nodalLoadVector(dataset), nodeId, !dataset.raw.model.elementLoads ? 'Load' : 'Nodal load');
    var reaction = compactVector(dataset, (dataset.raw.analysis.supportReactions || dataset.raw.analysis.reactions), nodeId, 'Support reaction');
    if (load) lines.push(load);
    if (reaction) lines.push(reaction);
    return lines;
  };

  S.nodeResultDetails = function (dataset, nodeId, mode) {
    if (mode === 'none' || mode === 'displacementMagnitude') return [];
    var connected = dataset.raw.model.elements.filter(function (element) { return element.nodeIds.indexOf(nodeId) >= 0; });
    if (!connected.length) return [];
    var labels = { sigmaMin:'normal stress σmin',sigmaMax:'normal stress σmax',axialForce: !dataset.raw.model.elementLoads ? 'axial force' : 'mean axial force', axialStress: !dataset.raw.model.elementLoads ? 'axial stress' : 'mean axial stress', N: 'axial force N', V: 'shear force V', M: 'bending moment M' };
    var suffix = (mode === 'axialStress' || M.sectionResults.isStress(mode)) ? unit(dataset, 'stress') : mode === 'M' ? unit(dataset, 'moment') : unit(dataset, 'force');
    var lines = ['Displayed ' + labels[mode] + ' at connected element' + (connected.length === 1 ? '' : 's') + ':'];
    connected.forEach(function (element) {
      var result = dataset.resultsByElementId.get(element.id), endIndex = element.nodeIds[0] === nodeId ? 0 : 1;
      if(M.sectionResults.isStress(mode)) {var data=M.sectionResults.get(dataset,element);lines.push('E'+element.id+' · end '+(endIndex+1)+' = '+(data ? withUnit(M.sectionResults.value(data,mode,endIndex),suffix) : 'unavailable (no section)'));return;}
      var value = mode === 'axialForce' || mode === 'axialStress' ? result[mode] : S.elementDiagram(dataset, element, mode, endIndex);
      lines.push('E' + element.id + ' · end ' + (endIndex + 1) + ' = ' + withUnit(value, suffix));
    });
    return lines;
  };

  S.frameDiagram = function (q, quantity, xi) {
    if (quantity === 'N') return -q[0];
    if (quantity === 'V') return -q[1];
    if (quantity === 'M') return (1 - xi) * (-q[2]) + xi * q[5];
    throw new Error('Unsupported frame diagram ' + quantity + '.');
  };

  function loadEnds(load) {
    return load.type===21 ? [{x:load.qx1,y:load.qy1},{x:load.qx2,y:load.qy2}] : [{x:load.qx,y:load.qy},{x:load.qx,y:load.qy}];
  }
  S.elementIntensities = function (dataset, element) {
    var a=dataset.nodesById.get(element.nodeIds[0]), b=dataset.nodesById.get(element.nodeIds[1]);
    var basis=M.geometry.elementBasis(a,b), ends=[{x:0,y:0},{x:0,y:0}];
    (dataset.raw.model.elementLoads || []).forEach(function (load) {
      if (load.elementId!==element.id) return;
      loadEnds(load).forEach(function (q,i) {
        ends[i].x+=load.coordinateSystem===1 ? q.x : basis.c*q.x+basis.s*q.y;
        ends[i].y+=load.coordinateSystem===1 ? q.y : -basis.s*q.x+basis.c*q.y;
      });
    });
    return {x:ends[0].x,y:ends[0].y,dx:ends[1].x-ends[0].x,dy:ends[1].y-ends[0].y,length:basis.length};
  };
  S.elementDiagram = function (dataset, element, quantity, xi) {
    var f=dataset.resultsByElementId.get(element.id).localEndForces;
    if (!dataset.raw.model.elementLoads) return S.frameDiagram(f,quantity,xi);
    var q=S.elementIntensities(dataset,element), x=xi*q.length;
    if (quantity==='N') return -f[0]-q.x*x-q.dx*x*xi/2;
    if (quantity==='V') return -f[1]-q.y*x-q.dy*x*xi/2;
    if (quantity==='M') return -f[2]+f[1]*x+q.y*x*x/2+q.dy*x*x*xi/6;
    throw new Error('Unsupported frame diagram '+quantity+'.');
  };
  // Scaled, cancellation-resistant roots; linear/uniform cases are exact limits.
  function quadraticRoots(a,b,c) {
    var scale=Math.max(Math.abs(a),Math.abs(b),Math.abs(c));
    if (!scale) return [];
    a/=scale; b/=scale; c/=scale;
    if (Math.abs(a)<1e-14) return b ? [-c/b] : [];
    var d=b*b-4*a*c;
    if (d < -1e-14) return [];
    var t=-.5*(b+(b>=0 ? 1 : -1)*Math.sqrt(Math.max(0,d)));
    return t ? [t/a,c/t] : [-b/(2*a)];
  }
  S.diagramSamples = function (dataset, element, quantity, count) {
    if(M.sectionResults.isStress(quantity))return M.sectionResults.samples(dataset,element,quantity,count);
    var positions=[], n=Math.max(5,count || M.config.defaultFrameSamples);
    for (var i=0;i<n;i++) positions.push(i/(n-1));
    if (dataset.raw.model.elementLoads) {
      var q=S.elementIntensities(dataset,element), f=dataset.resultsByElementId.get(element.id).localEndForces, roots=[];
      if (quantity==='M') roots=quadraticRoots(q.dy*q.length/2,q.y*q.length,f[1]);
      else if (quantity==='N' && q.dx) roots=[-q.x/q.dx];
      else if (quantity==='V' && q.dy) roots=[-q.y/q.dy];
      roots.forEach(function (xi) { if (xi>0 && xi<1) positions.push(xi); });
    }
    return Array.from(new Set(positions)).sort(function (a,b) { return a-b; }).map(function (xi) { return {xi:xi,value:S.elementDiagram(dataset,element,quantity,xi)}; });
  };

  S.drawElementLoads = function(renderer,dataset,scale,showOriginal,showDeformed) {renderer.setDistributedLoads(scale,showOriginal,showDeformed);};
  S.applyResult = function(renderer,dataset,mode,sampleCount) {return renderer.applyResult(dataset,mode,sampleCount);};
  S.format = format; S.withUnit = withUnit; S.viridis = viridis;
}(window.MKEFPost));
