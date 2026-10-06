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
  function diagonal(dataset) { var b = M.geometry.modelBounds(dataset.raw.model); return Math.max(Math.hypot(b.maxX - b.minX, b.maxY - b.minY), 1e-9); }
  function elementNormal(dataset, element) { var a = dataset.nodesById.get(element.nodeIds[0]), b = dataset.nodesById.get(element.nodeIds[1]), basis = M.geometry.elementBasis(a, b); return { x: -basis.s, y: basis.c }; }
  function appendPath(layer, attributes, titleText) { var path = M.svgElement('path', attributes), title = M.svgElement('title'); title.textContent = titleText; path.appendChild(title); layer.appendChild(path); return path; }
  function arrowPath(node, vector, length) {
    var magnitude = Math.hypot(vector.x, vector.y), ux = vector.x / magnitude, uy = vector.y / magnitude;
    var tail = { x: node.x - ux * length, y: node.y - uy * length }, normal = { x: -uy, y: ux }, head = length * .22;
    var a = M.geometry.svgPoint(tail), b = M.geometry.svgPoint(node);
    var h1 = M.geometry.svgPoint({ x: node.x - ux * head + normal.x * head * .45, y: node.y - uy * head + normal.y * head * .45 });
    var h2 = M.geometry.svgPoint({ x: node.x - ux * head - normal.x * head * .45, y: node.y - uy * head - normal.y * head * .45 });
    return 'M ' + a.x + ' ' + a.y + ' L ' + b.x + ' ' + b.y + ' M ' + h1.x + ' ' + h1.y + ' L ' + b.x + ' ' + b.y + ' L ' + h2.x + ' ' + h2.y;
  }
  function momentPath(node, value, maximum, size) {
    var radius = size * (.32 + .3 * Math.abs(value) / maximum), p = M.geometry.svgPoint(node);
    return window.MKEFMomentSymbol.path(p.x, p.y, radius, value);
  }
  function restrainedDOFs(dataset) {
    var restrained = new Set(), model = dataset.raw.model;
    model.supports.forEach(function (support) { var row = model.nodes.findIndex(function (node) { return node.id === support.nodeId; }); M.geometry.restrainedLocalDOFs(support.type, model.dofPerNode).forEach(function (local) { restrained.add(model.dofMap[row][local]); }); });
    return restrained;
  }
  function nodeComponents(dataset, vector, onlyRestrained) {
    var model = dataset.raw.model, restrained = onlyRestrained ? restrainedDOFs(dataset) : null;
    return model.nodes.map(function (node, row) { var dofs = model.dofMap[row]; function value(local) { return dofs[local] && (!restrained || restrained.has(dofs[local])) ? vector[dofs[local] - 1] : 0; } return { node: node, x: value(0), y: value(1), moment: dofs.length > 2 ? value(2) : 0 }; });
  }
  function tolerance(values) { return 1e-10 * Math.max(1, Math.max.apply(null, values.map(Math.abs).concat([0]))); }
  function textLegend(text) { return text ? { kind: 'text', text: text, summary: text } : { kind: 'none', summary: '' }; }
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

  S.mount = function (renderer, dataset) {
    renderer.clearLayer('loads'); renderer.clearLayer('reactions');
    var sets = [
      { layer: renderer.layers.loads, values: nodeComponents(dataset, S.nodalLoadVector(dataset), false), forceClass: 'load-symbol load-force-symbol', momentClass: 'load-symbol load-moment-symbol', prefix: 'Load' },
      { layer: renderer.layers.reactions, values: nodeComponents(dataset, (dataset.raw.analysis.supportReactions || dataset.raw.analysis.reactions), true), forceClass: 'reaction-symbol reaction-force-symbol', momentClass: 'reaction-symbol reaction-moment-symbol', prefix: 'Reaction' }
    ];
    var size = diagonal(dataset) * .13;
    sets.forEach(function (set) {
      var forceValues = [], momentValues = [];
      set.values.forEach(function (item) { forceValues.push(item.x, item.y); momentValues.push(item.moment); });
      var forceMaximum = Math.max.apply(null, forceValues.map(Math.abs).concat([0])), momentMaximum = Math.max.apply(null, momentValues.map(Math.abs).concat([0]));
      var forceTolerance = tolerance(forceValues), momentTolerance = tolerance(momentValues);
      set.values.forEach(function (item) {
        var vector = { x: Math.abs(item.x) <= forceTolerance ? 0 : item.x, y: Math.abs(item.y) <= forceTolerance ? 0 : item.y }, magnitude = Math.hypot(vector.x, vector.y);
        if (magnitude) {
          var length = size * (.35 + .65 * magnitude / (forceMaximum || 1));
          appendPath(set.layer, { d: arrowPath(item.node, vector, length), class: set.forceClass, 'data-node-id': item.node.id }, set.prefix + ' (' + withUnit(item.x, unit(dataset, 'force')) + ', ' + withUnit(item.y, unit(dataset, 'force')) + ')');
        }
        if (Math.abs(item.moment) > momentTolerance) appendPath(set.layer, { d: momentPath(item.node, item.moment, momentMaximum || 1, size), class: set.momentClass, 'data-node-id': item.node.id }, set.prefix + ' Mz ' + withUnit(item.moment, unit(dataset, 'moment')));
      });
    });
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
    var labels = { axialForce: !dataset.raw.model.elementLoads ? 'axial force' : 'mean axial force', axialStress: !dataset.raw.model.elementLoads ? 'axial stress' : 'mean axial stress', N: 'axial force N', V: 'shear force V', M: 'bending moment M' };
    var suffix = mode === 'axialStress' ? unit(dataset, 'stress') : mode === 'M' ? unit(dataset, 'moment') : unit(dataset, 'force');
    var lines = ['Displayed ' + labels[mode] + ' at connected element' + (connected.length === 1 ? '' : 's') + ':'];
    connected.forEach(function (element) {
      var result = dataset.resultsByElementId.get(element.id), endIndex = element.nodeIds[0] === nodeId ? 0 : 1;
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

  S.drawElementLoads = function (renderer,dataset,scale,showOriginal,showDeformed) {
    Array.prototype.forEach.call(renderer.layers.loads.querySelectorAll('.element-load-symbol'),function (item) { item.remove(); });
    (dataset.raw.model.elementLoads || []).forEach(function (load) {
      var element=dataset.elementsById.get(load.elementId), a=dataset.nodesById.get(element.nodeIds[0]), b=dataset.nodesById.get(element.nodeIds[1]), basis=M.geometry.elementBasis(a,b);
      var ends=loadEnds(load), vectors=ends.map(function (q) { return load.coordinateSystem===1 ? {x:basis.c*q.x-basis.s*q.y,y:basis.s*q.x+basis.c*q.y} : q; });
      var maximum=Math.max(Math.hypot(vectors[0].x,vectors[0].y),Math.hypot(vectors[1].x,vectors[1].y));
      if (!maximum) return;
      var size=Math.min(diagonal(dataset)*.08,basis.length*.3);
      var units=unit(dataset,'force') && unit(dataset,'length') ? unit(dataset,'force')+'/'+unit(dataset,'length') : 'force/length';
      function pair(q) { return '('+format(q.x)+', '+format(q.y)+')'; }
      var label=(load.coordinateSystem===1 ? 'Local' : 'Global')+(load.type===21 ? ' q1='+pair(ends[0])+' → q2='+pair(ends[1]) : ' q='+pair(ends[0]))+' '+units;
      [showOriginal ? 0 : null,showDeformed && (scale!==0 || !showOriginal) ? scale : null].forEach(function (factor) {
        if (factor===null) return;
        var points=M.geometry.sampleElement(dataset,element,dataset.raw.analysis.displacements,factor,9), tails=[];
        points.forEach(function (point,i) {
          var xi=i/(points.length-1), vector={x:vectors[0].x*(1-xi)+vectors[1].x*xi,y:vectors[0].y*(1-xi)+vectors[1].y*xi};
          var magnitude=Math.hypot(vector.x,vector.y);
          tails.push({x:point.x-vector.x/maximum*size,y:point.y-vector.y/maximum*size});
          if (magnitude) appendPath(renderer.layers.loads,{d:arrowPath(point,vector,size*magnitude/maximum),class:'element-load-symbol load-symbol load-force-symbol','data-element-id':element.id},label);
        });
        appendPath(renderer.layers.loads,{d:M.geometry.pathData(tails),class:'element-load-symbol element-load-envelope load-symbol load-force-symbol',fill:'none','data-element-id':element.id},label);
        var anchor=tails[Math.floor(tails.length/2)], position=M.geometry.svgPoint(anchor);
        position.y=Math.min.apply(null,tails.map(function (point) { return M.geometry.svgPoint(point).y; }));
        var text=M.svgElement('text',{x:position.x,y:position.y-diagonal(dataset)*.015,'font-size':diagonal(dataset)*.018,'text-anchor':'middle',fill:'#b63a3f',class:'element-load-symbol element-load-label'});
        text.textContent=label; renderer.layers.loads.appendChild(text);
      });
    });
  };
  function signedColor(value, maximum) { if (!maximum || value === 0) return '#66788a'; var amount = Math.min(Math.abs(value) / maximum, 1), light = Math.round(68 - amount * 25); return value > 0 ? 'hsl(4 58% ' + light + '%)' : 'hsl(218 55% ' + light + '%)'; }

  S.applyResult = function (renderer, dataset, mode, sampleCount) {
    renderer.clearLayer('diagrams'); renderer.resetElementColors(); renderer.layers.diagrams.style.display = mode === 'none' ? 'none' : '';
    if (mode === 'none') return textLegend('');
    if (mode === 'displacementMagnitude') return S.applyDisplacementMagnitude(renderer, dataset, sampleCount);
    if (mode === 'axialForce' || mode === 'axialStress') {
      var values = dataset.raw.model.elements.map(function (element) { return dataset.resultsByElementId.get(element.id)[mode]; });
      var maximum = Math.max.apply(null, values.map(Math.abs).concat([0]));
      dataset.raw.model.elements.forEach(function (element, index) { renderer.setElementColor(element.id, signedColor(values[index], maximum)); });
      return textLegend((mode === 'axialForce' ? 'Mean axial force' : 'Mean axial stress') + ' · max |value| ' + withUnit(maximum, unit(dataset, mode === 'axialForce' ? 'force' : 'stress')) + ' · red tension, blue compression');
    }
    var maximumDiagram = 0, samplesByElement = new Map();
    dataset.raw.model.elements.forEach(function (element) {
      var values = S.diagramSamples(dataset, element, mode, sampleCount); samplesByElement.set(element.id, values);
      values.forEach(function (item) { maximumDiagram = Math.max(maximumDiagram, Math.abs(item.value)); });
    });
    var diagramScale = maximumDiagram ? diagonal(dataset) * .14 / maximumDiagram : 0;
    dataset.raw.model.elements.forEach(function (element) {
      var first = dataset.nodesById.get(element.nodeIds[0]), second = dataset.nodesById.get(element.nodeIds[1]), normal = elementNormal(dataset, element);
      var points = samplesByElement.get(element.id).map(function (item) { return { x:first.x+(second.x-first.x)*item.xi+normal.x*item.value*diagramScale, y:first.y+(second.y-first.y)*item.xi+normal.y*item.value*diagramScale }; });
      appendPath(renderer.layers.diagrams, { class:'diagram-fill', 'data-element-id':element.id, d:M.geometry.pathData([first].concat(points,[second]))+' Z', tabindex:'0' }, mode+' diagram for element '+element.id);
      appendPath(renderer.layers.diagrams, { class:'diagram-line', 'data-element-id':element.id, d:M.geometry.pathData(points) }, mode+' diagram edge for element '+element.id);
    });
    return textLegend(mode + ' diagram · max |value| ' + withUnit(maximumDiagram, mode === 'M' ? unit(dataset, 'moment') : unit(dataset, 'force')) + ' · display scale ' + format(diagramScale));
  };
  S.format = format; S.withUnit = withUnit; S.viridis = viridis;
}(window.MKEFPost));
