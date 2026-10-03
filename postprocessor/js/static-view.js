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
    var radius = size * (.32 + .3 * Math.abs(value) / maximum), p = M.geometry.svgPoint(node), sweep = value > 0 ? 0 : 1;
    return 'M ' + (p.x + radius) + ' ' + p.y + ' A ' + radius + ' ' + radius + ' 0 1 ' + sweep + ' ' + (p.x - radius * .7) + ' ' + (p.y - radius * .7) +
      ' M ' + (p.x - radius * .7) + ' ' + (p.y - radius * .7) + ' l ' + (value > 0 ? radius * .05 : radius * .35) + ' ' + radius * .04 +
      ' M ' + (p.x - radius * .7) + ' ' + (p.y - radius * .7) + ' l ' + radius * .04 + ' ' + (value > 0 ? radius * .35 : radius * .05);
  }
  function restrainedDOFs(dataset) {
    var restrained = new Set(), model = dataset.raw.model;
    model.supports.forEach(function (support) { var row = model.nodes.findIndex(function (node) { return node.id === support.nodeId; }); M.geometry.restrainedLocalDOFs(support.type, model.dofPerNode).forEach(function (local) { restrained.add(model.dofMap[row][local]); }); });
    return restrained;
  }
  function nodeComponents(dataset, vector, onlyRestrained) {
    var model = dataset.raw.model, restrained = onlyRestrained ? restrainedDOFs(dataset) : null;
    return model.nodes.map(function (node, row) { var dofs = model.dofMap[row]; function value(local) { return !restrained || restrained.has(dofs[local]) ? vector[dofs[local] - 1] : 0; } return { node: node, x: value(0), y: value(1), moment: dofs.length > 2 ? value(2) : 0 }; });
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

  S.mount = function (renderer, dataset) {
    renderer.clearLayer('loads'); renderer.clearLayer('reactions');
    var sets = [
      { layer: renderer.layers.loads, values: nodeComponents(dataset, dataset.raw.analysis.loadVector, false), forceClass: 'load-symbol load-force-symbol', momentClass: 'load-symbol load-moment-symbol', prefix: 'Load' },
      { layer: renderer.layers.reactions, values: nodeComponents(dataset, dataset.raw.analysis.reactions, true), forceClass: 'reaction-symbol reaction-force-symbol', momentClass: 'reaction-symbol reaction-moment-symbol', prefix: 'Reaction' }
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
    if (values.length > 2) components.push((model.dofLabels[2] || 'rz') + '=' + withUnit(values[2], 'rad'));
    var lines = ['Node ' + nodeId, '|u| = ' + withUnit(Math.hypot(values[0], values[1]), lengthUnit), 'Components: ' + components.join(', ')];
    var load = compactVector(dataset, dataset.raw.analysis.loadVector, nodeId, 'Load');
    var reaction = compactVector(dataset, dataset.raw.analysis.reactions, nodeId, 'Reaction');
    if (load) lines.push(load);
    if (reaction) lines.push(reaction);
    return lines;
  };

  S.nodeResultDetails = function (dataset, nodeId, mode) {
    if (mode === 'none' || mode === 'displacementMagnitude') return [];
    var connected = dataset.raw.model.elements.filter(function (element) { return element.nodeIds.indexOf(nodeId) >= 0; });
    if (!connected.length) return [];
    var labels = { axialForce: 'axial force', axialStress: 'axial stress', N: 'axial force N', V: 'shear force V', M: 'bending moment M' };
    var suffix = mode === 'axialStress' ? unit(dataset, 'stress') : mode === 'M' ? unit(dataset, 'moment') : unit(dataset, 'force');
    var lines = ['Displayed ' + labels[mode] + ' at connected element' + (connected.length === 1 ? '' : 's') + ':'];
    connected.forEach(function (element) {
      var result = dataset.resultsByElementId.get(element.id), endIndex = element.nodeIds[0] === nodeId ? 0 : 1;
      var value = mode === 'axialForce' || mode === 'axialStress' ? result[mode] : S.frameDiagram(result.localEndForces, mode, endIndex);
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
  function signedColor(value, maximum) { if (!maximum || value === 0) return '#66788a'; var amount = Math.min(Math.abs(value) / maximum, 1), light = Math.round(68 - amount * 25); return value > 0 ? 'hsl(4 58% ' + light + '%)' : 'hsl(218 55% ' + light + '%)'; }

  S.applyResult = function (renderer, dataset, mode, sampleCount) {
    renderer.clearLayer('diagrams'); renderer.resetElementColors(); renderer.layers.diagrams.style.display = mode === 'none' ? 'none' : '';
    if (mode === 'none') return textLegend('');
    if (mode === 'displacementMagnitude') return S.applyDisplacementMagnitude(renderer, dataset, sampleCount);
    if (mode === 'axialForce' || mode === 'axialStress') {
      var values = dataset.raw.model.elements.map(function (element) { return dataset.resultsByElementId.get(element.id)[mode]; });
      var maximum = Math.max.apply(null, values.map(Math.abs).concat([0]));
      dataset.raw.model.elements.forEach(function (element, index) { renderer.setElementColor(element.id, signedColor(values[index], maximum)); });
      return textLegend((mode === 'axialForce' ? 'Axial force' : 'Axial stress') + ' · max |value| ' + withUnit(maximum, unit(dataset, mode === 'axialForce' ? 'force' : 'stress')) + ' · red tension, blue compression');
    }
    var maximumDiagram = 0;
    dataset.raw.model.elements.forEach(function (element) { var q = dataset.resultsByElementId.get(element.id).localEndForces; maximumDiagram = Math.max(maximumDiagram, Math.abs(S.frameDiagram(q, mode, 0)), Math.abs(S.frameDiagram(q, mode, 1))); });
    var diagramScale = maximumDiagram ? diagonal(dataset) * .14 / maximumDiagram : 0;
    dataset.raw.model.elements.forEach(function (element) {
      var first = dataset.nodesById.get(element.nodeIds[0]), second = dataset.nodesById.get(element.nodeIds[1]), normal = elementNormal(dataset, element), q = dataset.resultsByElementId.get(element.id).localEndForces;
      var v1 = S.frameDiagram(q, mode, 0), v2 = S.frameDiagram(q, mode, 1);
      var d1 = { x: first.x + normal.x * v1 * diagramScale, y: first.y + normal.y * v1 * diagramScale }, d2 = { x: second.x + normal.x * v2 * diagramScale, y: second.y + normal.y * v2 * diagramScale };
      var points = [first, d1, d2, second].map(M.geometry.svgPoint);
      var polygon = appendPath(renderer.layers.diagrams, { class: 'diagram-fill', 'data-element-id': element.id, d: 'M ' + points.map(function (p) { return p.x + ' ' + p.y; }).join(' L ') + ' Z', tabindex: '0' }, mode + ' diagram for element ' + element.id);
      appendPath(renderer.layers.diagrams, { class: 'diagram-line', 'data-element-id': element.id, d: M.geometry.pathData([d1, d2]) }, mode + ' diagram edge for element ' + element.id);
    });
    return textLegend(mode + ' diagram · max |value| ' + withUnit(maximumDiagram, mode === 'M' ? unit(dataset, 'moment') : unit(dataset, 'force')) + ' · display scale ' + format(diagramScale));
  };
  S.format = format; S.withUnit = withUnit; S.viridis = viridis;
}(window.MKEFPost));
