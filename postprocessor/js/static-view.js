(function (M) {
  'use strict';
  var S = M.staticResults = {};
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

  S.mount = function (renderer, dataset) {
    renderer.clearLayer('loads'); renderer.clearLayer('reactions');
    var sets = [
      { layer: renderer.layers.loads, values: nodeComponents(dataset, dataset.raw.analysis.loadVector, false), className: 'load-symbol', prefix: 'Load' },
      { layer: renderer.layers.reactions, values: nodeComponents(dataset, dataset.raw.analysis.reactions, true), className: 'reaction-symbol', prefix: 'Reaction' }
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
          appendPath(set.layer, { d: arrowPath(item.node, vector, length), class: set.className, 'data-node-id': item.node.id }, set.prefix + ' (' + withUnit(item.x, unit(dataset, 'force')) + ', ' + withUnit(item.y, unit(dataset, 'force')) + ')');
        }
        if (Math.abs(item.moment) > momentTolerance) appendPath(set.layer, { d: momentPath(item.node, item.moment, momentMaximum || 1, size), class: set.className, 'data-node-id': item.node.id }, set.prefix + ' Mz ' + withUnit(item.moment, unit(dataset, 'moment')));
      });
    });
  };

  S.frameDiagram = function (q, quantity, xi) {
    if (quantity === 'N') return -q[0];
    if (quantity === 'V') return -q[1];
    if (quantity === 'M') return (1 - xi) * (-q[2]) + xi * q[5];
    throw new Error('Unsupported frame diagram ' + quantity + '.');
  };
  function signedColor(value, maximum) { if (!maximum || value === 0) return '#66788a'; var amount = Math.min(Math.abs(value) / maximum, 1), light = Math.round(68 - amount * 25); return value > 0 ? 'hsl(4 58% ' + light + '%)' : 'hsl(218 55% ' + light + '%)'; }

  S.applyResult = function (renderer, dataset, mode) {
    renderer.clearLayer('diagrams'); renderer.resetElementColors(); renderer.layers.diagrams.style.display = mode === 'none' ? 'none' : '';
    if (mode === 'none') return '';
    if (mode === 'axialForce' || mode === 'axialStress') {
      var values = dataset.raw.model.elements.map(function (element) { return dataset.resultsByElementId.get(element.id)[mode]; });
      var maximum = Math.max.apply(null, values.map(Math.abs).concat([0]));
      dataset.raw.model.elements.forEach(function (element, index) { renderer.setElementColor(element.id, signedColor(values[index], maximum)); });
      return (mode === 'axialForce' ? 'Axial force' : 'Axial stress') + ' · max |value| ' + withUnit(maximum, unit(dataset, mode === 'axialForce' ? 'force' : 'stress')) + ' · red tension, blue compression';
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
    return mode + ' diagram · max |value| ' + withUnit(maximumDiagram, mode === 'M' ? unit(dataset, 'moment') : unit(dataset, 'force')) + ' · display scale ' + format(diagramScale);
  };
  S.format = format; S.withUnit = withUnit;
}(window.MKEFPost));
