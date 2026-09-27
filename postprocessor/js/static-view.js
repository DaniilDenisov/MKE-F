(function (M) {
  'use strict';
  var S = M.staticResults = {};

  function format(value) {
    if (value === 0) return '0';
    var absolute = Math.abs(value);
    return absolute >= 10000 || absolute < 0.001 ? value.toExponential(3) : Number(value.toPrecision(5)).toString();
  }
  function unit(dataset, name) { return dataset.raw.metadata.units[name] || ''; }
  function withUnit(value, suffix) { return format(value) + (suffix ? ' ' + suffix : ''); }
  function diagonal(dataset) {
    var bounds = M.geometry.modelBounds(dataset.raw.model);
    return Math.max(Math.hypot(bounds.maxX - bounds.minX, bounds.maxY - bounds.minY), 1e-9);
  }
  function elementNormal(dataset, element) {
    var first = dataset.nodesById.get(element.nodeIds[0]), second = dataset.nodesById.get(element.nodeIds[1]);
    var basis = M.geometry.elementBasis(first, second);
    return { x: -basis.s, y: basis.c };
  }
  function appendText(layer, point, text, className, anchor) {
    var svgPoint = M.geometry.svgPoint(point);
    var label = M.svgElement('text', { x: svgPoint.x, y: svgPoint.y, class: className || 'result-label', 'text-anchor': anchor || 'start' });
    label.textContent = text; layer.appendChild(label); return label;
  }
  function appendLine(layer, first, second, className) {
    var a = M.geometry.svgPoint(first), b = M.geometry.svgPoint(second);
    var line = M.svgElement('path', { d: 'M ' + a.x + ' ' + a.y + ' L ' + b.x + ' ' + b.y, class: className });
    layer.appendChild(line); return line;
  }

  S.frameDiagram = function (localEndForces, quantity, xi) {
    if (quantity === 'N') return -localEndForces[0];
    if (quantity === 'V') return -localEndForces[1];
    if (quantity === 'M') return (1 - xi) * (-localEndForces[2]) + xi * localEndForces[5];
    throw new Error('Unsupported frame diagram ' + quantity + '.');
  };

  function renderSupport(layer, dataset, support, size) {
    var node = dataset.nodesById.get(support.nodeId), p = M.geometry.svgPoint(node), path;
    if (support.type === 1) {
      path = 'M ' + (p.x - size) + ' ' + (p.y + size * 0.15) + ' L ' + (p.x + size) + ' ' + (p.y + size * 0.15) +
        ' M ' + (p.x - size) + ' ' + (p.y + size * 0.35) + ' L ' + (p.x - size * .55) + ' ' + (p.y + size * .75) +
        ' M ' + (p.x - size * .3) + ' ' + (p.y + size * .35) + ' L ' + (p.x + size * .15) + ' ' + (p.y + size * .75) +
        ' M ' + (p.x + size * .4) + ' ' + (p.y + size * .35) + ' L ' + (p.x + size * .85) + ' ' + (p.y + size * .75);
    } else {
      path = 'M ' + p.x + ' ' + p.y + ' L ' + (p.x - size) + ' ' + (p.y + size) + ' L ' + (p.x + size) + ' ' + (p.y + size) + ' Z';
      if (support.type === 2) path += ' M ' + (p.x - size * .7) + ' ' + (p.y + size * 1.25) + ' L ' + (p.x + size * .7) + ' ' + (p.y + size * 1.25);
      if (support.type === 3) path += ' M ' + (p.x - size * 1.25) + ' ' + (p.y - size * .7) + ' L ' + (p.x - size * 1.25) + ' ' + (p.y + size * .7);
      if (support.type === 4) path += ' M ' + (p.x - size * 1.2) + ' ' + (p.y + size * 1.2) + ' L ' + (p.x + size * 1.2) + ' ' + (p.y + size * 1.2);
    }
    var symbol = M.svgElement('path', { d: path, class: 'support-symbol', 'data-support-node-id': support.nodeId });
    symbol.appendChild(M.svgElement('title')); symbol.firstChild.textContent = 'Support type ' + support.type + ' at node ' + support.nodeId;
    layer.appendChild(symbol);
  }

  function arrowPath(node, vector, length) {
    var magnitude = Math.hypot(vector.x, vector.y), ux = vector.x / magnitude, uy = vector.y / magnitude;
    var tail = { x: node.x - ux * length, y: node.y - uy * length };
    var normal = { x: -uy, y: ux }, head = length * .22;
    var a = M.geometry.svgPoint(tail), b = M.geometry.svgPoint(node);
    var h1 = M.geometry.svgPoint({ x: node.x - ux * head + normal.x * head * .45, y: node.y - uy * head + normal.y * head * .45 });
    var h2 = M.geometry.svgPoint({ x: node.x - ux * head - normal.x * head * .45, y: node.y - uy * head - normal.y * head * .45 });
    return 'M ' + a.x + ' ' + a.y + ' L ' + b.x + ' ' + b.y + ' M ' + h1.x + ' ' + h1.y + ' L ' + b.x + ' ' + b.y + ' L ' + h2.x + ' ' + h2.y;
  }

  function renderVector(layer, node, vector, maximum, size, className, text, dataName) {
    var magnitude = Math.hypot(vector.x, vector.y);
    if (magnitude === 0) return;
    var length = size * (0.35 + 0.65 * magnitude / maximum);
    var path = M.svgElement('path', { d: arrowPath(node, vector, length), class: className });
    path.setAttribute(dataName, node.id); path.appendChild(M.svgElement('title')); path.firstChild.textContent = text;
    layer.appendChild(path);
    appendText(layer, { x: node.x - vector.x / magnitude * length, y: node.y - vector.y / magnitude * length }, text, 'result-label', 'middle');
  }

  function renderMoment(layer, node, value, maximum, size, className, text, dataName) {
    if (value === 0) return;
    var radius = size * (.32 + .3 * Math.abs(value) / maximum), p = M.geometry.svgPoint(node), sweep = value > 0 ? 0 : 1;
    var path = 'M ' + (p.x + radius) + ' ' + p.y + ' A ' + radius + ' ' + radius + ' 0 1 ' + sweep + ' ' + (p.x - radius * .7) + ' ' + (p.y - radius * .7) +
      ' M ' + (p.x - radius * .7) + ' ' + (p.y - radius * .7) + ' l ' + (value > 0 ? radius * .05 : radius * .35) + ' ' + radius * .04 +
      ' M ' + (p.x - radius * .7) + ' ' + (p.y - radius * .7) + ' l ' + radius * .04 + ' ' + (value > 0 ? radius * .35 : radius * .05);
    var symbol = M.svgElement('path', { d: path, class: className });
    symbol.setAttribute(dataName, node.id); symbol.appendChild(M.svgElement('title')); symbol.firstChild.textContent = text; layer.appendChild(symbol);
    appendText(layer, { x: node.x + radius, y: node.y + radius }, text, 'result-label');
  }

  function restrainedDOFs(dataset) {
    var restrained = new Set(), model = dataset.raw.model;
    model.supports.forEach(function (support) {
      var row = model.nodes.findIndex(function (node) { return node.id === support.nodeId; });
      M.geometry.restrainedLocalDOFs(support.type, model.dofPerNode).forEach(function (local) { restrained.add(model.dofMap[row][local]); });
    });
    return restrained;
  }

  function nodeComponents(dataset, vector, onlyRestrained) {
    var model = dataset.raw.model, restrained = onlyRestrained ? restrainedDOFs(dataset) : null;
    return model.nodes.map(function (node, row) {
      var dofs = model.dofMap[row], value = function (local) { return !restrained || restrained.has(dofs[local]) ? vector[dofs[local] - 1] : 0; };
      return { node: node, x: value(0), y: value(1), moment: dofs.length > 2 ? value(2) : 0 };
    });
  }

  function renderLoadsAndReactions(renderer, dataset, settings, resultInfo) {
    var forceUnit = unit(dataset, 'force'), momentUnit = unit(dataset, 'moment');
    var sets = [
      { show: settings.showLoads, layer: renderer.layers.loads, values: nodeComponents(dataset, dataset.raw.analysis.loadVector, false), className: 'load-symbol', prefix: 'Load', dataName: 'data-load-node-id' },
      { show: settings.showReactions, layer: renderer.layers.reactions, values: nodeComponents(dataset, dataset.raw.analysis.reactions, true), className: 'reaction-symbol', prefix: 'Reaction', dataName: 'data-reaction-node-id' }
    ];
    var forceMaximum = 0, momentMaximum = 0;
    sets.forEach(function (set) { set.values.forEach(function (item) { forceMaximum = Math.max(forceMaximum, Math.hypot(item.x, item.y)); momentMaximum = Math.max(momentMaximum, Math.abs(item.moment)); }); });
    var size = diagonal(dataset) * .13;
    sets.forEach(function (set) {
      set.layer.style.display = set.show ? '' : 'none';
      if (!set.show) return;
      set.values.forEach(function (item) {
        renderVector(set.layer, item.node, { x: item.x, y: item.y }, forceMaximum || 1, size, set.className,
          set.prefix + ' (' + withUnit(item.x, forceUnit) + ', ' + withUnit(item.y, forceUnit) + ')', set.dataName);
        renderMoment(set.layer, item.node, item.moment, momentMaximum || 1, size, set.className,
          set.prefix + ' Mz ' + withUnit(item.moment, momentUnit), set.dataName);
      });
    });
    resultInfo.forceMaximum = forceMaximum; resultInfo.momentMaximum = momentMaximum;
    resultInfo.forceDisplayLength = forceMaximum ? size : 0;
    resultInfo.momentDisplayRadius = momentMaximum ? size * .62 : 0;
  }

  function signedColor(value, maximum) {
    if (!maximum || value === 0) return '#66788a';
    var amount = Math.min(Math.abs(value) / maximum, 1), light = Math.round(68 - amount * 25);
    return value > 0 ? 'hsl(4 58% ' + light + '%)' : 'hsl(218 55% ' + light + '%)';
  }

  function renderTrussResults(renderer, dataset, settings) {
    if (settings.trussResult === 'none') return;
    var values = [];
    dataset.raw.model.elements.forEach(function (element) {
      if (element.type === 112) values.push(dataset.resultsByElementId.get(element.id)[settings.trussResult]);
    });
    if (!values.length) throw new Error('Truss force and stress coloring requires element type 112.');
    var maximum = Math.max.apply(null, values.map(Math.abs).concat([0]));
    dataset.raw.model.elements.forEach(function (element) {
      if (element.type !== 112) return;
      var result = dataset.resultsByElementId.get(element.id), value = result[settings.trussResult];
      var path = renderer.layers['deformed-geometry'].querySelector('[data-element-id="' + element.id + '"]');
      path.style.stroke = signedColor(value, maximum);
      var first = dataset.nodesById.get(element.nodeIds[0]), second = dataset.nodesById.get(element.nodeIds[1]);
      var suffix = settings.trussResult === 'axialStress' ? unit(dataset, 'stress') : unit(dataset, 'force');
      appendText(renderer.layers.labels, { x: (first.x + second.x) / 2, y: (first.y + second.y) / 2 }, withUnit(value, suffix), 'result-label', 'middle');
    });
  }

  function renderDiagrams(renderer, dataset, settings, resultInfo) {
    var quantity = settings.diagram;
    if (quantity === 'none') { renderer.layers.diagrams.style.display = 'none'; return; }
    var frames = dataset.raw.model.elements.filter(function (element) { return element.type === 113; });
    if (!frames.length) throw new Error('Frame diagrams require element type 113.');
    var maximum = 0;
    frames.forEach(function (element) {
      var q = dataset.resultsByElementId.get(element.id).localEndForces;
      maximum = Math.max(maximum, Math.abs(S.frameDiagram(q, quantity, 0)), Math.abs(S.frameDiagram(q, quantity, 1)));
    });
    var diagramScale = maximum ? diagonal(dataset) * .14 / maximum : 0;
    frames.forEach(function (element) {
      var first = dataset.nodesById.get(element.nodeIds[0]), second = dataset.nodesById.get(element.nodeIds[1]);
      var normal = elementNormal(dataset, element), q = dataset.resultsByElementId.get(element.id).localEndForces;
      var v1 = S.frameDiagram(q, quantity, 0), v2 = S.frameDiagram(q, quantity, 1);
      var d1 = { x: first.x + normal.x * v1 * diagramScale, y: first.y + normal.y * v1 * diagramScale };
      var d2 = { x: second.x + normal.x * v2 * diagramScale, y: second.y + normal.y * v2 * diagramScale };
      var points = [first, d1, d2, second].map(M.geometry.svgPoint);
      var polygon = M.svgElement('path', { class: 'diagram-fill', 'data-diagram-element-id': element.id,
        d: 'M ' + points.map(function (p) { return p.x + ' ' + p.y; }).join(' L ') + ' Z' });
      renderer.layers.diagrams.appendChild(polygon); appendLine(renderer.layers.diagrams, d1, d2, 'diagram-line');
      var suffix = quantity === 'M' ? unit(dataset, 'moment') : unit(dataset, 'force');
      appendText(renderer.layers.diagrams, d1, quantity + '1=' + withUnit(v1, suffix), 'result-label', 'middle');
      appendText(renderer.layers.diagrams, d2, quantity + '2=' + withUnit(v2, suffix), 'result-label', 'middle');
    });
    renderer.layers.diagrams.style.display = '';
    resultInfo.diagramScale = diagramScale; resultInfo.diagramMaximum = maximum;
  }

  function renderLegend(renderer, dataset, settings, info) {
    var bounds = M.geometry.modelBounds(dataset.raw.model), size = diagonal(dataset), x = bounds.minX - size * .13, y = bounds.maxY + size * .14;
    var lines = ['Static result', 'Deformation × ' + format(info.deformationScale)];
    if (info.forceMaximum) lines.push('Force: ' + withUnit(info.forceMaximum, unit(dataset, 'force')) + ' → ' + withUnit(info.forceDisplayLength, unit(dataset, 'length')));
    if (info.momentMaximum) lines.push('Moment: ' + withUnit(info.momentMaximum, unit(dataset, 'moment')) + ' → r ' + withUnit(info.momentDisplayRadius, unit(dataset, 'length')));
    if (settings.diagram !== 'none') lines.push(settings.diagram + ' diagram scale × ' + format(info.diagramScale || 0));
    if (settings.trussResult !== 'none') lines.push('Blue compression · red tension');
    var lineHeight = size * .035, width = size * .42, height = lineHeight * (lines.length + 1);
    var origin = M.geometry.svgPoint({ x: x, y: y });
    renderer.layers.labels.appendChild(M.svgElement('rect', { class: 'legend-background', x: origin.x, y: origin.y, width: width, height: height }));
    lines.forEach(function (line, index) {
      var label = M.svgElement('text', { class: 'legend-text', x: origin.x + lineHeight * .35, y: origin.y + lineHeight * (index + .85) });
      label.textContent = line; renderer.layers.labels.appendChild(label);
    });
  }

  S.render = function (renderer, dataset, settings, geometryInfo) {
    var info = { deformationScale: geometryInfo.scale, zeroDeformation: geometryInfo.zero };
    renderer.layers.supports.style.display = settings.showSupports ? '' : 'none';
    if (settings.showSupports) {
      var supportSize = diagonal(dataset) * .025;
      dataset.raw.model.supports.forEach(function (support) { renderSupport(renderer.layers.supports, dataset, support, supportSize); });
    }
    renderLoadsAndReactions(renderer, dataset, settings, info);
    renderTrussResults(renderer, dataset, settings);
    renderDiagrams(renderer, dataset, settings, info);
    renderLegend(renderer, dataset, settings, info);
    renderer.layers.labels.style.display = '';
    var textSize = diagonal(dataset) * .023;
    renderer.svg.querySelectorAll('.result-label, .legend-text').forEach(function (label) { label.setAttribute('font-size', textSize); });
    return info;
  };
}(window.MKEFPost));
