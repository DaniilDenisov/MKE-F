(function (M) {
  'use strict';

  var layerNames = ['original-geometry', 'diagrams', 'deformed-geometry', 'supports', 'loads', 'reactions', 'nodes', 'labels', 'selection-overlay'];
  var rendererSequence = 0;

  function clearElement(element) {
    while (element.firstChild) element.removeChild(element.firstChild);
  }

  function characteristicSize(model) {
    var bounds = M.geometry.modelBounds(model);
    return Math.max(Math.hypot(bounds.maxX - bounds.minX, bounds.maxY - bounds.minY), 1e-9);
  }

  function Renderer(svg, onSelection) {
    this.svg = svg;
    this.svg.classList.add('structural-viewport');
    this.onSelection = typeof onSelection === 'function' ? onSelection : function () {};
    this.dataset = null;
    this.viewBox = null;
    this.drag = null;
    this.layers = {};
    this.originalByElement = new Map();
    this.deformedByElement = new Map();
    this.nodeById = new Map();
    this.deformedNodeHitById = new Map();
    this.nodeLabels = new Map();
    this.elementLabels = new Map();
    this.selection = null;
    this.size = 1;
    this.gradientPrefix = 'mkef-result-' + (++rendererSequence) + '-';
    this.createLayers();
    this.bindNavigation();
  }

  Renderer.prototype.createLayers = function () {
    clearElement(this.svg);
    this.resultDefs = M.svgElement('defs', { 'data-result-defs': '' });
    this.svg.appendChild(this.resultDefs);
    this.svg.appendChild(M.svgElement('rect', { class: 'viewport-background', x: '-1000000', y: '-1000000', width: '2000000', height: '2000000' }));
    layerNames.forEach(function (name) {
      var group = M.svgElement('g', { 'data-layer': name });
      this.svg.appendChild(group);
      this.layers[name] = group;
    }, this);
  };

  Renderer.prototype.clearLayer = function (name) { clearElement(this.layers[name]); };

  Renderer.prototype.mount = function (dataset) {
    this.dataset = dataset;
    this.selection = null;
    this.originalByElement.clear();
    this.deformedByElement.clear();
    this.nodeById.clear();
    this.deformedNodeHitById.clear();
    this.nodeLabels.clear();
    this.elementLabels.clear();
    clearElement(this.resultDefs);
    layerNames.forEach(function (name) { this.clearLayer(name); }, this);
    this.size = characteristicSize(dataset.raw.model);
    var labelOffset = this.size * .018;
    var fontSize = this.size * .023;

    dataset.raw.model.elements.forEach(function (element) {
      var first = dataset.nodesById.get(element.nodeIds[0]);
      var second = dataset.nodesById.get(element.nodeIds[1]);
      var original = M.svgElement('path', { class: 'original-element', d: M.geometry.pathData([first, second]), 'data-element-id': element.id, tabindex: '0' });
      var originalTitle = M.svgElement('title');
      originalTitle.textContent = 'Element ' + element.id + ' · type ' + element.type;
      original.appendChild(originalTitle);
      this.layers['original-geometry'].appendChild(original);
      this.originalByElement.set(element.id, original);

      var deformed = M.svgElement('path', { class: 'deformed-element', d: M.geometry.pathData([first, second]), 'data-element-id': element.id, tabindex: '0' });
      var deformedTitle = M.svgElement('title');
      deformed.appendChild(deformedTitle);
      this.layers['deformed-geometry'].appendChild(deformed);
      this.deformedByElement.set(element.id, deformed);

      var middle = M.geometry.svgPoint({ x: (first.x + second.x) / 2, y: (first.y + second.y) / 2 });
      var elementLabel = M.svgElement('text', { class: 'label element-label', x: middle.x, y: middle.y - labelOffset, 'font-size': fontSize, 'stroke-width': fontSize * .18, 'text-anchor': 'middle', 'data-element-id': element.id });
      elementLabel.textContent = 'E' + element.id;
      this.layers.labels.appendChild(elementLabel);
      this.elementLabels.set(element.id, elementLabel);
    }, this);

    dataset.raw.model.nodes.forEach(function (node) {
      var point = M.geometry.svgPoint(node);
      var deformedHit = M.svgElement('circle', { class: 'node-hit deformed-node-hit', cx: point.x, cy: point.y, r: this.size * .008, 'data-node-id': node.id, 'aria-hidden': 'true' });
      this.layers['deformed-geometry'].appendChild(deformedHit);
      this.deformedNodeHitById.set(node.id, deformedHit);
      var hit = M.svgElement('circle', { class: 'node-hit', cx: point.x, cy: point.y, r: this.size * .008, 'data-node-id': node.id, 'aria-hidden': 'true' });
      this.layers.nodes.appendChild(hit);
      var circle = M.svgElement('circle', { class: 'node', cx: point.x, cy: point.y, r: this.size * .008, 'data-node-id': node.id, tabindex: '0' });
      var title = M.svgElement('title');
      title.textContent = 'Node ' + node.id + ' · (' + node.x + ', ' + node.y + ')';
      circle.appendChild(title);
      this.layers.nodes.appendChild(circle);
      this.nodeById.set(node.id, circle);
      var label = M.svgElement('text', { class: 'label node-label', x: point.x + labelOffset, y: point.y - labelOffset, 'font-size': fontSize, 'stroke-width': fontSize * .18, 'data-node-id': node.id });
      label.textContent = 'N' + node.id;
      this.layers.labels.appendChild(label);
      this.nodeLabels.set(node.id, label);
    }, this);

    dataset.raw.model.supports.forEach(function (support) {
      var point = M.geometry.svgPoint(dataset.nodesById.get(support.nodeId));
      var supportName = window.MKEFSupportMarkers.label(support.type, dataset.raw.model.dofPerNode);
      window.MKEFSupportMarkers.append(M.svgElement, this.layers.supports, point.x, point.y, dataset.raw.model.dofPerNode === 2 && support.type !== 1 && support.type !== 4 ? (support.type === 2 || support.type === 6 ? 6 : 5) : support.type, this.size * .035, { 'data-node-id': support.nodeId, 'data-support-node-id': support.nodeId }, 'Support ' + support.type + ' · ' + supportName + ' at node ' + support.nodeId);
    }, this);
    this.setVisibility({ showOriginal: true, showDeformed: true, showNodes: true, showNodeLabels: false, showElementLabels: false, showSupports: true, showLoads: true, showReactions: true });
    this.fit();
  };

  Renderer.prototype.updateDeformation = function (displacement, scale, samples) {
    if (!this.dataset) return;
    this.dataset.raw.model.elements.forEach(function (element) {
      var points = M.geometry.sampleElement(this.dataset, element, displacement, scale, samples);
      var path = this.deformedByElement.get(element.id);
      path.setAttribute('d', M.geometry.pathData(points));
      path.firstChild.textContent = 'Element ' + element.id + ' · deformation scale ' + Number(scale).toPrecision(5);
    }, this);
    this.dataset.raw.model.nodes.forEach(function (node) {
      var value = M.geometry.nodeDisplacement(this.dataset, node.id, displacement);
      var point = M.geometry.svgPoint({ x: node.x + scale * value.x, y: node.y + scale * value.y });
      var hit = this.deformedNodeHitById.get(node.id);
      hit.setAttribute('cx', point.x); hit.setAttribute('cy', point.y);
    }, this);
    this.reapplySelection();
  };

  Renderer.prototype.resetElementColors = function () {
    clearElement(this.resultDefs);
    this.deformedByElement.forEach(function (path) { path.style.removeProperty('stroke'); });
  };
  Renderer.prototype.setElementColor = function (elementId, color) {
    var path = this.deformedByElement.get(elementId);
    if (path) path.style.stroke = color;
  };
  Renderer.prototype.setElementGradient = function (elementId, first, second, stops) {
    var path = this.deformedByElement.get(elementId);
    if (!path) return;
    var id = this.gradientPrefix + elementId;
    var a = M.geometry.svgPoint(first), b = M.geometry.svgPoint(second);
    var gradient = M.svgElement('linearGradient', { id: id, gradientUnits: 'userSpaceOnUse', x1: a.x, y1: a.y, x2: b.x, y2: b.y });
    stops.forEach(function (stop) {
      gradient.appendChild(M.svgElement('stop', { offset: (stop.offset * 100) + '%', 'stop-color': stop.color }));
    });
    this.resultDefs.appendChild(gradient);
    path.style.stroke = 'url(#' + id + ')';
  };

  Renderer.prototype.setVisibility = function (settings) {
    function show(layer, visible) { layer.style.display = visible ? '' : 'none'; }
    show(this.layers['original-geometry'], settings.showOriginal);
    show(this.layers['deformed-geometry'], settings.showDeformed);
    show(this.layers.nodes, settings.showNodes);
    show(this.layers.supports, settings.showSupports);
    show(this.layers.loads, settings.showLoads);
    show(this.layers.reactions, settings.showReactions);
    this.nodeLabels.forEach(function (label) { label.style.display = settings.showNodeLabels ? '' : 'none'; });
    this.elementLabels.forEach(function (label) { label.style.display = settings.showElementLabels ? '' : 'none'; });
    show(this.layers.labels, settings.showNodeLabels || settings.showElementLabels);
  };

  Renderer.prototype.fit = function () {
    if (!this.dataset) return;
    var bounds = M.geometry.modelBounds(this.dataset.raw.model);
    var rawWidth = bounds.maxX - bounds.minX;
    var rawHeight = bounds.maxY - bounds.minY;
    var minimumSpan = this.size * .2;
    var width = Math.max(rawWidth, minimumSpan);
    var height = Math.max(rawHeight, minimumSpan);
    var centerX = (bounds.minX + bounds.maxX) / 2;
    var centerY = -(bounds.minY + bounds.maxY) / 2;
    var padding = this.size * .22;
    this.viewBox = { x: centerX - width / 2 - padding, y: centerY - height / 2 - padding, width: width + 2 * padding, height: height + 2 * padding };
    this.updateViewBox();
  };

  Renderer.prototype.updateViewBox = function () {
    if (this.viewBox) this.svg.setAttribute('viewBox', [this.viewBox.x, this.viewBox.y, this.viewBox.width, this.viewBox.height].join(' '));
  };

  Renderer.prototype.pointerInViewBox = function (event) {
    var matrix = this.svg.getScreenCTM();
    if (matrix && typeof DOMPoint === 'function') return new DOMPoint(event.clientX, event.clientY).matrixTransform(matrix.inverse());
    var rect = this.svg.getBoundingClientRect();
    return { x: this.viewBox.x + (event.clientX - rect.left) / rect.width * this.viewBox.width, y: this.viewBox.y + (event.clientY - rect.top) / rect.height * this.viewBox.height };
  };

  Renderer.prototype.nearestNodeId = function (point, deformed) {
    var bestId = null, bestDistance = Infinity;
    this.dataset.raw.model.nodes.forEach(function (node) {
      var candidate;
      if (deformed) {
        var hit = this.deformedNodeHitById.get(node.id);
        candidate = { x: Number(hit.getAttribute('cx')), y: Number(hit.getAttribute('cy')) };
      } else candidate = M.geometry.svgPoint(node);
      var distance = Math.pow(candidate.x - point.x, 2) + Math.pow(candidate.y - point.y, 2);
      if (distance < bestDistance) { bestDistance = distance; bestId = node.id; }
    }, this);
    return bestId;
  };

  Renderer.prototype.clickedNodeId = function (event, target) {
    var fallback = Number(target.getAttribute('data-node-id'));
    if (!target.classList.contains('node') && !target.classList.contains('node-hit')) return fallback;
    var rect = this.svg.getBoundingClientRect();
    var inside = event.clientX >= rect.left && event.clientX <= rect.right && event.clientY >= rect.top && event.clientY <= rect.bottom;
    if (!inside) return fallback;
    return this.nearestNodeId(this.pointerInViewBox(event), target.classList.contains('deformed-node-hit'));
  };

  Renderer.prototype.bindNavigation = function () {
    this.svg.addEventListener('wheel', function (event) {
      if (!this.viewBox) return;
      event.preventDefault();
      var anchor = this.pointerInViewBox(event);
      var factor = M.geometry.clamp(Math.exp(event.deltaY * .001), .5, 2);
      this.viewBox.x = anchor.x - (anchor.x - this.viewBox.x) * factor;
      this.viewBox.y = anchor.y - (anchor.y - this.viewBox.y) * factor;
      this.viewBox.width *= factor;
      this.viewBox.height *= factor;
      this.updateViewBox();
    }.bind(this), { passive: false });
    this.svg.addEventListener('pointerdown', function (event) {
      if (event.button !== 0 || !this.viewBox) return;
      if (event.target.closest('[data-node-id], [data-element-id]')) return;
      this.drag = { pointer: event.pointerId, start: this.pointerInViewBox(event), x: this.viewBox.x, y: this.viewBox.y, moved: false };
      this.svg.setPointerCapture(event.pointerId);
    }.bind(this));
    this.svg.addEventListener('pointermove', function (event) {
      if (!this.drag || this.drag.pointer !== event.pointerId) return;
      var point = this.pointerInViewBox(event);
      this.drag.moved = this.drag.moved || Math.abs(point.x - this.drag.start.x) + Math.abs(point.y - this.drag.start.y) > this.viewBox.width * .002;
      this.viewBox.x = this.drag.x - (point.x - this.drag.start.x);
      this.viewBox.y = this.drag.y - (point.y - this.drag.start.y);
      this.updateViewBox();
    }.bind(this));
    this.svg.addEventListener('pointerup', function () { this.drag = null; }.bind(this));
    this.svg.addEventListener('click', function (event) {
      var target = event.target.closest('[data-node-id], [data-element-id]');
      if (!target) this.setSelection(null);
      else if (target.hasAttribute('data-node-id')) this.setSelection({ kind: 'node', id: this.clickedNodeId(event, target) });
      else this.setSelection({ kind: 'element', id: Number(target.getAttribute('data-element-id')) });
    }.bind(this));
    this.svg.addEventListener('keydown', function (event) {
      if (event.key !== 'Enter' && event.key !== ' ') return;
      var target = event.target.closest('[data-node-id], [data-element-id]');
      if (!target) return;
      event.preventDefault();
      if (target.hasAttribute('data-node-id')) this.setSelection({ kind: 'node', id: Number(target.getAttribute('data-node-id')) });
      else this.setSelection({ kind: 'element', id: Number(target.getAttribute('data-element-id')) });
    }.bind(this));
  };

  Renderer.prototype.reapplySelection = function () {
    this.svg.querySelectorAll('.selected').forEach(function (item) { item.classList.remove('selected'); });
    if (!this.selection) return;
    if (this.selection.kind === 'node') {
      var node = this.nodeById.get(this.selection.id);
      if (node) node.classList.add('selected');
      var deformedHit = this.deformedNodeHitById.get(this.selection.id);
      if (deformedHit) deformedHit.classList.add('selected');
      return;
    }
    var attribute = 'data-element-id';
    this.svg.querySelectorAll('[' + attribute + '="' + this.selection.id + '"]').forEach(function (item) {
      if (item.tagName.toLowerCase() !== 'text' && !item.classList.contains('node-hit')) item.classList.add('selected');
    });
  };

  Renderer.prototype.setSelection = function (selection) {
    this.selection = selection;
    this.reapplySelection();
    this.onSelection(selection);
  };

  Renderer.prototype.sceneNodeCount = function () { return this.svg.querySelectorAll('*').length; };
  M.Renderer = Renderer;
}(window.MKEFPost));
