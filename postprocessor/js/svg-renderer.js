(function (M) {
  'use strict';

  function Renderer(svg, selectionOutput) {
    this.svg = svg;
    this.selectionOutput = selectionOutput;
    this.dataset = null;
    this.viewBox = null;
    this.drag = null;
    this.layers = {};
    this.createLayers();
    this.bindNavigation();
  }

  Renderer.prototype.createLayers = function () {
    while (this.svg.firstChild) this.svg.removeChild(this.svg.firstChild);
    var background = M.svgElement('rect', { class: 'viewport-background', x: '-1000000', y: '-1000000', width: '2000000', height: '2000000' });
    this.svg.appendChild(background);
    ['original-geometry', 'diagrams', 'deformed-geometry', 'supports', 'loads', 'reactions', 'nodes', 'labels', 'selection-overlay'].forEach(function (name) {
      var group = M.svgElement('g', { 'data-layer': name });
      this.svg.appendChild(group); this.layers[name] = group;
    }, this);
  };

  Renderer.prototype.clear = function () {
    Object.keys(this.layers).forEach(function (name) {
      var layer = this.layers[name]; while (layer.firstChild) layer.removeChild(layer.firstChild);
    }, this);
  };

  Renderer.prototype.render = function (dataset, displacement, settings) {
    this.dataset = dataset; this.clear();
    var model = dataset.raw.model;
    var scaleInfo = settings.scaleMode === 'auto' ? M.geometry.automaticScale(dataset, displacement, settings.samples) : { scale: settings.manualScale, zero: displacement.every(function (v) { return v === 0; }) };
    var scale = scaleInfo.scale;
    model.elements.forEach(function (element) {
      var first = dataset.nodesById.get(element.nodeIds[0]), second = dataset.nodesById.get(element.nodeIds[1]);
      var original = M.svgElement('path', { class: 'original-element', d: M.geometry.pathData([first, second]), 'data-element-id': element.id, tabindex: '0' });
      original.appendChild(M.svgElement('title'));
      original.firstChild.textContent = 'Element ' + element.id + ' · type ' + element.type;
      this.layers['original-geometry'].appendChild(original);
      var points = M.geometry.sampleElement(dataset, element, displacement, scale, settings.samples);
      var deformed = M.svgElement('path', { class: 'deformed-element', d: M.geometry.pathData(points), 'data-element-id': element.id, tabindex: '0' });
      deformed.appendChild(M.svgElement('title'));
      deformed.firstChild.textContent = 'Element ' + element.id + ' · deformation scale ' + scale.toPrecision(5);
      this.layers['deformed-geometry'].appendChild(deformed);
      if (settings.showElementLabels) {
        var middle = M.geometry.svgPoint(points[Math.floor(points.length / 2)]);
        var label = M.svgElement('text', { class: 'label', x: middle.x, y: middle.y - 7, 'text-anchor': 'middle' });
        label.textContent = 'E' + element.id; this.layers.labels.appendChild(label);
      }
    }, this);
    model.nodes.forEach(function (node) {
      var point = M.geometry.svgPoint(node);
      var circle = M.svgElement('circle', { class: 'node', cx: point.x, cy: point.y, r: 3.5, 'data-node-id': node.id, tabindex: '0' });
      circle.appendChild(M.svgElement('title'));
      circle.firstChild.textContent = 'Node ' + node.id + ' · (' + node.x + ', ' + node.y + ')';
      this.layers.nodes.appendChild(circle);
      if (settings.showNodeLabels) {
        var label = M.svgElement('text', { class: 'label', x: point.x + 7, y: point.y - 7 });
        label.textContent = 'N' + node.id; this.layers.labels.appendChild(label);
      }
    }, this);
    this.layers['original-geometry'].style.display = settings.showOriginal ? '' : 'none';
    this.layers['deformed-geometry'].style.display = settings.showDeformed ? '' : 'none';
    this.layers.nodes.style.display = settings.showNodes ? '' : 'none';
    this.layers.labels.style.display = (settings.showNodeLabels || settings.showElementLabels) ? '' : 'none';
    this.fit();
    return { scale: scale, zero: scaleInfo.zero };
  };

  Renderer.prototype.fit = function () {
    if (!this.dataset) return;
    var bounds = M.geometry.modelBounds(this.dataset.raw.model);
    var width = Math.max(bounds.maxX - bounds.minX, 1e-6), height = Math.max(bounds.maxY - bounds.minY, 1e-6);
    var padding = Math.max(width, height) * 0.18 + 0.05;
    this.viewBox = { x: bounds.minX - padding, y: -bounds.maxY - padding, width: width + 2 * padding, height: height + 2 * padding };
    this.updateViewBox();
  };
  Renderer.prototype.updateViewBox = function () {
    var v = this.viewBox; if (v) this.svg.setAttribute('viewBox', [v.x, v.y, v.width, v.height].join(' '));
  };
  Renderer.prototype.pointerInViewBox = function (event) {
    var rect = this.svg.getBoundingClientRect(), v = this.viewBox;
    return { x: v.x + (event.clientX - rect.left) / rect.width * v.width, y: v.y + (event.clientY - rect.top) / rect.height * v.height };
  };
  Renderer.prototype.bindNavigation = function () {
    this.svg.addEventListener('wheel', function (event) {
      if (!this.viewBox) return; event.preventDefault();
      var anchor = this.pointerInViewBox(event), factor = Math.exp(event.deltaY * 0.001);
      factor = M.geometry.clamp(factor, 0.5, 2);
      this.viewBox.x = anchor.x - (anchor.x - this.viewBox.x) * factor;
      this.viewBox.y = anchor.y - (anchor.y - this.viewBox.y) * factor;
      this.viewBox.width *= factor; this.viewBox.height *= factor; this.updateViewBox();
    }.bind(this), { passive: false });
    this.svg.addEventListener('pointerdown', function (event) {
      if (event.button !== 0 || !this.viewBox) return;
      this.drag = { pointer: event.pointerId, start: this.pointerInViewBox(event), x: this.viewBox.x, y: this.viewBox.y };
      this.svg.setPointerCapture(event.pointerId);
    }.bind(this));
    this.svg.addEventListener('pointermove', function (event) {
      if (!this.drag || this.drag.pointer !== event.pointerId) return;
      var point = this.pointerInViewBox(event);
      this.viewBox.x = this.drag.x - (point.x - this.drag.start.x);
      this.viewBox.y = this.drag.y - (point.y - this.drag.start.y); this.updateViewBox();
    }.bind(this));
    this.svg.addEventListener('pointerup', function () { this.drag = null; }.bind(this));
    this.svg.addEventListener('click', this.select.bind(this));
    this.svg.addEventListener('keydown', function (event) { if (event.key === 'Enter' || event.key === ' ') this.select(event); }.bind(this));
  };
  Renderer.prototype.select = function (event) {
    var target = event.target.closest('[data-node-id], [data-element-id]');
    this.svg.querySelectorAll('.selected').forEach(function (item) { item.classList.remove('selected'); });
    if (!target) { this.selectionOutput.textContent = 'Nothing selected'; return; }
    var selector, text;
    if (target.hasAttribute('data-node-id')) {
      var nodeId = target.getAttribute('data-node-id'); selector = '[data-node-id="' + nodeId + '"]'; text = 'Node ' + nodeId;
    } else {
      var elementId = target.getAttribute('data-element-id'); selector = '[data-element-id="' + elementId + '"]'; text = 'Element ' + elementId;
    }
    this.svg.querySelectorAll(selector).forEach(function (item) { item.classList.add('selected'); });
    this.selectionOutput.textContent = text;
  };
  M.Renderer = Renderer;
}(window.MKEFPost));
