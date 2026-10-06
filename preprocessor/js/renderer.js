(function (global) {
  'use strict';
  var M = global.MKEFPre;

  function Renderer(svg, callbacks) {
    this.svg = svg;
    this.callbacks = callbacks || {};
    this.model = null;
    this.tool = 'select';
    this.selection = null;
    this.pendingNode = null;
    this.view = { x: -1, y: -1, width: 2, height: 2 };
    this.gridEnabled = true;
    this.gridSpacing = 1;
    this.gridLayer = null;
    this.axisLayer = null;
    this.previewLayer = null;
    this.previewElement = null;
    this.previewPointerPoint = null;
    this.previewPoint = null;
    this.drag = null;
    this.bindNavigation();
    this.bindResize();
  }

  Renderer.prototype.screenPoint = function (event) {
    var rect = this.svg.getBoundingClientRect();
    return {
      x: this.view.x + (event.clientX - rect.left) / Math.max(rect.width, 1) * this.view.width,
      y: -(this.view.y + (event.clientY - rect.top) / Math.max(rect.height, 1) * this.view.height)
    };
  };

  Renderer.prototype.bindNavigation = function () {
    var self = this;
    function clearNativeSelection() {
      var selection = global.getSelection ? global.getSelection() : null;
      if (selection && selection.removeAllRanges) selection.removeAllRanges();
    }
    this.svg.style.userSelect = 'none';
    this.svg.style.webkitUserSelect = 'none';
    this.svg.setAttribute('unselectable', 'on');
    this.svg.addEventListener('selectstart', function (event) { event.preventDefault(); });
    this.svg.addEventListener('dragstart', function (event) { event.preventDefault(); });
    this.svg.addEventListener('mousedown', function (event) { if (event.button === 0) event.preventDefault(); });
    document.addEventListener('selectionchange', function () { if (self.drag) clearNativeSelection(); });
    this.svg.addEventListener('wheel', function (event) {
      event.preventDefault();
      var rect = self.svg.getBoundingClientRect(), factor = event.deltaY < 0 ? 0.86 : 1.16;
      var px = self.view.x + (event.clientX - rect.left) / Math.max(rect.width, 1) * self.view.width;
      var py = self.view.y + (event.clientY - rect.top) / Math.max(rect.height, 1) * self.view.height;
      self.view.x = px - (px - self.view.x) * factor;
      self.view.y = py - (py - self.view.y) * factor;
      self.view.width *= factor; self.view.height *= factor; self.applyView();
    }, { passive: false });
    this.svg.addEventListener('pointerdown', function (event) {
      if (event.button !== 0) return;
      event.preventDefault();
      clearNativeSelection();
      if (event.target !== self.svg) return;
      if (self.tool === 'node') {
        if (self.callbacks.addNode) self.callbacks.addNode(self.screenPoint(event));
        return;
      }
      self.drag = { x: event.clientX, y: event.clientY, viewX: self.view.x, viewY: self.view.y };
      self.svg.setPointerCapture(event.pointerId);
    });
    this.svg.addEventListener('pointermove', function (event) {
      if (self.tool === 'node') self.updateNodePreview(self.screenPoint(event));
      if (!self.drag) return;
      event.preventDefault();
      clearNativeSelection();
      var rect = self.svg.getBoundingClientRect();
      self.view.x = self.drag.viewX - (event.clientX - self.drag.x) / Math.max(rect.width, 1) * self.view.width;
      self.view.y = self.drag.viewY - (event.clientY - self.drag.y) / Math.max(rect.height, 1) * self.view.height;
      self.applyView();
    });
    this.svg.addEventListener('pointerup', function () { clearNativeSelection(); self.drag = null; });
    this.svg.addEventListener('pointercancel', function () { clearNativeSelection(); self.drag = null; });
    this.svg.addEventListener('pointerleave', function () { self.previewPointerPoint = null; self.previewPoint = null; self.updateNodePreview(); });
  };

  Renderer.prototype.bindResize = function () {
    var self = this;
    function synchronize() { self.synchronizeAspectRatio(); }
    if (typeof ResizeObserver !== 'undefined') {
      this.resizeObserver = new ResizeObserver(synchronize);
      this.resizeObserver.observe(this.svg);
    } else global.addEventListener('resize', synchronize);
  };

  Renderer.prototype.synchronizeAspectRatio = function () {
    var rect = this.svg.getBoundingClientRect();
    if (!(rect.width > 0) || !(rect.height > 0)) return;
    var targetAspect = rect.width / rect.height;
    var currentAspect = this.view.width / this.view.height;
    if (Math.abs(targetAspect - currentAspect) <= targetAspect * 1e-6) return;
    var centerX = this.view.x + this.view.width / 2;
    var centerY = this.view.y + this.view.height / 2;
    if (currentAspect < targetAspect) this.view.width = this.view.height * targetAspect;
    else this.view.height = this.view.width / targetAspect;
    this.view.x = centerX - this.view.width / 2;
    this.view.y = centerY - this.view.height / 2;
    this.applyView();
  };

  Renderer.prototype.applyView = function () {
    this.svg.setAttribute('viewBox', [this.view.x, this.view.y, this.view.width, this.view.height].join(' '));
    this.drawCoordinateSystem();
  };

  Renderer.prototype.setGrid = function (enabled, spacing) {
    this.gridEnabled = !!enabled;
    if (Number.isFinite(spacing) && spacing > 0) this.gridSpacing = spacing;
    this.drawCoordinateSystem();
  };

  Renderer.prototype.drawCoordinateSystem = function () {
    if (!this.gridLayer || !this.axisLayer) return;
    while (this.gridLayer.firstChild) this.gridLayer.removeChild(this.gridLayer.firstChild);
    while (this.axisLayer.firstChild) this.axisLayer.removeChild(this.axisLayer.firstChild);
    var left = this.view.x, right = left + this.view.width, top = this.view.y, bottom = top + this.view.height;
    var spacing = this.gridSpacing, displaySpacing = spacing;
    while (Math.max(this.view.width, this.view.height) / displaySpacing > 50) displaySpacing *= 2;
    var epsilon = displaySpacing * 1e-9, index = 0, value;

    if (this.gridEnabled) {
      for (value = Math.ceil((left - epsilon) / displaySpacing) * displaySpacing; value <= right + epsilon; value += displaySpacing) {
        if (Math.abs(value) <= epsilon) continue;
        index = Math.round(value / displaySpacing);
        this.gridLayer.appendChild(M.svgElement('line', { x1: value, y1: top, x2: value, y2: bottom, class: index % 5 === 0 ? 'major' : '' }));
      }
      for (value = Math.ceil((top - epsilon) / displaySpacing) * displaySpacing; value <= bottom + epsilon; value += displaySpacing) {
        if (Math.abs(value) <= epsilon) continue;
        index = Math.round(value / displaySpacing);
        this.gridLayer.appendChild(M.svgElement('line', { x1: left, y1: value, x2: right, y2: value, class: index % 5 === 0 ? 'major' : '' }));
      }
    }

    var fontSize = Math.min(this.view.width, this.view.height) * 0.027;
    var tickOffset = fontSize * 0.85;
    function label(layer, x, y, text, className) {
      var item = M.svgElement('text', { x: x, y: y, 'font-size': fontSize, class: className || '' }); item.textContent = text; layer.appendChild(item);
    }
    function formatted(number) {
      var precision = Math.max(0, Math.min(8, Math.ceil(-Math.log10(displaySpacing)) + 1));
      return Number(number.toFixed(precision)).toString();
    }

    if (top <= 0 && bottom >= 0) {
      this.axisLayer.appendChild(M.svgElement('line', { x1: left, y1: 0, x2: right, y2: 0 }));
      for (value = Math.ceil((left - epsilon) / displaySpacing) * displaySpacing; value <= right + epsilon; value += displaySpacing) {
        if (Math.abs(value) <= epsilon) continue;
        this.axisLayer.appendChild(M.svgElement('line', { x1: value, y1: -tickOffset * 0.2, x2: value, y2: tickOffset * 0.2 }));
        label(this.axisLayer, value, tickOffset, formatted(value));
      }
      label(this.axisLayer, right - fontSize, -tickOffset, 'X', 'axis-name');
    }
    if (left <= 0 && right >= 0) {
      this.axisLayer.appendChild(M.svgElement('line', { x1: 0, y1: top, x2: 0, y2: bottom }));
      for (value = Math.ceil((top - epsilon) / displaySpacing) * displaySpacing; value <= bottom + epsilon; value += displaySpacing) {
        if (Math.abs(value) <= epsilon) continue;
        this.axisLayer.appendChild(M.svgElement('line', { x1: -tickOffset * 0.2, y1: value, x2: tickOffset * 0.2, y2: value }));
        label(this.axisLayer, -tickOffset * 1.3, value, formatted(-value));
      }
      label(this.axisLayer, tickOffset, top + fontSize, 'Y', 'axis-name');
    }
    if (left <= 0 && right >= 0 && top <= 0 && bottom >= 0) label(this.axisLayer, -tickOffset, tickOffset, '0');
    label(this.axisLayer, left + fontSize * 3.1, bottom - fontSize, 'X →  ·  Y ↑', 'axis-name');
  };

  Renderer.prototype.fit = function () {
    if (!this.model || !this.model.nodes.length) this.view = { x: -1, y: -1, width: 2, height: 2 };
    else {
      var xs = this.model.nodes.map(function (node) { return node.x; });
      var ys = this.model.nodes.map(function (node) { return -node.y; });
      var minX = Math.min.apply(null, xs), maxX = Math.max.apply(null, xs), minY = Math.min.apply(null, ys), maxY = Math.max.apply(null, ys);
      var span = Math.max(maxX - minX, maxY - minY, 1), pad = span * 0.2;
      var rect = this.svg.getBoundingClientRect(), aspect = Math.max(rect.width, 1) / Math.max(rect.height, 1);
      var width = Math.max(maxX - minX + 2 * pad, span * 0.5), height = Math.max(maxY - minY + 2 * pad, span * 0.5);
      if (width / height < aspect) width = height * aspect; else height = width / aspect;
      this.view = { x: (minX + maxX - width) / 2, y: (minY + maxY - height) / 2, width: width, height: height };
    }
    this.applyView();
  };

  Renderer.prototype.setTool = function (tool, pendingNode) { this.tool = tool; this.pendingNode = Number.isInteger(pendingNode) ? pendingNode : null; this.previewPointerPoint = null; this.previewPoint = null; this.draw(); };
  Renderer.prototype.setSelection = function (selection) { this.selection = selection; this.draw(); };

  Renderer.prototype.updateNodePreview = function (point) {
    if (point) this.previewPointerPoint = point;
    if (this.previewPointerPoint) this.previewPoint = this.callbacks.placementPoint ? this.callbacks.placementPoint(this.previewPointerPoint) : this.previewPointerPoint;
    if (!this.previewElement) return;
    var visible = this.tool === 'node' && this.previewPoint;
    this.previewElement.setAttribute('visibility', visible ? 'visible' : 'hidden');
    if (!visible) return;
    this.previewElement.setAttribute('cx', this.previewPoint.x);
    this.previewElement.setAttribute('cy', -this.previewPoint.y);
  };

  Renderer.prototype.refreshNodePreview = function () {
    if (this.previewPointerPoint) this.updateNodePreview();
  };

  Renderer.prototype.nodeClicked = function (index, event) {
    event.stopPropagation();
    if (this.tool === 'member') {
      if (this.pendingNode === null) {
        this.pendingNode = index;
        if (this.callbacks.elementPending) this.callbacks.elementPending(index);
      }
      else if (this.pendingNode !== index) {
        if (!this.callbacks.addElement || this.callbacks.addElement(this.pendingNode, index) !== false) this.pendingNode = null;
      }
      this.draw();
    } else if (this.callbacks.select) this.callbacks.select({ kind: 'node', index: index });
  };

  Renderer.prototype.drawMPCs = function (layer,scale) {
    var self=this;
    (this.model.mpcs || []).forEach(function (m,i) {
      var p=self.model.nodes[m.depNode-1]; if (!p) return;
      var group=M.svgElement('g',{'class':'mpc-constraint','data-mpc-id':i+1,stroke:'#7c3aed','stroke-width':scale*.003,fill:'none'});
      var title=M.svgElement('title'); title.textContent='MPC '+(i+1)+': '+global.MKEFMPC.equation(m); group.appendChild(title);
      m.masters.forEach(function (a) { var q=self.model.nodes[a.node-1]; if (q) group.appendChild(M.svgElement('line',{x1:p.x,y1:-p.y,x2:q.x,y2:-q.y,'stroke-dasharray':scale*.02+' '+scale*.01})); });
      group.appendChild(M.svgElement('circle',{cx:p.x,cy:-p.y,r:scale*.025})); layer.appendChild(group);
    });
  };

  Renderer.prototype.drawSupport = function (layer, support, scale) {
    var node = this.model.nodes[support.node - 1];
    if (!node) return;
    var dofs = this.model.elementType === 112 ? 2 : 3, type;
    try { type = global.MKEFReleases.effectiveSupport(support,global.MKEFReleases.inactive(this.model,support.node),dofs); } catch (_) { return; }
    if (!type) return;
    if (dofs === 2 && type !== 1) type = type === 2 ? 6 : 5;
    global.MKEFSupportMarkers.append(M.svgElement, layer, node.x, -node.y, type, scale * 0.035, { 'aria-hidden': 'true' });
  };

  Renderer.prototype.drawLoad = function (layer, load, scale) {
    var node = this.model.nodes[load.node - 1];
    if (!node) return;
    var magnitude = Math.hypot(load.fx, load.fy), size = scale * 0.12;
    if (magnitude > 0) {
      var dx = load.fx / magnitude * size, dy = -load.fy / magnitude * size;
      layer.appendChild(M.svgElement('line', { x1: node.x - dx, y1: -node.y - dy, x2: node.x, y2: -node.y, class: 'load-symbol', 'marker-end': 'url(#arrow)' }));
    }
    if (load.mz) layer.appendChild(M.svgElement('circle', { cx: node.x, cy: -node.y, r: size * 0.35, class: 'moment-symbol' }));
  };

  Renderer.prototype.draw = function (model) {
    if (model) this.model = model;
    var self = this;
    while (this.svg.firstChild) this.svg.removeChild(this.svg.firstChild);
    var defs = M.svgElement('defs'), marker = M.svgElement('marker', { id: 'arrow', markerWidth: 8, markerHeight: 8, refX: 7, refY: 3, orient: 'auto', markerUnits: 'strokeWidth' });
    marker.appendChild(M.svgElement('path', { d: 'M0,0 L0,6 L8,3 z', class: 'arrow-head' })); defs.appendChild(marker); this.svg.appendChild(defs);
    if (!this.model) return;
    this.gridLayer = M.svgElement('g', { class: 'coordinate-grid', 'aria-hidden': 'true' });
    this.axisLayer = M.svgElement('g', { class: 'coordinate-axes', 'aria-hidden': 'true' });
    this.svg.appendChild(this.gridLayer); this.svg.appendChild(this.axisLayer); this.drawCoordinateSystem();
    var elementLayer = M.svgElement('g', { class: 'elements-layer' }), supportLayer = M.svgElement('g', { class: 'supports-layer' }), nodeLayer = M.svgElement('g', { class: 'nodes-layer' }), labelLayer = M.svgElement('g', { class: 'labels-layer' }), loadLayer = M.svgElement('g', { class: 'loads-layer' });
    this.previewLayer = M.svgElement('g', { class: 'node-preview-layer', 'aria-hidden': 'true' });
    var scale = Math.max(this.view.width, this.view.height, 1);
    var mpcLayer=M.svgElement('g',{'class':'mpcs-layer'}); this.svg.appendChild(mpcLayer); this.drawMPCs(mpcLayer,scale);
    this.svg.appendChild(elementLayer); this.svg.appendChild(supportLayer); this.svg.appendChild(nodeLayer); this.svg.appendChild(labelLayer); this.svg.appendChild(this.previewLayer); this.svg.appendChild(loadLayer);
    this.previewElement = M.svgElement('circle', { r: scale * 0.012, class: 'node-placement-preview', visibility: 'hidden' });
    this.previewLayer.appendChild(this.previewElement); this.updateNodePreview();
    this.model.elements.forEach(function (element, index) {
      var first = self.model.nodes[element.node1 - 1], second = self.model.nodes[element.node2 - 1];
      if (!first || !second) return;
      var line = M.svgElement('line', { x1: first.x, y1: -first.y, x2: second.x, y2: -second.y, class: 'model-element' + (self.selection && self.selection.kind === 'element' && self.selection.index === index ? ' selected' : ''), 'data-element-id': index + 1, tabindex: 0, role: 'button', 'aria-label': 'Element ' + (index + 1) });
      line.addEventListener('click', function (event) { event.stopPropagation(); if (self.callbacks.select) self.callbacks.select({ kind: 'element', index: index }); });
      line.addEventListener('keydown', function (event) { if (event.key === 'Enter' || event.key === ' ') { event.preventDefault(); if (self.callbacks.select) self.callbacks.select({ kind: 'element', index: index }); } });
      elementLayer.appendChild(line);
      [1,2].forEach(function (end) { if (global.MKEFReleases.has(self.model,index+1,end)) global.MKEFReleases.append(M.svgElement,elementLayer,{x:first.x,y:-first.y},{x:second.x,y:-second.y},end,scale*.035,{'data-release-element':index+1,'data-release-end':end}); });
      var text = M.svgElement('text', { x: (first.x + second.x) / 2, y: -(first.y + second.y) / 2, class: 'element-label' }); text.textContent = 'E' + (index + 1); labelLayer.appendChild(text);
    });
    this.model.supports.forEach(function (support) { self.drawSupport(supportLayer, support, scale); });
    this.model.loads.forEach(function (load) { self.drawLoad(loadLayer, load, scale); });
    this.drawElementLoads(loadLayer, scale);
    this.model.nodes.forEach(function (node, index) {
      var circle = M.svgElement('circle', { cx: node.x, cy: -node.y, r: scale * 0.012, class: 'model-node' + (self.selection && self.selection.kind === 'node' && self.selection.index === index ? ' selected' : '') + (self.pendingNode === index ? ' pending' : ''), 'data-node-id': index + 1, tabindex: 0, role: 'button', 'aria-label': 'Node ' + (index + 1) });
      circle.addEventListener('click', function (event) { self.nodeClicked(index, event); }); nodeLayer.appendChild(circle);
      circle.addEventListener('keydown', function (event) { if (event.key === 'Enter' || event.key === ' ') { event.preventDefault(); self.nodeClicked(index, event); } });
      var text = M.svgElement('text', { x: node.x + scale * 0.018, y: -node.y - scale * 0.018, class: 'node-label' }); text.textContent = String(index + 1); labelLayer.appendChild(text);
    });
  };

  Renderer.prototype.drawElementLoads = function (layer, scale) {
    var self = this;
    (this.model.elementLoads || []).forEach(function (load) {
      var element = self.model.elements[load.elementId - 1], a = element && self.model.nodes[element.node1 - 1], b = element && self.model.nodes[element.node2 - 1];
      if (!a || !b) return;
      var L = Math.hypot(b.x-a.x, b.y-a.y), c = (b.x-a.x)/L, s = (b.y-a.y)/L;
      var ends = load.type===21 ? [[load.qx1,load.qy1],[load.qx2,load.qy2]] : [[load.qx,load.qy],[load.qx,load.qy]];
      var vectors=ends.map(function (q) { return load.coordinateSystem===1 ? [c*q[0]-s*q[1],s*q[0]+c*q[1]] : q; });
      var maximum=Math.max(Math.hypot.apply(null,vectors[0]),Math.hypot.apply(null,vectors[1]));
      if (!Number.isFinite(maximum) || !maximum || !L) return;
      var size=Math.min(scale*.08,L*.3), group=M.svgElement('g',{class:'element-load-symbol','data-element-id':load.elementId}), tails=[];
      for (var i=0;i<=8;i++) {
        var xi=i/8, px=a.x+(b.x-a.x)*xi, py=a.y+(b.y-a.y)*xi;
        var x=vectors[0][0]*(1-xi)+vectors[1][0]*xi, y=vectors[0][1]*(1-xi)+vectors[1][1]*xi;
        var tx=px-x/maximum*size, ty=-py+y/maximum*size;
        tails.push(tx+','+ty);
        if (Math.hypot(x,y)>0) group.appendChild(M.svgElement('line',{x1:tx,y1:ty,x2:px,y2:-py,class:'load-symbol','marker-end':'url(#arrow)'}));
      }
      group.appendChild(M.svgElement('polyline',{points:tails.join(' '),class:'load-symbol element-load-envelope',fill:'none'}));
      var midX=(vectors[0][0]+vectors[1][0])/2, labelY=Math.min.apply(null,tails.map(function (point) { return Number(point.split(',')[1]); }));
      var label=M.svgElement('text',{x:(a.x+b.x)/2-midX/maximum*size,y:labelY-scale*.018,'text-anchor':'middle','font-size':scale*.018,class:'element-load-label'});
      function pair(q) { return '('+q.map(function (v) { return Number(v.toPrecision(5)); }).join(', ')+')'; }
      label.textContent=(load.coordinateSystem===1 ? 'Local' : 'Global') + (load.type===21 ? ' q1='+pair(ends[0])+' → q2='+pair(ends[1]) : ' q='+pair(ends[0]));
      group.appendChild(label); layer.appendChild(group);
    });
    if (this.selection && this.selection.kind === 'element') {
      var e = this.model.elements[this.selection.index], first = e && this.model.nodes[e.node1-1], second = e && this.model.nodes[e.node2-1];
      if (!first || !second) return;
      var length = Math.hypot(second.x-first.x, second.y-first.y); if (!Number.isFinite(length) || !length) return;
      var c = (second.x-first.x)/length, s = (second.y-first.y)/length, size = Math.min(scale*.07, length*.3);
      [[c, s, 'x local'], [-s, c, 'y local']].forEach(function (axis) {
        var x = first.x+axis[0]*size, y = first.y+axis[1]*size;
        layer.appendChild(M.svgElement('line', { x1:first.x, y1:-first.y, x2:x, y2:-y, class:'local-axis', 'marker-end':'url(#arrow)' }));
        var text = M.svgElement('text', { x:x, y:-y-scale*.012, 'font-size':scale*.018, class:'element-load-label' }); text.textContent=axis[2]; layer.appendChild(text);
      });
    }
  };

  M.Renderer = Renderer;
}(window));
