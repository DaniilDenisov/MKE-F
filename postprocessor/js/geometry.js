(function (M) {
  'use strict';
  var G = M.geometry = {};

  G.clamp = function (value, minimum, maximum) { return Math.max(minimum, Math.min(maximum, value)); };
  G.svgPoint = function (point) { return { x: point.x, y: -point.y }; };
  G.nodeDisplacement = function (dataset, nodeId, displacement) {
    var row = dataset.raw.model.nodes.findIndex(function (node) { return node.id === nodeId; });
    var dofs = dataset.raw.model.dofMap[row];
    return { x: displacement[dofs[0] - 1], y: displacement[dofs[1] - 1], r: dofs.length > 2 ? displacement[dofs[2] - 1] : 0 };
  };
  G.elementBasis = function (first, second) {
    var dx = second.x - first.x, dy = second.y - first.y;
    var length = Math.hypot(dx, dy);
    if (!(length > 0)) throw new Error('Element length must be positive.');
    return { c: dx / length, s: dy / length, length: length };
  };
  G.sample112 = function (first, second, firstD, secondD, scale) {
    return [
      { x: first.x + scale * firstD.x, y: first.y + scale * firstD.y },
      { x: second.x + scale * secondD.x, y: second.y + scale * secondD.y }
    ];
  };
  G.sample112Displacements = function (firstD, secondD, sampleCount) {
    var values = [];
    for (var i = 0; i < sampleCount; i += 1) {
      var xi = i / (sampleCount - 1);
      values.push({
        x: (1 - xi) * firstD.x + xi * secondD.x,
        y: (1 - xi) * firstD.y + xi * secondD.y
      });
    }
    return values;
  };
  G.sample113Displacements = function (first, second, firstD, secondD, sampleCount) {
    var basis = G.elementBasis(first, second), c = basis.c, s = basis.s, L = basis.length;
    var u1 = c * firstD.x + s * firstD.y;
    var v1 = -s * firstD.x + c * firstD.y;
    var u2 = c * secondD.x + s * secondD.y;
    var v2 = -s * secondD.x + c * secondD.y;
    var values = [];
    for (var i = 0; i < sampleCount; i += 1) {
      var xi = i / (sampleCount - 1), xi2 = xi * xi, xi3 = xi2 * xi;
      var axial = (1 - xi) * u1 + xi * u2;
      var n1 = 1 - 3 * xi2 + 2 * xi3;
      var n2 = L * (xi - 2 * xi2 + xi3);
      var n3 = 3 * xi2 - 2 * xi3;
      var n4 = L * (-xi2 + xi3);
      var transverse = n1 * v1 + n2 * firstD.r + n3 * v2 + n4 * secondD.r;
      values.push({ x: c * axial - s * transverse, y: s * axial + c * transverse });
    }
    return values;
  };
  G.sample113 = function (first, second, firstD, secondD, scale, sampleCount) {
    var basis = G.elementBasis(first, second), displacements = G.sample113Displacements(first, second, firstD, secondD, sampleCount), points = [];
    for (var i = 0; i < sampleCount; i += 1) {
      var xi = i / (sampleCount - 1);
      points.push({
        x: first.x + basis.c * xi * basis.length + scale * displacements[i].x,
        y: first.y + basis.s * xi * basis.length + scale * displacements[i].y
      });
    }
    return points;
  };
  G.sampleElementDisplacements = function (dataset, element, displacement, sampleCount) {
    var first = dataset.nodesById.get(element.nodeIds[0]);
    var second = dataset.nodesById.get(element.nodeIds[1]);
    var firstD = G.nodeDisplacement(dataset, first.id, displacement);
    var secondD = G.nodeDisplacement(dataset, second.id, displacement);
    return element.type === 112 ? G.sample112Displacements(firstD, secondD, sampleCount) : G.sample113Displacements(first, second, firstD, secondD, sampleCount);
  };
  G.sampleElement = function (dataset, element, displacement, scale, samples) {
    var first = dataset.nodesById.get(element.nodeIds[0]);
    var second = dataset.nodesById.get(element.nodeIds[1]);
    var firstD = G.nodeDisplacement(dataset, first.id, displacement);
    var secondD = G.nodeDisplacement(dataset, second.id, displacement);
    return element.type === 112 ? G.sample112(first, second, firstD, secondD, scale) : G.sample113(first, second, firstD, secondD, scale, samples);
  };
  G.modelBounds = function (model) {
    var xs = model.nodes.map(function (node) { return node.x; });
    var ys = model.nodes.map(function (node) { return node.y; });
    return { minX: Math.min.apply(null, xs), maxX: Math.max.apply(null, xs), minY: Math.min.apply(null, ys), maxY: Math.max.apply(null, ys) };
  };
  G.automaticScale = function (dataset, displacement, samples) {
    var bounds = G.modelBounds(dataset.raw.model);
    var diagonal = Math.hypot(bounds.maxX - bounds.minX, bounds.maxY - bounds.minY);
    if (!(diagonal > 0)) diagonal = 1;
    var maximum = 0;
    dataset.raw.model.elements.forEach(function (element) {
      var deformed = G.sampleElement(dataset, element, displacement, 1, samples);
      var first = dataset.nodesById.get(element.nodeIds[0]), second = dataset.nodesById.get(element.nodeIds[1]);
      deformed.forEach(function (point, index) {
        var xi = index / (deformed.length - 1);
        var x = first.x + xi * (second.x - first.x), y = first.y + xi * (second.y - first.y);
        maximum = Math.max(maximum, Math.hypot(point.x - x, point.y - y));
      });
    });
    return maximum <= 1e-12 * Math.max(diagonal, 1) ? { scale: 1, zero: true } : { scale: 0.1 * diagonal / maximum, zero: false };
  };
  G.pathData = function (points) {
    return points.map(function (point, index) {
      var svg = G.svgPoint(point);
      return (index ? 'L' : 'M') + ' ' + svg.x.toPrecision(12) + ' ' + svg.y.toPrecision(12);
    }).join(' ');
  };
  G.restrainedLocalDOFs = function (type, dofPerNode) {
    if (type === 1) return Array.from({ length: dofPerNode }, function (_, i) { return i; });
    if (type === 2) return Array.from({ length: dofPerNode - 1 }, function (_, i) { return i + 1; });
    if (type === 3) return [0].concat(dofPerNode === 3 ? [2] : []);
    return [0, 1];
  };
}(window.MKEFPost));
