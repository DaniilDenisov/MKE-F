(function (M) {
  'use strict';
  var V = M.modalView = {};
  V.displacement = function (dataset, modeIndex) {
    return dataset.raw.analysis.modeShapes.map(function (row) { return row[modeIndex]; });
  };
  V.display = function (dataset, modeIndex, multiplier, samples) {
    var vector = V.displacement(dataset, modeIndex);
    var bounds = M.geometry.modelBounds(dataset.raw.model);
    var diagonal = Math.max(Math.hypot(bounds.maxX - bounds.minX, bounds.maxY - bounds.minY), 1e-9);
    var maximum = 0;
    dataset.raw.model.dofMap.forEach(function (dofs) {
      maximum = Math.max(maximum, Math.hypot(vector[dofs[0] - 1], vector[dofs[1] - 1]));
    });
    var zero = maximum <= 1e-12 * Math.max(diagonal, 1);
    var scale = zero ? 0 : .1 * diagonal / maximum;
    return { vector: vector, scale: scale * multiplier, zero: zero,
      frequencyHz: dataset.raw.analysis.frequenciesHz[modeIndex],
      angularFrequency: dataset.raw.analysis.angularFrequenciesRadPerSec[modeIndex] };
  };
}(window.MKEFPost));
