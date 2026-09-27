(function (M) {
  'use strict';
  var T = M.transientView = {};

  T.hasFullDisplacements = function (dataset) {
    var analysis = dataset.raw.analysis, count = dataset.raw.model.nodes.length * dataset.raw.model.dofPerNode;
    if (!analysis.displacements || analysis.globalDOFIds.length !== count) return false;
    var ids = analysis.globalDOFIds.slice().sort(function (a, b) { return a - b; });
    return ids.every(function (id, index) { return id === index + 1; });
  };

  T.vectorAt = function (dataset, field, timeIndex) {
    var analysis = dataset.raw.analysis, count = dataset.raw.model.nodes.length * dataset.raw.model.dofPerNode;
    if (!analysis[field]) return null;
    var vector = new Array(count);
    analysis.globalDOFIds.forEach(function (id, row) { vector[id - 1] = analysis[field][row][timeIndex]; });
    return vector;
  };

  T.displayScale = function (dataset) {
    var analysis = dataset.raw.analysis;
    if (!analysis.displacements) return 1;
    var bounds = M.geometry.modelBounds(dataset.raw.model);
    var size = Math.max(Math.hypot(bounds.maxX - bounds.minX, bounds.maxY - bounds.minY), 1e-9);
    var maximum = 0;
    analysis.globalDOFIds.forEach(function (id, row) {
      var location = dataset.dofById.get(id), factor = location && location.localIndex === 2 ? size : 1;
      analysis.displacements[row].forEach(function (value) { maximum = Math.max(maximum, Math.abs(value) * factor); });
    });
    return maximum <= 1e-12 * Math.max(size, 1) ? 1 : .1 * size / maximum;
  };

  T.availableQuantities = function (dataset) {
    var analysis = dataset.raw.analysis, values = [];
    ['displacements', 'velocities', 'accelerations', 'loadHistory', 'reactions'].forEach(function (name) { if (analysis[name]) values.push(name); });
    if (analysis.spectrumFrequencyHz && analysis.displacementAmplitudeSpectrum) values.push('spectrum');
    return values;
  };

  T.series = function (dataset, dofId, quantity) {
    var analysis = dataset.raw.analysis, row = analysis.globalDOFIds.indexOf(dofId);
    if (row < 0) return null;
    if (quantity === 'spectrum') return { x: analysis.spectrumFrequencyHz, y: analysis.displacementAmplitudeSpectrum[row], xLabel: 'Frequency (Hz)', yLabel: 'Amplitude' };
    if (!analysis[quantity]) return null;
    return { x: analysis.time, y: analysis[quantity][row], xLabel: 'Time' + (dataset.raw.metadata.units.time ? ' (' + dataset.raw.metadata.units.time + ')' : ''), yLabel: quantity };
  };

  T.downsample = function (series, maximumPoints) {
    if (!series || series.x.length <= maximumPoints) return series;
    var x = [], y = [], lastIndex = -1, bucketCount = Math.max(1, Math.floor(maximumPoints / 2));
    for (var bucket = 0; bucket < bucketCount; bucket += 1) {
      var start = Math.floor(bucket * series.x.length / bucketCount), end = Math.floor((bucket + 1) * series.x.length / bucketCount);
      if (end <= start) continue;
      var minimum = start, maximum = start;
      for (var i = start + 1; i < end; i += 1) { if (series.y[i] < series.y[minimum]) minimum = i; if (series.y[i] > series.y[maximum]) maximum = i; }
      [minimum, maximum].sort(function (a, b) { return a - b; }).forEach(function (index) { if (index !== lastIndex) { x.push(series.x[index]); y.push(series.y[index]); lastIndex = index; } });
    }
    return { x: x, y: y, xLabel: series.xLabel, yLabel: series.yLabel };
  };
}(window.MKEFPost));
