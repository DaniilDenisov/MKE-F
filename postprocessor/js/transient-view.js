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
    ['displacements', 'velocities', 'accelerations', 'loadHistory', 'reactions', 'supportReactions', 'mpcForces'].forEach(function (name) { if (analysis[name]) values.push(name); });
    if (analysis.spectrumFrequencyHz && analysis.displacementAmplitudeSpectrum) values.push('spectrum');
    return values;
  };

  T.quantityDescriptor = function (dataset, dofId, quantity) {
    var units = dataset.raw.metadata.units || {}, location = dataset.dofById.get(dofId), rotational = location && location.localIndex === 2;
    var timeUnit = units.time || 'time unit', displacementUnit = rotational ? 'rad' : (units.length || ''), actionUnit = rotational ? (units.moment || '') : (units.force || '');
    var descriptors = {
      displacements: { label: 'Displacement', unit: displacementUnit },
      velocities: { label: 'Velocity', unit: displacementUnit ? displacementUnit + '/' + timeUnit : '' },
      accelerations: { label: 'Acceleration', unit: displacementUnit ? displacementUnit + '/' + timeUnit + '²' : '' },
      loadHistory: { label: 'Applied load', unit: actionUnit },
      supportReactions: { label: 'Support reaction', unit: actionUnit },
      mpcForces: { label: 'MPC force', unit: actionUnit },
      reactions: { label: 'Dynamic residual / reaction', unit: actionUnit },
      spectrum: { label: 'Displacement amplitude', unit: displacementUnit }
    };
    return descriptors[quantity] || { label: quantity, unit: '' };
  };

  T.series = function (dataset, dofId, quantity) {
    var analysis = dataset.raw.analysis, row = analysis.globalDOFIds.indexOf(dofId);
    if (row < 0) return null;
    var descriptor = T.quantityDescriptor(dataset, dofId, quantity), location = dataset.dofById.get(dofId);
    var yLabel = descriptor.label + ' · ' + location.label + (descriptor.unit ? ' (' + descriptor.unit + ')' : '');
    if (quantity === 'spectrum') return { x: analysis.spectrumFrequencyHz, y: analysis.displacementAmplitudeSpectrum[row], xLabel: 'Frequency (Hz)', yLabel: yLabel, quantity: quantity, dofId: dofId, unit: descriptor.unit };
    if (!analysis[quantity]) return null;
    return { x: analysis.time, y: analysis[quantity][row], xLabel: 'Time' + (dataset.raw.metadata.units.time ? ' (' + dataset.raw.metadata.units.time + ')' : ''), yLabel: yLabel, quantity: quantity, dofId: dofId, unit: descriptor.unit };
  };

  T.fieldStatus = function (dataset) {
    var analysis = dataset.raw.analysis, names = ['displacements', 'velocities', 'accelerations', 'loadHistory', 'reactions', 'supportReactions', 'mpcForces'], present = [], omitted = [];
    names.forEach(function (name) { (analysis[name] ? present : omitted).push(name); });
    if (analysis.spectrumFrequencyHz && analysis.displacementAmplitudeSpectrum) present.push('spectrum'); else omitted.push('spectrum');
    return { present: present, omitted: omitted };
  };

  T.samplingSummary = function (dataset) {
    var sampling = dataset.raw.analysis.sampling;
    return { original: sampling.originalSampleCount, exported: sampling.exportedSampleCount, stride: sampling.timeStride,
      decimated: sampling.timeStride > 1 || sampling.originalSampleCount > sampling.exportedSampleCount };
  };

  T.clampFrameIndex = function (index, count) {
    return Math.max(0, Math.min(count - 1, Math.round(Number(index) || 0)));
  };

  T.stepFrameIndex = function (index, direction, count, stride) {
    var step = Math.max(1, Math.round(Number(stride) || 1));
    return T.clampFrameIndex(Number(index) + (direction < 0 ? -step : step), count);
  };

  T.playbackFrames = function (count, startIndex, stride) {
    var start = T.clampFrameIndex(startIndex, count), step = Math.max(1, Math.round(Number(stride) || 1)), indices = [start];
    for (var index = start + step; index < count - 1; index += step) indices.push(index);
    if (indices[indices.length - 1] !== count - 1) indices.push(count - 1);
    return indices;
  };

  T.playbackFrame = function (time, startIndex, elapsedMilliseconds, speed, stride, durationMilliseconds) {
    var frames = T.playbackFrames(time.length, startIndex, stride), duration = durationMilliseconds || 5000;
    var totalSpan = time[time.length - 1] - time[0], remainingSpan = time[frames[frames.length - 1]] - time[frames[0]];
    var remainingDuration = totalSpan > 0 ? duration * remainingSpan / totalSpan : 0;
    var progress = remainingDuration > 0 ? Math.max(0, Math.min(1, elapsedMilliseconds * Math.max(Number(speed) || 1, 0) / remainingDuration)) : 1;
    if (progress >= 1 || frames.length === 1) return { index: frames[frames.length - 1], complete: true, progress: 1 };
    var firstTime = time[frames[0]], lastTime = time[frames[frames.length - 1]], target = firstTime + progress * (lastTime - firstTime), selected = frames[0];
    for (var i = 1; i < frames.length && time[frames[i]] <= target; i += 1) selected = frames[i];
    return { index: selected, complete: false, progress: progress };
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
    return { x: x, y: y, xLabel: series.xLabel, yLabel: series.yLabel, quantity: series.quantity, dofId: series.dofId, unit: series.unit };
  };
}(window.MKEFPost));
