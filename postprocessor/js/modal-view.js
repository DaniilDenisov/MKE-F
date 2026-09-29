(function (M) {
  'use strict';
  var V = M.modalView = {};
  V.visualCycleMilliseconds = 2000;

  function modeColumn(dataset, modeIndex) {
    return dataset.raw.analysis.modeShapes.map(function (row) { return row[modeIndex]; });
  }

  function largestComponent(vector, ids) {
    var index = -1, maximum = 0;
    ids.forEach(function (id) {
      var candidate = Math.abs(vector[id - 1]);
      if (candidate > maximum) { maximum = candidate; index = id - 1; }
    });
    return { index: index, maximum: maximum };
  }

  V.displacement = function (dataset, modeIndex) {
    return modeColumn(dataset, modeIndex);
  };

  V.canonicalMode = function (dataset, modeIndex) {
    var source = modeColumn(dataset, modeIndex), translational = [];
    dataset.raw.model.dofMap.forEach(function (row) {
      translational.push(row[0], row[1]);
    });
    var pivot = largestComponent(source, translational);
    if (pivot.index < 0) pivot = largestComponent(source, source.map(function (_, index) { return index + 1; }));
    if (pivot.index < 0) return { vector: source.slice(), zero: true, pivotDOF: null, normalization: 0, sign: 1 };
    var sign = source[pivot.index] < 0 ? -1 : 1;
    return {
      vector: source.map(function (value) { return sign * value / pivot.maximum; }),
      zero: false,
      pivotDOF: pivot.index + 1,
      normalization: pivot.maximum,
      sign: sign
    };
  };

  V.phaseVector = function (vector, factor) {
    return vector.map(function (value) { return value * factor; });
  };

  V.phaseAngleForFactor = function (factor) {
    return Math.acos(Math.max(-1, Math.min(1, Number(factor))));
  };

  V.phaseAt = function (startAngle, elapsedMilliseconds, speed) {
    var angle = startAngle + 2 * Math.PI * elapsedMilliseconds * speed / V.visualCycleMilliseconds;
    angle %= 2 * Math.PI;
    if (angle < 0) angle += 2 * Math.PI;
    return { angle: angle, factor: Math.cos(angle) };
  };

  V.clampModeIndex = function (index, count) {
    return Math.max(0, Math.min(count - 1, Number(index)));
  };

  V.display = function (dataset, modeIndex, multiplier, samples) {
    var canonical = V.canonicalMode(dataset, modeIndex);
    var automatic = canonical.zero ? { scale: 1, zero: true } : M.geometry.automaticScale(dataset, canonical.vector, samples);
    var zero = canonical.zero || automatic.zero;
    var peakScale = zero ? 0 : automatic.scale * multiplier;
    return { modeIndex: modeIndex, vector: canonical.vector, scale: peakScale, peakScale: peakScale, displayMultiplier: multiplier, zero: zero,
      pivotDOF: canonical.pivotDOF, normalization: canonical.normalization, sign: canonical.sign,
      frequencyHz: dataset.raw.analysis.frequenciesHz[modeIndex],
      angularFrequency: dataset.raw.analysis.angularFrequenciesRadPerSec[modeIndex] };
  };
}(window.MKEFPost));
