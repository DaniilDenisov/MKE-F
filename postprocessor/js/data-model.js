(function (M) {
  'use strict';

  function DataError(path, message) {
    this.name = 'DataError';
    this.path = path;
    this.message = path + ': ' + message;
  }
  DataError.prototype = Object.create(Error.prototype);
  DataError.prototype.constructor = DataError;

  function fail(path, message) { throw new DataError(path, message); }
  function objectAt(value, path) {
    if (!value || typeof value !== 'object' || Array.isArray(value)) fail(path, 'expected an object');
    return value;
  }
  function arrayAt(value, path) {
    if (!Array.isArray(value)) fail(path, 'expected an array');
    return value;
  }
  function numberAt(value, path) {
    if (typeof value !== 'number' || !Number.isFinite(value)) fail(path, 'expected a finite number');
    return value;
  }
  function positiveInteger(value, path) {
    numberAt(value, path);
    if (!Number.isInteger(value) || value < 1) fail(path, 'expected a positive integer');
    return value;
  }
  function vector(value, length, path) {
    arrayAt(value, path);
    if (value.length !== length) fail(path, 'expected ' + length + ' values');
    value.forEach(function (item, index) { numberAt(item, path + '[' + index + ']'); });
    return value;
  }
  function matrix(value, rows, columns, path) {
    arrayAt(value, path);
    if (value.length !== rows) fail(path, 'expected ' + rows + ' rows');
    value.forEach(function (row, index) { vector(row, columns, path + '[' + index + ']'); });
    return value;
  }

  function validateModel(model) {
    objectAt(model, '$.model');
    if (model.dimension !== 2) fail('$.model.dimension', 'only two-dimensional models are supported');
    if (model.dofPerNode !== 2 && model.dofPerNode !== 3) fail('$.model.dofPerNode', 'expected 2 or 3');
    var labels = arrayAt(model.dofLabels, '$.model.dofLabels');
    if (labels.length !== model.dofPerNode) fail('$.model.dofLabels', 'length does not match dofPerNode');

    var nodes = arrayAt(model.nodes, '$.model.nodes');
    if (!nodes.length || nodes.length > M.config.maxNodes) fail('$.model.nodes', 'invalid or unreasonable node count');
    var nodeIds = new Set();
    nodes.forEach(function (node, index) {
      var path = '$.model.nodes[' + index + ']';
      objectAt(node, path);
      positiveInteger(node.id, path + '.id');
      if (nodeIds.has(node.id)) fail(path + '.id', 'duplicate node ID');
      nodeIds.add(node.id);
      numberAt(node.x, path + '.x'); numberAt(node.y, path + '.y');
    });

    var dofCount = nodes.length * model.dofPerNode;
    matrix(model.dofMap, nodes.length, model.dofPerNode, '$.model.dofMap');
    var dofs = new Set();
    model.dofMap.forEach(function (row, rowIndex) {
      row.forEach(function (value, columnIndex) {
        var path = '$.model.dofMap[' + rowIndex + '][' + columnIndex + ']';
        positiveInteger(value, path);
        if (value > dofCount) fail(path, 'global DOF exceeds model size');
        if (dofs.has(value)) fail(path, 'duplicate global DOF');
        dofs.add(value);
      });
    });
    if (dofs.size !== dofCount) fail('$.model.dofMap', 'does not contain every global DOF');

    var elements = arrayAt(model.elements, '$.model.elements');
    if (!elements.length || elements.length > M.config.maxElements) fail('$.model.elements', 'invalid or unreasonable element count');
    var elementIds = new Set();
    elements.forEach(function (element, index) {
      var path = '$.model.elements[' + index + ']';
      objectAt(element, path);
      positiveInteger(element.id, path + '.id');
      if (elementIds.has(element.id)) fail(path + '.id', 'duplicate element ID');
      elementIds.add(element.id);
      if (element.type !== 112 && element.type !== 113) fail(path + '.type', 'unsupported element type');
      if ((element.type === 112 ? 2 : 3) !== model.dofPerNode) fail(path + '.type', 'does not match model DOFs');
      vector(element.nodeIds, 2, path + '.nodeIds');
      element.nodeIds.forEach(function (id, j) {
        positiveInteger(id, path + '.nodeIds[' + j + ']');
        if (!nodeIds.has(id)) fail(path + '.nodeIds[' + j + ']', 'unknown node ID');
      });
      if (element.nodeIds[0] === element.nodeIds[1]) fail(path + '.nodeIds', 'element endpoints must differ');
      var properties = objectAt(element.properties, path + '.properties');
      ['area', 'youngsModulus', 'density'].concat(element.type === 113 ? ['momentOfInertia'] : []).forEach(function (name) {
        numberAt(properties[name], path + '.properties.' + name);
        if (properties[name] <= 0) fail(path + '.properties.' + name, 'expected a positive value');
      });
    });

    var supports = arrayAt(model.supports, '$.model.supports');
    supports.forEach(function (support, index) {
      var path = '$.model.supports[' + index + ']';
      objectAt(support, path);
      positiveInteger(support.type, path + '.type');
      try { window.MKEFSupports.mask(support, model.dofPerNode); } catch (error) { fail(path + '.type', error.message); }
      positiveInteger(support.nodeId, path + '.nodeId');
      if (!nodeIds.has(support.nodeId)) fail(path + '.nodeId', 'unknown node ID');
    });
    return { dofCount: dofCount, nodeIds: nodeIds, elementIds: elementIds };
  }

  function validateStatic(analysis, context) {
    var count = context.dofCount;
    vector(analysis.displacements, count, '$.analysis.displacements');
    vector(analysis.loadVector, count, '$.analysis.loadVector');
    vector(analysis.reactions, count, '$.analysis.reactions');
    vector(analysis.equilibriumResidual, 3, '$.analysis.equilibriumResidual');
    var results = arrayAt(analysis.elementResults, '$.analysis.elementResults');
    if (results.length !== context.elementIds.size) fail('$.analysis.elementResults', 'expected one record per element');
    var resultIds = new Set();
    results.forEach(function (result, index) {
      var path = '$.analysis.elementResults[' + index + ']';
      objectAt(result, path);
      positiveInteger(result.elementId, path + '.elementId');
      if (!context.elementIds.has(result.elementId)) fail(path + '.elementId', 'unknown element ID');
      if (resultIds.has(result.elementId)) fail(path + '.elementId', 'duplicate element result');
      resultIds.add(result.elementId);
      if (result.type !== 112 && result.type !== 113) fail(path + '.type', 'unsupported result type');
      vector(result.localEndForces, result.type === 112 ? 4 : 6, path + '.localEndForces');
      ['axialStrain', 'axialStress', 'axialForce'].forEach(function (name) { numberAt(result[name], path + '.' + name); });
    });
  }

  function validateModal(analysis, context) {
    var frequencies = arrayAt(analysis.frequenciesHz, '$.analysis.frequenciesHz');
    if (!frequencies.length) fail('$.analysis.frequenciesHz', 'expected at least one mode');
    frequencies.forEach(function (v, i) { numberAt(v, '$.analysis.frequenciesHz[' + i + ']'); });
    vector(analysis.angularFrequenciesRadPerSec, frequencies.length, '$.analysis.angularFrequenciesRadPerSec');
    matrix(analysis.modeShapes, context.dofCount, frequencies.length, '$.analysis.modeShapes');
  }

  function validateTransient(analysis, context) {
    var time = arrayAt(analysis.time, '$.analysis.time');
    if (!time.length) fail('$.analysis.time', 'expected at least one time sample');
    time.forEach(function (value, index) {
      numberAt(value, '$.analysis.time[' + index + ']');
      if (index && value <= time[index - 1]) fail('$.analysis.time[' + index + ']', 'time must increase');
    });
    var ids = arrayAt(analysis.globalDOFIds, '$.analysis.globalDOFIds');
    if (!ids.length) fail('$.analysis.globalDOFIds', 'expected at least one global DOF');
    var unique = new Set();
    ids.forEach(function (id, index) {
      positiveInteger(id, '$.analysis.globalDOFIds[' + index + ']');
      if (id > context.dofCount || unique.has(id)) fail('$.analysis.globalDOFIds[' + index + ']', 'invalid or duplicate global DOF');
      unique.add(id);
    });
    ['displacements', 'velocities', 'accelerations', 'loadHistory', 'reactions'].forEach(function (name) {
      if (Object.prototype.hasOwnProperty.call(analysis, name)) matrix(analysis[name], ids.length, time.length, '$.analysis.' + name);
    });
    var hasSpectrumFrequency = Object.prototype.hasOwnProperty.call(analysis, 'spectrumFrequencyHz'), hasSpectrumAmplitude = Object.prototype.hasOwnProperty.call(analysis, 'displacementAmplitudeSpectrum');
    if (hasSpectrumFrequency !== hasSpectrumAmplitude) fail('$.analysis', 'spectrumFrequencyHz and displacementAmplitudeSpectrum must be provided together');
    if (hasSpectrumFrequency) {
      var frequency = arrayAt(analysis.spectrumFrequencyHz, '$.analysis.spectrumFrequencyHz');
      frequency.forEach(function (v, i) { numberAt(v, '$.analysis.spectrumFrequencyHz[' + i + ']'); });
      matrix(analysis.displacementAmplitudeSpectrum, ids.length, frequency.length, '$.analysis.displacementAmplitudeSpectrum');
    }
    objectAt(analysis.sampling, '$.analysis.sampling');
    positiveInteger(analysis.sampling.originalSampleCount, '$.analysis.sampling.originalSampleCount');
    positiveInteger(analysis.sampling.exportedSampleCount, '$.analysis.sampling.exportedSampleCount');
    positiveInteger(analysis.sampling.timeStride, '$.analysis.sampling.timeStride');
    if (analysis.sampling.exportedSampleCount !== time.length) fail('$.analysis.sampling.exportedSampleCount', 'does not match time length');
    if (analysis.sampling.originalSampleCount < analysis.sampling.exportedSampleCount) fail('$.analysis.sampling.originalSampleCount', 'cannot be smaller than exportedSampleCount');
    if (analysis.sampling.timeStride === 1 && analysis.sampling.originalSampleCount !== analysis.sampling.exportedSampleCount) fail('$.analysis.sampling', 'timeStride 1 requires every original sample to be exported');
  }

  M.validateDataset = function (data) {
    objectAt(data, '$');
    if (data.format !== 'mkef-postprocessor') fail('$.format', 'expected "mkef-postprocessor"');
    if (data.version !== 1 && data.version !== 2) fail('$.version', 'unsupported format version ' + String(data.version));
    var metadata = objectAt(data.metadata, '$.metadata');
    if (typeof metadata.title !== 'string') fail('$.metadata.title', 'expected a string');
    objectAt(metadata.units, '$.metadata.units');
    var context = validateModel(data.model);
    var analysis = objectAt(data.analysis, '$.analysis');
    if (analysis.type === 'static') validateStatic(analysis, context);
    else if (analysis.type === 'modal') validateModal(analysis, context);
    else if (analysis.type === 'transient') validateTransient(analysis, context);
    else fail('$.analysis.type', 'unsupported analysis type');
    if (data.version === 2) {
      if (analysis.type !== 'static') fail('$.analysis.type', 'version 2 element loads require static analysis');
      var memberMap = new Map(data.model.elements.map(function (item) { return [item.id, item]; }));
      arrayAt(data.model.nodalLoads, '$.model.nodalLoads').forEach(function (load, index) {
        var path = '$.model.nodalLoads[' + index + ']'; objectAt(load, path);
        if (load.type !== 10 || !context.nodeIds.has(load.nodeId)) fail(path, 'expected a static load on an existing node');
        ['fx', 'fy', 'mz'].forEach(function (name) { numberAt(load[name], path + '.' + name); });
      });
      arrayAt(data.model.elementLoads, '$.model.elementLoads').forEach(function (load, index) {
        var path = '$.model.elementLoads[' + index + ']'; objectAt(load, path);
        var member = memberMap.get(load.elementId);
        if (load.type !== 20 || !member || member.type !== 113 || [1, 2].indexOf(load.coordinateSystem) < 0) fail(path, 'expected a uniform load on a frame in Local/Global coordinates');
        numberAt(load.qx, path + '.qx'); numberAt(load.qy, path + '.qy');
        if (!load.qx && !load.qy) fail(path, 'intensity cannot be zero');
      });
      analysis.elementResults.forEach(function (result, index) { vector(result.equivalentLocalLoadVector, result.type === 113 ? 6 : 4, '$.analysis.elementResults[' + index + '].equivalentLocalLoadVector'); });
    }

    var nodesById = new Map(data.model.nodes.map(function (node) { return [node.id, node]; }));
    var elementsById = new Map(data.model.elements.map(function (element) { return [element.id, element]; }));
    var nodeIndexById = new Map(data.model.nodes.map(function (node, index) { return [node.id, index]; }));
    var dofById = new Map();
    data.model.dofMap.forEach(function (row, nodeIndex) {
      row.forEach(function (dofId, localIndex) {
        dofById.set(dofId, { id: dofId, nodeId: data.model.nodes[nodeIndex].id, nodeIndex: nodeIndex,
          localIndex: localIndex, label: data.model.dofLabels[localIndex] });
      });
    });
    var resultsByElementId = new Map();
    if (analysis.type === 'static') analysis.elementResults.forEach(function (result) { resultsByElementId.set(result.elementId, result); });
    return { raw: data, nodesById: nodesById, elementsById: elementsById,
      nodeIndexById: nodeIndexById, dofById: dofById, resultsByElementId: resultsByElementId };
  };
  M.DataError = DataError;
}(window.MKEFPost));
