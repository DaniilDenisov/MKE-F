(function (global) {
  'use strict';
  var M = global.MKEFPre;

  function fail(line, message) {
    throw new Error('Line ' + line + ': ' + message);
  }

  function finite(value) {
    return typeof value === 'number' && Number.isFinite(value);
  }

  function integer(value) {
    return finite(value) && Math.floor(value) === value;
  }

  function numberText(value) {
    if (!finite(value)) throw new Error('Cannot serialize a non-finite number.');
    return String(value);
  }

  function newModel(analysisType, elementType) {
    return {
      name: 'Case',
      analysis: {
        type: analysisType || '',
        timeStep: 0.0001,
        duration: 0.01,
        monitorNode: 1,
        monitorDOF: 1
      },
      elementType: elementType === undefined ? 113 : elementType,
      nodes: [],
      elements: [],
      supports: [],
      mpcs: [],
      releases: [],
      loads: [],
      elementLoads: []
    };
  }

  function parseRecord(text, count, line, description) {
    var fields = text.split(',');
    if (fields.length !== count) fail(line, description + ' requires exactly ' + count + ' fields.');
    return fields.map(function (field) {
      var value = field.trim() === '' ? NaN : Number(field.trim());
      if (!Number.isFinite(value)) fail(line, description + ' contains a nonnumeric or non-finite value.');
      return value;
    });
  }

  function parseCase(text) {
    if (typeof text !== 'string') throw new Error('Case input must be text.');
    var lines = text.replace(/^\uFEFF/, '').split(/\r?\n/), cursor = 0;
    var model = newModel('', 113), hasNodes = false, hasElements = false, hasAnalysis = false;

    function nextData(description) {
      while (cursor < lines.length) {
        var lineNumber = cursor + 1, value = lines[cursor++].trim();
        if (value && value.charAt(0) !== '#') return { text: value, line: lineNumber };
      }
      fail(lines.length + 1, 'Unexpected end of file while reading ' + description + '.');
    }

    function count(description, requirePositive) {
      var item = nextData(description + ' count'), value = Number(item.text);
      if (!integer(value) || value < (requirePositive ? 1 : 0)) {
        fail(item.line, description + ' requires an integer count >= ' + (requirePositive ? 1 : 0) + '.');
      }
      return value;
    }

    function requireNodes(line) {
      if (!hasNodes) fail(line, 'The nodes section must appear before element sections.');
    }

    function requireElements(line) {
      if (!hasElements) fail(line, 'An element section must appear before boundary conditions.');
    }

    while (cursor < lines.length) {
      var raw = lines[cursor++], marker = raw.trim(), markerLine = cursor;
      if (!marker || marker.charAt(0) === '#') continue;

      if (marker === 'analysis') {
        if (hasAnalysis) fail(markerLine, 'The analysis section may appear only once.');
        var analysisRecord = nextData('analysis configuration');
        var analysisFields = analysisRecord.text.split(',').map(function (field) { return field.trim(); });
        var type = analysisFields[0].toLowerCase();
        if (type === 'static' || type === 'modal') {
          if (analysisFields.length !== 1) fail(analysisRecord.line, 'Analysis ' + type + ' requires exactly one field.');
          model.analysis.type = type;
        } else if (type === 'transient') {
          if (analysisFields.length !== 5) fail(analysisRecord.line, 'Transient analysis requires five fields.');
          var settings = analysisFields.slice(1).map(Number);
          if (!settings.every(Number.isFinite)) fail(analysisRecord.line, 'Transient settings must be finite numbers.');
          model.analysis = { type: type, timeStep: settings[0], duration: settings[1], monitorNode: settings[2], monitorDOF: settings[3] };
        } else fail(analysisRecord.line, 'Unsupported analysis type "' + type + '".');
        hasAnalysis = true;
      } else if (marker === 'nodes') {
        if (hasNodes) fail(markerLine, 'The nodes section may appear only once.');
        var nodeCount = count('nodes', true);
        if (nodeCount > M.config.maxNodes) fail(markerLine, 'Node count exceeds the browser limit.');
        for (var n = 0; n < nodeCount; n += 1) {
          var nodeRecord = nextData('node coordinates');
          var nodeValues = parseRecord(nodeRecord.text, 3, nodeRecord.line, 'Node coordinates');
          model.nodes.push({ x: nodeValues[0], y: nodeValues[1] });
        }
        hasNodes = true;
      } else if (marker === 'elems_112' || marker === 'elems_113') {
        requireNodes(markerLine);
        var elementType = Number(marker.slice(-3));
        if (hasElements && model.elementType !== elementType) fail(markerLine, 'Mixed truss/frame meshes are not supported.');
        model.elementType = elementType;
        var elementCount = count(marker, true);
        if (model.elements.length + elementCount > M.config.maxElements) fail(markerLine, 'Element count exceeds the browser limit.');
        for (var e = 0; e < elementCount; e += 1) {
          var elementRecord = nextData('element');
          var expected = elementType === 112 ? 6 : 7;
          var elementValues = parseRecord(elementRecord.text, expected, elementRecord.line, 'Element');
          if (elementValues[0] !== elementType) fail(elementRecord.line, 'Element type does not match its section.');
          model.elements.push({
            node1: elementValues[1], node2: elementValues[2], area: elementValues[3],
            youngsModulus: elementValues[4], density: elementValues[5],
            momentOfInertia: elementType === 113 ? elementValues[6] : null
          });
        }
        hasElements = true;
      } else if (marker === 'releases') {
        requireElements(markerLine);
        var releaseCount = count('releases',false);
        for (var ri=0;ri<releaseCount;ri++) {
          var releaseRecord=nextData('release Mz'), releaseFields=releaseRecord.text.split(',').map(function (s) { return s.trim(); });
          if (releaseFields.length !== 3) fail(releaseRecord.line,'Expected elementId,end,Mz.');
          model.releases.push({elementId:Number(releaseFields[0]),end:Number(releaseFields[1]),component:releaseFields[2],sourceLine:releaseRecord.line});
        }
      } else if (marker === 'mpc') {
        requireElements(markerLine);
        var mpcCount = count('mpc', false);
        for (var mi=0; mi<mpcCount; mi++) {
          var record=nextData('MPC'), parts=record.text.split(','), masters=Number(parts[3]);
          if (!integer(masters) || masters < 1 || parts.length !== 4+3*masters) fail(record.line,'MPC masterCount does not match its fields.');
          var values=parseRecord(record.text,4+3*masters,record.line,'MPC');
          var mpc={depNode:values[0],depDOF:values[1],rhs:values[2],masters:[],sourceLine:record.line};
          for (var mt=0;mt<masters;mt++) mpc.masters.push({node:values[4+3*mt],dof:values[5+3*mt],coefficient:values[6+3*mt]});
          model.mpcs.push(mpc);
        }
      } else if (marker === 'bcfix') {
        requireElements(markerLine);
        var supportCount = count('bcfix', false);
        for (var s = 0; s < supportCount; s += 1) {
          var supportRecord = nextData('fixed condition');
          var supportValues = parseRecord(supportRecord.text, 5, supportRecord.line, 'Fixed condition');
          try { model.supports.push(global.MKEFSupports.fromType(supportValues[0], supportValues[1], model.elementType === 112 ? 2 : 3)); } catch (error) { fail(supportRecord.line, error.message); }
          if (supportValues.slice(2).some(function (value) { return value !== 0; })) fail(supportRecord.line, 'Only zero prescribed displacements are supported.');
        }
      } else if (marker === 'eload_uniform') {
        requireElements(markerLine);
        var elementLoadCount = count(marker, false);
        for (var el = 0; el < elementLoadCount; el += 1) {
          var elementLoadRecord = nextData('uniform element load');
          var v = parseRecord(elementLoadRecord.text, 5, elementLoadRecord.line, 'Uniform element load');
          if (v[0] !== 20 || !integer(v[1]) || v[1] < 1 || v[1] > model.elements.length || model.elementType !== 113 || [1, 2].indexOf(v[2]) < 0 || (v[3] === 0 && v[4] === 0)) fail(elementLoadRecord.line, 'Invalid uniform element load; requires frame 113, existing element, Local/Global and nonzero intensity.');
          model.elementLoads.push({ type: 20, elementId: v[1], coordinateSystem: v[2], qx: v[3], qy: v[4] });
        }
      } else if (['bcforce_stat', 'bcforce_harm', 'bcforce_pulse', 'bcforce_step'].indexOf(marker) >= 0) {
        requireElements(markerLine);
        var loadTypes = { bcforce_stat: 10, bcforce_harm: 11, bcforce_pulse: 12, bcforce_step: 13 };
        var loadType = loadTypes[marker], loadCount = count(marker, false);
        for (var l = 0; l < loadCount; l += 1) {
          var loadRecord = nextData('nodal load');
          var loadValues = parseRecord(loadRecord.text, loadType === 11 ? 6 : 5, loadRecord.line, 'Nodal load');
          if (loadValues[0] !== loadType) fail(loadRecord.line, 'Load type does not match its section.');
          model.loads.push({ type: loadType, node: loadValues[1], fx: loadValues[2], fy: loadValues[3], mz: loadValues[4], frequency: loadType === 11 ? loadValues[5] : null });
        }
      } else fail(markerLine, 'Unknown section marker "' + marker + '".');
    }

    if (!hasNodes) throw new Error('The input file contains no nodes section.');
    if (!hasElements) throw new Error('The input file contains no element section.');
    var errors = validateModel(model, { requireAnalysis: false });
    if (errors.length) throw new Error(errors[0]);
    return model;
  }

  function validateModel(model, options) {
    options = options || {};
    var errors = [], nodeCount = model.nodes.length, dofPerNode = model.elementType === 112 ? 2 : 3;
    function add(message) { errors.push(message); }
    if (options.requireAnalysis !== false && ['static', 'modal', 'transient'].indexOf(model.analysis.type) < 0) add('Select static, modal, or transient analysis.');
    if ([112, 113].indexOf(model.elementType) < 0) add('Select truss element 112 or frame element 113.');
    if (!nodeCount) add('Add at least one node.');
    if (!model.elements.length) add('Add at least one element.');
    model.nodes.forEach(function (node, index) {
      if (!finite(node.x) || !finite(node.y)) add('Node ' + (index + 1) + ' coordinates must be finite.');
    });
    model.elements.forEach(function (element, index) {
      var label = 'Element ' + (index + 1);
      if (!integer(element.node1) || !integer(element.node2) || element.node1 < 1 || element.node1 > nodeCount || element.node2 < 1 || element.node2 > nodeCount) add(label + ' has an invalid node reference.');
      else if (element.node1 === element.node2) add(label + ' must connect two distinct nodes.');
      else {
        var first = model.nodes[element.node1 - 1], second = model.nodes[element.node2 - 1];
        if (first.x === second.x && first.y === second.y) add(label + ' has zero length.');
      }
      ['area', 'youngsModulus', 'density'].forEach(function (field) { if (!finite(element[field]) || element[field] <= 0) add(label + ' ' + field + ' must be positive and finite.'); });
      if (model.elementType === 113 && (!finite(element.momentOfInertia) || element.momentOfInertia <= 0)) add(label + ' momentOfInertia must be positive and finite.');
    });
    var constrained = {};
    model.supports.forEach(function (support, index) {
      if (!integer(support.node) || support.node < 1 || support.node > nodeCount) add('Support ' + (index + 1) + ' has an invalid node.');
      try {
        global.MKEFSupports.dofs(support, dofPerNode).forEach(function (dof) {
          if (dof === 3 && global.MKEFReleases.inactive(model,support.node)) return;
          var key = support.node + ':' + dof;
          if (constrained[key]) add('Support ' + (index + 1) + ' duplicates node ' + support.node + ' DOF ' + ['ux','uy','thetaZ'][dof - 1] + '.');
          constrained[key] = true;
        });
      } catch (error) { add(error.message); }
    });
    errors = errors.concat(global.MKEFMPC.validate(model),global.MKEFReleases.validate(model));
    model.loads.forEach(function (load, index) {
      var label = 'Load ' + (index + 1);
      if ([10, 11, 12, 13].indexOf(load.type) < 0) add(label + ' has an unsupported type.');
      if (!integer(load.node) || load.node < 1 || load.node > nodeCount) add(label + ' has an invalid node.');
      if (![load.fx, load.fy, load.mz].every(finite)) add(label + ' components must be finite.');
      if (model.elementType === 112 && load.mz !== 0) add(label + ' gives Mz to a truss node.');
      if (load.type === 11 && (!finite(load.frequency) || load.frequency <= 0)) add(label + ' harmonic frequency must be positive and finite.');
      if (model.analysis.type === 'static' && load.type !== 10) add(label + ' is not valid for static analysis.');
      if (model.analysis.type === 'modal') add(label + ' is not valid for modal analysis.');
    });
    (model.elementLoads || []).forEach(function (load, index) {
      var label = 'Element load ' + (index + 1);
      if (load.type !== 20 || model.elementType !== 113) add(label + ' requires uniform load type 20 and frame element 113.');
      if (!integer(load.elementId) || load.elementId < 1 || load.elementId > model.elements.length) add(label + ' has an invalid element.');
      if ([1, 2].indexOf(load.coordinateSystem) < 0) add(label + ' requires Local or Global coordinates.');
      if (![load.qx, load.qy].every(finite) || (load.qx === 0 && load.qy === 0)) add(label + ' requires finite, nonzero intensity.');
      if (model.analysis.type && model.analysis.type !== 'static') add(label + ' is valid only for static analysis.');
    });
    if (model.analysis.type === 'transient') {
      if (!finite(model.analysis.timeStep) || model.analysis.timeStep <= 0) add('Transient time step must be positive and finite.');
      if (!finite(model.analysis.duration) || model.analysis.duration < model.analysis.timeStep) add('Transient duration must cover at least one time step.');
      if (!integer(model.analysis.monitorNode) || model.analysis.monitorNode < 1 || model.analysis.monitorNode > nodeCount) add('Transient monitor node is outside the model.');
      if (!integer(model.analysis.monitorDOF) || model.analysis.monitorDOF < 1 || model.analysis.monitorDOF > dofPerNode) add('Transient monitor DOF is outside the element type.');
    }
    return errors;
  }

  function serializeCase(model) {
    var errors = validateModel(model);
    if (errors.length) throw new Error(errors[0]);
    var out = ['analysis'];
    if (model.analysis.type === 'transient') out.push(['transient', numberText(model.analysis.timeStep), numberText(model.analysis.duration), numberText(model.analysis.monitorNode), numberText(model.analysis.monitorDOF)].join(','));
    else out.push(model.analysis.type);
    out.push('nodes', String(model.nodes.length));
    model.nodes.forEach(function (node) { out.push([numberText(node.x), numberText(node.y), '0'].join(',')); });
    out.push('elems_' + model.elementType, String(model.elements.length));
    model.elements.forEach(function (element) {
      var values = [model.elementType, element.node1, element.node2, element.area, element.youngsModulus, element.density];
      if (model.elementType === 113) values.push(element.momentOfInertia);
      out.push(values.map(numberText).join(','));
    });
    if ((model.releases || []).length) {
      var releases=model.releases.slice().sort(function (a,b) { return a.elementId-b.elementId || a.end-b.end; });
      out.push('releases',String(releases.length));
      releases.forEach(function (r) { out.push([r.elementId,r.end,'Mz'].join(',')); });
    }
    if ((model.mpcs || []).length) {
      out.push('mpc', String(model.mpcs.length));
      model.mpcs.forEach(function (m) {
        var values=[m.depNode,m.depDOF,m.rhs,m.masters.length];
        m.masters.forEach(function (a) { values.push(a.node,a.dof,a.coefficient); });
        out.push(values.map(numberText).join(','));
      });
    }
    out.push('bcfix', String(model.supports.length));
    model.supports.forEach(function (support) { out.push([global.MKEFSupports.type(support, model.elementType === 112 ? 2 : 3), support.node, 0, 0, 0].join(',')); });
    [{ type: 10, marker: 'bcforce_stat' }, { type: 11, marker: 'bcforce_harm' }, { type: 12, marker: 'bcforce_pulse' }, { type: 13, marker: 'bcforce_step' }].forEach(function (group) {
      var loads = model.loads.filter(function (load) { return load.type === group.type; });
      if (!loads.length) return;
      out.push(group.marker, String(loads.length));
      loads.forEach(function (load) {
        var values = [load.type, load.node, load.fx, load.fy, load.mz];
        if (load.type === 11) values.push(load.frequency);
        out.push(values.map(numberText).join(','));
      });
    });
    if ((model.elementLoads || []).length) {
      out.push('eload_uniform', String(model.elementLoads.length));
      model.elementLoads.forEach(function (load) { out.push([20, load.elementId, load.coordinateSystem, load.qx, load.qy].map(numberText).join(',')); });
    }
    return out.join('\n') + '\n';
  }

  M.caseFormat = { newModel: newModel, parse: parseCase, validate: validateModel, serialize: serializeCase };
}(window));
