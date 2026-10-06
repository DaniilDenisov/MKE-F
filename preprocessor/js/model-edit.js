(function (M) {
  'use strict';
  var E = M.modelEdit = {};
  function counts(model) { return [model.nodes.length, model.elements.length, model.supports.length, model.loads.length, (model.elementLoads || []).length, (model.mpcs || []).length]; }
  E.describe = function (before, model) {
    var after = counts(model), labels = ['nodes', 'elements', 'supports', 'nodal loads', 'element loads', 'MPCs'];
    var removed = before.map(function (value, index) { return value > after[index] ? (value - after[index]) + ' ' + labels[index] : ''; }).filter(Boolean);
    return removed.length ? 'Removed: ' + removed.join(', ') + '. Undo restores the entire operation.' : '';
  };
  E.counts = counts;
  function removeElements(model, ids) {
    var newIds = {}, kept = [];
    model.elements.forEach(function (element, index) { if (!ids.has(index + 1)) { kept.push(element); newIds[index + 1] = kept.length; } });
    model.elements = kept;
    model.releases = (model.releases || []).filter(function (r) { return !ids.has(r.elementId); });
    model.releases.forEach(function (r) { r.elementId = newIds[r.elementId]; });
    model.elementLoads = (model.elementLoads || []).filter(function (load) { return !ids.has(load.elementId); });
    model.elementLoads.forEach(function (load) { load.elementId = newIds[load.elementId]; });
  }
  E.deleteElement = function (model, index) { removeElements(model, new Set([index + 1])); };
  E.setRelease = function (model, index, end, released) {
    model.releases = (model.releases || []).filter(function (r) { return r.elementId !== index+1 || r.end !== end; });
    if (released) model.releases.push({elementId:index+1,end:end,component:'Mz'});
  };
  E.deleteNode = function (model, index) {
    var id = index + 1, ids = new Set();
    model.elements.forEach(function (element, i) { if (element.node1 === id || element.node2 === id) ids.add(i + 1); });
    removeElements(model, ids);
    model.nodes.splice(index, 1);
    model.elements.forEach(function (element) { if (element.node1 > id) element.node1 -= 1; if (element.node2 > id) element.node2 -= 1; });
    ['supports', 'loads'].forEach(function (key) {
      model[key] = model[key].filter(function (item) { return item.node !== id; });
      model[key].forEach(function (item) { if (item.node > id) item.node -= 1; });
    });
    model.mpcs = (model.mpcs || []).filter(function (m) { return m.depNode !== id && !m.masters.some(function (a) { return a.node === id; }); });
    model.mpcs.forEach(function (m) { if (m.depNode > id) m.depNode--; m.masters.forEach(function (a) { if (a.node > id) a.node--; }); });
    if (model.analysis.monitorNode === id) model.analysis.monitorNode = 0;
    else if (model.analysis.monitorNode > id) model.analysis.monitorNode -= 1;
  };
  E.updateElement = function (model, index, field, value) {
    if (Object.is(model.elements[index][field], value)) return false;
    model.elements[index][field] = value;
    model.elementLoads = (model.elementLoads || []).filter(function (load) { return load.elementId !== index + 1; });
    return true;
  };
  E.updateNode = function (model, index, field, value) {
    if (Object.is(model.nodes[index][field], value)) return false;
    model.nodes[index][field] = value;
    var affected = new Set();
    model.elements.forEach(function (element, i) { if (element.node1 === index + 1 || element.node2 === index + 1) affected.add(i + 1); });
    model.elementLoads = (model.elementLoads || []).filter(function (load) { return !affected.has(load.elementId); });
    return true;
  };
  E.changeLoadType = function (model, index, type) {
    var load=model.elementLoads[index];
    if (load.type===type) return;
    if (type===21) {
      load.qx1=load.qx2=load.qx; load.qy1=load.qy2=load.qy;
      delete load.qx; delete load.qy;
    } else {
      load.qx=(load.qx1+load.qx2)/2; load.qy=(load.qy1+load.qy2)/2;
      ['qx1','qy1','qx2','qy2'].forEach(function (key) { delete load[key]; });
    }
    load.type=type;
  };
  E.changeCoordinates = function (model, index, system) {
    var load = model.elementLoads[index];
    if (load.coordinateSystem === system) return;
    var element = model.elements[load.elementId - 1];
    var a = element && model.nodes[element.node1 - 1], b = element && model.nodes[element.node2 - 1];
    var length = a && b ? Math.hypot(b.x - a.x, b.y - a.y) : NaN;
    var fields = load.type === 21 ? [['qx1','qy1'],['qx2','qy2']] : [['qx','qy']];
    if (!Number.isFinite(length) || length <= 0 || !fields.every(function (pair) { return pair.every(function (key) { return Number.isFinite(load[key]); }); })) throw new Error('Valid element geometry and finite components are required to change coordinates.');
    var c = (b.x - a.x)/length, s = (b.y - a.y)/length;
    fields.forEach(function (pair) {
      var x=load[pair[0]], y=load[pair[1]];
      load[pair[0]] = system === 1 ? c*x + s*y : c*x - s*y;
      load[pair[1]] = system === 1 ? -s*x + c*y : s*x + c*y;
    });
    load.coordinateSystem = system;
  };
}(window.MKEFPre));
