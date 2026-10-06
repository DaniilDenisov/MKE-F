(function (global) {
  'use strict';
  var R = global.MKEFReleases = {};
  R.has = function (model, id, end) { return (model.releases || []).some(function (r) { return r.elementId === id && r.end === end && r.component === 'Mz'; }); };
  R.inactive = function (model, node) {
    if (model.elementType !== 113) return false;
    var attached = false, rigid = false;
    model.elements.forEach(function (e, i) {
      [e.node1, e.node2].forEach(function (n, j) { if (n === node) { attached = true; if (!R.has(model, i + 1, j + 1)) rigid = true; } });
    });
    return attached && !rigid;
  };
  R.validate = function (model) {
    var errors = [], seen = new Set();
    (model.releases || []).forEach(function (r) {
      var label = 'Release ' + r.component + (r.sourceLine ? ' at line ' + r.sourceLine : '') + ': ';
      if (model.elementType !== 113 || !Number.isInteger(r.elementId) || r.elementId < 1 || r.elementId > model.elements.length || [1,2].indexOf(r.end) < 0 || r.component !== 'Mz') errors.push(label + 'expected frame element, end 1 or 2, Mz.');
      var key = r.elementId + ':' + r.end;
      if (seen.has(key)) errors.push(label + 'duplicate element end.');
      seen.add(key);
    });
    (model.loads || []).forEach(function (l) { if (l.mz !== 0 && R.inactive(model,l.node)) errors.push('Node ' + l.node + ': cannot apply Mz to absent thetaZ.'); });
    (model.mpcs || []).forEach(function (m,i) {
      [{node:m.depNode,dof:m.depDOF}].concat(m.masters).forEach(function (a) { if (a.dof === 3 && R.inactive(model,a.node)) errors.push('MPC ' + (i+1) + ': node ' + a.node + ' has no shared thetaZ.'); });
    });
    if (model.analysis && model.analysis.type === 'transient' && model.analysis.monitorDOF === 3 && R.inactive(model,model.analysis.monitorNode)) errors.push('Monitor node ' + model.analysis.monitorNode + ' has no shared thetaZ.');
    return errors;
  };
  R.warnings = function (model) {
    return model.supports.filter(function (s) { return R.inactive(model,s.node) && global.MKEFSupports.mask(s,3)[2]; }).map(function (s) { return 'Node ' + s.node + ': thetaZ restraint ignored because all connected ends release Mz.'; });
  };
  R.effectiveSupport = function (support, inactive, dofs) {
    var mask = global.MKEFSupports.mask(support,dofs);
    if (inactive) mask[2] = false;
    if (!mask.some(Boolean)) return null;
    return global.MKEFSupports.type({node:support.node,fixUx:mask[0],fixUy:mask[1],fixThetaZ:mask[2]},dofs);
  };
  R.append = function (create, layer, first, second, end, size, attrs) {
    var dx=second.x-first.x, dy=second.y-first.y, length=Math.hypot(dx,dy);
    if (!length) return;
    var offset=Math.min(size*1.1,length*.15), fraction=end === 1 ? offset/length : 1-offset/length;
    var circle=create('circle',Object.assign({cx:first.x+fraction*dx,cy:first.y+fraction*dy,r:Math.min(size*.22,length*.04),fill:'white',stroke:'#92400e','stroke-width':size*.08,class:'end-release','pointer-events':'none'},attrs || {}));
    var title=create('title'); title.textContent='Mz'+end+' released'; circle.appendChild(title); layer.appendChild(circle);
    return circle;
  };
}(window));
