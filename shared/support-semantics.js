(function (global) {
  'use strict';
  var S = global.MKEFSupports = {}, names = ['ux', 'uy', 'thetaZ'];
  S.fromType = function (type, node, dofs) {
    var rows = global.MKEFSupportDefinitions[dofs === 2 ? 'truss' : 'frame'];
    if (!Number.isInteger(type) || !rows[type - 1]) throw new Error('Unsupported support type ' + type + ' for ' + dofs + ' DOFs.');
    var m = rows[type - 1];
    return {node: node, fixUx: !!m[0], fixUy: !!m[1], fixThetaZ: !!m[2]};
  };
  S.mask = function (support, dofs) {
    var s = Object.prototype.hasOwnProperty.call(support, 'fixUx') ? support : S.fromType(support.type, support.node, dofs);
    var m = [s.fixUx, s.fixUy, s.fixThetaZ];
    if (!m.every(function (v) { return typeof v === 'boolean'; }) || !m.some(Boolean) || (dofs === 2 && m[2])) throw new Error('Invalid fixed DOFs at node ' + s.node + '.');
    return m;
  };
  S.type = function (support, dofs) {
    var mask = S.mask(support, dofs), rows = global.MKEFSupportDefinitions[dofs === 2 ? 'truss' : 'frame'];
    return rows.findIndex(function (row) { return row.every(function (v, i) { return !!v === mask[i]; }); }) + 1;
  };
  S.dofs = function (support, count) { return S.mask(support, count).map(function (v, i) { return v ? i + 1 : 0; }).filter(Boolean); };
  S.label = function (support, count) { return 'Fixed: ' + S.dofs(support, count).map(function (i) { return names[i - 1]; }).join(', '); };
  S.changeFamily = function (model, family) {
    var oldCount = model.elementType === 112 ? 2 : 3;
    var supports = model.supports.map(function (s) { var m = S.mask(s, oldCount); return {node:s.node, fixUx:m[0], fixUy:m[1], fixThetaZ:m[2]}; });
    if (family === 112 && supports.some(function (s) { return s.fixThetaZ; })) throw new Error('Remove thetaZ restraints before changing to a truss.');
    model.supports = supports; model.elementType = family;
  };
}(window));
