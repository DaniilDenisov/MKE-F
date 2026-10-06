(function (global) {
  'use strict';
  var Q = global.MKEFMPC = {}, labels = ['ux','uy','thetaZ'];
  function key(node,dof) { return node + ':' + dof; }
  Q.validate = function (model) {
    var equations = model.mpcs || [], dofs = model.elementType === 112 ? 2 : 3;
    var fixed = new Set(), owner = new Map(), children = [], degree = [], errors = [];
    (model.supports || []).forEach(function (s) {
      try { global.MKEFSupports.dofs(s,dofs).forEach(function (d) { fixed.add(key(s.node,d)); }); } catch (_) { /* SPC validator reports this. */ }
    });
    var nodeIds = new Set(model.nodes.map(function (node,i) { return node.id === undefined ? i+1 : node.id; }));
    function valid(node,dof) { return Number.isInteger(node) && nodeIds.has(node) && Number.isInteger(dof) && dof >= 1 && dof <= dofs; }
    equations.forEach(function (m,i) {
      function fail(message) { errors.push('MPC ' + (i+1) + (m.sourceLine ? ' (source line ' + m.sourceLine + ')' : '') + ': ' + message); }
      var dep = key(m.depNode,m.depDOF);
      if (!valid(m.depNode,m.depDOF)) fail('Invalid dependent node or DOF.');
      if (m.rhs !== 0) fail('Only rhs = 0 is supported.');
      if (fixed.has(dep)) fail('Dependent DOF is also restrained by an SPC.');
      if (owner.has(dep)) fail('Dependent DOF occurs in more than one MPC.');
      owner.set(dep,i); children[i] = []; degree[i] = 0;
      if (!Array.isArray(m.masters) || !m.masters.length) { fail('Expected nonzero masters.'); return; }
      var seen = new Set(), nonzero = false;
      m.masters.forEach(function (a) {
        if (!valid(a.node,a.dof) || !Number.isFinite(a.coefficient)) fail('Invalid master node, DOF or coefficient.');
        var k = key(a.node,a.dof);
        if (k === dep) fail('Dependent DOF cannot be its own master.');
        if (seen.has(k)) fail('Repeated master DOF.');
        seen.add(k); if (Number.isFinite(a.coefficient) && a.coefficient !== 0) nonzero = true;
      });
      if (!nonzero) fail('Expected nonzero masters.');
    });
    if (errors.length) return errors;
    equations.forEach(function (m,i) { m.masters.forEach(function (a) {
      var parent = owner.get(key(a.node,a.dof));
      if (parent !== undefined && a.coefficient !== 0) { children[parent].push(i); degree[i]++; }
    }); });
    var queue = [], head = 0;
    degree.forEach(function (d,i) { if (!d) queue.push(i); });
    while (head < queue.length) children[queue[head++]].forEach(function (i) { if (--degree[i] === 0) queue.push(i); });
    if (head !== equations.length) errors.push('MPC ' + (degree.findIndex(function (d) { return d > 0; })+1) + ': Cycle in dependency graph.');
    return errors;
  };
  Q.equation = function (m) {
    return labels[m.depDOF-1] + '(' + m.depNode + ') = ' + m.masters.map(function (a) { return a.coefficient + '·' + labels[a.dof-1] + '(' + a.node + ')'; }).join(' + ');
  };
  Q.onAxis = function (model,i,j,k) {
    if (![i,j,k].every(function (n) { return Number.isInteger(n) && n >= 1 && n <= model.nodes.length; }) || new Set([i,j,k]).size !== 3) throw new Error('Choose three distinct existing nodes.');
    var a=model.nodes[i-1], b=model.nodes[j-1], p=model.nodes[k-1], dx=b.x-a.x, dy=b.y-a.y, L=Math.hypot(dx,dy);
    if (!Number.isFinite(L) || L <= 0) throw new Error('Axis length must be positive.');
    var nx=-dy/L, ny=dx/L, xi=((p.x-a.x)*dx+(p.y-a.y)*dy)/(L*L);
    if (!Number.isFinite(xi) || xi <= 0 || xi >= 1 || Math.abs(nx*(p.x-a.x)+ny*(p.y-a.y)) > 1e-9*L) throw new Error('Intermediate node must lie inside the axis segment (tolerance 1e-9 L).');
    var normal=[nx,ny], candidates=Math.abs(nx)>=Math.abs(ny)?[1,2]:[2,1];
    for (var c=0;c<candidates.length;c++) {
      var d=candidates[c], divisor=normal[d-1]; if (divisor === 0) continue;
      var terms=[];
      [1,2].forEach(function (component) {
        var value=normal[component-1]/divisor;
        if (value !== 0) {
          terms.push({node:i,dof:component,coefficient:(1-xi)*value});
          terms.push({node:j,dof:component,coefficient:xi*value});
          if (component !== d) terms.push({node:k,dof:component,coefficient:-value});
        }
      });
      var m={depNode:k,depDOF:d,rhs:0,masters:terms};
      var trial=Object.assign({},model,{mpcs:(model.mpcs||[]).concat([m])});
      if (!Q.validate(trial).length) return m;
    }
    throw new Error('No available dependent DOF: check SPCs and MPC dependencies.');
  };
}(window));
