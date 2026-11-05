(function (global) {
  'use strict';
  var C = global.MKEFSpatialConditions = {}, labels = ['ux','uy','uz','thetaX','thetaY','thetaZ'];
  C.keys = ['supports','loads','elementLoads','mpcs','releases']; C.labels = labels;
  var dot = function (a,b) { return a.reduce(function (s,v,i) { return s+v*b[i]; },0); };
  var norm = function (a) { return Math.hypot.apply(null,a); };
  var scale = function (a,s) { return a.map(function (v) { return v*s; }); };
  var add = function (a,b,s) { return a.map(function (v,i) { return v+s*b[i]; }); };
  var unit = function (a) { return scale(a,1/norm(a)); };
  var cross = function (a,b) { return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]; };
  C.basis = function (m,e) {
    var a=m.nodes[e.node1-1],b=m.nodes[e.node2-1],x=unit([b.x-a.x,b.y-a.y,b.z-a.z]),ref=e.referenceY;
    if (!ref) { var axes=[1,2,0].sort(function (i,j) { return Math.abs(x[i])-Math.abs(x[j]); }); ref=[0,0,0];ref[axes[0]]=1; }
    ref=unit(ref);var y=unit(add(ref,x,-dot(ref,x))),z=cross(x,y),angle=(e.angleDegrees || 0)%360*Math.PI/180;
    return [x,add(scale(y,Math.cos(angle)),z,Math.sin(angle)),add(scale(y,-Math.sin(angle)),z,Math.cos(angle))];
  };
  C.toLocal = function (R,v) { return R.map(function (row) { return dot(row,v); }); };
  C.toGlobal = function (R,v) { return [0,1,2].map(function (j) { return R.reduce(function (s,row,i) { return s+row[j]*v[i]; },0); }); };
  C.transformLoadComponents = function (components,basis,target) {
    var result=[],transform=target===1 ? C.toLocal : C.toGlobal;
    for(var i=0;i<components.length;i+=3) result.push.apply(result,transform(basis,components.slice(i,i+3)));
    return result;
  };
  // One-sided Jacobi SVD of a tall matrix with three columns. Do not use A'A:
  // its squared condition number loses the solver's 1e-10 singular-value cutoff.
  C.activeBasis = function (rows) {
    if (!rows.length) return [];
    var b=rows.map(function (r) { return r.slice(); }),V=[[1,0,0],[0,1,0],[0,0,1]];
    for (var sweep=0;sweep<30;sweep++) {
      var changed=false;
      [[0,1],[0,2],[1,2]].forEach(function (pair) {
        var p=pair[0],q=pair[1],aa=0,bb=0,ab=0;
        b.forEach(function (r) { aa+=r[p]*r[p];bb+=r[q]*r[q];ab+=r[p]*r[q]; });
        if (Math.abs(ab)<=Number.EPSILON*Math.sqrt(aa*bb)*4 || ab===0) return;
        var tau=(bb-aa)/(2*ab),t=(tau>=0 ? 1 : -1)/(Math.abs(tau)+Math.hypot(1,tau)),c=1/Math.hypot(1,t),s=t*c;
        b.concat(V).forEach(function (r) { var x=r[p],y=r[q];r[p]=c*x-s*y;r[q]=s*x+c*y; }); changed=true;
      });
      if (!changed) break;
    }
    var singular=[0,1,2].map(function (j) { return Math.sqrt(b.reduce(function (s,r) { return s+r[j]*r[j]; },0)); }),maximum=Math.max.apply(null,singular);
    var ids=[0,1,2].filter(function (j) { return singular[j]>1e-10*maximum; }).sort(function (i,j) { return singular[j]-singular[i]; });
    return ids.length===3 ? [[1,0,0],[0,1,0],[0,0,1]] : ids.map(function (j) { return V.map(function (r) { return r[j]; }); });
  };
  C.validate = function (m) {
    var errors=[],warnings=[],inactiveSupports=[],info=global.MKEFSpatial.registry.types.find(function (t) { return t.code===m.elementType; });
    var frame=info.family==='frame',dofs=info.dofLabels.length,rotation=m.nodes.map(function () { return []; }),bases=[],released=new Set();
    var validNode=function (n) { return Number.isInteger(n)&&n>=1&&n<=m.nodes.length; };
    var validElement=function (n) { return Number.isInteger(n)&&n>=1&&n<=m.elements.length; };
    var validDOF=function (d) { return Number.isInteger(d)&&d>=1&&d<=dofs; };
    var finiteArray=function (a,n) { return Array.isArray(a)&&a.length===n&&a.every(Number.isFinite); };
    function fail(key,i,text) { errors.push(key+' '+(i+1)+': '+text); }
    if (C.keys.reduce(function (s,k) { return s+m[k].length; },0)>100000) errors.push('Limit: 100000 condition records.');
    m.releases.forEach(function (r,i) {
      var key=r.elementId+':'+r.end+':'+r.component;
      if (!frame||!validElement(r.elementId)||[1,2].indexOf(r.end)<0||['Mx','My','Mz'].indexOf(r.component)<0) fail('releases',i,'Invalid frame, end or component.');
      if (released.has(key)) fail('releases',i,'Duplicate release.');released.add(key);
    });
    m.elements.forEach(function (e,i) {
      if (!frame) return;
      var R=C.basis(m,e); bases.push(R);
      if (frame) [1,2].forEach(function (end) { var node=end===1 ? e.node1 : e.node2; ['Mx','My','Mz'].forEach(function (axis,j) { if (!released.has((i+1)+':'+end+':'+axis)) rotation[node-1].push(R[j]); }); });
    });
    var active=rotation.map(function (rows) { return frame ? C.activeBasis(rows) : []; });
    function physical(node,dof) {
      var row=new Map(); if (!validNode(node)||!validDOF(dof)) return row;
      if (dof<=3) row.set((node-1)*6+dof-1,1);
      else active[node-1].forEach(function (v,j) { if (v[dof-4]!==0) row.set((node-1)*6+3+j,v[dof-4]); });
      return row;
    }
    function rowNorm(row) { return Math.hypot.apply(null,Array.from(row.values())); }
    var fixed=new Set(),supportRows=[];
    m.supports.forEach(function (r,i) {
      if (!validNode(r.node)||!finiteArray(r.mask,6)||r.mask.some(function (v) { return v!==0&&v!==1; })||!r.mask.some(Boolean)||!frame&&r.mask.slice(3).some(Boolean)) { fail('supports',i,'Invalid node or nonempty component mask.');return; }
      r.mask.forEach(function (v,j) { if (!v) return; var key=r.node+':'+(j+1);
        if (fixed.has(key)) fail('supports',i,'Duplicate restraint '+labels[j]+' at node '+r.node); fixed.add(key);
        var row=physical(r.node,j+1);
        if (rowNorm(row)<1e-12) { warnings.push('supports '+(i+1)+': '+labels[j]+' at node '+r.node+' is inactive; restraint is ignored by the solver.');inactiveSupports.push({index:i,component:j}); }
        else supportRows.push({id:(r.node-1)*6+j,row:row});
      });
    });
    m.loads.forEach(function (r,i) {
      if (!validNode(r.node)||!finiteArray(r.components,6)||[10,11,12,13].indexOf(r.type)<0||r.type===11&&(!Number.isFinite(r.frequency)||r.frequency<=0)) { fail('loads',i,'Invalid node, components, law or frequency.');return; }
      if (m.analysis.type==='modal'||m.analysis.type==='static'&&r.type!==10||m.analysis.type==='transient'&&r.type===10) fail('loads',i,'Load law conflicts with analysis '+m.analysis.type+'.');
      var moment=r.components.slice(3),projected=[0,0,0]; active[r.node-1].forEach(function (q) { projected=add(projected,q,dot(moment,q)); });
      if (!frame&&moment.some(Boolean)||norm(add(moment,projected,-1))>1e-10*norm(moment)) fail('loads',i,'Moment at node '+r.node+' acts in an inactive rotation direction.');
    });
    m.elementLoads.forEach(function (r,i) {
      if (!frame||!validElement(r.elementId)||[20,21].indexOf(r.type)<0||[1,2].indexOf(r.coordinateSystem)<0||!finiteArray(r.components,r.type===20?3:6)||!r.components.some(Boolean)) fail('elementLoads',i,'Expected a frame, axes and nonzero finite intensities.');
      if (m.analysis.type&&m.analysis.type!=='static') fail('elementLoads',i,'Distributed loads are static only.');
    });
    var owner=new Map(),rows=[],children=m.mpcs.map(function () { return []; }),degree=m.mpcs.map(function () { return 0; });
    m.mpcs.forEach(function (r,i) {
      var dep=r.depNode+':'+r.depDOF;
      if (!validNode(r.depNode)||!validDOF(r.depDOF)||r.rhs!==0) fail('mpcs',i,'Invalid dependent component or nonzero rhs.');
      if (fixed.has(dep)) fail('mpcs',i,'Dependent component is restrained by an SPC.');
      if (owner.has(dep)) fail('mpcs',i,'Duplicate dependent component.'); owner.set(dep,i);
      if (!Array.isArray(r.masters)||!r.masters.length) { fail('mpcs',i,'Expected nonzero masters.'); return; }
      var seen=new Set(),row=physical(r.depNode,r.depDOF),nonzero=false,referenceScale=0;row.forEach(function(v){referenceScale=Math.max(referenceScale,Math.abs(v));});
      if (rowNorm(row)<1e-12) fail('mpcs',i,'Dependent component is inactive.');
      r.masters.forEach(function (a) {
        var key=a.node+':'+a.dof;
        if (!validNode(a.node)||!validDOF(a.dof)||!Number.isFinite(a.coefficient)) fail('mpcs',i,'Invalid master component or coefficient.');
        if (key===dep||seen.has(key)) fail('mpcs',i,'Self-reference or repeated master component.'); seen.add(key);
        var master=physical(a.node,a.dof); if (rowNorm(master)<1e-12) fail('mpcs',i,'Master component is inactive.');
        master.forEach(function(v){referenceScale=Math.max(referenceScale,Math.abs(a.coefficient*v));});
        if (a.coefficient!==0) nonzero=true;
        master.forEach(function (v,k) { row.set(k,(row.get(k)||0)-a.coefficient*v); });
      });
      var projectedPeak=0;row.forEach(function(v){projectedPeak=Math.max(projectedPeak,Math.abs(v));});
      if(projectedPeak<1e-10*referenceScale)fail('mpcs',i,'Redundant after projection onto active rotations.');
      if (!nonzero) fail('mpcs',i,'Expected nonzero masters.'); rows[i]=row;
    });
    if (!errors.length) {
      m.mpcs.forEach(function (r,i) { r.masters.forEach(function (a) { var parent=owner.get(a.node+':'+a.dof); if (parent!==undefined&&a.coefficient!==0) { children[parent].push(i);degree[i]++; } }); });
      var queue=[],head=0;degree.forEach(function (v,i) { if (!v) queue.push(i); });
      while(head<queue.length) children[queue[head++]].forEach(function (i) { if (--degree[i]===0) queue.push(i); });
      if (head!==m.mpcs.length) fail('mpcs',degree.findIndex(function (v) { return v>0; }),'Cycle in dependency graph.');
    }
    if (!errors.length&&rows.length) {
      // Sparse row elimination, in the solver's SPC-then-MPC order. Only visit
      // pivots present in a row; independent restraints do not cause quadratic scans.
      var pivots=new Map(),order=0;
      function eliminate(row,index) {
        row=new Map(row);var maximum=0;row.forEach(function (v) { maximum=Math.max(maximum,Math.abs(v)); });
        if (maximum) row.forEach(function (v,k) { row.set(k,v/maximum); });
        while(true) {
          var selected=null;row.forEach(function (v,k) { var pivot=pivots.get(k); if(v!==0&&pivot&&(!selected||pivot.order<selected.order)) selected=pivot; });
          if (!selected) break;var coefficient=row.get(selected.id);row.delete(selected.id);
          selected.row.forEach(function (v,k) { if(k!==selected.id) { var value=(row.get(k)||0)-coefficient*v; if(value===0) row.delete(k);else row.set(k,value); } });
        }
        var peak=0,id=-1;row.forEach(function (v,k) { if(Math.abs(v)>peak||Math.abs(v)===peak&&k<id) {peak=Math.abs(v);id=k;} });
        if (!Number.isFinite(peak)) { fail('mpcs',Math.max(0,index),'Constraint elimination overflow.');return; }
        if (peak<1e-10) { if(index>=0) fail('mpcs',index,'Redundant after projection onto active rotations.');return; }
        var divisor=row.get(id);row.forEach(function (v,k) { if(Math.abs(v)<1e-13*peak) row.delete(k);else row.set(k,v/divisor); });
        pivots.set(id,{id:id,row:row,order:order++});
      }
      supportRows.sort(function (a,b) { return a.id-b.id; }).forEach(function (r) { eliminate(r.row,-1); }); rows.forEach(eliminate);
    }
    var analysis=m.analysis;
    if (['','static','modal','transient'].indexOf(analysis.type)<0) errors.push('analysis: unsupported type.');
    if (analysis.type==='transient') {
      if (!Number.isFinite(analysis.timeStep)||analysis.timeStep<=0||!Number.isFinite(analysis.duration)||analysis.duration<analysis.timeStep||!validNode(analysis.monitorNode)||!validDOF(analysis.monitorDOF)) errors.push('analysis: invalid time settings or monitor component.');
      else if (rowNorm(physical(analysis.monitorNode,analysis.monitorDOF))<1e-12) errors.push('analysis: monitor rotation is inactive.');
    }
    return {errors:errors,warnings:warnings,bases:bases,active:active,inactiveSupports:inactiveSupports};
  };
  C.equalComponents=function(m,dependent,master,components){
    if(!components.length)throw new Error('Choose at least one component.');
    var rows=components.map(function(d){return {depNode:dependent,depDOF:d,rhs:0,masters:[{node:master,dof:d,coefficient:1}]};});
    var errors=C.validate(Object.assign({},m,{mpcs:m.mpcs.concat(rows)})).errors;if(errors.length)throw new Error(errors.join('\n'));return rows;
  };
  // Two transverse interpolation equations. Eliminate two physical columns together,
  // so the generated dependencies never reference each other. Prefer the middle node.
  C.axisHelper=function(m,start,end,middle){
    var ids=[middle,start,end];if(new Set(ids).size!==3||ids.some(function(id){return !Number.isInteger(id)||!m.nodes[id-1];}))throw new Error('Choose three distinct existing nodes.');
    var xyz=function(id){var n=m.nodes[id-1];return [n.x,n.y,n.z];},a=xyz(start),d=add(xyz(end),a,-1),L=norm(d),x=unit(d),relative=add(xyz(middle),a,-1),t=dot(relative,x)/L;
    if(!Number.isFinite(L)||!L||!(t>0&&t<1)||norm(add(relative,x,-t*L))>1e-9*L)throw new Error('The intermediate node must lie inside the axis segment (tolerance 1e-9 L).');
    var order=[0,1,2].sort(function(i,j){return Math.abs(x[i])-Math.abs(x[j]);}),unitAxis=[0,0,0];unitAxis[order[0]]=1;
    var n=unit(cross(x,unitAxis)),q=cross(x,n),weights=[1,-(1-t),-t],columns=[];
    ids.forEach(function(id,k){for(var j=0;j<3;j++)columns.push({node:id,dof:j+1,a:n[j]*weights[k],b:q[j]*weights[k]});});
    var blocked=new Set();m.supports.forEach(function(r){r.mask.forEach(function(v,j){if(v)blocked.add(r.node+':'+(j+1));});});m.mpcs.forEach(function(r){blocked.add(r.depNode+':'+r.depDOF);});
    var pairs=[];for(var i=0;i<9;i++)for(var j=i+1;j<9;j++){var u=columns[i],v=columns[j],det=u.a*v.b-v.a*u.b;if(!blocked.has(u.node+':'+u.dof)&&!blocked.has(v.node+':'+v.dof)&&Math.abs(det)>1e-12)pairs.push({i:i,j:j,det:det,priority:(i<3?0:1)+(j<3?0:1)});}
    pairs.sort(function(a,b){return a.priority-b.priority||Math.abs(b.det)-Math.abs(a.det)||a.i-b.i||a.j-b.j;});
    for(var p of pairs){var u=columns[p.i],v=columns[p.j];var rows=[u,v].map(function(dep,k){return {depNode:dep.node,depDOF:dep.dof,rhs:0,masters:columns.map(function(c,index){if(index===p.i||index===p.j)return null;var coef=k===0?(-v.b*c.a+v.a*c.b)/p.det:(u.b*c.a-u.a*c.b)/p.det;return coef===0?null:{node:c.node,dof:c.dof,coefficient:coef};}).filter(Boolean)};});
      if(!C.validate(Object.assign({},m,{mpcs:m.mpcs.concat(rows)})).errors.length)return rows;
    }
    throw new Error('No compatible dependent components for the axis helper: check supports and existing MPC records.');
  };
  C.equation=function(r){return labels[r.depDOF-1]+'('+r.depNode+') = '+r.masters.map(function(a){return a.coefficient+'·'+labels[a.dof-1]+'('+a.node+')';}).join(' + ');};
  C.convertLoadAxes=function(m,row,target){
    var out=Object.assign({},row,{components:row.components.slice(),coordinateSystem:target});if(row.coordinateSystem===target)return out;
    if(!m.elements[row.elementId-1])throw new Error('Choose an existing member before converting axes.');
    out.components=C.transformLoadComponents(row.components,C.basis(m,m.elements[row.elementId-1]),target);
    return out;
  };
  C.convertLoadType=function(row,target,confirmed){
    if(row.type===target)return Object.assign({},row,{components:row.components.slice()});
    if(target===20&&row.components.slice(0,3).some(function(x,i){return x!==row.components[i+3];})&&!confirmed)throw new Error('Confirm averaging the two different end intensities.');
    return Object.assign({},row,{type:target,components:target===21?row.components.concat(row.components):row.components.slice(0,3).map(function(x,i){return x/2+row.components[i+3]/2;})});
  };
}(window));
