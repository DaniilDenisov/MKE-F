const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
const ctx={window:{MKEFSpatial:{registry:require('../shared/element-types.json')}},TextEncoder};vm.createContext(ctx);
for(const file of ['preprocessor/js/namespace.js','shared/spatial-conditions.js','preprocessor/js/case-format.js','preprocessor/js/spatial-format.js','preprocessor/js/spatial-edit.js']) vm.runInContext(fs.readFileSync(file,'utf8'),ctx);
const C=ctx.window.MKEFSpatialConditions,F=ctx.window.MKEFPre.spatialFormat,E=ctx.window.MKEFPre.spatialEdit;
const frame=()=>F.parse(fs.readFileSync('examples/cases/CaseSpatialFrame.txt','utf8'));
const m=frame(),h=new E.History(m),before=F.serialize(m);
assert.throws(()=>h.condition('supports',null,{node:1,mask:[1,0,0,0,0,0]}),/Duplicate/);
assert.throws(()=>h.condition('releases',null,{elementId:1,end:2,component:'Mx'}),/inactive/);
assert.equal(F.serialize(m),before);assert.equal(h.undo.length,0);
h.condition('releases',null,{elementId:1,end:2,component:'Mz'});assert.equal(m.releases.length,1);h.travel(true);assert.equal(F.serialize(m),before);
assert.throws(()=>h.setAnalysis({...m.analysis,type:'modal'}),/loads/);
assert.throws(()=>h.setAnalysis({...m.analysis,type:'transient'}),/Confirm/);
h.setAnalysis({...m.analysis,type:'transient'},true);h.travel(true);assert.equal(F.serialize(m),before);
const mp={depNode:2,depDOF:2,rhs:0,masters:[{node:2,dof:3,coefficient:2}]};h.condition('mpcs',null,mp);
assert.throws(()=>h.condition('mpcs',null,{depNode:2,depDOF:3,rhs:0,masters:[{node:2,dof:2,coefficient:1}]}),/Cycle/);
assert.throws(()=>h.condition('supports',null,{node:2,mask:[0,1,0,0,0,0]}),/restrained/);
// Compare the recovered rotation subspace, not SVD-vector signs.
for(const epsilon of [0,1e-12,1e-8]) {
  const q=C.activeBasis([[1,0,0],[1,epsilon,0]]);assert.equal(q.length,epsilon>1e-10?2:1);
  for(const a of q) assert(Math.abs(a.reduce((s,v)=>s+v*v,0)-1)<1e-12);
}
const r=frame();r.loads=[];r.releases=['Mx','My','Mz'].map(component=>({elementId:1,end:1,component}));
assert.equal(C.validate(r).warnings.length,3);
for(const type of [10,11,12,13]) {
  const record={type,frequency:2};const actual=[0,.1,.2].map(t=>C.law(record,'transient',t,.1));
  const expected=type===11?[0,Math.sin(.4*Math.PI),Math.sin(.8*Math.PI)]:type===13?[1,1,1]:[0,1,0];
  assert(actual.every((v,i)=>Math.abs(v-expected[i])<1e-12));
}
const rotated=frame();rotated.nodes[1]={x:2,y:3,z:4};rotated.elements[0].angleDegrees=37;
const q={type:21,elementId:1,coordinateSystem:1,components:[1,2,-3,-4,5,6]},globalQ=C.convertLoadAxes(rotated,q,2),backQ=C.convertLoadAxes(rotated,globalQ,1);
assert(backQ.components.every((x,i)=>Math.abs(x-q.components[i])<1e-12));assert.throws(()=>C.convertLoadType(q,20,false),/Confirm/);
assert.deepEqual(Array.from(C.convertLoadType(q,20,true).components),[-1.5,3.5,1.5]);
const axis=F.empty(1133);axis.nodes=[{x:1,y:2,z:3},{x:4,y:8,z:12},{x:2,y:4,z:6}];
const axisRows=C.axisHelper(axis,1,2,3);assert.equal(axisRows.length,2);assert.equal(new Set(axisRows.map(r=>r.depNode+':'+r.depDOF)).size,2);
// Arbitrary endpoint displacement plus an independent axial slide satisfies both equations.
const u1=[2,-3,5],u2=[7,1,-2],middle=u1.map((x,i)=>x*2/3+u2[i]/3+[3,6,9][i]*.73),disp=[u1,u2,middle];
for(const r of axisRows)assert(Math.abs(disp[r.depNode-1][r.depDOF-1]-r.masters.reduce((s,a)=>s+a.coefficient*disp[a.node-1][a.dof-1],0))<1e-12);
const ah=new E.History(axis);ah.addMPCs(axisRows);assert.equal(ah.undo.length,1);ah.travel(true);assert.equal(axis.mpcs.length,0);ah.travel(false);assert.equal(axis.mpcs.length,2);
axis.nodes[2].z+=1e-5;assert.throws(()=>C.axisHelper(axis,1,2,3),/inside/);
const released=frame(),rh=new E.History(released);rh.setReleaseMask(1,2,['My','Mz']);assert.equal(released.releases.length,2);assert.equal(rh.undo.length,1);rh.travel(true);assert.equal(released.releases.length,0);rh.travel(false);assert.equal(released.releases.length,2);
assert.throws(()=>rh.setReleaseMask(1,2,['Mx','My','Mz']),/inactive/);assert.equal(released.releases.length,2);rh.setReleaseMask(1,2,[]);assert.equal(released.releases.length,0);rh.travel(true);assert.equal(released.releases.length,2);
rh.setReleaseMask(1,1,['Mx','My','Mz']);assert.equal(C.validate(released).warnings.length,3);assert.equal(released.releases.length,5);rh.travel(true);assert.equal(released.releases.length,2);
assert.throws(()=>rh.setAnalysis({...released.analysis,type:'transient',monitorNode:2,monitorDOF:5},true),/inactive/);const unchanged=F.serialize(released);assert.throws(()=>rh.moveNode(1,{x:2,y:1,z:0}),/inactive/);assert.equal(F.serialize(released),unchanged);
const large=F.empty(1133);large.nodes=Array.from({length:100000},(_,i)=>({x:i,y:0,z:0}));large.mpcs=large.nodes.slice(1).map((_,i)=>({depNode:i+2,depDOF:1,rhs:0,masters:[{node:i+1,dof:1,coefficient:1}]}));large.supports=[{node:1,mask:[1,0,0,0,0,0]}];
const start=Date.now();assert.equal(C.validate(large).errors.length,0);large.loads=[{type:10,node:2,components:[1,0,0,0,0,0],frequency:null}];assert(C.validate(large).errors.some(x=>x.includes('100000')));console.log('100000-condition sparse chain validation: '+(Date.now()-start)+' ms.');
console.log('Spatial conditions: validation, history, subspaces, laws, transforms, MPC helpers and end masks passed.');
