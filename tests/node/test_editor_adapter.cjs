const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
const ctx={window:{MKEFSpatial:{registry:require('../../shared/element-types.json')}},TextEncoder};vm.createContext(ctx);
for(const f of ['sections','spatial-sections','support-catalog','support-semantics','mpc','releases','spatial-conditions'])vm.runInContext(fs.readFileSync('shared/'+f+'.js','utf8'),ctx);
for(const f of ['namespace','case-format','spatial-format','editor-adapter','spatial-edit'])vm.runInContext(fs.readFileSync('preprocessor/js/'+f+'.js','utf8'),ctx);
const M=ctx.window.MKEFPre,A=M.editorAdapter,json=x=>JSON.stringify(x,(_,v)=>v&&typeof v==='object'&&!Array.isArray(v)?Object.fromEntries(Object.keys(v).sort().map(k=>[k,v[k]])):v);
let count=0;
for(const file of fs.readdirSync('examples/cases').filter(f=>f.endsWith('.txt'))){
 const text=fs.readFileSync('examples/cases/'+file,'utf8');let native;try{native=M.caseFormat.parse(text);}catch(e){continue;}
 const m=A.fromNative(native,text);if(m.dimension!==2)continue;
 const round=A.toNative(m);for(const k of ['nodes','elements','supports','loads','elementLoads','mpcs','releases','sections','analysis']){
  if(['supports','loads'].includes(k))continue;assert.equal(json(round[k]),json(native[k]),file+' '+k);
 }
 if(!A.validate(m,true).length){const reread=A.parse(A.serialize(m));assert.equal(json(A.toNative(reread)),json(A.toNative(A.parse(M.caseFormat.serialize(native)))),file);}
 const frame=A.sceneConditions(m);m.mpcs.forEach((r,i)=>assert.equal(frame.mpcs[i].depDOF,r.depDOF===3?6:r.depDOF));
 m.loads.forEach((r,i)=>{assert.equal(r.components[5],native.loads[i].mz);assert.equal(r.components[2],0);});
 count++;
}
for(const code of [112,113]){
 const m=A.empty(code),h=new M.spatialEdit.History(m);h.addNode({x:0,y:0,z:0});h.addNode({x:2,y:0,z:0});
 h.member(null,{node1:1,node2:2,area:.01,youngsModulus:200e9,density:7800,...(code===113?{momentOfInertia:.0001}:{})});
 h.setAnalysis({...m.analysis,type:'static'});h.condition('supports',null,{node:1,mask:[1,1,0,0,0,code===113?1:0]});
 h.condition('loads',null,{node:2,type:10,components:[1,2,0,0,0,code===113?3:0]});
 assert.equal(A.validate(m,true).length,0);const before=json(m);assert.throws(()=>h.moveNode(1,{x:2,y:0,z:1}),/Z = 0/);assert.equal(json(m),before);
 if(code===113){h.condition('elementLoads',null,{type:20,elementId:1,coordinateSystem:1,components:[0,-1,0]});h.moveNode(1,{x:3,y:0,z:0});assert.equal(m.elementLoads.length,0);h.travel(true);assert.equal(m.elementLoads.length,1);}
}
console.log('Planar adapter: '+count+' native cases, physical Mz/MPC mapping, atomic edits and load invalidation PASS');
// An inclined planar basis differs from the spatial automatic reference axis.
// Check known components as well as reversibility through both public callers.
const C=ctx.window.MKEFSpatialConditions;
for(const dimension of [2,3]) for(const type of [20,21]) {
 const model={dimension,nodes:[{x:0,y:0,z:0},{x:3,y:4,z:0}],elements:[{node1:1,node2:2,referenceY:[0,0,1],angleDegrees:37}]};
 const source={elementId:1,type,coordinateSystem:2,components:type===20?[0,-7,0]:[0,-7,0,2,3,4]};
 const saved=JSON.stringify(source),basis=A.basis(model,model.elements[0]);
 const expected=[];
 for(let offset=0;offset<source.components.length;offset+=3) {
  const [x,y,z]=source.components.slice(offset,offset+3);
  for(const row of basis)expected.push(row[0]*x+row[1]*y+row[2]*z);
 }
 const local=M.editorConditions.convertLoadAxes(model,source,1);
 assert(local.components.every((v,i)=>Math.abs(v-expected[i])<1e-12));
 if(dimension===2) { assert(Math.abs(local.components[0]+5.6)<1e-12);assert(Math.abs(local.components[1]+4.2)<1e-12); }
 if(dimension===3)assert.equal(json(local),json(C.convertLoadAxes(model,source,1)));
 const back=M.editorConditions.convertLoadAxes(model,local,2);
 assert(back.components.every((v,i)=>Math.abs(v-source.components[i])<1e-12));
 assert.equal(JSON.stringify(source),saved);
 for(const [convert,message] of [[M.editorConditions.convertLoadAxes,'Choose an existing element.'],[C.convertLoadAxes,'Choose an existing member before converting axes.']]) {
  const bad={...source,elementId:99};
  const copied=convert(model,bad,2);assert.notEqual(copied,bad);assert.notEqual(copied.components,bad.components);assert.equal(json(copied),json(bad));
  assert.throws(()=>convert(model,bad,1),e=>e.message===message);
 }
}
console.log('Load axes: known planar components, spatial reference/angle, both endpoints, copies and errors PASS');
