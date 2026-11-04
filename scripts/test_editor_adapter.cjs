const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
const ctx={window:{MKEFSpatial:{registry:require('../shared/element-types.json')}},TextEncoder};vm.createContext(ctx);
for(const f of ['sections','spatial-sections','support-catalog','support-semantics','mpc','releases','spatial-conditions'])vm.runInContext(fs.readFileSync('shared/'+f+'.js','utf8'),ctx);
for(const f of ['namespace','case-format','spatial-format','model-edit','editor-adapter','spatial-edit'])vm.runInContext(fs.readFileSync('preprocessor/js/'+f+'.js','utf8'),ctx);
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
