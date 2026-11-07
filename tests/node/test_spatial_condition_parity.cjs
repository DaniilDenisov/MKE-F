// Generates solver inputs through the same public operations as the editor.
const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
const ctx={window:{MKEFSpatial:{registry:require('../../shared/element-types.json')}},TextEncoder};vm.createContext(ctx);
for(const name of ['preprocessor/js/namespace.js','shared/spatial-conditions.js','preprocessor/js/case-format.js','preprocessor/js/spatial-format.js','preprocessor/js/spatial-edit.js'])vm.runInContext(fs.readFileSync(name,'utf8'),ctx);
const {spatialFormat:F,spatialEdit:E}=ctx.window.MKEFPre,C=ctx.window.MKEFSpatialConditions;
const cases=[],frame=()=>F.parse(fs.readFileSync('examples/cases/CaseSpatialFrame.txt','utf8'));
function record(name,m,solve=false){const d=C.validate(m);assert.equal(d.errors.length,0,d.errors.join('\n'));const text=F.serialize(m);assert.equal(F.serialize(F.parse(text)),text);cases.push({name,text,valid:true,solve,warnings:d.warnings.length,projectors:d.active.map(q=>[0,1,2].map(i=>[0,1,2].map(j=>q.reduce((s,v)=>s+v[i]*v[j],0))))});}
const simple=frame(),sh=new E.History(simple);sh.condition('loads',0,{node:2,type:10,components:[1000,0,0,0,0,0],frequency:null});record('axial-static',simple,true);
const coupled=frame(),h=new E.History(coupled);h.condition('loads',0,{node:2,type:10,components:[0,1000,0,0,0,0],frequency:null});h.condition('mpcs',null,{depNode:2,depDOF:2,rhs:0,masters:[{node:2,dof:3,coefficient:2}]});h.setReleaseMask(1,2,['Mx']);h.condition('elementLoads',null,{type:21,elementId:1,coordinateSystem:1,components:[0,100,0,0,-50,0]});record('coupled-static',coupled,true);
h.condition('elementLoads',0,null);h.condition('loads',0,null);h.setAnalysis({...coupled.analysis,type:'modal'});record('coupled-modal',coupled,true);
h.setAnalysis({...coupled.analysis,type:'transient',timeStep:.0001,duration:.01,monitorNode:2,monitorDOF:2});
for(const type of [11,12,13]){h.condition('loads',coupled.loads.length?0:null,{type,node:2,components:[0,1000,0,0,0,0],frequency:type===11?10:null});record('coupled-transient-'+type,coupled,true);}
const tilted=frame();tilted.nodes[1]={x:2,y:3,z:4};tilted.elements[0].angleDegrees=37;tilted.loads=[];tilted.releases=[{elementId:1,end:2,component:'My'},{elementId:1,end:2,component:'Mz'}];record('tilted-rank-one',tilted);
const R=C.basis(tilted,tilted.elements[0]);tilted.loads=[{type:10,node:2,components:[0,0,0,...R[0].map(v=>v*400)],frequency:null}];record('tilted-moment-in-span',tilted);
const rankTwo=frame();rankTwo.loads=[];rankTwo.nodes=[{x:0,y:0,z:0},{x:2,y:3,z:4},{x:-1,y:2,z:1}];rankTwo.elements.push({...rankTwo.elements[0],node1:3,node2:2,referenceY:[1,0,0]});rankTwo.releases=[1,2].flatMap(elementId=>['Mx','My'].map(component=>({elementId,end:2,component})));record('two-oriented-ends',rankTwo);
const inactive=frame();inactive.loads=[];inactive.releases=['Mx','My','Mz'].map(component=>({elementId:1,end:1,component}));record('inactive-support-warning',inactive);
const small=frame();small.loads=[];small.nodes[1]={x:1e-11,y:1,z:0};small.elements[0].referenceY=[0,0,1];small.releases=['My','Mz'].map(component=>({elementId:1,end:2,component}));small.mpcs=[{depNode:2,depDOF:4,rhs:0,masters:[{node:2,dof:1,coefficient:1e-11}]}];record('small-active-projection',small);
function bad(name,base,suffix){const text=F.serialize(base)+suffix;assert.throws(()=>F.parse(text),undefined,name);cases.push({name,text,valid:false,solve:false});}
bad('duplicate-spc',simple,'bcfix_3d\n1\n1,1,0,0,0,0,0\n');
bad('duplicate-release',tilted,'releases\n1\n1,2,My\n');
bad('moment-out-of-span',tilted,'bcforce_stat_3d\n1\n10,2,0,0,0,0,1,0\n');
bad('projected-redundant-mpc',tilted,'mpc\n1\n2,5,0,1,2,4,1.5\n');
bad('spc-mpc-conflict',simple,'mpc\n1\n1,1,0,1,2,1,1\n');
bad('mpc-cycle',simple,'mpc\n2\n2,2,0,1,2,3,1\n2,3,0,1,2,2,1\n');
bad('duplicate-master',simple,'mpc\n1\n2,2,0,2,2,3,1,2,3,2\n');
bad('inactive-mpc',inactive,'mpc\n1\n2,4,0,1,1,4,1\n');
// Near-parallel transferred axes exercise the same singular-value cutoff as Octave.
for(const epsilon of [1e-12,1e-8]){const m=frame();m.loads=[];m.nodes=[{x:-1,y:0,z:0},{x:0,y:0,z:0},{x:-1,y:epsilon,z:0}];m.elements.push({...m.elements[0],node1:3,node2:2});m.releases=[1,2].flatMap(elementId=>['My','Mz'].map(component=>({elementId,end:2,component})));record('svd-'+epsilon,m);}
fs.mkdirSync('output',{recursive:true});fs.writeFileSync('output/spatial-condition-parity.json',JSON.stringify(cases));
console.log('Spatial browser/Octave fixture set: '+cases.length+' cases prepared.');
