const fs=require('fs'),vm=require('vm'),assert=require('node:assert/strict');
const window={};window.window=window;const ctx=vm.createContext({window,console,TextEncoder,TextDecoder,Map,Set});
for(const f of ['sections','spatial-sections','support-catalog','support-semantics','releases','mpc'])vm.runInContext(fs.readFileSync('shared/'+f+'.js','utf8'),ctx);
for(const f of ['namespace','geometry','data-model','static-view','modal-view','transient-view','section-view','planar-presentation'])vm.runInContext(fs.readFileSync('postprocessor/js/'+f+'.js','utf8'),ctx);
const M=window.MKEFPost,P=M.planarPresentation;let count=0;
for(const file of fs.readdirSync('output').filter(f=>f.endsWith('.json'))){let raw;try{raw=JSON.parse(fs.readFileSync('output/'+file));}catch(_){continue;}if(!raw.model||raw.model.dimension!==2)continue;let ds;try{ds=M.validateDataset(raw);}catch(_){continue;}if(raw.analysis.type!=='static')continue;
 const frame=P.geometry(ds,raw.analysis.displacements,3,31,{showOriginal:true,showDeformed:true,showNodes:true});let offset=0;
 raw.model.elements.forEach(e=>{const points=M.geometry.sampleElement(ds,e,raw.analysis.displacements,3,e.type===112?2:31);for(let j=1;j<points.length;j++)for(const p of [points[j-1],points[j]]){assert.equal(frame.deformed.positions[offset++],p.x);assert.equal(frame.deformed.positions[offset++],p.y);assert.equal(frame.deformed.positions[offset++],0);}});
 const mode=raw.model.elements[0].type===113?'M':'axialForce',result=P.result(ds,mode,31);assert(Number.isFinite(result.maximum));assert(result.diagrams.positions.every(Number.isFinite));assert(result.diagramFill.positions.every(Number.isFinite));
 if(mode==='M'){let max=0;for(const e of raw.model.elements)for(const r of M.staticResults.diagramSamples(ds,e,mode,31))max=Math.max(max,Math.abs(r.value));assert.equal(max,result.maximum);}
 const reactions=P.forces(ds,raw.analysis.reactions,'#c56a00',true);assert(reactions.positions.every(Number.isFinite));count++;
}
assert(count>=5,'Expected generated planar fixtures');console.log('Prepared planar geometry, diagrams, fills and reactions match existing numerics: '+count+' static fixtures PASS');
// Independent cubic coefficients check the two sampling paths, including
// private end rotations that differ from the shared nodal rotations.
const close=(a,b)=>assert(Math.abs(a-b)<=1e-12*Math.max(1,Math.abs(a),Math.abs(b)),`${a} != ${b}`);
for(const direction of [[1,0],[0,1],[.6,.8],[-.8,.6]]) for(const released of [false,true]) {
 const [c,s]=direction,L=5,a={id:1,x:2,y:-3},b={id:2,x:2+L*c,y:-3+L*s};
 const u=[.1,.2,.3],v=[-.4,.5,-.6];
 const first={x:c*u[0]-s*u[1],y:s*u[0]+c*u[1],r:u[2]},second={x:c*v[0]-s*v[1],y:s*v[0]+c*v[1],r:v[2]};
 const element={type:113,nodeIds:[1,2],...(released?{globalDOFs:[1,2,7,4,5,8]}:{})};
 const vector=[first.x,first.y,released?99:first.r,second.x,second.y,released?-99:second.r,first.r,second.r];
 const dataset={nodesById:new Map([[1,a],[2,b]]),raw:{model:{nodes:[a,b],dofMap:[[1,2,3],[4,5,6]]}}};
 const transverse=[u[1],L*u[2],-3*u[1]-2*L*u[2]+3*v[1]-L*v[2],2*u[1]+L*u[2]-2*v[1]+L*v[2]];
 function expected(t) {
  const axial=u[0]+(v[0]-u[0])*t,bend=transverse.reduce((sum,value,i)=>sum+value*t**i,0);
  return {x:c*axial-s*bend,y:s*axial+c*bend};
 }
 for(const samples of [5,31,201])M.geometry.sampleElementDisplacements(dataset,element,vector,samples).forEach((d,i)=>{
  const reference=expected(i/(samples-1));close(d.x,reference.x);close(d.y,reference.y);
 });
 for(const scale of [0,-3,1e-12,1,1e4]) {
  const at=M.geometry.elementInterpolator(dataset,element,vector,scale);
  for(const t of [0,1e-12,.173,.5,.827,1-1e-12,1]) {
   const reference=expected(t),point=at(t);close(point.x,a.x+(b.x-a.x)*t+scale*reference.x);close(point.y,a.y+(b.y-a.y)*t+scale*reference.y);
  }
 }
 const zero=M.geometry.elementInterpolator(dataset,element,vector.map(()=>0),1);
 close(zero(.173).x,a.x+(b.x-a.x)*.173);close(zero(.173).y,a.y+(b.y-a.y)*.173);
}
console.log('Hermite interpolation: independent cubic, axes, private rotations, nonuniform stations and scales PASS');
