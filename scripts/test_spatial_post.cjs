const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
const ctx={window:{MKEFPost:{},MKEFSpatial:{registry:require('../shared/element-types.json')}}};vm.createContext(ctx);
for(const f of ['shared/sections.js','shared/spatial-sections.js','postprocessor/js/data-model.js','postprocessor/js/spatial-data.js','postprocessor/js/spatial-results.js','postprocessor/js/spatial-stress.js'])vm.runInContext(fs.readFileSync(f,'utf8'),ctx);
const M=ctx.window.MKEFPost;
for(const name of fs.readdirSync('output').filter(n=>/^CaseSpatial.*\.json$/.test(n))){const data=JSON.parse(fs.readFileSync('output/'+name));M.spatialData.validate(data);}
const data=JSON.parse(fs.readFileSync('output/CaseSpatialFrame.json'));const bad=structuredClone(data);bad.model.nodeKinematics.entries[0][2]=2;assert.throws(()=>M.spatialData.validate(bad),/registry mismatch/);
const bad2=structuredClone(data);bad2.model.elements[0].localBasis[0][0]=-1;assert.throws(()=>M.spatialData.validate(bad2),/direction/);
console.log('Spatial v8 contract: fixtures, singleton/partial results and invalid kinematics passed.');
const ds=M.spatialData.validate(data),S=M.spatialResults,e=data.model.elements[0],u=S.vector(ds,'displacements'),r=S.recover(ds,e,u),op=JSON.parse(fs.readFileSync('output/spatial-frame-operators.json')),stored=data.analysis.elementResults[0];
for(let i=0;i<12;i++)assert(Math.abs(op.K[i].reduce((sum,x,j)=>sum+x*r.u[j],0)-stored.localEndForces[i])<1e-8);
const tip=S.sample(e,r.u,1);assert(Math.abs(tip.displacement[1]-2000*8/(3*2e11*2e-6))<1e-12);assert(Math.abs(tip.displacement[2]-3000*8/(3*2e11*8e-6))<1e-12);assert(Math.abs(tip.rotation[0]-400*2/(8e10*1e-6))<1e-12);
const torsion=Array(12).fill(0);torsion[9]=.1;const middle=S.sample(e,torsion,.5);assert.equal(middle.displacement.every(x=>x===0),true);assert.equal(middle.rotation[0],.05);
for(const mat of [op.K,op.M])for(let i=0;i<12;i++)for(let j=0;j<12;j++)assert.equal(mat[i][j],mat[j][i]);
for(const name of fs.readdirSync('output').filter(n=>/^CaseSpatial.*\.json$/.test(n))){const ds=M.spatialData.validate(JSON.parse(fs.readFileSync('output/'+name))),a=ds.raw.analysis;if(a.type!=='static')continue;
for(const e of ds.raw.model.elements){const r=ds.resultsByElementId.get(e.id),u=S.local(e,S.vector(ds,'displacements'));assert(u.every((x,i)=>Math.abs(x-r.localDisplacements[i])<1e-10));}}
console.log('Spatial recovery: local stiffness, force parity, two bending planes, torsion and symmetry passed.');
const tilted=M.spatialData.validate(JSON.parse(fs.readFileSync('output/CaseSpatialTiltedTransient.json'))),reference=JSON.parse(fs.readFileSync('output/spatial-post-parity.json')),te=tilted.raw.model.elements[0];
function near(a,b,tolerance=1e-9){assert(Math.abs(a-b)<=tolerance*Math.max(1,Math.abs(a),Math.abs(b)),`${a} != ${b}`);}
const operators=S.operators(te,tilted.raw.version,tilted.raw.model);for(let i=0;i<12;i++)for(let j=0;j<12;j++){near(operators.K[i][j],reference.K[i][j],1e-12);near(operators.M[i][j],reference.M[i][j],1e-12);}
for(const [k,frame]of reference.frames.entries()){
  const u=S.vector(tilted,'displacements',frame),a=S.vector(tilted,'accelerations',frame),r=S.recover(tilted,te,u,a);
  for(let i=0;i<12;i++){near(r.u[i],reference.u[i][k]);near(r.end[i],reference.elastic[i][k]);near(r.inertia[i],reference.inertia[i][k]);near(S.physical(tilted,u,i),reference.physical[i][frame]);}
  assert(Math.abs(r.total[11])<1e-8*Math.max(...r.end.map(Math.abs)));assert(Math.abs(r.end[11])>1e-6);
  const lambda=tilted.raw.analysis.mpcMultipliers[0][frame];near(S.reaction(tilted,'mpcForces',7,frame),lambda);near(S.reaction(tilted,'mpcForces',8,frame),-2*lambda);
}
const mixed=tilted.physicalRows.findIndex(row=>row.length>1),newSpectrum=S.spectrum(tilted,mixed);assert(newSpectrum);
const old=structuredClone(tilted.raw);delete old.analysis.nodalDisplacementAmplitudeSpectrum;delete old.analysis.physicalDOFIds;delete old.analysis.nodalSupportReactions;delete old.analysis.nodalMPCForces;delete old.analysis.timeStep;
const legacy=M.spatialData.validate(old);assert.equal(S.spectrum(legacy,mixed),null);assert(S.spectrum(legacy,0));assert.equal(S.reaction(legacy,'supportReactions',0),null);
const partial=M.spatialData.validate(JSON.parse(fs.readFileSync('output/CaseSpatialTiltedPartial.json')));assert.equal(S.full(partial),false);assert.equal(S.physical(partial,S.vector(partial,'displacements'),7),null);assert.equal(S.spectrum(partial,mixed),null);
for(const mutate of [x=>x.model.elements[0].localFromSystem[0][0]+=.01,x=>x.model.dofRegistry.find(r=>r.kind==='elementEnd').direction[0]+=.01]){const bad=structuredClone(tilted.raw);mutate(bad);assert.throws(()=>M.spatialData.validate(bad));}
console.log('Tilted/released frame: Octave K/M, local and physical recovery, elastic/inertial balance, legacy and partial spectra passed.');
const sectionDS=M.spatialData.validate(JSON.parse(fs.readFileSync('output/CaseSpatialSectionReleasedTransient.json'))),se=sectionDS.raw.model.elements[0];
for(const row of JSON.parse(fs.readFileSync('output/spatial-section-dynamic-parity.json'))){const r=S.recover(sectionDS,se,S.vector(sectionDS,'displacements',row.frame),S.vector(sectionDS,'accelerations',row.frame)),stress=S.sectionStress(se,r);
 for(let i=0;i<6;i++)for(let j=0;j<4;j++)near(r.coeff[i][j],row.coefficients[i][j]);
 for(const bound of ['minimum','maximum'])near(stress[bound].value,row.stress[bound].value);
 assert(Math.abs(r.total[11])<1e-9*Math.max(1,...r.end.map(Math.abs),...r.inertia.map(Math.abs)),'Released end elastic/inertial balance');
}
const sectionPartial=M.spatialData.validate(JSON.parse(fs.readFileSync('output/CaseSpatialSectionReleasedPartial.json')));assert.equal(S.recover(sectionPartial,sectionPartial.raw.model.elements[0],S.vector(sectionPartial,'displacements')),null);
console.log('Released catalog frame: Octave normal stress parity, dynamic balance and missing partial data PASS');
