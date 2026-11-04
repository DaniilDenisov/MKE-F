const assert=require('node:assert/strict'),fs=require('node:fs');
const context=require('./recovery-test-context.cjs'),M=context(),S=M.spatialResults;
const read=file=>JSON.parse(fs.readFileSync(file,'utf8'));
const manifest=read('output/recovery/manifest.json');assert.equal(manifest.length,12);
assert.equal(new Set(manifest.map(row=>row.name)).size,12);
function near(a,b,rtol=1e-10,atol=1e-10){
  assert(Number.isFinite(a)&&Number.isFinite(b));
  assert(Math.abs(a-b)<=atol+rtol*Math.max(Math.abs(a),Math.abs(b)),`${a} != ${b}`);
}
function same(a,b){
  if(a===null||b===null){assert.equal(a,b);return;}
  assert.equal(a.length,b.length);
  a.forEach((v,i)=>Array.isArray(v)?same(v,b[i]):near(v,b[i]));
}
function sameProduct(actual,expected,operator,vector,label){
  assert.equal(actual.length,expected.length);
  actual.forEach((value,i)=>{
    // A near-zero resultant can subtract large terms. Bound roundoff using
    // the dot-product contributions, not the residual's relative magnitude.
    const scale=operator[i].reduce((sum,v,j)=>sum+Math.abs(v*vector[j]),0);
    assert(Number.isFinite(value)&&Number.isFinite(expected[i]));
    assert(Math.abs(value-expected[i])<=1e-10+1e-12*scale,
      `${label}, row ${i}: ${value} != ${expected[i]}, contribution scale ${scale}`);
  });
}
let frames=0;
for(const row of manifest){
  const stem='output/recovery/'+row.name;
  const raw=read(stem+'-v8.json'),ds=M.validateDataset(raw),references=read(stem+'-reference.json');
  assert.equal(raw.model.elements.length,row.elements);
  for(const [i,e] of raw.model.elements.entries()){
    if(row.analysis==='static')assert.equal(e.recoveryOperators,undefined);
    else{same(e.recoveryOperators.K,references[i].K);same(e.recoveryOperators.M,references[i].M);}
    const count=row.analysis==='modal'?raw.analysis.frequenciesHz.length:row.analysis==='transient'?raw.analysis.time.length:1;
    for(let frame=0;frame<count;frame++){
      const r=S.recover(ds,e,S.vector(ds,'displacements',frame),S.vector(ds,'accelerations',frame));
      if(!r){assert(/partial|no-displacements/.test(row.name));continue;}
      if(row.analysis==='static')same(r.end,references[i].elastic);
      if(row.analysis!=='static'){
        sameProduct(r.end,references[i].elastic.map(values=>values[frame]),e.recoveryOperators.K,r.u,row.name+' elastic '+frame);
        if(r.inertia)sameProduct(r.inertia,references[i].inertia.map(values=>values[frame]),e.recoveryOperators.M,
          S.local(e,S.vector(ds,'accelerations',frame)),row.name+' inertia '+frame);
      }
      frames++;
    }
  }
}
assert(frames>100,'No dynamic recovery coverage');
// Prove the exported matrix is applied (not just transported/validated).
const source=read('output/recovery/CaseSpatialFrameModal-v8.json');
const baseline=M.validateDataset(structuredClone(source)),baseElement=baseline.raw.model.elements[0];
const base=S.recover(baseline,baseElement,S.vector(baseline,'displacements'),S.vector(baseline,'accelerations'));
const altered=structuredClone(source);for(const key of ['K','M'])altered.model.elements[0].recoveryOperators[key]=altered.model.elements[0].recoveryOperators[key].map(row=>row.map(v=>2*v));
const changed=M.validateDataset(altered),e=changed.raw.model.elements[0];
const actual=S.recover(changed,e,S.vector(changed,'displacements'),S.vector(changed,'accelerations'));
same(actual.end,base.end.map(v=>2*v));same(actual.inertia,base.inertia.map(v=>2*v));
const invalid=[
  x=>delete x.model.elements[0].recoveryOperators,
  x=>delete x.model.elements[0].recoveryOperators.M,
  x=>x.model.elements[0].recoveryOperators.K.pop(),
  x=>x.model.elements[0].recoveryOperators.M[0].pop(),
  x=>x.model.elements[0].recoveryOperators.K[0][0]=NaN,
  x=>x.model.elements[0].recoveryOperators.M[0][0]=Infinity,
  x=>x.model.elements[0].recoveryOperators.K[0][1]=1e9,
  x=>x.model.recovery.coordinates='system',
  x=>x.model.recovery.stiffnessFormulation='unknown',
  x=>x.model.recovery.massFormulation='lumped',
  x=>x.model.recovery.version=2,
  x=>x.metadata.exporter.contractVersion=2,
  x=>x.model.dimension=2
];
for(const mutate of invalid){const bad=structuredClone(source);mutate(bad);assert.throws(()=>M.validateDataset(bad));}
const missing=structuredClone(source);delete missing.model.elements[0].recoveryOperators;
assert.throws(()=>S.operators(missing.model.elements[0],8),/Missing v8/);
const staticBad=read('output/recovery/CaseSpatialFrame-v8.json');staticBad.model.elements[0].recoveryOperators=source.model.elements[0].recoveryOperators;
assert.throws(()=>M.validateDataset(staticBad),/static results/);
// K has rigid modes; it must not be rejected as non-positive-definite.
const K=source.model.elements[0].recoveryOperators.K,rigid=[1,0,0,0,0,0,1,0,0,0,0,0];
assert(K.every(row=>row.reduce((sum,v,i)=>sum+v*rigid[i],0)===0));
console.log(`v8: 12 solver pairs, ${frames} recovered element frames, matrix provenance, malformed contracts and rigid modes PASS`);
