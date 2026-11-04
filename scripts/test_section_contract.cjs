const planarStress=require('./planar-stress-reference.cjs');
const fs=require('node:fs'),assert=require('node:assert/strict'),path=require('node:path');
const root=require('./recovery-test-context.cjs')(true),M=root.MKEFPost,S=root.MKEFSections,P=root.MKEFSpatialSections;
const directory='output/section-contract',files=fs.readdirSync(directory).filter(n=>n.endsWith('.json'));
assert.equal(S.stress,undefined,'Planar extrema oracle must not be shipped in the application.');
assert.equal(files.length,42);
const near=(a,b)=>assert(Math.abs(a-b)<=1e-8*Math.max(1,Math.abs(a),Math.abs(b)),`${a} != ${b}`);
const cases=files.map(name=>JSON.parse(fs.readFileSync(path.join(directory,name))));
// Independent formula implementation remains a test oracle, outside the result path.
for(const raw of cases){
 const e=raw.model.elements[0],s=raw.model.sections[0],r=raw.analysis.elementResults[0],spatial=raw.model.dimension===3;
 const p=spatial?P.properties(s):S.properties(s,e.sectionAngle),c=r.sectionForceCoefficients;
 for(const k of ['area','Iy','Iz','momentOfInertia','centroidY','centroidZ','etaMin','etaMax'])near(e.sectionProperties[k],p[k]);
 if(spatial)near(e.sectionProperties.torsionConstant,p.torsionConstant);
 const expected=spatial?P.stress(p,c):planarStress(S,p,e.type,e.length,r.localEndForces,[0,-1000,0,-1000]);
 for(const k of ['minimum','maximum'])near(expected[k].value,r.sectionStress[k].value);
}
// Loading, coloring and point evaluation must not call property formulae or root search.
for(const [o,k]of [[S,'properties'],[P,'properties'],[S,'stress'],[P,'stress'],[S,'roots']])o[k]=()=>{throw Error('Unexpected browser calculation: '+k);};
for(const raw of cases){
 raw.analysis.elementResults[0].sectionStress.at=null; // Unknown JSON keys cannot replace evaluators.
 const ds=M.validateDataset(raw),e=raw.model.elements[0],stored=raw.analysis.elementResults[0],r=M.normalStressResults.get(ds,e);
 assert.equal(r.status,'ok',r.reason);near(r.sigmaMin,stored.sectionStress.minimum.value);near(r.sigmaMax,stored.sectionStress.maximum.value);
 for(let i=0;i<=200;i++)assert(r.at(i/200).sigmaAbsMax<=r.sigmaAbsMax+1e-8*Math.max(1,r.sigmaAbsMax));
 const profile=raw.model.dimension===3?M.spatialResults.sectionProfile(e):M.sectionResults.get(ds,e).properties;
 for(const k of Object.keys(e.sectionProperties))assert.equal(profile[k],e.sectionProperties[k]);
}
let invalid=0;
for(const dimension of [2,3]){
 const original=cases.find(d=>d.model.dimension===dimension&&d.model.sections[0].type==='rectangle');
 for(const change of [
  d=>delete d.model.elements[0].sectionProperties,
  d=>d.model.elements[0].sectionProperties.area*=2,
  d=>d.model.elements[0].sectionProperties.Iy=-1,
  d=>d.model.elements[0].sectionProperties.centroidY=Infinity,
  d=>d.model.elements[0].sectionProperties.etaMin=100,
  d=>delete d.analysis.elementResults[0].sectionStress,
  d=>delete d.analysis.elementResults[0].sectionStress.positions,
  d=>d.analysis.elementResults[0].sectionStress.minimum.value+=1e6,
  d=>d.analysis.elementResults[0].sectionStress.minimum.xi=2,
  d=>d.analysis.elementResults[0].sectionStress.positions=[0,.7,.6,1],
  d=>d.analysis.elementResults[0].sectionStress.minimum[dimension===2?'eta':'y']=1e6,
  d=>d.model.elements[0].sectionProperties[dimension===2?'momentOfInertia':'torsionConstant']*=2
 ]){const bad=structuredClone(original);change(bad);assert.throws(()=>M.validateDataset(bad));invalid++;}
}
console.log(`Solver section contract: ${cases.length} cases, independent extrema/properties, no browser formula/root calls, ${invalid} malformed contracts PASS`);
