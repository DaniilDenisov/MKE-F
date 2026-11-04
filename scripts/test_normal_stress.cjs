const planarStress=require('./planar-stress-reference.cjs');
const assert=require('node:assert/strict'),fs=require('node:fs'),vm=require('node:vm');
const ctx={window:{MKEFSpatial:{registry:require('../shared/element-types.json')}},TextEncoder};vm.createContext(ctx);
for(const file of ['shared/sections.js','shared/spatial-sections.js','shared/support-catalog.js','shared/support-semantics.js','shared/mpc.js','shared/releases.js',...['namespace','data-model','spatial-data','spatial-results','spatial-stress','geometry','static-view','section-view','normal-stress-results'].map(s=>'postprocessor/js/'+s+'.js')])vm.runInContext(fs.readFileSync(file,'utf8'),ctx);
const M=ctx.window.MKEFPost,V=M.normalStressResults,S=ctx.window.MKEFSections,P=ctx.window.MKEFSpatialSections;
const near=(a,b)=>assert(Math.abs(a-b)<=1e-9*Math.max(1,Math.abs(b)),`${a} != ${b}`);
const read=name=>JSON.parse(fs.readFileSync('output/'+name+'.json'));
const load=name=>M.validateDataset(read(name));
let cases=0;
for(const name of fs.readdirSync('output').filter(s=>/^(SectionCase|Mises[23]D_|CaseSectionBeam|CaseSpatialFrameSection|CaseSpatialTrussSection).*\.json$/.test(s))){
 const raw=JSON.parse(fs.readFileSync('output/'+name));if(raw.analysis?.type!=='static')continue;
 const ds=M.validateDataset(raw);
 for(const e of raw.model.elements){
  if(!e.sectionId&&![112,1133].includes(e.type)){assert.equal(V.get(ds,e).status,'unavailable');continue;}
  const r=V.get(ds,e);assert.equal(r.status,'ok',name+': '+r.reason);assert.strictEqual(V.get(ds,e),r);
  const expected=ds.resultsByElementId.get(e.id).sectionStress;
  if(expected?.minimum){near(r.sigmaMin,expected.minimum.value);near(r.sigmaMax,expected.maximum.value);}
  near(r.sigmaAbsMax,Math.max(Math.abs(r.sigmaMin),Math.abs(r.sigmaMax)));
  assert(r.positions.includes(0)&&r.positions.includes(1)&&r.positions.includes(r.peak.xi));
  for(let i=0;i<=200;i++){const p=r.at(i/200);near(p.sigmaAbsMax,Math.max(Math.abs(p.sigmaMin),Math.abs(p.sigmaMax)));assert(p.sigmaAbsMax<=r.sigmaAbsMax+1e-8*Math.max(1,r.sigmaAbsMax));}
  cases++;
 }
}
// Independent analytic signs, both bending axes, combined loading and pure torsion.
for(const forces of [[100,0,0],[-100,0,0],[0,30,0],[0,0,40],[-100,30,40],[0,0,0]]){
 const ds=load('Mises3D_rectangle'),e=ds.raw.model.elements[0],p=M.spatialResults.sectionProfile(e),stored=ds.resultsByElementId.get(e.id);
 stored.sectionForceCoefficients=Array.from({length:6},()=>[0,0,0,0]);[0,4,5].forEach((k,i)=>stored.sectionForceCoefficients[k][0]=forces[i]);stored.sectionForceCoefficients[3][0]=999;
 stored.sectionStress=P.stress(p,stored.sectionForceCoefficients);
 const r=V.get(ds,e),bend=Math.abs(forces[1])*Math.max(...p.vertices.map(v=>Math.abs(v.z)))/p.Iy+Math.abs(forces[2])*Math.max(...p.vertices.map(v=>Math.abs(v.y)))/p.Iz;
 near(r.sigmaAbsMax,Math.abs(forces[0])/p.area+bend);near(r.peak.sigmaX,forces[0]<0?r.sigmaMin:r.sigmaMax);
 // Changing orientation/global geometry cannot apply the section angle twice.
 e.angleDegrees=127;e.localBasis=[[0,1,0],[-1,0,0],[0,0,1]];V.clear(ds);near(V.get(ds,e).sigmaAbsMax,r.sigmaAbsMax);
}
// Manual truss with a nonconstant axial force: internal peak, not exported mean.
for(const spatial of [false,true]){
 const ds=load(spatial?'CaseSpatialTrussSection':'CaseSectionBeam'),e=ds.raw.model.elements[0],stored=ds.resultsByElementId.get(e.id);e.type=spatial?1133:112;e.sectionId=0;e.properties.area=2;
 stored.axialStress=999;stored.localEndForces[0]=2;
 stored.sectionForceCoefficients=Array.from({length:spatial?6:3},()=>[0,0,0,0]);stored.sectionForceCoefficients[0]=[-2,-6,6,0];
 stored.sectionStress={minimum:{value:-1.75,xi:.5},maximum:{value:-1,xi:0},positions:[0,.5,1]};
 const r=V.get(ds,e);assert.equal(r.status,'ok');near(r.peak.xi,.5);near(r.sigmaAbsMax,1.75);assert.equal(r.peak.eta,null);assert.equal(r.peak.y,null);
 V.clear(ds);stored.localEndForces[0]=NaN;stored.sectionForceCoefficients[0][1]=NaN;assert.equal(V.get(ds,e).status,'unavailable');
}
// Equivalent rotated planar and spatial states for all eight shapes.
for(const type of Object.keys(S.types))for(const angle of [0,90,180,270]){
 const ds=load('CaseSectionBeam'),e=ds.raw.model.elements[0],shape={id:1,name:type,type,dimensions:{b:.1,h:.2,d:.1,t:.005,tw:.008,tf:.012}};
 ds.raw.model.sections=[shape];e.sectionAngle=angle;const p=S.properties(shape,angle);e.properties.area=p.area;e.properties.momentOfInertia=p.momentOfInertia;e.sectionProperties=p;
 const stored=ds.resultsByElementId.get(e.id),length=e.length;
 // This synthetic comparison changes the section after loading: supply its
 // matching extrema explicitly, as a new solver export would do.
 stored.sectionStress=planarStress(S,p,e.type,length,stored.localEndForces,[0,-1000,0,-1000]);
 const r=V.get(ds,e),data=M.sectionResults.get(ds,e),a=angle*Math.PI/180,c=Math.round(Math.cos(a)),s=Math.round(Math.sin(a)),coeff=Array.from({length:6},()=>[0,0,0,0]);
 coeff[0]=data.stress.axial;coeff[4]=data.stress.moment.map(v=>v*s);coeff[5]=data.stress.moment.map(v=>v*c);
 const spatial=P.stress(P.properties(shape),coeff);near(r.sigmaAbsMax,Math.max(Math.abs(spatial.minimum.value),Math.abs(spatial.maximum.value)));
}
for(const name of ['CaseSectionBeam','CaseSpatialFrameSection']){
 const ds=load(name),e=ds.raw.model.elements[0];e.sectionId=0;
 if(ds.raw.model.dimension===3)M.spatialResults.sectionStress=()=>null;
 assert.match(V.get(ds,e).reason,/geometry/);
}
const fresh=load('CaseSectionBeam'),e=fresh.raw.model.elements[0],r=V.get(fresh,e),other=load('CaseSectionBeam');assert.notStrictEqual(V.get(other,other.raw.model.elements[0]),r);
assert.equal(V.collect(fresh).maximum,r.sigmaAbsMax);
console.log(`Normal stress adapter: ${cases} exported members; analytic signs, torsion, interior truss peak, eight shapes/orientations, 2D/3D parity, unavailable data and dataset cache PASS`);

