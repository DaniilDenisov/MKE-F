// Run after tests/test_sections.m to compare the two independent runtimes.
const fs=require('fs'),path=require('path'),vm=require('vm'),assert=require('assert/strict');
const root=path.resolve(__dirname,'..'),S=require('../shared/sections.js');
function near(a,b) {assert(Math.abs(a-b)<=1e-8*Math.max(Math.abs(b),1e-14),`${a} != ${b}`);}
for(const row of JSON.parse(fs.readFileSync(path.join(root,'output/section-parity.json'),'utf8'))) {
  const p=S.properties(row.section,row.angle),stress=S.stress(p,113,2,[-40,150,-30,0,0,0],[20,-300,100,200]);
  for(const key of Object.keys(row.properties))near(p[key],row.properties[key]);
  for(const bound of ['minimum','maximum'])for(const key of ['value','xi','eta'])near(stress[bound][key],row.stress[bound][key]);
}
const rect={type:'rectangle',dimensions:{b:2,h:4}},p=S.properties(rect);
for(const root of [1e-12,1-1e-12,0,1]) {
  const endpointStress=S.stress(p,113,1,[0,1,root,0,0,0],[0,0,0,0]);
  assert.deepEqual(endpointStress.positions,[0,1],'Near-endpoint roots must preserve exact endpoints');
}
assert.deepEqual(S.unique([0,1,1-1e-12,.5,1e-12]),[0,.5,1]);
assert.deepEqual(S.unique([0,1,1-1e-8,1e-8]),[0,1e-8,1-1e-8,1]);
near(p.area,8);near(p.Iz,2*4**3/12);near(S.properties(rect,90).momentOfInertia,4*2**3/12);
for(const type of ['channel','tee','ibeam'])assert.throws(()=>S.properties({type,dimensions:{b:1,h:1,tw:1,tf:.1}}));
for(const type of ['square_tube','rectangular_tube','circular_tube']) {
  const section={type,dimensions:{b:2,h:4,d:2,t:.25}},p=S.properties(section);
  const expectedArea=type==='circular_tube' ? Math.PI*(2**2-1.5**2)/4 : type==='square_tube' ? 4-1.5**2 : 8-1.5*3.5;
  const expectedIz=type==='circular_tube' ? Math.PI*(2**4-1.5**4)/64 : type==='square_tube' ? (2**4-1.5**4)/12 : (2*4**3-1.5*3.5**3)/12;
  near(p.area,expectedArea);near(p.Iz,expectedIz);
  for(const angle of [0,90,180,270])near(S.properties(section,angle).momentOfInertia,angle%180 ? p.Iy : p.Iz);
  for(const t of [undefined,0,-.1,1,2,NaN,Infinity])assert.throws(()=>S.properties({type,dimensions:{...section.dimensions,t}}));
}
assert.throws(()=>S.properties({type:'rectangular_tube',dimensions:{b:4,h:2,t:1}}));
assert.throws(()=>S.properties(rect,45));assert.throws(()=>S.properties({type:'angle',dimensions:{}}));
// Pure axial and pure bending signs; no sampled approximation of interior extrema.
let stress=S.stress(p,112,2,[-80,0,80,0]);near(stress.minimum.value,10);near(stress.maximum.value,10);
stress=S.stress(p,113,2,[0,0,-32,0,0,32]);near(S.evaluate(stress.fibers[1],0),-6);near(S.evaluate(stress.fibers[0],0),6);
const tee={type:'tee',dimensions:{b:2,h:4,tw:.2,tf:.4}},t0=S.properties(tee),t180=S.properties(tee,180);
near(t180.etaMin,-t0.etaMax);near(t180.etaMax,-t0.etaMin);near(t180.momentOfInertia,t0.momentOfInertia);
const context={window:{},console,Set,Map};vm.createContext(context);
for(const file of ['preprocessor/js/namespace.js','postprocessor/js/namespace.js','shared/sections.js','shared/support-catalog.js','shared/support-semantics.js','shared/mpc.js','shared/releases.js','preprocessor/js/case-format.js','preprocessor/js/model-edit.js','postprocessor/js/data-model.js','postprocessor/js/geometry.js','postprocessor/js/static-view.js','postprocessor/js/section-view.js'])vm.runInContext(fs.readFileSync(path.join(root,file),'utf8'),context,{filename:file});
const F=context.window.MKEFPre.caseFormat;
const text=fs.readFileSync(path.join(root,'examples/cases/CaseSectionBeam.txt'),'utf8'),model=F.parse(text);
assert.equal(model.sections.length,1);near(model.elements[0].area,.02);assert.equal(F.serialize(F.parse(F.serialize(model))),F.serialize(model));
model.sections[0].name='Сечение, «Т»';assert.equal(F.parse(F.serialize(model)).sections[0].name,model.sections[0].name);
for(const type of ['square_tube','rectangular_tube','circular_tube']) {
  model.sections[0]={id:1,name:type,type,dimensions:{b:.1,h:.2,d:.1,t:.005}};
  for(const angle of [0,90,180,270]) {
    model.elements[0].sectionAngle=angle;
    const restored=F.parse(F.serialize(model)),p=S.properties(model.sections[0],angle);
    assert.equal(restored.sections[0].type,type);near(restored.sections[0].dimensions.t,.005);
    near(restored.elements[0].area,p.area);near(restored.elements[0].momentOfInertia,p.momentOfInertia);
  }
}
assert.throws(()=>F.parse(text+'\nelement_sections\n1\n1,1,90\n'));
assert.throws(()=>F.parse(text.replace('1,1,0','1,2,0')));
for(const name of fs.readdirSync(path.join(root,'output')).filter(n=>/^SectionCase.*\.json$/.test(n)||n==='CaseSectionBeam.json')) {
  const data=JSON.parse(fs.readFileSync(path.join(root,'output',name),'utf8'));
  context.window.MKEFPost.validateDataset(data);
  if(name==='CaseSectionBeam.json') {
    const bad=structuredClone(data);bad.model.elements[0].properties.area*=2;assert.throws(()=>context.window.MKEFPost.validateDataset(bad));
    const badAngle=structuredClone(data);badAngle.model.elements[0].sectionAngle=45;assert.throws(()=>context.window.MKEFPost.validateDataset(badAngle));
  }
}
console.log('Sections: geometry, stress parity, signs, orientations, TXT roundtrip, v6 combinations PASS');
