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
near(p.area,8);near(p.Iz,2*4**3/12);near(S.properties(rect,90).momentOfInertia,4*2**3/12);
for(const type of ['channel','tee','ibeam'])assert.throws(()=>S.properties({type,dimensions:{b:1,h:1,tw:1,tf:.1}}));
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
