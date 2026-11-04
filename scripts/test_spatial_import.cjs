const fs = require('node:fs'), vm = require('node:vm'), assert = require('node:assert/strict');
const ctx = {window:{MKEFSpatial:{registry:require('../shared/element-types.json')}},TextEncoder}; vm.createContext(ctx);
vm.runInContext(fs.readFileSync('shared/spatial-conditions.js','utf8'),ctx);
for (const f of ['sections','spatial-sections']) vm.runInContext(fs.readFileSync('shared/'+f+'.js','utf8'),ctx);
for (const f of ['namespace','case-format','spatial-format']) vm.runInContext(fs.readFileSync('preprocessor/js/'+f+'.js','utf8'),ctx);
const F = ctx.window.MKEFPre.caseFormat;
for (const name of fs.readdirSync('examples/cases').filter(f=>/^CaseSpatial.*\.txt$/.test(f))) {
  const source = fs.readFileSync('examples/cases/'+name,'utf8'), model = F.parse(source);
  assert.equal(model.sourceText,source); assert.equal(model.dimension,3); assert(model.nodes.length); assert(model.elements.length);
}
const sample = fs.readFileSync('examples/cases/CaseSpatialFrame.txt','utf8');
for (const text of [sample+'\nunknown\n0',sample.replace('0.01,2e11','-0.01,2e11'),sample.replace('0,1,0,0','1,0,0,0'),sample+'\nelems_1133\n0']) assert.throws(()=>F.parse(text),/Line \d+/);
assert.throws(()=>F.parse(sample.replace('2,0,0','NaN,0,0')),/finite/);
const conditions = '\nmpc\n1\n2,2,0,1,1,2,1\nreleases\n3\n1,2,Mx\n1,2,My\n1,2,Mz\neload_uniform_3d\n1\n20,1,1,0,-2,3\neload_linear_3d\n1\n21,1,2,0,1,2,0,3,4\n';
const withConditions = F.parse(sample+conditions.replaceAll('1,2,M','1,1,M'));
assert.equal(withConditions.mpcs.length,1); assert.equal(withConditions.releases.length,3); assert.equal(withConditions.elementLoads.length,2);
const transient = sample.replace('static','transient,0.01,1,2,6').replace('bcforce_stat_3d','bcforce_pulse_3d').replace('10,2,','12,2,')+'\nbcforce_pulse_3d\n1\n12,2,1,2,3,4,5,6\nbcforce_step_3d\n1\n13,2,1,2,3,4,5,6';
assert.equal(F.parse(transient).loads.length,3);
for (const tail of ['mpc\n1\n2,7,0,1,1,1,1','mpc\n1\n2,2,0,2,1,2,1','releases\n1\n2,1,Mz','releases\n1\n1,3,Mx','bcfix_3d\n1\n1,2,0,0,0,0,0','bcforce_stat_3d\n1\n10,3,0,0,0,0,0,0','eload_uniform_3d\n1\n20,1,3,0,1,0','nodes\n1\n0,0,0']) assert.throws(()=>F.parse(sample+'\n'+tail),/Line \d+/);
assert.throws(()=>F.parse(sample.replace('nodes\n2','nodes\n100001')),/limit/);
// Non-zero Z must not switch the display dimension.
const nonzeroZ = sample.replace('elems_1436','elems_113').split('bcfix_3d')[0].replace(/1436,1,2,.*/, '113,1,2,0.01,2e11,7850,8e-6').replace('2,0,0','2,0,9');
assert.equal(ctx.window.MKEFPre.spatialFormat.detects(nonzeroZ),false);
console.log('Spatial imports, source preservation and malformed records passed.');
