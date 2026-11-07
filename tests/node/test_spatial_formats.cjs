// Spatial TXT is viewable; spatial v8 results require the spatial viewer module.
const fs = require('fs'), path = require('path'), vm = require('vm'), assert = require('assert/strict');
const root = path.resolve(__dirname, '../..');
const context = {window:{MKEFSpatial:{registry:require('../../shared/element-types.json')}}, console, Set, Map, TextEncoder}; vm.createContext(context);
vm.runInContext(fs.readFileSync(path.join(root,'shared/spatial-conditions.js'),'utf8'),context);
for (const file of ['preprocessor/js/namespace.js','postprocessor/js/namespace.js','preprocessor/js/case-format.js','preprocessor/js/spatial-format.js','postprocessor/js/data-model.js']) {
  vm.runInContext(fs.readFileSync(path.join(root,file),'utf8'), context, {filename:file});
}
const registry = require('../../shared/element-types.json');
assert.deepEqual(registry.types.map(t=>t.code), [112,113,1133,1436]);
assert.equal(registry.types.find(t=>t.code===113).dimension, 2);
assert.equal(registry.types.find(t=>t.code===1133).dimension, 3);
for (const type of registry.types) assert.equal(type.dofLabels.length, type.nodalLoads.length);
for (const name of ['CaseSpatialTruss','CaseSpatialFrame','CaseSpatialTrussModal','CaseSpatialFrameModal','CaseSpatialTrussTransient','CaseSpatialFrameTransient']) {
  const text = fs.readFileSync(path.join(root,'examples/cases',name+'.txt'),'utf8');
  const model = context.window.MKEFPre.caseFormat.parse(text);
  assert.equal(model.dimension,3); assert.equal(model.sourceText,text);
  const data=JSON.parse(fs.readFileSync(path.join(root,'output',name+'.json'),'utf8'));
  assert.throws(()=>context.window.MKEFPost.validateDataset(data), /Spatial viewer is not loaded/);
  assert(Array.isArray(data.model.elements)); assert(Array.isArray(data.model.supports));
  if (data.analysis.type==='modal') assert(Array.isArray(data.analysis.modeShapes[0]));
}
console.log('Spatial format and UI boundary checks passed.');
