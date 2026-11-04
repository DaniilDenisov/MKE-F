const fs=require('node:fs'),assert=require('node:assert/strict');
const read=name=>JSON.parse(fs.readFileSync('output/result-numbers/'+name+'.json','utf8')),reference=read('reference');
const hex=x=>{const b=Buffer.alloc(8);b.writeDoubleBE(x);return b.toString('hex');};
let count=0;
function bits(values,expected){assert.equal(values.length,expected.length);values.forEach((x,i)=>{if(Array.isArray(x))bits(x,expected[i]);else{assert.equal(hex(x),expected[i]==='8000000000000000'?'0000000000000000':expected[i]);count++;}});}
for(const name of ['compact','pretty']){
 const d=read(name);bits(d.values,reference.values);assert.deepEqual(d.singleton,[[1e-40]]);assert.deepEqual(d.empty,[]);assert.deepEqual(d.flags,[true,false]);assert.equal(d.text,'Quote " slash \\ newline\nUnicode: тест');
}
const modal=read('tiny-modal'),op=modal.model.elements[0].recoveryOperators;
bits(op.K,reference.K);bits(op.M,reference.M);bits(modal.analysis.modeShapes,reference.modeShapes);
const transient=read('tiny-transient');
for(const key of ['displacements','velocities','accelerations','loadHistory','reactions'])bits(transient.analysis[key],reference[key]);
console.log(`JSON numerical transport: ${count} IEEE bit checks, subnormals, normalized signed zero, operators and histories PASS`);
