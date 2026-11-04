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
