const assert=require('node:assert/strict'),fs=require('node:fs');
const S=require('../shared/sections.js'),P=require('../shared/spatial-sections.js');
const near=(a,b)=>assert(Math.abs(a-b)<=1e-11*Math.max(Math.abs(b),1e-15),`${a} != ${b}`);
for(const row of JSON.parse(fs.readFileSync('output/spatial-section-parity.json','utf8'))){
 const p=P.properties(row.section);
 for(const k of ['area','centroidY','centroidZ','Iy','Iz','torsionConstant'])near(p[k],row.properties[k]);
 const scaled=structuredClone(row.section);for(const k in scaled.dimensions)scaled.dimensions[k]*=3;
 near(P.properties(scaled).torsionConstant,p.torsionConstant*81);
 assert(p.loops.length);for(const loop of p.loops)for(const v of loop)assert(P.contains(p,v.y,v.z));
 if(row.section.type.includes('tube'))assert(!P.contains(p,0,0));
 if(row.section.dimensions.inner_radius){
   assert(!P.contains(p,p.roundedTube.h/2,p.roundedTube.b/2));
   assert(P.contains(p,p.roundedTube.h/2,0));
   const d=p.roundedTube,r=d.ro,coeff=Array.from({length:6},()=>[0,0,0,0]);coeff[4][0]=p.Iy;coeff[5][0]=-p.Iz;
   const peak=P.stress(p,coeff).maximum;
   near(peak.value,d.h/2+d.b/2-2*r+Math.SQRT2*r);assert(P.contains(p,peak.y,peak.z));
   for(const inner_radius of [-.01,NaN,Infinity,.045])assert.throws(()=>P.properties({...row.section,dimensions:{...row.section.dimensions,inner_radius}}));
 }
 assert.throws(()=>P.properties({...row.section,torsion:{mode:'manual',value:0}}));
 near(P.properties({...row.section,torsion:{mode:'manual',value:.123}}).torsionConstant,.123);
 const stress=P.stress(p,row.coefficients);for(const k of ['minimum','maximum'])near(stress[k].value,row.stress[k].value);
 for(let i=0;i<=500;i++){const at=stress.at(i/500);assert(at.minimum.value>=stress.minimum.value-1e-7);assert(at.maximum.value<=stress.maximum.value+1e-7);}
 const axial=Array.from({length:6},(_,i)=>[i===0?100:0,0,0,0]);near(P.stress(p,axial).maximum.value,100/p.area);
 if(p.roundedTube){const peak=P.stress(p,axial).maximum;assert(P.contains(p,peak.y,peak.z));}
 axial[0][0]=0;axial[3][0]=1000;near(P.stress(p,axial).maximum.value,0);
}
near(P.properties({type:'rectangle',dimensions:{b:1,h:1}}).torsionConstant,.140577014955);
{
 const p=P.properties({type:'circle',dimensions:{d:2}}),coeff=Array.from({length:6},()=>[0,0,0,0]);coeff[0][0]=2;coeff[4][0]=3;coeff[5][0]=4;
 const stress=P.stress(p,coeff);near(stress.maximum.value,2/Math.PI+20/Math.PI);near(stress.maximum.y,-.8);near(stress.maximum.z,.6);near(p.torsionConstant,Math.PI/2);
 const rect=P.properties({type:'rectangle',dimensions:{b:2,h:4}});coeff[0][0]=0;coeff[4][0]=0;coeff[5]=[-.123456,1,0,0];const envelope=P.stress(rect,coeff);assert(envelope.positions.some(x=>Math.abs(x-.123456)<1e-12));assert(Math.abs(envelope.at(.123456).maximum.value)<1e-14);
 const poses=P.positions(Array(100000).fill({}),99999,10,.123);assert.equal(poses.length,600);assert.equal(poses[0].index,99999);assert.equal(poses.filter(p=>p.selected).length,1);assert(poses.some(p=>p.selected&&p.xi===.123));
}
console.log('Spatial profiles, J parity, material boundaries PASS');
const vm=require('node:vm'),ctx={window:{MKEFSpatial:{registry:require('../shared/element-types.json')}},TextEncoder};vm.createContext(ctx);
for(const f of ['shared/sections.js','shared/spatial-sections.js','shared/spatial-conditions.js','preprocessor/js/namespace.js','preprocessor/js/case-format.js','preprocessor/js/spatial-format.js','preprocessor/js/spatial-edit.js'])vm.runInContext(fs.readFileSync(f,'utf8'),ctx);
const M=ctx.window.MKEFPre,F=M.spatialFormat;
for(const code of [1133,1436]){
 const m=F.empty(code);m.nodes=[{x:0,y:0,z:0},{x:2,y:0,z:0}];
 m.sections=[{id:1,name:'Rectangle',type:'rectangle',dimensions:{b:.1,h:.2}}];
 const h=new M.spatialEdit.History(m),e={code,node1:1,node2:2,sectionId:1,youngsModulus:2e11,density:7850,shearModulus:8e10,referenceY:[0,1,0],angleDegrees:-37.5};
 const before=F.serialize(m);
 assert.throws(()=>h.member(null,{...e,sectionId:2}),/Unknown section/);
 assert.throws(()=>h.member(null,{...e,density:NaN}),/positive and finite/);
 assert.equal(h.undo.length,0);assert.equal(F.serialize(m),before);
 h.member(null,e);assert.equal(h.undo.length,1);near(m.elements[0].area,.02);
 if(code===1436){near(m.elements[0].Iy,P.properties(m.sections[0]).Iy);near(m.elements[0].torsionConstant,P.properties(m.sections[0]).torsionConstant);assert.equal(m.elements[0].angleDegrees,-37.5);}
 const saved=F.serialize(m);assert.equal(F.parse(saved).elements[0].sectionId,1);
 assert.throws(()=>h.member(null,e),/already/);assert.equal(h.undo.length,1);
 h.travel(true);assert.equal(F.serialize(m),before);h.travel(false);assert.equal(F.serialize(m),saved);
}
console.log('Creation with catalog properties, atomic rejection, TXT and undo/redo PASS');
for(const name of ['CaseSpatialFrame','CaseSpatialFrameModal','CaseSpatialFrameTransient','CaseSpatialTruss']){
 const source=fs.readFileSync('output/'+name+'Section.txt','utf8'),m=F.parse(source),h=new M.spatialEdit.History(m),original=F.serialize(m);
 assert.equal(h.text(),source);assert.equal(F.serialize(F.parse(original)),original);
 const s=structuredClone(m.sections[0]),oldA=m.elements[0].area;s.dimensions.b*=1.2;h.section(0,s);
 assert.notEqual(m.elements[0].area,oldA);assert.equal(h.undo[0].patches.length,1);assert.equal(h.undo[0].patches[0].key,'sections');
 h.travel(true);assert.equal(F.serialize(m),original);h.travel(false);assert.notEqual(m.elements[0].area,oldA);
 assert.throws(()=>h.section(0,null),/used/);const area=m.elements[0].area;h.assignSection(0,0,-12.5);assert.equal(m.elements[0].area,area);
 const detached=F.parse(F.serialize(m));assert.equal(detached.elements[0].sectionId||0,0);assert.equal(detached.elements[0][m.elementType===1436?'angleDegrees':'sectionAngle'],-12.5);
 h.travel(true);assert.equal(m.elements[0].sectionId,s.id);
}
console.log('Spatial catalog assignment, detachment, TXT and delta history PASS');
for(const type of ['channel','tee','ibeam']){
 const source=fs.readFileSync('output/LegacyFrame-'+type+'.txt','utf8'),m=F.parse(source),h=new M.spatialEdit.History(m),before=F.serialize(m);
 assert.equal(F.validate(m,false).length,0);assert(F.validate(m,true).some(s=>s.includes('1 ('+type+')')));
 assert.equal(F.serialize(F.parse(before)),before);assert.equal(h.text(),source);
 const sid=m.sections[0].id,section=structuredClone(m.sections[0]);
 assert.throws(()=>h.assignSection(0,sid,0),/Frames 1436/);
 assert.throws(()=>h.section(null,{...section,id:99}),/Frames 1436/);
 m.nodes.push({x:4,y:0,z:0});
 assert.throws(()=>h.member(null,{...m.elements[0],node1:2,node2:3}),/Frames 1436/);
 assert.equal(h.undo.length,0);
 // Editing a legacy member and profile is permitted without changing its assignment.
 h.member(0,{...m.elements[0],youngsModulus:210e9});h.travel(true);
 h.section(0,{...section,name:'Legacy renamed'});h.travel(true);
 h.section(0,{...section,type:'rectangle',dimensions:{b:.1,h:.2}});
 assert.equal(F.validate(m,true).length,0);
 h.travel(true);assert(F.validate(m,true).some(s=>s.includes('Frames 1436')));
 h.travel(false);assert.equal(F.validate(m,true).length,0);
}
// One shared catalog policy leaves every shape usable in trusses and 2D.
for(const code of [112,113,1133])for(const type of Object.keys(S.types))assert(P.allowedForElement(code,type));
for(const type of Object.keys(S.types)){
 const m=F.parse(fs.readFileSync('examples/cases/CaseSpatialTruss.txt','utf8')),h=new M.spatialEdit.History(m);
 h.section(null,{id:1,name:type,type,dimensions:{b:.1,h:.2,d:.1,t:.005,tw:.008,tf:.012}});
 h.assignSection(0,1,0);assert.equal(F.validate(m,true).length,0);
}
console.log('Frame restrictions, legacy editing, round-trip, atomic rejection, repair and undo PASS');
{
 const m=F.parse(fs.readFileSync('output/CaseSpatialFrameSection.txt','utf8'));m.elements[0].angleDegrees=0;m.loads=[];
 m.releases=[{elementId:1,end:2,component:'My'}];m.mpcs=[{depNode:2,depDOF:6,rhs:0,masters:[{node:1,dof:6,coefficient:1}]}];
 const h=new M.spatialEdit.History(m),before=F.serialize(m);assert.throws(()=>h.assignSection(0,m.sections[0].id,90),/inactive|absent|MPC/i);assert.equal(F.serialize(m),before);assert.equal(h.undo.length,0);
 console.log('Section rotation rejects incompatible projected MPC atomically PASS');
}
