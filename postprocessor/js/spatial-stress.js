(function(M,P){
 'use strict';
 const S=M.spatialResults,sections=window.MKEFSections,profiles=new WeakMap(),stored=new WeakMap(),validate=M.spatialData.validate;
 M.spatialData.validate=function(raw){
  const ds=validate(raw),catalog=new Map((raw.model.sections||[]).map(s=>[s.id,s]));
  for(const e of raw.model.elements){
   if(!e.sectionId)continue;
   const p=P.decorate(catalog.get(e.sectionId),sections.resultProperties(catalog.get(e.sectionId),e));profiles.set(e,p);
   if(raw.analysis.type==='static'){
    const r=ds.resultsByElementId.get(e.id),coeff=r.sectionForceCoefficients||Array.from({length:6},(_,i)=>[i===0?-r.localEndForces[0]:0,0,0,0]),field=P.stressField(p,coeff);
    sections.validateStress(r.sectionStress,v=>field.point(v.xi,v.y,v.z).value,v=>P.contains(p,v.y,v.z));
    stored.set(e,r);
   }
  }
  return ds;
 };
 S.sectionProfile=e=>profiles.get(e);
 S.sectionStress=(e,r)=>{
  const p=profiles.get(e);if(!p||!r)return null;
  const result=stored.get(e);
  if(!result)return P.stress(p,r.coeff);
  const {minimum,maximum,positions}=result.sectionStress;
  return {...P.stressField(p,r.coeff),minimum,maximum,positions};
 };
 const coefficient=S.coefficient,value=S.value,extrema=S.extrema;
 S.coefficient=function(e,r,q){if(!['sigmaMin','sigmaMax'].includes(q))return coefficient(e,r,q);const stress=S.sectionStress(e,r);return stress?{sectionStress:stress,branch:q==='sigmaMin'?'minimum':'maximum'}:null;};
 S.value=(c,t)=>c&&c.sectionStress?c.sectionStress.at(t)[c.branch].value:c?value(c,t):NaN;
 S.extrema=c=>c&&c.sectionStress?c.sectionStress.positions:c?extrema(c):[];
})(window.MKEFPost,window.MKEFSpatialSections);
