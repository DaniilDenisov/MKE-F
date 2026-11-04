(function(M,P){
 'use strict';const S=M.spatialResults,profiles=new WeakMap(),validate=M.spatialData.validate;
 M.spatialData.validate=function(raw){const ds=validate(raw),catalog=P.catalog(raw.model);raw.model.elements.forEach(e=>{if(e.sectionId)profiles.set(e,catalog.get(e.sectionId).properties);});
  if(raw.analysis.type==='static')for(const e of raw.model.elements){const r=ds.resultsByElementId.get(e.id);if(!r.sectionStress||Array.isArray(r.sectionStress)&&!r.sectionStress.length)continue;const p=profiles.get(e);if(!p)throw Error('Stress summary without catalog.');const coeff=r.sectionForceCoefficients||Array.from({length:6},(_,i)=>[i===0?-r.localEndForces[0]:0,0,0,0]),stress=P.stress(p,coeff);
   for(const kind of ['minimum','maximum']){const v=r.sectionStress[kind];if(!v||!['value','xi','y','z'].every(k=>Number.isFinite(v[k]))||v.xi<0||v.xi>1||!P.contains(p,v.y,v.z)||Math.abs(v.value-stress[kind].value)>1e-8*Math.max(1,Math.abs(stress[kind].value))||Math.abs(stress.point(v.xi,v.y,v.z).value-v.value)>1e-8*Math.max(1,Math.abs(v.value)))throw Error('Inconsistent section stress summary.');}
  }return ds;
 };
 S.sectionProfile=e=>profiles.get(e);
 S.sectionStress=(e,r)=>{const p=profiles.get(e);return p&&r?P.stress(p,r.coeff):null;};
 const coefficient=S.coefficient,value=S.value,extrema=S.extrema;
 S.coefficient=function(e,r,q){if(!['sigmaMin','sigmaMax'].includes(q))return coefficient(e,r,q);const stress=S.sectionStress(e,r);return stress?{sectionStress:stress,branch:q==='sigmaMin'?'minimum':'maximum'}:null;};
 S.value=(c,t)=>c&&c.sectionStress?c.sectionStress.at(t)[c.branch].value:c?value(c,t):NaN;
 S.extrema=c=>c&&c.sectionStress?c.sectionStress.positions:c?extrema(c):[];
})(window.MKEFPost,window.MKEFSpatialSections);
