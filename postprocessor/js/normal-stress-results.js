(function(M,S,P){
  'use strict';
  const V=M.normalStressResults={},datasets=new WeakMap();
  const finiteRow=row=>Array.isArray(row)&&row.length===4&&row.every(Number.isFinite);
  const unavailable=reason=>({status:'unavailable',reason});
  const governing=(min,max)=>Math.abs(min.value)>Math.abs(max.value)?min:max;
  function section(min,max){
    const p=governing(min,max);
    return {xi:p.xi,sigmaMin:min.value,sigmaMax:max.value,sigmaAbsMax:Math.abs(p.value),sigmaX:p.value,eta:p.eta??null,y:p.y??null,z:p.z??null};
  }
  function recover(ds,e){
    if(ds.raw.analysis.type!=='static')return unavailable('Normal stress coloring is available only for static results.');
    const spatial=ds.raw.model.dimension===3,truss=[112,1133].includes(e.type),r=ds.resultsByElementId.get(e.id);
    if(!r)return unavailable('Element force results are missing.');
    if(!Number.isFinite(e.properties?.area)||e.properties.area<=0)return unavailable('A finite positive area is required.');
    let stress,at;
    if(truss){
      let n;
        if(r.sectionForceCoefficients){if(!finiteRow(r.sectionForceCoefficients[0]))return unavailable('Axial force coefficients are missing or non-finite.');n=r.sectionForceCoefficients[0];}
        else {if(!Number.isFinite(r.localEndForces?.[0]))return unavailable('Axial end force is missing or non-finite.');n=[-r.localEndForces[0],0,0,0];}
      const poly=n.map(v=>v/e.properties.area);
      // Uniform axial stress has no unique critical fiber. Geometry is optional.
      const point=xi=>({xi,value:S.evaluate(poly,xi)});
      stress=r.sectionStress;
      if(!stress?.minimum||!stress?.maximum||!Array.isArray(stress.positions))return unavailable('Static normal stress summary was not exported.');
      at=xi=>{const p=point(xi);return section(p,p);};
    }else if(spatial){
      const c=r.sectionForceCoefficients;
      if(!c||![0,4,5].every(i=>finiteRow(c[i])))return unavailable('Normal force and bending coefficients were not exported or are non-finite.');
      stress=M.spatialResults.sectionStress(e,{coeff:c});
      if(!stress)return unavailable('Section geometry is required for axial and bending normal stresses.');
      at=xi=>{const p=stress.at(xi);return section(p.minimum,p.maximum);};
    }else{
      if(!e.sectionId)return unavailable('Section geometry is required for axial and bending normal stresses.');
      const data=M.sectionResults.get(ds,e);stress=data.stress;
      at=xi=>{const points=stress.fibers.map((poly,i)=>({xi,value:S.evaluate(poly,xi),eta:i?data.properties.etaMax:data.properties.etaMin}));return section(points.reduce((a,b)=>b.value<a.value?b:a),points.reduce((a,b)=>b.value>a.value?b:a));};
    }
    if(![stress.minimum.value,stress.maximum.value].every(Number.isFinite))return unavailable('Recovered normal stresses are non-finite.');
    const critical=governing(stress.minimum,stress.maximum),positions=S.unique([0,1,...stress.positions,...Array.from({length:41},(_,i)=>i/40)]),samples=positions.map(at);
    if(samples.some(p=>![p.sigmaMin,p.sigmaMax,p.sigmaAbsMax].every(Number.isFinite)))return unavailable('Recovered normal stresses are non-finite.');
    return {status:'ok',reason:'',sigmaMin:stress.minimum.value,sigmaMax:stress.maximum.value,sigmaAbsMax:Math.abs(critical.value),minimum:stress.minimum,maximum:stress.maximum,peak:{...at(critical.xi),sigmaX:critical.value,eta:critical.eta??null,y:critical.y??null,z:critical.z??null},positions,samples,at};
  }
  V.get=function(ds,e){
    let cache=datasets.get(ds);if(!cache){cache=new Map();datasets.set(ds,cache);}
    if(!cache.has(e.id)){
      let result;try{result=recover(ds,e);}catch(error){result=unavailable('Normal stress recovery unavailable: '+error.message);}
      cache.set(e.id,result);
    }
    return cache.get(e.id);
  };
  V.clear=ds=>datasets.delete(ds);
  V.collect=ds=>{
    const results=new Map();let maximum=0;
    for(const e of ds.raw.model.elements){const r=V.get(ds,e);results.set(e.id,r);if(r.status==='ok')maximum=Math.max(maximum,r.sigmaAbsMax);}
    return {results,maximum};
  };
  V.at=(r,xi)=>r?.status==='ok'?r.at(xi).sigmaAbsMax:null;
}(window.MKEFPost,window.MKEFSections,window.MKEFSpatialSections));
