(function (M) {
  'use strict';
  var S=window.MKEFSections, V=M.sectionResults={};
  V.isStress=function (mode) {return mode==='sigmaMin' || mode==='sigmaMax';};
  V.get=function (dataset,element) {
    if(!element.sectionId) return null;
    if(!dataset.sectionStressCache) dataset.sectionStressCache=new Map();
    if(dataset.sectionStressCache.has(element.id)) return dataset.sectionStressCache.get(element.id);
    var section=dataset.raw.model.sections.find(function (s) {return s.id===element.sectionId;}), p=S.resultProperties(section,element);
    var a=dataset.nodesById.get(element.nodeIds[0]),b=dataset.nodesById.get(element.nodeIds[1]),L=Math.hypot(b.x-a.x,b.y-a.y),stress;
    var r=dataset.resultsByElementId.get(element.id),c=r.sectionForceCoefficients;
    stress={...r.sectionStress,axial:c[0],moment:c[2],fibers:[p.etaMin,p.etaMax].map(eta=>c[0].map((v,i)=>v/p.area-c[2][i]*eta/p.momentOfInertia))};
    var result={section:section,properties:p,stress:stress,length:L};dataset.sectionStressCache.set(element.id,result);return result;
  };
  V.value=function (data,mode,xi) {var v=data.stress.fibers.map(function (p) {return S.evaluate(p,xi);});return mode==='sigmaMin' ? Math.min.apply(null,v) : Math.max.apply(null,v);};
  V.samples=function (dataset,element,mode,count) {
    var data=V.get(dataset,element);if(!data)return [];
    var positions=data.stress.positions.slice(),n=Math.max(5,count || M.config.defaultFrameSamples);
    for(var i=0;i<n;i++) positions.push(i/(n-1));
    return S.unique(positions).map(function (xi) {return {xi:xi,value:V.value(data,mode,xi)};});
  };
  V.details=function (dataset,element) {
    var data=V.get(dataset,element);if(!data)return ['Section stresses: unavailable — assign a section and rerun.'];
    var unit=dataset.raw.metadata.units.stress || '',s=data.stress;
    function line(label,p) {return label+' = '+M.staticResults.withUnit(p.value,unit)+' at x/L='+M.staticResults.format(p.xi)+' (x='+M.staticResults.format(p.xi*data.length)+'), η='+M.staticResults.format(p.eta);}
    return ['Section: '+data.section.name+' · '+element.sectionAngle+'°',line('σmin',s.minimum),line('σmax',s.maximum)];
  };
}(window.MKEFPost));
