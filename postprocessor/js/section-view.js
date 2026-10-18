(function (M) {
  'use strict';
  var S=window.MKEFSections, V=M.sectionResults={};
  V.isStress=function (mode) {return mode==='sigmaMin' || mode==='sigmaMax';};
  V.get=function (dataset,element) {
    if(!element.sectionId) return null;
    if(!dataset.sectionStressCache) dataset.sectionStressCache=new Map();
    if(dataset.sectionStressCache.has(element.id)) return dataset.sectionStressCache.get(element.id);
    var section=dataset.raw.model.sections.find(function (s) {return s.id===element.sectionId;}), p=S.properties(section,element.sectionAngle);
    var a=dataset.nodesById.get(element.nodeIds[0]),b=dataset.nodesById.get(element.nodeIds[1]),L=Math.hypot(b.x-a.x,b.y-a.y),q=M.staticResults.elementIntensities(dataset,element);
    var stress=S.stress(p,element.type,L,dataset.resultsByElementId.get(element.id).localEndForces,[q.x,q.y,q.x+q.dx,q.y+q.dy]);
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
  var lastElement=null,lastDataset=null,position=.5,resizeObserver=null;
  V.show=function (dataset,selection) {
    if(resizeObserver) {resizeObserver.disconnect();resizeObserver=null;}
    var panel=document.getElementById('section-result-panel');if(!panel)return;
    var element=selection && selection.kind==='element' && dataset && dataset.elementsById.get(selection.id);
    panel.hidden=!element || dataset.raw.analysis.type!=='static' || !element.sectionId;
    if(panel.hidden) {
      lastElement=null;
      document.getElementById('section-result-svg').replaceChildren();
      var exportButton=document.getElementById('export-section-svg');
      exportButton.disabled=true;exportButton.onclick=null;
      return;
    }
    if(lastDataset!==dataset || lastElement!==element.id)position=.5;
    lastDataset=dataset;lastElement=element.id;
    var data=V.get(dataset,element),svg=document.getElementById('section-result-svg'),info=document.getElementById('section-result-info'),range=document.getElementById('section-position'),input=document.getElementById('section-position-number');
    range.disabled=input.disabled=!data;svg.toggleAttribute('hidden',!data);
    document.getElementById('export-section-svg').disabled=!data;
    if(!data) {info.textContent='Assign a section in the preprocessor and rerun to view bending stresses.';return;}
    function draw() {
      range.value=input.value=String(position);
      window.MKEFSectionDrawing.preview(svg,data.section,element.sectionAngle,data.stress,position,dataset.raw.metadata.units.stress);
      var a=V.value(data,'sigmaMin',position),b=V.value(data,'sigmaMax',position),u=dataset.raw.metadata.units.stress || '';
      var limit=Math.max(Math.abs(data.stress.minimum.value),Math.abs(data.stress.maximum.value));
      info.textContent=data.section.name+' · '+element.sectionAngle+'° · x/L='+M.staticResults.format(position)+'\nσmin='+M.staticResults.withUnit(a,u)+' · σmax='+M.staticResults.withUnit(b,u)+'\nScale: '+M.staticResults.withUnit(-limit,u)+' … '+M.staticResults.withUnit(limit,u)+'\nRed: tension (+) · blue: compression (−). Scale fixed over this element.';
    }
    function change(event) {var v=Number(event.target.value);if(event.target.value.trim()==='' || !Number.isFinite(v) || v<0 || v>1) {event.target.setCustomValidity('Use a position between 0 and 1.');return;} event.target.setCustomValidity('');position=v;draw();}
    range.oninput=change;input.oninput=change;draw();
    if(typeof ResizeObserver!=='undefined') {resizeObserver=new ResizeObserver(function () {draw();});resizeObserver.observe(svg);}
    document.getElementById('export-section-svg').onclick=function () {
      var text=new XMLSerializer().serializeToString(svg),blob=new Blob([text],{type:'image/svg+xml'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download='section-E'+element.id+'.svg';a.click();setTimeout(function () {URL.revokeObjectURL(url);},1000);
    };
  };
}(window.MKEFPost));
