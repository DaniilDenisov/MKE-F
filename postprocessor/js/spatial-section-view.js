(function(M,P){
 'use strict';const proto=M.SpatialPost.prototype,S=M.spatialResults,U=M.spatialUI;
 const controls=proto.makeControls;proto.makeControls=function(){controls.call(this);this.sectionCatalog=P.resultCatalog(this.model);
  this.sectionContours=U.check(this.settings,'section-contours','Section outlines',false);this.sectionCount=U.field(this.settings,'section-count','Positions per member (1–10)',1);this.sectionCount.step=1;
  this.sectionScale=U.field(this.settings,'section-scale','Outline scale',1);this.sectionNote=U.el('p','',this.settings);
  for(const e of [this.sectionContours,this.sectionCount,this.sectionScale])e.onchange=()=>this.update();
  for(const [key,label]of [['sigmaMin','σmin'],['sigmaMax','σmax']]){const o=U.el('option',label,this.quantity);o.value=key;}
  this.cut=U.el('section',undefined,this.panel);this.cut.id='spatial-section-cut';U.el('h3','Normal stress σx · N, My, Mz',this.cut);
  this.cutSVG=M.svgElement('svg',{id:'spatial-section-svg',viewBox:'0 0 500 420',width:500,height:420,role:'img','aria-label':'Section: Y up, Z right, viewed from −X'});this.cutSVG.style.maxWidth='100%';this.cut.append(this.cutSVG);
  this.pointY=U.field(this.cut,'point-y','Material point · Y',0);this.pointZ=U.field(this.cut,'point-z','Material point · Z',0);this.pointNote=U.el('p','',this.cut);this.cutNote=U.el('p','',this.cut);
  const apply=U.el('button','Show point stress',this.cut);apply.id='spatial-section-point';apply.onclick=()=>{this.drawSection();this.drawChart();};
  this.cutSVG.onclick=event=>{if(!this.cutTransform)return;const point=this.cutSVG.createSVGPoint();point.x=event.clientX;point.y=event.clientY;const local=point.matrixTransform(this.cutSVG.getScreenCTM().inverse()),t=this.cutTransform;this.pointY.value=(t.cy-local.y)/t.scale;this.pointZ.value=(local.x-t.cx)/t.scale;this.drawSection();this.drawChart();};
  for(const kind of ['svg','csv']){const b=U.el('button',kind.toUpperCase()+' section view',this.cut);b.id='spatial-section-export-'+kind;b.onclick=()=>this.exportSection(kind);}
  if(this.historyKind){const o=U.el('option','Material point σx',this.historyKind);o.value='pointStress';}
 };
 const update=proto.update;proto.update=function(){update.call(this);this.updateSections();this.drawSection();};
 const unit=proto.unit;proto.unit=function(q){return q==='sigmaMin'||q==='sigmaMax'?this.ds.raw.metadata.units.stress||'':unit.call(this,q);};
 proto.drawSection=function(){
  const selected=this.selection&&this.selection.kind==='element',e=selected?this.model.elements[this.selection.index]:null,p=e&&S.sectionProfile(e),r=p&&S.recover(this.ds,e,this.vector,this.acceleration),stress=r&&S.sectionStress(e,r),xi=Number(this.position.value);
  this.cutSVG.replaceChildren();this.cutTransform=null;this.boundaryRows=null;
  if(!stress||!Number.isFinite(xi)||xi<0||xi>1){this.cutNote.textContent=!e?'Select a member.':!p?'No catalog section is available. N/A is available; the profile is not inferred from A/I/J.':'Required data was not exported, or x/L is outside [0,1].';this.pointNote.textContent='';return;}
  if(this.analysis.type==='transient'&&!this.prepared){this.cutNote.textContent='The stress scale requires completed preparation of the entire saved history.';this.pointNote.textContent='';return;}
  const at=stress.at(xi),f=U.fmt,maximum=this.stressMaximum===undefined?Math.max(Math.abs(stress.minimum.value),Math.abs(stress.maximum.value)):this.stressMaximum,points=p.loops.flat();
  const minY=Math.min(...points.map(v=>v.y)),maxY=Math.max(...points.map(v=>v.y)),minZ=Math.min(...points.map(v=>v.z)),maxZ=Math.max(...points.map(v=>v.z)),scale=Math.min(360/(maxZ-minZ),270/(maxY-minY)),cx=250-scale*(minZ+maxZ)/2,cy=190+scale*(minY+maxY)/2;
  this.cutTransform={scale,cx,cy};const project=(y,z)=>[cx+scale*z,cy-scale*y],svg=this.cutSVG,append=(tag,attrs,text,parent=svg)=>{const item=M.svgElement(tag,attrs);if(text!==undefined)item.textContent=text;parent.append(item);return item;};
  append('rect',{width:500,height:420,fill:'white'});const defs=append('defs',{}),gradient=append('linearGradient',{id:'spatial-stress-gradient',gradientUnits:'userSpaceOnUse'},undefined,defs);
  const a=window.MKEFSections.evaluate(stress.a,xi),b=window.MKEFSections.evaluate(stress.b,xi),c=window.MKEFSections.evaluate(stress.c,xi),norm=b*b+c*c;
  let fill='#ffffff';if(norm&&maximum){const start=project((-maximum-a)*b/norm,(-maximum-a)*c/norm),end=project((maximum-a)*b/norm,(maximum-a)*c/norm);['x1','y1','x2','y2'].forEach((key,i)=>gradient.setAttribute(key,start.concat(end)[i]));append('stop',{offset:0,'stop-color':'#2166ac'},undefined,gradient);append('stop',{offset:.5,'stop-color':'#ffffff'},undefined,gradient);append('stop',{offset:1,'stop-color':'#c61b32'},undefined,gradient);fill='url(#spatial-stress-gradient)';}
  else if(maximum){const t=Math.min(1,Math.abs(a)/maximum),base=a>=0?[198,27,50]:[33,102,172];fill='rgb('+base.map(v=>Math.round(255+(v-255)*t)).join(',')+')';}
  const path=p.loops.map(loop=>loop.map((v,i)=>(i?'L':'M')+project(v.y,v.z).join(',')).join(' ')+'Z').join(' ');
  append('path',{d:path,'fill-rule':'evenodd',fill,stroke:'#34495e','stroke-width':1.2});const clip=append('clipPath',{id:'spatial-stress-clip'},undefined,defs);append('path',{d:path,'clip-rule':'evenodd'},undefined,clip);
  if(norm){const center=[-a*b/norm,-a*c/norm],length=2*Math.hypot(maxY-minY,maxZ-minZ)/Math.sqrt(norm),v1=project(center[0]-c*length,center[1]+b*length),v2=project(center[0]+c*length,center[1]-b*length);append('line',{x1:v1[0],y1:v1[1],x2:v2[0],y2:v2[1],stroke:'#333','stroke-dasharray':'5 4','clip-path':'url(#spatial-stress-clip)'});}
  for(const key of ['minimum','maximum']){const v=at[key],q=project(v.y,v.z);append('circle',{cx:q[0],cy:q[1],r:4,fill:key==='maximum'?'#c61b32':'#2166ac',stroke:'white'});}
  append('text',{x:12,y:20,fill:'#263847','font-size':14},`Y ↑ · Z → · viewed from −X · E${e.id} · x/L=${f(xi)}`);
  append('text',{x:12,y:355,fill:'#2166ac','font-size':13},`σmin=${f(at.minimum.value)} · y=${f(at.minimum.y)}, z=${f(at.minimum.z)}`);
  append('text',{x:12,y:377,fill:'#c61b32','font-size':13},`σmax=${f(at.maximum.value)} · y=${f(at.maximum.y)}, z=${f(at.maximum.z)}`);
  append('text',{x:12,y:402,fill:'#263847','font-size':13},`Scale ±${f(maximum)} ${this.unit('sigmaMax')} · red: tension`);
  const y=Number(this.pointY.value),z=Number(this.pointZ.value),inside=this.pointY.value.trim()&&this.pointZ.value.trim()&&P.contains(p,y,z);
  this.pointNote.textContent=inside?`σx(${f(y)}, ${f(z)}) = ${f(stress.point(xi,y,z).value)} ${this.unit('sigmaMax')}`:'Point outside the material or inside a hole: stress is undefined.';
  if(inside){const q=project(y,z);append('circle',{cx:q[0],cy:q[1],r:4,fill:'#111',stroke:'white'});}
  this.cutNote.textContent=`Along member length: σmin=${f(stress.minimum.value)} @ x/L=${f(stress.minimum.xi)}, y=${f(stress.minimum.y)}, z=${f(stress.minimum.z)}; σmax=${f(stress.maximum.value)} @ x/L=${f(stress.maximum.xi)}, y=${f(stress.maximum.y)}, z=${f(stress.maximum.z)}.`+(p.open?' Open profile: restrained warping and automatic shear-center offsets are not included.':'')+(this.analysis.type==='modal'?' Modal stresses at the saved amplitude × phase.':'');
  this.boundaryRows=p.loops.flatMap((loop,boundary)=>loop.map(v=>({boundary:boundary+1,...stress.point(xi,v.y,v.z)})));this.boundaryRows.push({boundary:0,...at.minimum},{boundary:0,...at.maximum});
 };
 proto.exportSection=function(kind){if(!this.boundaryRows){this.status.textContent='Select an available section view.';return;}let text,type;
  if(kind==='csv'){text='boundary,x/L,y,z,sigmaX\r\n'+this.boundaryRows.map(v=>[v.boundary,v.xi,v.y,v.z,v.value].join(',')).join('\r\n');type='text/csv;charset=utf-8';}
  else{const svg=this.cutSVG.cloneNode(true);svg.setAttribute('xmlns',M.svgNS);const meta=M.svgElement('metadata');meta.textContent=this.caption().join('\n');svg.prepend(meta);text=new XMLSerializer().serializeToString(svg);type='image/svg+xml;charset=utf-8';}
  const url=URL.createObjectURL(new Blob([text],{type})),a=document.createElement('a');a.href=url;a.download='MKE-F-3D-section.'+kind;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
 };
 const prepare=proto.prepare;proto.prepare=function(){
  if(this.analysis.type==='modal'){let maximum=0;const vector=S.vector(this.ds,'displacements',this.index,1);for(const e of this.model.elements){const stress=S.sectionStress(e,S.recover(this.ds,e,vector));if(stress)maximum=Math.max(maximum,Math.abs(stress.minimum.value),Math.abs(stress.maximum.value));}this.stressMaximum=maximum;}
  return prepare.call(this);
 };
 const draw=proto.drawChart;proto.drawChart=function(){
  if(!this.historyKind||this.historyKind.value!=='pointStress')return draw.call(this);
  const e=this.model.elements[Number(this.historyElement.value)-1],p=e&&S.sectionProfile(e),y=Number(this.pointY.value),z=Number(this.pointZ.value),xi=Number(this.position.value),a=this.analysis;
  if(!p||!P.contains(p,y,z)||!Number.isFinite(xi)||xi<0||xi>1||!a.displacements||!e.globalDOFs.every(id=>this.ds.systemRows.has(id-1))){this.chart.hidden=true;this.series=null;this.legend.textContent='σx history is unavailable: select a material point and an exported member.';return;}
  const key=JSON.stringify(['pointStress',e.id,xi,y,z]);let series=this.historyCache&&this.historyCache.key===key?this.historyCache.series:null;
  if(!series){const values=a.time.map((_,j)=>{const vector={};e.globalDOFs.forEach(id=>vector[id-1]=a.displacements[this.ds.systemRows.get(id-1)][j]);const c=S.recover(this.ds,e,vector).coeff,at=window.MKEFSections.evaluate;return at(c[0],xi)/p.area+at(c[4],xi)*z/p.Iy-at(c[5],xi)*y/p.Iz;});series={x:a.time,y:values,xLabel:'Analysis time',yLabel:`σx · E${e.id} · x/L=${xi} · Y=${y}, Z=${z}`};this.historyCache={key,series};}
  this.series=series;this.chart.hidden=false;M.charts.render(this.chart,series,a.time[this.index],i=>this.setFrame(i));
 };
 const caption=proto.caption;proto.caption=function(){return caption.call(this).concat(this.sectionNote?this.sectionNote.textContent:'');};
 proto.updateSections=function(){
  const count=Number(this.sectionCount.value),scale=Number(this.sectionScale.value);if(!Number.isInteger(count)||count<1||count>10||!(scale>0)||!Number.isFinite(scale)){this.sectionNote.textContent='Position count must be 1–10; scale must be positive.';return;}
  const selected=this.selection&&this.selection.kind==='element'?this.selection.index:null,xi=this.position?Number(this.position.value):.5,full=S.full(this.ds),deformed=this.layers.deformed&&full;
  const poses=P.positions(this.model.elements,selected,count,xi).flatMap(row=>{
   const e=this.model.elements[row.index],a=this.model.nodes[e.nodeIds[0]-1],b=this.model.nodes[e.nodeIds[1]-1],av=[a.x,a.y,a.z],bv=[b.x,b.y,b.z],entry=this.sectionCatalog.get(e.sectionId),u=deformed?S.local(e,this.vector):null;
   let center=av.map((x,i)=>x+(bv[i]-x)*row.xi),basis=e.localBasis;
   if(deformed&&!u)return [];
   if(u){const p=S.sample(e,u,row.xi),disp=p.displacement;center=center.map((x,i)=>x+this.currentScale*disp[i]);
    if(e.type===1436)basis=[0,1,2].map(k=>S.global(e.localBasis,P.rotate([0,1,2].map(i=>i===k?1:0),p.rotation.map(v=>v*this.currentScale))));
    else{const ends=[0,1].map(t=>{const d=S.sample(e,u,t).displacement;return av.map((x,i)=>x+(bv[i]-x)*t+d[i]*this.currentScale);});basis=P.transport(e.localBasis,ends[1].map((v,i)=>v-ends[0][i]));}
   }
   if(e.type===1133)basis=P.roll(basis,e.sectionAngle||0);
   return [{element:row.index,center,basis,loops:entry&&entry.properties.loops,selected:row.selected}];
  });
  if(this.scene)this.scene.setSectionFrame({poses,contours:this.sectionContours.checked,axes:this.layers.axes,scale});
  this.sectionNote.textContent=`Section positions ${poses.length}`+(this.model.elements.length*count>600?' · limited to 600':'')+' · outline scale ×'+scale;
 };
})(window.MKEFPost,window.MKEFSpatialSections);
