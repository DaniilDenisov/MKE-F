(function(M){
 'use strict';
 const P=M.planarPresentation={},G=M.geometry;
 const rgb=hex=>[1,3,5].map(i=>parseInt(hex.slice(i,i+2),16)/255);
 const group=(color,kind)=>({positions:[],owners:[],colors:[],color,kind:kind||'lines'});
 function segment(out,a,b,owner,color){out.positions.push(...a,...b);out.owners.push(owner);if(color)out.colors.push(...rgb(color),...rgb(color));}
 function line(out,points,owner,color){for(let i=1;i<points.length;i++)segment(out,points[i-1],points[i],owner,color);}
 const xyz=p=>[p.x,p.y,0];
 function arrow(out,p,v,length,owner){const mag=Math.hypot(v[0],v[1]);if(!mag)return;const x=v[0]/mag,y=v[1]/mag,h=length*.22;
  segment(out,[p[0]-x*length,p[1]-y*length,0],p,owner);
  line(out,[[p[0]-x*h-y*h*.45,p[1]-y*h+x*h*.45,0],p,[p[0]-x*h+y*h*.45,p[1]-y*h-x*h*.45,0]],owner);
 }
 function moment(out,p,value,radius,owner){if(!value)return;const sign=Math.sign(value),points=Array.from({length:33},(_,i)=>{const a=sign*(Math.PI/6+i/32*Math.PI*1.65);return [p[0]+Math.cos(a)*radius,p[1]+Math.sin(a)*radius,0];});line(out,points,owner);const tip=points[32],a=sign*(Math.PI/6+Math.PI*1.65);arrow(out,tip,[-Math.sin(a)*sign,Math.cos(a)*sign],radius*.35,owner);}
 P.size=ds=>{const b=G.modelBounds(ds.raw.model);return Math.max(Math.hypot(b.maxX-b.minX,b.maxY-b.minY),1e-9);};
 P.conditions=ds=>{const m=ds.raw.model,physical=d=>d===3?6:d,node=id=>ds.nodeIndexById.get(id)+1,element=id=>m.elements.findIndex(e=>e.id===id)+1;
  return {supports:m.supports.map(s=>{const a=window.MKEFSupports.mask(s,m.dofPerNode);return {node:node(s.nodeId),mask:[+a[0],+a[1],0,0,0,+a[2]]};}),loads:[],elementLoads:[],
   mpcs:(m.mpcs||[]).map(r=>({...r,depNode:node(r.depNode),depDOF:physical(r.depDOF),masters:r.masters.map(a=>({...a,node:node(a.node),dof:physical(a.dof)}))})),releases:(m.releases||[]).map(r=>({...r,elementId:element(r.elementId)})),
   bases:m.elements.map(e=>{const b=G.elementBasis(ds.nodesById.get(e.nodeIds[0]),ds.nodesById.get(e.nodeIds[1]));return [[b.c,b.s,0],[-b.s,b.c,0],[0,0,1]];}),analysisType:ds.raw.analysis.type,
   inactiveSupports:m.supports.flatMap((s,i)=>m.dofMap[ds.nodeIndexById.get(s.nodeId)][2]===0?[{index:i,component:5}]:[])};
 };
 P.geometry=function(ds,vector,scale,count,visibility,colors){
  const original=group('#8797a7'),deformed=group('#087f9c'),nodes=group('#087f9c','points');
  original.visible=visibility.showOriginal!==false;deformed.visible=!!vector&&visibility.showDeformed!==false;nodes.visible=deformed.visible&&visibility.showNodes!==false;
  ds.raw.model.elements.forEach((e,index)=>{
   const a=ds.nodesById.get(e.nodeIds[0]),b=ds.nodesById.get(e.nodeIds[1]),n=e.type===112?2:count,points=vector?G.sampleElement(ds,e,vector,scale,n):[],style=colors&&colors.get(e.id);
   function color(t){if(!style)return null;if(typeof style==='string')return style;const right=style.findIndex(s=>s.offset>=t);if(right<=0)return style[0].color;const u=style[right-1],v=style[right],f=(t-u.offset)/(v.offset-u.offset),x=rgb(u.color),y=rgb(v.color);return '#'+x.map((c,i)=>Math.round(255*(c+(y[i]-c)*f)).toString(16).padStart(2,'0')).join('');}
   const steps=style&&typeof style!=='string'?Math.max(count,2):2;
   const normal=ds.normalStress?.active;
   for(let j=1;j<steps;j++)segment(original,[a.x+(b.x-a.x)*(j-1)/(steps-1),a.y+(b.y-a.y)*(j-1)/(steps-1),0],[a.x+(b.x-a.x)*j/(steps-1),a.y+(b.y-a.y)*j/(steps-1),0],index,normal?'#8797a7':color((j-.5)/(steps-1))||'#8797a7');
   if(normal&&vector){
    const r=ds.normalStress.results.get(e.id),ts=[...new Set(Array.from({length:Math.max(2,n)},(_,i)=>i/(Math.max(2,n)-1)).concat(r?.status==='ok'?r.positions:[]))].sort((a,b)=>a-b),at=G.elementInterpolator(ds,e,vector,scale);
    for(let j=1;j<ts.length;j++){deformed.positions.push(...xyz(at(ts[j-1])),...xyz(at(ts[j])));deformed.owners.push(index);for(const t of [ts[j-1],ts[j]]){const v=M.normalStressResults.at(r,t);deformed.colors.push(...(v===null?rgb('#8797a7'):M.normalStressResults.rgb(v,ds.normalStress.maximum)));}}
   }else for(let j=1;j<points.length;j++)segment(deformed,xyz(points[j-1]),xyz(points[j]),index,color((j-.5)/(points.length-1))||'#087f9c');
  });
  if(vector)ds.raw.model.nodes.forEach((node,i)=>{const u=G.nodeDisplacement(ds,node.id,vector);nodes.positions.push(node.x+scale*u.x,node.y+scale*u.y,0);nodes.owners.push(i);});
  return {original,deformed,resultNodes:nodes};
 };
 P.result=function(ds,mode,count){
  const colors=new Map(),diagrams=group('#3274ac'),diagramFill=group('#9fc4e2','triangles');diagramFill.opacity=.25;
  let legend={kind:'none',summary:''},maximum=0;
  if(mode==='displacementMagnitude')legend=M.staticResults.applyDisplacementMagnitude({setElementColor:(id,c)=>colors.set(id,c),setElementGradient:(id,a,b,s)=>colors.set(id,s)},ds,count);
  else if(['axialForce','axialStress'].includes(mode)){
   const values=ds.raw.model.elements.map(e=>ds.resultsByElementId.get(e.id)[mode]);maximum=Math.max(0,...values.map(Math.abs));ds.raw.model.elements.forEach((e,i)=>colors.set(e.id,!values[i]?'#66788a':values[i]>0?'#bd493e':'#3971b9'));
   legend={kind:'text',text:mode+' · max |value| '+M.staticResults.withUnit(maximum,ds.raw.metadata.units[mode==='axialForce'?'force':'stress'])+' · red tension, blue compression'};
  }else if(mode!=='none'){
   const values=ds.raw.model.elements.map(e=>M.staticResults.diagramSamples(ds,e,mode,count));for(const rows of values)for(const r of rows)maximum=Math.max(maximum,Math.abs(r.value));const scale=maximum?P.size(ds)*.14/maximum:0;
   ds.raw.model.elements.forEach((e,index)=>{const a=ds.nodesById.get(e.nodeIds[0]),b=ds.nodesById.get(e.nodeIds[1]),basis=G.elementBasis(a,b),rows=values[index];
    const base=r=>[a.x+(b.x-a.x)*r.xi,a.y+(b.y-a.y)*r.xi,0],tip=r=>{const p=base(r);return [p[0]-basis.s*r.value*scale,p[1]+basis.c*r.value*scale,0];};line(diagrams,rows.map(tip),index);
    for(let i=1;i<rows.length;i++){const a0=base(rows[i-1]),a1=tip(rows[i-1]),b0=base(rows[i]),b1=tip(rows[i]);diagramFill.positions.push(...a0,...a1,...b1,...a0,...b1,...b0);diagramFill.owners.push(index,index);}
   });legend={kind:'text',text:mode+' diagram · max |value| '+M.staticResults.withUnit(maximum,ds.raw.metadata.units[M.sectionResults.isStress(mode)?'stress':mode==='M'?'moment':'force'])+' · display scale '+M.staticResults.format(scale)};
  }
  if(legend.kind==='text')legend.summary=legend.text;
  return {colors,diagrams,diagramFill,legend,maximum};
 };
 P.forces=function(ds,vector,color,restrained){
  const out=group(color);out.ownerKind='node';if(!vector)return out;const rows=ds.raw.model.nodes.map((node,i)=>{const ids=ds.raw.model.dofMap[i],support=ds.raw.model.supports.find(s=>s.nodeId===node.id),mask=restrained?(support?window.MKEFSupports.mask(support,ds.raw.model.dofPerNode):[]):[true,true,true];return {node,values:ids.map((id,j)=>mask[j]&&id?vector[id-1]||0:0)};});
  let maxF=0,maxM=0;rows.forEach(r=>{maxF=Math.max(maxF,Math.abs(r.values[0]),Math.abs(r.values[1]));maxM=Math.max(maxM,Math.abs(r.values[2]||0));});const size=P.size(ds)*.13;
  rows.forEach((r,i)=>{const [x,y,m]=r.values,p=xyz(r.node),mag=Math.hypot(x,y);if(mag>Math.max(1,maxF)*1e-10)arrow(out,p,[x,y],size*(.35+.65*mag/maxF),i);if(Math.abs(m||0)>Math.max(1,maxM)*1e-10)moment(out,p,m,size*(.32+.3*Math.abs(m)/maxM),i);});return out;
 };
 P.distributed=function(ds,scale,original,deformed){const out=group('#d44a3a');out.labels=[];
  const indices=new Map(ds.raw.model.elements.map((e,i)=>[e.id,i]));
  for(const load of ds.raw.model.elementLoads||[]){const e=ds.elementsById.get(load.elementId),a=ds.nodesById.get(e.nodeIds[0]),b=ds.nodesById.get(e.nodeIds[1]),R=G.elementBasis(a,b),ends=load.type===21?[[load.qx1,load.qy1],[load.qx2,load.qy2]]:[[load.qx,load.qy],[load.qx,load.qy]],v=ends.map(q=>load.coordinateSystem===1?[R.c*q[0]-R.s*q[1],R.s*q[0]+R.c*q[1]]:q),max=Math.max(...v.map(q=>Math.hypot(...q)));if(!max)continue;
   for(const factor of [original?0:null,deformed&&(scale!==0||!original)?scale:null]){if(factor===null)continue;const points=G.sampleElement(ds,e,ds.raw.analysis.displacements,factor,9),tails=[],size=Math.min(P.size(ds)*.08,R.length*.3);points.forEach((point,i)=>{const t=i/(points.length-1),q=v[0].map((x,j)=>x*(1-t)+v[1][j]*t),p=xyz(point);arrow(out,p,q,size*Math.hypot(...q)/max,indices.get(e.id));tails.push([p[0]-q[0]/max*size,p[1]-q[1]/max*size,0]);});line(out,tails,indices.get(e.id));if(out.labels.length<200){const units=ds.raw.metadata.units,pair=q=>'('+q.map(M.staticResults.format).join(', ')+')';out.labels.push({position:tails[4],text:(load.coordinateSystem===1?'Local':'Global')+(load.type===21?' q₁='+pair(ends[0])+' → q₂='+pair(ends[1]):' q='+pair(ends[0]))+' '+(units.force||'force')+'/'+(units.length||'length'),color:'#b63a3f'});}}
  }return out;
 };
}(window.MKEFPost));
