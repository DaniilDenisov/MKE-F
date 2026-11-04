(function(M){
 'use strict';
 const E=M.exporting;
 // Export already prepared buffers. This creates a standalone artifact, not a second viewport.
 E.sceneSource=function(scene,raster){
  scene.capture();const width=scene.width,height=scene.pixels,svg=M.svgElement('svg',{xmlns:M.svgNS,width,height,viewBox:`0 0 ${width} ${height}`});
  svg.getBoundingClientRect=()=>({width,height});svg.append(M.svgElement('rect',{width,height,fill:'#f4f7fb'}));
  if(raster){const canvas=document.createElement('canvas');canvas.width=scene.canvas.width;canvas.height=scene.canvas.height;const ctx=canvas.getContext('2d');ctx.drawImage(scene.canvas,0,0);ctx.drawImage(scene.overlay,0,0,canvas.width,canvas.height);svg.append(M.svgElement('image',{width,height,href:canvas.toDataURL('image/png')}));return svg;}
  const names=new Map(Object.entries(scene.results.objects).map(([key,obj])=>[obj,key]));names.set(scene.nodes,'original-nodes');names.set(scene.members,'original-geometry');
  const objects=scene.scene.children.filter(o=>o.visible&&o.geometry&&o.material&&o.material.visible!==false).sort((a,b)=>(a.isMesh?-1:0)-(b.isMesh?-1:0)||(a.renderOrder||0)-(b.renderOrder||0));
  for(const obj of objects){const attr=obj.geometry.attributes.position;if(!attr)continue;const group=M.svgElement('g',{'data-layer':names.get(obj)||'conditions',opacity:obj.material.opacity??1}),points=[];
   for(let i=0;i<attr.count;i++)points.push(scene.projectPoint(attr.getX(i),attr.getY(i),attr.getZ(i)));
   const colors=obj.geometry.attributes.color,color=i=>colors&&obj.material.vertexColors?obj.material.color.clone().setRGB(colors.getX(i),colors.getY(i),colors.getZ(i)).getStyle():obj.material.color.getStyle();
   if(obj.isPoints){points.forEach((p,i)=>{if(Math.abs(p.z)>1)return;group.append(M.svgElement('circle',{cx:p.x,cy:p.y,r:obj.material.size/2,fill:color(i)}));});}
   else {const step=obj.isMesh?3:2;for(let i=0;i+step<=points.length;i+=step){const row=points.slice(i,i+step);if(row.every(p=>Math.abs(p.z)>1))continue;const d=row.map((p,j)=>(j?'L':'M')+p.x+' '+p.y).join(' ')+(obj.isMesh?' Z':'');group.append(M.svgElement('path',{d,fill:obj.isMesh?color(i):'none',stroke:obj.isMesh?'none':color(i),'stroke-width':obj===scene.highlight?3:1.5}));}}
   svg.append(group);
  }
  for(const label of scene.labelFrame||[]){const text=M.svgElement('text',{x:label.x,y:label.y,fill:label.color,style:'font:'+label.font});text.textContent=label.text;svg.append(text);}
  return svg;
 };
 E.downloadSceneSvg=function(scene,context,chart,legend){return E.downloadSvg(E.sceneSource(scene,false),context,chart,legend);};
 E.downloadScenePng=function(scene,context,resolution,chart,legend){return E.downloadPng(E.sceneSource(scene,true),context,resolution,chart,legend);};
}(window.MKEFPost));
