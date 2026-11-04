(function(M){
 'use strict';
 const P=M.planarPresentation;
 class Renderer {
  constructor(host,onSelection){this.host=host;this.onSelection=onSelection;this.webgl=true;this.visibility={};this.layers={};this.extra={};this.colors=new Map();
   for(const key of ['loads','reactions','diagrams','mpcs','spc-forces','mpc-forces']){let display='';this.layers[key]={key,clear:()=>{delete this.extra[this.layerKey(key)];this.request();},style:{get display(){return display;},set display(value){display=value;this.changed();},changed:()=>this.request()}};}
  }
  layerKey(key){return {'spc-forces':'spcForces','mpc-forces':'mpcForces'}[key]||key;}
  mount(ds){this.dispose();this.dataset=ds;this.elementIndexById=new Map(ds.raw.model.elements.map((e,i)=>[e.id,i]));this.vector=null;this.selection=null;this.colors=new Map();this.extra={};this.visibility={};for(const layer of Object.values(this.layers))layer.style.display='';this.size=P.size(ds);this.samples=31;this.scale=1;
   this.host.replaceChildren();
   try{this.scene=new window.MKEFSpatial.SpatialScene(this.host,{select:s=>{if(s&&!['node','element'].includes(s.kind))return;this.setSelection(s?{kind:s.kind,id:(s.kind==='node'?this.dataset.raw.model.nodes:this.dataset.raw.model.elements)[s.index].id}:null);},contextLost:()=>{if(this.onContextLost)this.onContextLost();}},{dimension:2});this.scene.setModel({nodes:ds.raw.model.nodes.map(n=>({...n,z:0})),elements:ds.raw.model.elements.map(e=>({id:e.id,node1:ds.nodeIndexById.get(e.nodeIds[0])+1,node2:ds.nodeIndexById.get(e.nodeIds[1])+1}))});this.scene.setConditions(P.conditions(ds));}
   catch(e){this.scene=null;const note=document.createElement('p');note.textContent='WebGL 2 unavailable. Tables and charts remain available. '+e.message;this.host.append(note);}
   this.makeObjectList();this.updateSectionAxes();this.request();
  }
  makeObjectList(){let box=document.getElementById('planar-result-objects');if(!box){box=document.createElement('section');box.id='planar-result-objects';document.getElementById('result-sidebar').prepend(box);}this.list=box;let page=0,kind='node';
   const draw=()=>{box.replaceChildren();const title=document.createElement('h2');title.textContent='Scene objects';box.append(title);const select=document.createElement('select');select.setAttribute('aria-label','Scene object type');for(const value of ['node','element']){const o=document.createElement('option');o.value=value;o.textContent=value==='node'?'Nodes':'Elements';select.append(o);}select.value=kind;select.onchange=()=>{kind=select.value;page=0;draw();};box.append(select);const items=kind==='node'?this.dataset.raw.model.nodes:this.dataset.raw.model.elements;
    for(const [label,delta] of [['Previous',-1],['Next',1]]){const b=document.createElement('button');b.textContent=label;b.type='button';b.disabled=delta<0?page===0:(page+1)*200>=items.length;b.onclick=()=>{page+=delta;draw();};box.append(b);}
    const table=document.createElement('table');table.id='planar-result-rows';items.slice(page*200,(page+1)*200).forEach(item=>{const tr=document.createElement('tr');tr.tabIndex=0;for(const text of [item.id,kind==='node'?[item.x,item.y].join(', '):item.nodeIds.join(' → ')]){const td=document.createElement('td');td.textContent=text;tr.append(td);}tr.onclick=()=>this.setSelection({kind,id:item.id});tr.onkeydown=e=>{if(e.key==='Enter')tr.click();};table.append(tr);});box.append(table);};draw();
  }
  dispose(){if(this.frameRequest)cancelAnimationFrame(this.frameRequest);this.frameRequest=null;if(this.scene)this.scene.dispose();this.scene=null;if(this.list)this.list.remove();this.list=null;this.dataset=null;this.frame=null;this.vector=null;this.extra={};this.colors.clear();}
  request(){if(this.frameRequest||!this.dataset)return;this.frameRequest=requestAnimationFrame(()=>{this.frameRequest=null;this.flush();});}
  flush(){if(!this.dataset)return;this.frame={...P.geometry(this.dataset,this.vector,this.scale,this.samples,this.visibility,this.colors),...this.extra};
   for(const [key,layer] of Object.entries(this.layers)){const data=this.frame[this.layerKey(key)];if(data)data.visible=layer.style.display!=='none';}
   if(this.frame.diagramFill)this.frame.diagramFill.visible=this.layers.diagrams.style.display!=='none';
   for(const key of ['loads','distributed'])if(this.frame[key])this.frame[key].visible=this.visibility.showLoads!==false;
   if(this.frame.reactions)this.frame.reactions.visible=this.visibility.showReactions!==false;
   this.frame.labels=this.visibility.showLoads!==false?this.extra.distributed?.labels||[]:[];
   if(this.scene){this.scene.members.visible=false;this.scene.nodes.visible=this.visibility.showOriginal!==false&&this.visibility.showNodes!==false;this.scene.labelVisibility={node:this.visibility.showNodeLabels!==false,element:this.visibility.showElementLabels!==false};this.scene.setConditionLayers({loads:false,elementLoads:false,supports:this.visibility.showSupports!==false,mpcs:this.layers.mpcs.style.display!=='none',releases:true,axes:false});this.scene.setResultFrame(this.frame);}
  }
  updateDeformation(vector,scale,samples){this.vector=vector;this.scale=scale;this.samples=samples;this.request();}
  setVisibility(settings){this.visibility={...settings};this.request();}
  setSelection(s){this.selection=s;if(this.scene)this.scene.setSelection(s?{kind:s.kind,index:(s.kind==='node'?this.dataset.nodeIndexById:this.elementIndexById).get(s.id)}:null);if(this.onSelection)this.onSelection(s);}
  reapplySelection(){if(this.scene)this.scene.setSelection(this.selection?{kind:this.selection.kind,index:(this.selection.kind==='node'?this.dataset.nodeIndexById:this.elementIndexById).get(this.selection.id)}:null);}
  clearLayer(key){this.layers[key].clear();}
  resetElementColors(){this.colors.clear();this.request();}
  setElementColor(id,color){this.colors.set(id,color);this.request();}
  setElementGradient(id,a,b,stops){this.colors.set(id,stops);this.request();}
  applyResult(ds,mode,count){const r=P.result(ds,mode,count);this.colors=r.colors;this.extra.diagrams=r.diagrams;this.extra.diagramFill=r.diagramFill;this.layers.diagrams.style.display=mode==='none'?'none':'';this.request();return r.legend;}
  setForces(key,vector,color,restrained){this.extra[this.layerKey(key)]=P.forces(this.dataset,vector,color,restrained);this.request();}
  setDistributedLoads(scale,original,deformed){this.extra.distributed=P.distributed(this.dataset,scale,original,deformed);this.request();}
  updateSectionAxes(){if(!this.scene||!this.dataset)return;const m=this.dataset.raw.model,catalog=window.MKEFSpatialSections.catalog(m),bases=P.conditions(this.dataset).bases;
   this.scene.setSectionFrame({contours:false,axes:document.getElementById('show-section-axes')?.checked!==false,count:1,scale:1,poses:m.elements.slice(0,600).map((e,i)=>{const a=this.dataset.nodesById.get(e.nodeIds[0]),b=this.dataset.nodesById.get(e.nodeIds[1]),section=catalog.get(e.sectionId);return {element:i,center:[(a.x+b.x)/2,(a.y+b.y)/2,0],basis:window.MKEFSpatialSections.roll(bases[i],e.sectionAngle||0),loops:section&&section.properties.loops};})});}
  fit(){this.flush();if(this.scene)this.scene.fitModel();}
  sceneNodeCount(){return this.scene?this.scene.scene.children.length:0;}
 }
 M.WebGLRenderer=Renderer;
}(window.MKEFPost));
