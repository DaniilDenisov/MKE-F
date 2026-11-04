(function(M){
  'use strict';
  const S=M.spatialResults,el=(tag,text,parent)=>{const e=document.createElement(tag);if(text!==undefined)e.textContent=text;if(parent)parent.append(e);return e;},fmt=x=>x===null||x===undefined?'not exported':Number(x.toPrecision(7)).toString();
  const labels=['ux','uy','uz','thetaX','thetaY','thetaZ'],forceNames=['N','Vy','Vz','Mx','My','Mz'];
  function field(parent,key,text,value,choices){const label=el('label',text,parent),input=el(choices?'select':'input',undefined,label);input.id='spatial-result-'+key;if(choices)choices.forEach(([v,t])=>{const o=el('option',t,input);o.value=v;});else {input.type='number';input.step='any';}input.value=value;return input;}
  function check(parent,key,text,value){const label=el('label',text,parent),input=el('input',undefined,label);input.id='spatial-result-'+key;input.type='checkbox';input.checked=value;return input;}
  class SpatialPost {
    constructor(ds){
      this.ds=ds;this.model=ds.raw.model;this.analysis=ds.raw.analysis;this.index=0;this.phase=1;this.selection=null;this.page=0;this.token=0;
      this.layers={original:true,deformed:true,nodes:true,axes:true,loads:true,supports:true,mpcs:true,releases:true};
      document.body.classList.add('spatial-result-mode');this.settings=el('div',undefined,document.querySelector('.settings'));this.settings.id='spatial-result-settings';
      this.host=document.getElementById('spatial-result-viewport');this.host.replaceChildren();
      this.panel=document.getElementById('spatial-result-panel');this.panel.replaceChildren();this.panel.hidden=false;this.status=el('p','',this.panel);this.status.id='spatial-result-status';
      this.details=el('pre','Select a node or member.',this.panel);this.details.id='spatial-result-details';this.tables=el('div',undefined,this.panel);
      this.makeControls();
      const availability=[];if(this.model.dofPerNode===6&&!this.analysis.nodalSupportReactions)availability.push('Physical support reactions were not exported.');if(this.model.dofPerNode===6&&!this.analysis.nodalMPCForces)availability.push('Physical MPC forces were not exported.');if(this.analysis.type==='transient'&&!this.analysis.nodalDisplacementAmplitudeSpectrum)availability.push('Spectra of mixed physical rotations were not exported.');
      el('p',availability.concat(this.model.warnings).join('\n'),this.panel);
      try{this.scene=new window.MKEFSpatial.SpatialScene(this.host,{select:s=>this.select(s),choose:hits=>this.choose(hits),contextLost:()=>this.stop&&this.stop(),change:s=>{this.viewNote.textContent=s.workPlane.plane+' · offset '+fmt(s.workPlane.offset)+' · grid spacing '+fmt(s.displayStep);}});this.scene.setModel({nodes:this.model.nodes,elements:this.model.elements.map(e=>({node1:e.nodeIds[0],node2:e.nodeIds[1]}))});this.scene.members.material.color.set('#8797a7');this.scene.setConditions(this.conditions());}
      catch(error){this.host.replaceChildren();el('p','WebGL 2 is unavailable. Tables and charts remain available. '+error.message,this.host);}
      const min=[Infinity,Infinity,Infinity],max=[-Infinity,-Infinity,-Infinity];this.model.nodes.forEach(n=>[n.x,n.y,n.z].forEach((v,i)=>{min[i]=Math.min(min[i],v);max[i]=Math.max(max[i],v);}));this.size=Math.max(Math.hypot(...max.map((v,i)=>v-min[i])),1e-9);
      this.prepare();this.table();this.visibility=()=>{if(document.hidden&&this.stop)this.stop();};document.addEventListener('visibilitychange',this.visibility);
    }
    makeControls(){
      const controls=this.settings;el('h2','Spatial results',controls);
      for(const [key,title]of [['original','Original geometry'],['deformed','Deformation'],['nodes','Nodes'],['axes','Local axes'],['loads','Loads'],['supports','Supports'],['mpcs','MPC'],['releases','Releases']]){const input=check(controls,'show-'+key,title,this.layers[key]);input.onchange=()=>{this.layers[key]=input.checked;this.update();};}
      this.scaleMode=field(controls,'scale-mode','Deformation scale','auto',[['auto','Automatic'],['manual','Manual']]);this.manualScale=field(controls,'scale','Factor',1);this.manualScale.min='0';this.scaleMode.onchange=this.manualScale.onchange=()=>this.update();
      this.samples=field(controls,'samples','Interpolation points',31);this.samples.min=5;this.samples.max=201;this.samples.step=1;this.samples.onchange=()=>this.update();
      const nav=el('div',undefined,controls);for(const [view,title]of [['fit','Fit all'],['origin','Reset to origin'],['iso','Isometric'],['plane','View work plane']]){const b=el('button',title,nav);b.dataset.resultView=view;b.onclick=()=>{if(!this.scene)return;if(view==='fit')this.scene.fitModel();else if(view==='origin'){this.scene.resetOrigin();this.plane.value='XY';this.offset.value=0;}else this.scene.setView(view);};}
      this.plane=field(controls,'plane','Work plane','XY',[['XY','XY'],['XZ','XZ'],['YZ','YZ']]);this.offset=field(controls,'offset','Offset',0);this.gridStep=field(controls,'step','Grid spacing',1);this.grid=check(controls,'grid','Grid',true);
      for(const input of [this.plane,this.offset,this.gridStep,this.grid])input.onchange=()=>{try{if(this.scene)this.scene.setWorkPlane({plane:this.plane.value,offset:Number(this.offset.value),step:Number(this.gridStep.value),visible:this.grid.checked});}catch(e){this.status.textContent=e.message;}};
      this.viewNote=el('p','',controls);this.scaleNote=el('p','',controls);this.scaleNote.className='control-note';
    }
    conditions(){const m=this.model;return {supports:m.supports.map(r=>({node:r.node,mask:['fixUx','fixUy','fixUz','fixThetaX','fixThetaY','fixThetaZ'].map(k=>r[k]||0)})),loads:m.nodalLoads.map(r=>({type:r[0],node:r[1],components:r.slice(2,8),frequency:r[8]})),elementLoads:m.elementLoads.map(r=>({type:r.type,elementId:r.elementId,coordinateSystem:r.coordinateSystem,components:r.type===20?r.ends.map(v=>v[0]):[0,1].flatMap(j=>r.ends.map(v=>v[j]))})),mpcs:m.mpcs.map(r=>({...r,masters:r.masters.map(a=>({node:a[0],dof:a[1],coefficient:a[2]}))})),releases:m.releases,bases:m.elements.map(e=>e.localBasis),analysisType:this.analysis.type};}
    prepare(){this.maximum=0;this.maxRotation=0;const vector=S.vector(this.ds,'displacements',this.index);if(vector)for(const e of this.model.elements){const u=S.local(e,vector);if(!u||u.some(v=>v===null))continue;const b=S.bounds(e,u);this.maximum=Math.max(this.maximum,b.translation);this.maxRotation=Math.max(this.maxRotation,b.rotation);}this.update();if(this.scene&&!this.initialFit){this.scene.fitModel();this.initialFit=true;}}
    scale(){const gain=Math.min(this.maximum>0?.1*this.size/this.maximum:Infinity,this.maxRotation>0?.25/this.maxRotation:Infinity);return this.scaleMode.value==='manual'?Number(this.manualScale.value):Number.isFinite(gain)?gain:1;}
    update(){
      const scale=this.scale(),requested=Number(this.samples.value);if(!(scale>0)||!Number.isFinite(scale)||!Number.isInteger(requested)||requested<5||requested>201){this.status.textContent='Scale must be positive; point count must be between 5 and 201.';return;}
      this.vector=S.vector(this.ds,'displacements',this.index,this.phase);this.acceleration=S.vector(this.ds,'accelerations',this.index,this.phase);this.currentScale=scale;
      const count=this.model.elements.length,samples=Math.max(5,Math.min(requested,Math.floor(1000000/(2*count))+1)),full=S.full(this.ds),positions=[],owners=[],parameters=[],nodes=[],nodeOwners=[];
      const put=(array,p)=>array.push(...p);
      const order=this.selection&&this.selection.kind==='element'?[this.model.elements[this.selection.index],...this.model.elements.filter((_,i)=>i!==this.selection.index)]:this.model.elements;
      for(const e of order){const u=S.local(e,this.vector);if(!u||u.some(v=>v===null))continue;const a=this.model.nodes[e.nodeIds[0]-1],b=this.model.nodes[e.nodeIds[1]-1],p=S.polynomials(e,u),at=t=>{const disp=S.global(e.localBasis,p.u.map(c=>S.value(c,t)));return [a.x+(b.x-a.x)*t+scale*disp[0],a.y+(b.y-a.y)*t+scale*disp[1],a.z+(b.z-a.z)*t+scale*disp[2]];};
        if(full&&this.layers.deformed){const n=e.type===1133?2:samples,ts=this.memberSamplePositions?this.memberSamplePositions(e,n):Array.from({length:n},(_,i)=>i/(n-1));let previous=at(0);for(let j=1;j<ts.length;j++){const next=at(ts[j]);put(positions,previous);put(positions,next);owners.push(e.id-1);parameters.push(ts[j-1],ts[j]);previous=next;}}
        // Section poses and screen-offset triads are supplied by spatial-section-view.
      }
      if(full&&this.layers.deformed&&this.layers.nodes)this.model.nodes.forEach((n,i)=>{const u=[0,1,2].map(j=>S.physical(this.ds,this.vector,i*this.model.dofPerNode+j));if(u.some(x=>x===null))return;nodes.push(n.x+u[0]*scale,n.y+u[1]*scale,n.z+u[2]*scale);nodeOwners.push(i);});
      this.frame={deformed:{positions,owners,parameters,visible:full&&this.layers.deformed},resultNodes:{positions:nodes,owners:nodeOwners,visible:this.layers.nodes&&full&&this.layers.deformed}};
      if(this.decorateFrame)this.decorateFrame(this.frame);
      if(this.scene){this.scene.members.visible=this.layers.original;this.scene.nodes.visible=this.layers.nodes&&this.layers.original;this.scene.setConditionLayers({loads:this.layers.loads,supports:this.layers.supports,mpcs:this.layers.mpcs,releases:this.layers.releases,axes:false});this.scene.setResultFrame(this.frame);}
      this.status.textContent=this.analysis.type+' · '+(full?'All displacements are available':'Partial export: available charts only, no full-model animation')+' · '+samples+' points per frame · labels ≤ 200';
      this.scaleNote.textContent='Display scale ×'+fmt(scale)+'. Numerical values are not scaled.';this.showDetails();if(this.afterUpdate)this.afterUpdate();
    }
    select(s){this.selection=s;if(this.scene)this.scene.setSelection(s);this.update();this.table();}
    choose(hits){if(!this.dialog){this.dialog=el('dialog',undefined,document.body);}this.dialog.replaceChildren();el('h3','Objects at pick location',this.dialog);hits.slice(0,200).forEach(hit=>{const b=el('button',hit.kind+' '+(hit.index+1),this.dialog);b.onclick=()=>{this.dialog.close();this.select(hit);};});const close=el('button','Close',this.dialog);close.onclick=()=>this.dialog.close();this.dialog.showModal();}
    showDetails(){const s=this.selection;if(!s){this.details.textContent='Select a node or member.';return;}const lines=[];
      if(s.kind==='node'){const n=this.model.nodes[s.index];lines.push('Node '+n.id+' · XYZ '+[n.x,n.y,n.z].join(', '));labels.slice(0,this.model.dofPerNode).forEach((name,j)=>lines.push(name+' = '+(this.ds.physicalRows[s.index*this.model.dofPerNode+j].length?fmt(S.physical(this.ds,this.vector,s.index*this.model.dofPerNode+j)):'inactive')));const u=[0,1,2].map(j=>S.physical(this.ds,this.vector,s.index*this.model.dofPerNode+j));lines.push('|u| = '+(u.some(x=>x===null)?'not exported':fmt(Math.hypot(...u))));
        if(this.analysis.type==='transient')for(const field of ['velocities','accelerations']){const vector=S.vector(this.ds,field,this.index);lines.push(field+': '+labels.slice(0,this.model.dofPerNode).map((name,j)=>name+'='+fmt(S.physical(this.ds,vector,s.index*this.model.dofPerNode+j))).join(', '));}}
      else if(s.kind==='element'){const e=this.model.elements[s.index],r=S.recover(this.ds,e,this.vector,this.acceleration);lines.push('Member '+e.id+' · L='+fmt(e.length));if(r){lines.push('Local end displacements: '+r.u.map(fmt).join(', '));lines.push('End actions: '+r.end.map(fmt).join(', '));if(r.inertia)lines.push('Inertial: '+r.inertia.map(fmt).join(', '),'Total: '+r.total.map(fmt).join(', '));}else lines.push('Required displacements were not exported.');}
      else{lines.push(s.kind+' '+(s.index+1));const record=this.model[s.kind]&&this.model[s.kind][s.index];if(record)lines.push(JSON.stringify(record,null,2));}this.details.textContent=lines.join('\n');
    }
    table(){
      this.tables.replaceChildren();const kind=this.tableKind||'node',choices=[['node','Nodes'],['element','Members'],['supports','Supports'],['loads','Nodal loads'],['elementLoads','Distributed loads'],['mpcs','MPC'],['releases','Releases']];
      const items=kind==='node'?this.model.nodes:kind==='element'?this.model.elements:kind==='loads'?this.model.nodalLoads:this.model[kind];
      const select=field(this.tables,'table','Table',kind,choices);select.onchange=()=>{this.tableKind=select.value;this.page=0;this.table();};const prev=el('button','Previous',this.tables),next=el('button','Next',this.tables);prev.disabled=this.page===0;next.disabled=(this.page+1)*200>=items.length;prev.onclick=()=>{this.page--;this.table();};next.onclick=()=>{this.page++;this.table();};el('span',' Page '+(this.page+1)+' / '+Math.max(1,Math.ceil(items.length/200)),this.tables);
      const body=el('tbody',undefined,el('table',undefined,this.tables));body.id='spatial-result-rows';items.slice(this.page*200,(this.page+1)*200).forEach((r,j)=>{const index=this.page*200+j,tr=el('tr',undefined,body);tr.tabIndex=0;tr.classList.toggle('selected-row',!!this.selection&&this.selection.kind===kind&&this.selection.index===index);el('td',index+1,tr);el('td',kind==='node'?[r.x,r.y,r.z].map(fmt).join(', '):kind==='element'?r.nodeIds.join(' → '):JSON.stringify(r),tr);tr.onclick=()=>this.select({kind,index});tr.onkeydown=e=>{if(e.key==='Enter')tr.click();};});
    }
    dispose(){this.disposed=true;this.token++;if(this.stop)this.stop();if(this.scene)this.scene.dispose();document.removeEventListener('visibilitychange',this.visibility);this.settings.remove();this.host.replaceChildren();this.panel.replaceChildren();this.panel.hidden=true;if(this.dialog)this.dialog.remove();document.body.classList.remove('spatial-result-mode');}
  }
  M.SpatialPost=SpatialPost;M.spatialUI={el,field,check,fmt,labels,forceNames};
}(window.MKEFPost));
