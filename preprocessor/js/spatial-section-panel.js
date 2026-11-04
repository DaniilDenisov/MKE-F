(function(M,P,S){
 'use strict';
 const proto=M.SpatialPanel.prototype,el=(t,text,parent)=>{const e=document.createElement(t);if(text!==undefined)e.textContent=text;if(parent)parent.append(e);return e;};
 function input(parent,label,value,id){const l=el('label',label,parent),e=el('input',undefined,l);e.id=id;e.type='number';e.step='any';e.required=true;e.value=value;return e;}
 const steel={youngsModulus:210e9,density:7850,shearModulus:210e9/2.6};
 const render=proto.render;
 proto.render=function(m,s){if(this.model!==m){this.sectionEditor=null;this.assignSequential=false;this.sectionPage=0;this.materialDraft=null;this.materialValues={...steel};this.memberDefaults={112:{youngsModulus:200e9,density:7800,area:.0001},113:{youngsModulus:200e9,density:7800,area:.0001,momentOfInertia:3.33e-9},1133:{...steel},1436:{...steel}};this.materialValues={...steel,...this.memberDefaults[m.elementType]};if(this.catalogDialog)this.catalogDialog.close();}render.call(this,m,s);this.renderCatalog();this.renderSectionView();};
 proto.openCatalog=function(){this.catalogReturnFocus=document.activeElement;this.renderCatalog();this.catalogDialog.showModal();};
 proto.closeCatalog=function(){this.catalogDialog.close();if(this.catalogReturnFocus&&this.catalogReturnFocus.isConnected)this.catalogReturnFocus.focus();};
 proto.ensureCatalog=function(){
  if(this.catalogDialog)return;
  const dialog=el('dialog',undefined,document.body);dialog.id='spatial-catalog-tool';dialog.setAttribute('aria-labelledby','spatial-catalog-title');this.catalogDialog=dialog;
  const header=el('header',undefined,dialog);el('h2','Materials and Sections',header).id='spatial-catalog-title';
  for(const [key,title,back]of [['undo','Undo',true],['redo','Redo',false]]){const b=el('button',title,header);b.id='spatial-catalog-'+key;b.onclick=()=>this.travel(back);}
  const close=el('button','Close ×',header);close.id='spatial-catalog-close';close.type='button';close.onclick=()=>this.closeCatalog();
  dialog.addEventListener('cancel',event=>{event.preventDefault();this.closeCatalog();});
  dialog.addEventListener('keydown',event=>{if((event.ctrlKey||event.metaKey)&&!event.isComposing&&!event.target.closest('input,textarea,select,[contenteditable="true"]')){const key=event.key.toLowerCase();if(key==='z'||key==='y'){event.preventDefault();this.travel(key==='z'&&!event.shiftKey);}}});
  const tabs=el('nav',undefined,dialog);tabs.setAttribute('aria-label','Material and section settings');
  for(const [key,title]of [['sections','Sections'],['material','Material defaults for new members']]){const b=el('button',title,tabs);b.id='spatial-catalog-tab-'+key;b.onclick=()=>{this.catalogTab=key;this.showCatalogTab();};}
  this.catalogMessage=el('p','',dialog);this.catalogMessage.id='spatial-catalog-message';this.catalogMessage.setAttribute('role','status');
  this.catalogBody=el('div',undefined,dialog);this.catalogBody.className='spatial-catalog-body';
  document.getElementById('spatial-catalog-toggle').onclick=()=>this.openCatalog();
  document.getElementById('assign-section').onclick=()=>{if(this.assignSequential){this.setTool('select');return;}this.catalogTab='sections';this.openCatalog();};
 };
 proto.showCatalogTab=function(){const key=this.catalogTab||'sections';for(const tab of ['sections','material']){document.getElementById('spatial-catalog-'+tab).hidden=tab!==key;document.getElementById('spatial-catalog-tab-'+tab).setAttribute('aria-pressed',String(tab===key));}};
 proto.renderMaterial=function(){
  const box=el('section',undefined,this.catalogBody);box.id='spatial-catalog-material';el('h3','Material defaults for new members',box);
  el('p','Default steel: E = 210 GPa, ρ = 7850 kg/m³, G = E / [2(1 + 0.3)]. Enter consistent SI units: Pa and kg/m³.',box);
  const form=el('form',undefined,box),fields={},values=this.materialDraft||this.materialValues;
  for(const [key,label]of [['youngsModulus','E, Pa'],['density','ρ, kg/m³'],['shearModulus','G, Pa (frames)']]){fields[key]=input(form,label,values[key],'spatial-material-'+key);fields[key].min='0';if(key==='shearModulus'&&this.model.elementType!==1436)fields[key].parentElement.hidden=true;}
  const read=()=>Object.fromEntries(Object.entries(fields).map(([key,e])=>[key,e.value]));form.oninput=()=>{this.materialDraft=read();};
  const apply=el('button','Apply to new members',form);apply.id='spatial-material-apply';apply.type='submit';
  const reset=el('button','Default steel',form);reset.type='button';reset.id='spatial-material-steel';reset.onclick=()=>{for(const key in fields)fields[key].value=steel[key];this.materialDraft=read();};
  el('p','Existing members retain their properties. These defaults apply to new members and can be overridden in the member form.',box);
  form.onsubmit=event=>{event.preventDefault();const next=Object.fromEntries(Object.entries(fields).map(([key,e])=>[key,e.value.trim()?Number(e.value):NaN]));if(Object.values(next).some(v=>!Number.isFinite(v)||v<=0)){this.catalogMessage.textContent='E, ρ and G must be positive finite numbers.';return;}
   this.materialValues=next;this.materialDraft=null;for(const code of [112,113,1133,1436])Object.assign(this.memberDefaults[code],next);
   if(this.tool==='member'){this.memberDraft={...this.readMemberDraft(),...next};this.form.replaceChildren();this.renderMemberForm();}
   this.catalogMessage.textContent='Material defaults apply to subsequently created members.';
  };
 };
 proto.renderSectionView=function(){
  const settings=this.sectionView||(this.sectionView={contours:false,axes:true,count:1,scale:1,xi:.5}),details=el('details',undefined,this.root);details.id='spatial-section-view';details.open=!!this.sectionViewOpen;details.ontoggle=()=>{if(details.isConnected)this.sectionViewOpen=details.open;};el('summary','Section display',details);const box=el('div',undefined,details);
  for(const [key,label]of [['contours','Section outlines'],['axes','Local axes']]){const check=el('input',undefined,el('label',label,box));check.type='checkbox';check.checked=settings[key];check.id='spatial-sections-'+key;check.onchange=()=>{settings[key]=check.checked;this.updateSections();};}
  for(const [key,label]of [['count','Positions per member (1–10)'],['scale','Outline scale']]){const e=input(box,label,settings[key],'spatial-sections-'+key);e.onchange=()=>{const v=Number(e.value);if(!(v>0)||!Number.isFinite(v)||key==='count'&&(!Number.isInteger(v)||v>10)){this.status('Position count must be 1–10; scale must be positive.');return;}settings[key]=v;this.updateSections();};}
  const xi=input(box,'Selected member section · x/L',settings.xi,'spatial-sections-position');xi.min=0;xi.max=1;xi.onchange=()=>{const v=Number(xi.value);if(!xi.value.trim()||!Number.isFinite(v)||v<0||v>1){this.status('x/L must be between 0 and 1.');return;}settings.xi=v;this.updateSections();};
  this.sectionViewNote=el('p','',box);this.updateSections();
 };
 proto.updateSections=function(){
  const settings=this.sectionView,m=this.model,catalog=P.catalog(m),selected=this.selection&&this.selection.kind==='element'?this.selection.index:null,poses=P.positions(m.elements,selected,settings.count,settings.xi).map(row=>{
   const e=m.elements[row.index],a=m.nodes[e.node1-1],b=m.nodes[e.node2-1],av=[a.x,a.y,a.z],bv=[b.x,b.y,b.z],entry=catalog.get(e.sectionId);
   return {element:row.index,center:av.map((x,i)=>x+(bv[i]-x)*row.xi),basis:m.dimension===2?P.roll(M.editorAdapter.basis(m,e),e.sectionAngle||0):m.elementType===1436?window.MKEFSpatialConditions.basis(m,e):P.roll(P.trussBasis(av,bv),e.sectionAngle||0),loops:entry&&entry.properties.loops,selected:row.selected};
  });
  if(this.scene)this.scene.setSectionFrame({...settings,poses});this.sectionViewNote.textContent=`Positions ${poses.length}`+(m.elements.length*settings.count>600?' · limited to 600; the selected member has priority.':'');
 };
 proto.renderCatalog=function(){
  this.ensureCatalog();const tool=document.getElementById('assign-section');tool.disabled=!this.model.elements.length;tool.textContent=this.assignSequential?'Finish assignment':'Assign section';tool.setAttribute('aria-pressed',String(!!this.assignSequential));this.catalogBody.replaceChildren();this.renderMaterial();
  const m=this.model,panel=el('section',undefined,this.catalogBody);panel.id='spatial-catalog-sections';
  const box=el('section',undefined,panel);box.id='spatial-catalog';
  const editor=el('section',undefined,panel);editor.id='spatial-catalog-editor';el('p','Create a section or choose Edit in the catalog.',editor);
  if(m.elementType===1436){el('p',P.frameRestriction+' '+P.manualNote,box);const error=P.calculationError(m);if(error)el('p',error,box).setAttribute('role','status');}
  const add=el('button','Create section',box);add.id='spatial-section-new';add.onclick=()=>this.editSection(null);
  const page=this.sectionPage||0,total=Math.max(1,Math.ceil(m.sections.length/200)),prev=el('button','Previous',box),next=el('button','Next',box);prev.disabled=!page;next.disabled=page+1>=total;
  prev.onclick=()=>{this.sectionPage=page-1;this.render(m,this.selection);};next.onclick=()=>{this.sectionPage=page+1;this.render(m,this.selection);};el('span',` ${page+1} / ${total}`,box);
  const usage=new Map();m.elements.forEach((e,i)=>{if(e.sectionId){if(!usage.has(e.sectionId))usage.set(e.sectionId,[]);usage.get(e.sectionId).push(i+1);}});
  m.sections.slice(page*200,(page+1)*200).forEach((s,j)=>{const row=el('div',undefined,box),index=page*200+j,ids=usage.get(s.id)||[];el('span',`${s.id} · ${s.name} · ${s.type} · assignments ${ids.length} `,row);
   const edit=el('button','Edit',row),copy=el('button','Copy',row),remove=el('button','Delete',row);edit.onclick=()=>this.editSection(index);copy.onclick=()=>this.editSection(null,s);remove.onclick=()=>this.execute(h=>h.section(index,null));
   const show=el('button','Assignments',row);show.onclick=()=>{this.assignmentPage=0;this.showAssignments(s.id);};
  });
  if(this.sectionEditor){const draft=this.sectionEditor;this.editSection(draft.index,draft.value);}
  const assign=el('form',undefined,box);el('h4','Assign to selected member',assign);
  const selected=this.selection&&this.selection.kind==='element'?m.elements[this.selection.index]:null;
  const id=input(assign,'Section ID (0 — remove assignment)',selected&&selected.sectionId||this.assignId||0,'spatial-section-id');id.min=0;id.step=1;
  const angle=input(assign,'Final angle, °',selected?(m.elementType===1436?selected.angleDegrees:selected.sectionAngle||0):0,'spatial-section-angle');
  const apply=el('button','Apply assignment',assign);apply.id='spatial-section-assign';apply.type='submit';apply.disabled=!selected;
  assign.onsubmit=e=>{e.preventDefault();this.assignId=Number(id.value);this.execute(h=>h.assignSection(this.selection.index,Number(id.value),Number(angle.value)));};
  const seq=el('button',this.assignSequential?'Finish click assignment':'Assign by clicking',box);seq.id='spatial-section-sequence';seq.onclick=()=>{const sid=Number(id.value),section=m.sections.find(s=>s.id===sid);if(!Number.isInteger(sid)||!section){this.status('Choose an existing section.');return;}if(!P.allowedForElement(m.elementType,section.type)){this.catalogMessage.textContent=P.frameRestriction;return;}const active=!this.assignSequential;this.assignId=sid;this.assignAngle=Number(angle.value);this.closeCatalog();this.setTool('select');this.assignSequential=active;this.render(this.model,this.selection);this.status(active?'Click elements to assign section '+sid+'. Escape finishes assignment.':'');};
  el('p','Rotation changes local loads and releases with the local basis. Global loads and MPC coefficients remain unchanged. Removing an assignment retains numerical properties; the full σx map becomes unavailable.',box);
  this.showCatalogTab();
  document.getElementById('spatial-catalog-undo').disabled=!this.history.undo.length;document.getElementById('spatial-catalog-redo').disabled=!this.history.redo.length;
 };
 proto.showAssignments=function(id){
  const box=document.getElementById('spatial-section-assignments')||el('div',undefined,document.getElementById('spatial-catalog'));box.id='spatial-section-assignments';box.replaceChildren();
  const ids=this.model.elements.flatMap((e,i)=>e.sectionId===id?[i+1]:[]),page=this.assignmentPage||0;
  el('p',`Section ${id}: ${ids.length} members · page ${page+1}`,box);ids.slice(page*200,(page+1)*200).forEach(id=>{const b=el('button',id,box);b.onclick=()=>this.select({kind:'element',index:id-1});});
  for(const delta of [-1,1]){const b=el('button',delta<0?'Previous':'Next',box);b.disabled=delta<0?!page:(page+1)*200>=ids.length;b.onclick=()=>{this.assignmentPage=page+delta;this.showAssignments(id);};}
 };
 proto.editSection=function(index,copy){
  const existing=this.model.sections[index],value=copy||existing||{name:'New section',type:'rectangle',dimensions:{b:.1,h:.2,d:.1,t:.005,tw:.008,tf:.012},torsion:{mode:'auto'}};
  const draft=JSON.parse(JSON.stringify(value));if(index===null){let id=1;const ids=new Set(this.model.sections.map(s=>s.id));while(ids.has(id))id++;draft.id=id;}
  const parent=document.getElementById('spatial-catalog-editor');parent.replaceChildren();const box=el('form',undefined,parent);box.id='spatial-section-editor';
  this.sectionEditor={index,value:draft};el('h4',(index===null?'New section':'Section '+draft.id)+' · preview',box);
  const name=el('input',undefined,el('label','Name',box));name.id='spatial-section-name';name.value=draft.name;name.required=true;
  const type=el('select',undefined,el('label','Shape',box));type.id='spatial-section-type';Object.entries(S.types).forEach(([v,t])=>{const allowed=P.allowedForElement(this.model.elementType,v);if(!allowed&&v!==draft.type)return;const o=el('option',t+(allowed?'':' — legacy, calculation unavailable'),type);o.value=v;o.disabled=!allowed;});type.value=draft.type;
  if(!P.allowedForElement(this.model.elementType,draft.type))el('p',P.frameRestriction,box);
  const dimensions=el('div',undefined,box),inputs={};function fields(){dimensions.replaceChildren();Object.keys(inputs).forEach(k=>delete inputs[k]);S.fields(type.value).forEach(k=>{inputs[k]=input(dimensions,k==='inner_radius'?'Inside corner radius':k,draft.dimensions[k]??({b:.1,h:.2,d:.1,t:.005,tw:.008,tf:.012,inner_radius:0})[k],'spatial-section-'+k);if(k==='inner_radius')inputs[k].min='0';});if(inputs.inner_radius)el('p','Use model length units. 0 keeps sharp corners. For a positive inside radius, outside radius = inside radius + wall thickness. This geometry updates assigned members; rerun the analysis after editing.',dimensions);}fields();
  const mode=el('select',undefined,el('label','J',box));mode.id='spatial-section-j-mode';['auto','manual'].forEach(v=>{const o=el('option',v,mode);o.value=v;});mode.value=(draft.torsion||{}).mode||'auto';
  const manual=input(box,'Manual J (not recalculated when dimensions change)',draft.torsion&&draft.torsion.value||'','spatial-section-j');
  const preview=document.createElementNS('http://www.w3.org/2000/svg','svg');preview.id='spatial-section-preview';preview.setAttribute('width','260');preview.setAttribute('height','220');box.append(preview);
  const note=el('p','',box);note.id='spatial-section-properties';
  function read(){draft.name=name.value;draft.type=type.value;draft.dimensions=Object.fromEntries(Object.entries(inputs).map(([k,e])=>[k,e.value.trim()?Number(e.value):NaN]));draft.torsion=mode.value==='auto'?{mode:'auto'}:{mode:'manual',value:Number(manual.value)};return draft;}
  function update(){manual.disabled=mode.value!=='manual';manual.required=!manual.disabled;try{const p=P.properties(read());note.textContent=`A=${p.area} · Iy=${p.Iy} · Iz=${p.Iz} · J=${p.torsionConstant} · ${p.torsionMethod}`+(p.roundedTube?' · Inside radius='+p.roundedTube.ri+' · Outside radius='+Number(p.roundedTube.ro.toPrecision(12)):'')+(p.approximate?' · Approximate thin-walled J. For thick walls, specify a justified manual J.':'')+(p.open?' · Restrained warping and automatic shear-center offsets are not included; manual J does not remove these limitations.':'')+(mode.value==='manual'?' · Manual J remains unchanged when dimensions change.':'');
    const points=p.loops.flat(),ymin=Math.min(...points.map(v=>v.y)),ymax=Math.max(...points.map(v=>v.y)),zmin=Math.min(...points.map(v=>v.z)),zmax=Math.max(...points.map(v=>v.z)),margin=Math.max(ymax-ymin,zmax-zmin)*.1;
    preview.setAttribute('viewBox',[zmin-margin,-ymax-margin,zmax-zmin+2*margin,ymax-ymin+2*margin].join(' '));preview.replaceChildren();const path=document.createElementNS(preview.namespaceURI,'path');path.setAttribute('d',p.loops.map(loop=>loop.map((v,i)=>(i?'L':'M')+v.z+','+(-v.y)).join(' ')+'Z').join(' '));path.setAttribute('fill-rule','evenodd');path.setAttribute('fill','#9cc6df');preview.append(path);
   }catch(error){note.textContent=error.message;preview.replaceChildren();}}
  type.onchange=()=>{fields();update();};box.addEventListener('input',update);mode.onchange=update;update();
  const apply=el('button','Apply section',box);apply.id='spatial-section-apply';apply.type='submit';const cancel=el('button','Cancel',box);cancel.type='button';cancel.onclick=()=>{this.sectionEditor=null;box.remove();};
  box.onsubmit=e=>{e.preventDefault();const next=read();if(this.execute(h=>{h.section(index,next);this.sectionEditor=null;}))this.status('Section applied as one operation.');};
 };
 const pick=proto.handlePick;proto.handlePick=function(s){if(this.assignSequential&&s&&s.kind==='element'){const e=this.model.elements[s.index];this.execute(h=>h.assignSection(s.index,this.assignId,this.assignAngle??(this.model.elementType===1436?e.angleDegrees:e.sectionAngle||0)));return;}pick.call(this,s);};
 const member=proto.applyMember;proto.applyMember=function(draft,index){this.editingMember=index;try{const ok=member.call(this,draft,index);if(ok&&index===null){const defaults=this.memberDefaults[this.model.elementType];for(const key of Object.keys(steel))if(defaults[key]!==undefined)this.materialValues[key]=defaults[key];}return ok;}finally{this.editingMember=null;}};
 const value=proto.memberValue;proto.memberValue=function(draft){
  const e=value.call(this,draft),old=this.editingMember!==null&&this.model.elements[this.editingMember];
  e.sectionId=old?old.sectionId||0:Number(draft.sectionId||0);
  if(!Number.isInteger(e.sectionId)||e.sectionId<0)throw Error('Choose an existing catalog section.');
  if(this.model.elementType!==1436)e.sectionAngle=old?old.sectionAngle||0:0;
  return M.editorAdapter.resolve({...this.model,elements:[e]})[0];
 };
 const form=proto.renderMemberForm;proto.renderMemberForm=function(){
  form.call(this);
  if(this.tool!=='member'){
   const e=this.selection&&this.model.elements[this.selection.index];if(e&&e.sectionId){['area','momentOfInertia','Iy','Iz','torsionConstant'].forEach(k=>{const input=document.getElementById('spatial-member-'+k);if(input)input.readOnly=true;});el('p','Section properties come from the catalog. Remove the assignment to edit them independently.',this.form);}else if(e&&this.model.elementType===1436)el('p',P.manualNote,this.form);return;
  }
  const body=this.form.querySelector('form'),box=el('div',undefined,body);box.className='spatial-member-section-picker';body.prepend(box);
  const search=el('input',undefined,el('label','Search sections by name or ID',box));search.type='search';search.id='spatial-member-section-search';
  const select=el('select',undefined,el('label','Member section',box));select.id='spatial-member-sectionId';select.dataset.memberField='sectionId';
  let chosen=Number(this.memberDraft.sectionId??(this.memberDefaults[this.model.elementType]||{}).sectionId??0);
  const note=el('p','',box),catalog=this.model.sections;
  if(this.model.elementType===1436)el('p',P.manualNote,box);
  const populate=()=>{
   select.replaceChildren();el('option','No catalog — numerical properties',select).value='0';
   const query=search.value.trim().toLocaleLowerCase(),matches=catalog.filter(s=>P.allowedForElement(this.model.elementType,s.type)&&`${s.id} ${s.name}`.toLocaleLowerCase().includes(query)),rows=matches.slice(0,200),current=catalog.find(s=>s.id===chosen);
   if(current&&!rows.includes(current)){if(rows.length===200)rows.pop();rows.unshift(current);}
   rows.forEach(s=>{const allowed=P.allowedForElement(this.model.elementType,s.type),o=el('option',`${s.id} · ${s.name} · ${S.types[s.type]||s.type}`+(allowed?'':' — unavailable for 1436'),select);o.value=s.id;o.disabled=!allowed;});
   if(chosen&&!current)el('option',`Section ${chosen} was deleted — choose another`,select).value=chosen;
   select.value=String(chosen);note.textContent=matches.length>200?'Showing up to 200 sections. Refine the search.':this.model.elementType===1436?'A / Iy / Iz / J come from the catalog; E, ρ and G are specified separately.':'A comes from the catalog; E and ρ are specified separately.';
  };
  const properties=()=>{
   const section=catalog.find(s=>s.id===chosen),p=section?P.properties(section):null;
   for(const key of ['area','momentOfInertia','Iy','Iz','torsionConstant']){const input=document.getElementById('spatial-member-'+key);if(!input)continue;input.readOnly=!!chosen;if(p)input.value=key==='momentOfInertia'?S.properties(section,0).momentOfInertia:p[key];}
  };
  search.oninput=populate;select.onchange=()=>{chosen=Number(select.value);properties();this.memberDraft=this.readMemberDraft();};populate();properties();
  // Preserve unfinished material inputs when catalog operations rebuild the panel.
  const save=()=>{this.memberDraft=this.readMemberDraft();};body.addEventListener('input',save);body.addEventListener('change',save);
  const manage=el('button','Materials and Sections…',box);manage.type='button';manage.id='spatial-member-catalog';manage.onclick=()=>this.openCatalog();
 };
 const leave=proto.leave;proto.leave=function(){this.sectionEditor=null;this.assignSequential=false;if(this.catalogDialog)this.catalogDialog.close();leave.call(this);};
})(window.MKEFPre,window.MKEFSpatialSections,window.MKEFSections);
