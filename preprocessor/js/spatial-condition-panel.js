(function(global){
  'use strict';
  var M=global.MKEFPre,C=M.editorConditions,P=M.SpatialPanel.prototype;
  var names={supports:'Supports',loads:'Nodal loads',elementLoads:'Distributed loads',mpcs:'MPC',releases:'Releases'};
  function el(tag,text,parent){var e=document.createElement(tag);if(text!==undefined)e.textContent=text;if(parent)parent.appendChild(e);return e;}
  function field(form,key,title,value,options){
    var label=el('label',title,form),input=el(options?'select':'input',undefined,label);input.id='condition-'+key;
    if(options)options.forEach(function(o){var option=el('option',o[1],input);option.value=o[0];});
    else {input.type='number';input.step='any';input.required=true;}
    input.value=value===undefined?'':value;return input;
  }
  function number(input){return input.value.trim()?Number(input.value):NaN;}
  function mask(form,labels,values){return labels.map(function(name,i){var label=el('label',name,form),input=el('input',undefined,label);input.type='checkbox';input.id='condition-'+name;input.checked=!!values[i];return input;});}
  function describe(kind,r,analysis,model){
    if(kind==='supports')return 'Node '+r.node+' · '+C.labels.filter(function(_,i){return r.mask[i];}).join(', ');
    if(kind==='mpcs')return C.equation(r,model);
    if(kind==='releases')return 'Member '+r.elementId+' · end '+r.end+' · '+r.component;
    if(kind==='elementLoads')return 'Member '+r.elementId+' · '+(r.coordinateSystem===1?'local':'global')+' · '+(r.type===20?'Uniform':'Linear')+' · q₁=('+r.components.slice(0,3).join(', ')+')'+(r.type===21?' · q₂=('+r.components.slice(3).join(', ')+')':'');
    var law=r.type===11?'sin · '+r.frequency+' Hz':r.type===13?'step':r.type===12||analysis==='transient'?'pulse':'static';
    return 'Node '+r.node+' · '+law+' · '+r.components.flatMap(function(x,i){return x?['Fx','Fy','Fz','Mx','My','Mz'][i]+'='+x:[];}).join(' · ');
  }
  P.conditionTools=['supports','loads','elementLoads','mpcs','releases'];
  P.conditionEditors={};
  P.conditionUI={el:el,field:field,number:number,mask:mask};
  var oldPick=P.handlePick;
  P.handlePick=function(next){
    if(next&&C.keys.includes(next.kind)){this.tool='select';if(this.scene)this.scene.setTool('select');this.select(next);return;}
    if(C.keys.includes(this.tool)){
      var target=this.tool==='supports'||this.tool==='loads'||this.tool==='mpcs'?'node':'element';
      if(next&&next.kind===target){this.conditionTarget=next.index+1;var input=document.getElementById('condition-'+(target==='node'?'node':'elementId'));if(input){input.value=this.conditionTarget;input.dispatchEvent(new Event('change'));}if(this.scene)this.scene.setSelection(next);this.updateConditionFocus();}
      return;
    }
    oldPick.call(this,next);
  };
  P.startCondition=function(kind){
    this.assignSequential=false;
    this.conditionTarget='';this.tool=kind;this.status('');if(this.scene)this.scene.setTool(kind);this.select(null);
  };
  P.chooseCondition=function(hits){
    var self=this,dialog=this.deleteDialog,page=0;
    function render(){dialog.replaceChildren();el('h3','Objects at pick location ('+hits.length+')',dialog);
      hits.slice(page*200,(page+1)*200).forEach(function(hit){var b=el('button',(names[hit.kind]||hit.kind)+' '+(hit.index+1),dialog);b.onclick=function(){dialog.close();self.handlePick(hit);};});
      var prev=el('button','Previous',dialog),next=el('button','Next',dialog),cancel=el('button','Cancel',dialog);prev.disabled=!page;next.disabled=(page+1)*200>=hits.length;prev.onclick=function(){page--;render();};next.onclick=function(){page++;render();};cancel.onclick=function(){dialog.close();};
    }render();dialog.showModal();
  };
  P.renderConditions=function(updated){
    var self=this,root=el('section',undefined,this.root);root.id='spatial-conditions';
    if(!this.conditionDefaults)this.conditionDefaults={};if(!this.conditionPages)this.conditionPages={};
    if(!this.conditionLayers)this.conditionLayers={loads:true,supports:true,mpcs:true,releases:true,axes:false};
    var tools=document.getElementById('spatial-condition-tools');if(!tools){tools=el('div');tools.id='spatial-condition-tools';tools.className='spatial-edit-tools';document.getElementById('spatial-controls').insertBefore(tools,document.getElementById('spatial-edit-status'));}tools.replaceChildren();
    this.conditionTools.forEach(function(kind){var b=el('button',names[kind],tools);b.dataset.conditionTool=kind;b.classList.toggle('active',self.tool===kind);b.onclick=function(){self.startCondition(kind);};b.disabled=!self.model.elements.length||(kind==='releases'||kind==='elementLoads')&&!M.editorAdapter.info(self.model.elementType).releases.length||(kind==='loads'&&self.model.analysis.type==='modal')||(kind==='elementLoads'&&self.model.analysis.type!=='static');});
    var layers=el('div');layers.className='spatial-condition-layers';self.root.insertBefore(layers,self.root.children[1]);el('span','Show: ',layers);
    Object.keys(this.conditionLayers).forEach(function(key){var label=el('label',names[key]||'Local axes',layers),input=el('input',undefined,label);input.type='checkbox';input.id='condition-layer-'+key;input.checked=self.conditionLayers[key];input.onchange=function(){self.conditionLayers[key]=input.checked;if(self.scene)self.scene.setConditionLayers(self.conditionLayers);};});
    if(updated||!this.conditionDiagnostics)this.conditionDiagnostics=C.validate(this.model);
    var warning=el('p',this.conditionDiagnostics.warnings.join('\n'),root);warning.id='condition-warnings';warning.setAttribute('role','status');
    C.keys.forEach(function(kind){
      var items=self.model[kind],section=el('section',undefined,root);section.className='spatial-list';el('h3',names[kind]+' ('+items.length+')',section);
      var page=self.conditionPages[kind]||0;if(self.selection&&self.selection.kind===kind&&self.conditionLastSelection!==self.selection)page=Math.floor(self.selection.index/200);
      page=Math.max(0,Math.min(page,Math.ceil(items.length/200)-1));self.conditionPages[kind]=page;
      var prev=el('button','Previous',section),label=el('span',' '+(page+1)+' / '+Math.max(1,Math.ceil(items.length/200))+' ',section),next=el('button','Next',section);prev.disabled=!page;next.disabled=(page+1)*200>=items.length;
      prev.onclick=function(){self.conditionPages[kind]--;self.render(self.model,self.selection);};next.onclick=function(){self.conditionPages[kind]++;self.render(self.model,self.selection);};
      var scroll=el('div',undefined,section);scroll.className='table-scroll';var body=el('tbody',undefined,el('table',undefined,scroll));body.id='condition-'+kind+'-rows';
      items.slice(page*200,(page+1)*200).forEach(function(row,j){var i=page*200+j,tr=el('tr',undefined,body);tr.tabIndex=0;el('td',i+1,tr);el('td',describe(kind,row,self.model.analysis.type,self.model),tr);if(self.selection&&self.selection.kind===kind&&self.selection.index===i)tr.className='selected-row';tr.onclick=function(){self.handlePick({kind:kind,index:i});};tr.onkeydown=function(e){if(e.key==='Enter')tr.onclick();};});
    });this.conditionLastSelection=this.selection;
    if(this.scene){this.scene.setConditionLayers(this.conditionLayers);if(updated)this.scene.setConditions(M.editorAdapter.sceneConditions(this.model,this.conditionDiagnostics));}
  };
  P.renderConditionForm=function(){
    var creating=C.keys.includes(this.tool),kind=creating?this.tool:this.selection&&this.selection.kind;
    if(!C.keys.includes(kind))return false;
    var self=this,index=creating?null:this.selection.index;
    var value=creating?JSON.parse(JSON.stringify(this.conditionDefaults[kind]||{})):this.model[kind][index];
    if(creating){value.node=this.conditionTarget;value.elementId=this.conditionTarget;}
    el('h3',(creating?'Create: ':'Edit: ')+names[kind]+(creating?'':' '+(index+1)),this.form);
    var form=el('form',undefined,this.form);form.id='condition-form';
    if(!this.conditionEditors[kind]){el('p','Select a record in the list.',form);return true;}
    var read=this.conditionEditors[kind].call(this,form,value,creating);
    var apply=el('button',creating?'Create condition':'Apply condition',form);apply.type='submit';apply.id='condition-apply';
    form.onsubmit=function(event){event.preventDefault();self.execute(function(h){var row=read();if(row===null)return;
      if(row.helperRows)h.addMPCs(row.helperRows);else if(row.releaseMask)h.setReleaseMask(row.elementId,row.end,row.releaseMask);else h.condition(kind,index,row);self.conditionDefaults[kind]=row.helperRows?row.defaults:row;if(creating){self.conditionTarget='';h.selection=null;}
    });};
    if(!creating){var remove=el('button',kind==='releases'?'Clear end release mask':'Delete condition',form);remove.type='button';remove.id='condition-delete';remove.onclick=function(){self.execute(function(h){if(kind==='releases')h.setReleaseMask(value.elementId,value.end,[]);else h.condition(kind,index,null);});};}
    el('p',this.model.dimension===2?'Editing member geometry removes its distributed loads. MPC coefficients are not recalculated automatically.':'Local load and release components remain unchanged when geometry changes. MPC coefficients are not recalculated automatically.',this.form);
    this.updateConditionFocus();return true;
  };
  P.updateConditionFocus=function(){
    if(!this.scene)return;
    var kind=C.keys.includes(this.tool)?this.tool:this.selection&&this.selection.kind,member=document.getElementById('condition-elementId'),end=document.getElementById('condition-end'),axes=document.getElementById('condition-axes');
    var focus=(kind==='releases'||kind==='elementLoads'&&axes&&number(axes)===1)&&member&&this.model.elements[number(member)-1]?{elementId:number(member),end:kind==='releases'?number(end):null}:null;
    this.scene.setConditionFocus(focus);
  };
  P.conditionEditors.supports=function(form,value){
    var node=field(form,'node','Node',value.node),ids=M.editorAdapter.physicalDOFs(this.model),checks=mask(form,M.editorAdapter.info(this.model.elementType).dofLabels,ids.map(i=>(value.mask||[])[i]));
    return function(){return {node:number(node),mask:Array.from({length:6},function(_,i){return checks[ids.indexOf(i)]&&checks[ids.indexOf(i)].checked?1:0;})};};
  };
  P.conditionEditors.loads=function(form,value,creating){
    var analysis=this.model.analysis.type;
    var options=analysis==='static'?[[10,'10 · static']]:analysis==='transient'?[[11,'11 · harmonic'],[12,'12 · pulse'],[13,'13 · step']]:[];
    if(!creating&&value.type===10&&analysis==='transient')options.unshift([10,'10 · imported pulse']);
    if(!options.length)el('p','Choose static or transient analysis first.',form);
    var node=field(form,'node','Node',value.node),type=field(form,'type','Law',options.some(function(o){return o[0]===value.type;})?value.type:options[0]?.[0],options),freq=field(form,'frequency','Frequency, Hz',value.frequency||1);
    var ids=M.editorAdapter.physicalDOFs(this.model),inputs=ids.map(function(i){var name=['Fx','Fy','Fz','Mx','My','Mz'][i];return field(form,name,name,(value.components||[])[i]||0);});
    function law(){freq.parentElement.hidden=Number(type.value)!==11;freq.required=Number(type.value)===11;}type.onchange=law;law();
    el('p','11: F₀ sin(2πft); 12: first step; 13: constant from the start of the history. The scene shows amplitudes without animation.',form);
    return function(){if(!options.length)throw new Error('Choose an analysis before creating a load.');return {node:number(node),type:number(type),frequency:Number(type.value)===11?number(freq):null,components:Array.from({length:6},function(_,i){return inputs[ids.indexOf(i)]?number(inputs[ids.indexOf(i)]):0;})};};
  };
  P.conditionEditors.elementLoads=function(form,value){
    var self=this,planar=this.model.dimension===2,oldType=value.type||20,oldAxes=value.coordinateSystem||1;
    var member=field(form,'elementId','Member',value.elementId),type=field(form,'type','Distribution',oldType,[[20,'Uniform'],[21,'Linear']]),axes=field(form,'axes','Axes',oldAxes,[[1,'Local'],[2,'Global']]);
    var inputs=Array.from({length:6},function(_,i){return field(form,'q'+i,(i<3?'End 1 · ':'End 2 · ')+['qx','qy','qz'][i%3],(value.components||[])[i]||0);});
    function visibility(){inputs.forEach(function(input,i){input.parentElement.hidden=(planar&&i%3===2)||i>=3&&oldType===20;input.required=!input.parentElement.hidden;});}
    function row(){return {type:oldType,coordinateSystem:oldAxes,elementId:number(member),components:inputs.slice(0,oldType===20?3:6).map(number)};}
    axes.onchange=function(){try{var next=C.convertLoadAxes(self.model,row(),number(axes));next.components.forEach(function(x,i){inputs[i].value=x;});oldAxes=next.coordinateSystem;self.updateConditionFocus();}catch(e){axes.value=oldAxes;self.status(e.message);}};member.onchange=function(){self.updateConditionFocus();};
    type.onchange=function(){var before=row(),target=number(type),confirmed=target!==20||!before.components.slice(0,3).some(function(x,i){return x!==before.components[i+3];})||global.confirm('Replace different end intensities with their average?');
      try{var next=C.convertLoadType(before,target,confirmed);next.components.forEach(function(x,i){inputs[i].value=x;});oldType=next.type;visibility();}catch(e){type.value=oldType;self.status(e.message);}
    };visibility();el('p','Full member length. Changing the coordinate system preserves the physical load vector. Available for static frames only.',form);
    return function(){if(self.model.analysis.type!=='static')throw new Error('Distributed loads require static analysis.');return row();};
  };
  P.conditionEditors.mpcs=function(form,value,creating){
    var self=this,options=M.editorAdapter.info(this.model.elementType).dofLabels.map(function(name,i){return [i+1,name];});
    var mode=field(form,'mode','Mode',value.mode||'manual',creating?[['manual','Equation'],['equal','Equal components'],['axis','Node on axis']]:[['manual','Equation']]);
    var node=field(form,'node','Dependent / intermediate node',creating?value.node:value.depNode),dependent=field(form,'depDOF','Dependent component',value.depDOF||1,options);
    var master=field(form,'masterNode','Master node / axis start',value.masterNode||''),end=field(form,'endNode','Axis end',value.endNode||'');
    var equalBox=el('div',undefined,form),checks=mask(equalBox,options.map(function(o){return o[1];}),value.mask||[]);
    var termsBox=el('div',undefined,form),terms=[];
    function addTerm(value){var box=el('div',undefined,termsBox),i=terms.length,term={box:box,node:field(box,'master-'+i,'Master node',value.node),dof:field(box,'dof-'+i,'Component',value.dof||1,options),coefficient:field(box,'coefficient-'+i,'Coefficient',value.coefficient===undefined?1:value.coefficient)};terms.push(term);var remove=el('button','Remove term',box);remove.type='button';remove.onclick=function(){term.removed=true;box.remove();};}
    (value.masters||[{}]).forEach(addTerm);var add=el('button','Add term',form);add.type='button';add.id='condition-add-master';add.onclick=function(){addTerm({});};
    function visibility(){var kind=mode.value;dependent.parentElement.hidden=kind!=='manual';termsBox.hidden=kind!=='manual';add.hidden=kind!=='manual';master.parentElement.hidden=kind==='manual';master.required=kind!=='manual';end.parentElement.hidden=kind!=='axis';end.required=kind==='axis';equalBox.hidden=kind!=='equal';terms.forEach(function(t){t.node.required=t.coefficient.required=kind==='manual';});}
    mode.onchange=visibility;visibility();
    var preview=el('pre','Right-hand side: 0. Review the equations before applying.',form);preview.id='condition-equations';var show=el('button','Show equations',form);show.type='button';show.id='condition-preview';
    var accepted=null;
    function read(){if(mode.value==='equal')return C.equalComponents(self.model,number(node),number(master),checks.flatMap(function(c,i){return c.checked?[i+1]:[];}));if(mode.value==='axis')return C.axisHelper(self.model,number(master),number(end),number(node));return [{depNode:number(node),depDOF:number(dependent),rhs:0,masters:terms.filter(function(t){return !t.removed;}).map(function(t){return {node:number(t.node),dof:number(t.dof),coefficient:number(t.coefficient)};})}];}
    show.onclick=function(){try{var rows=read();preview.textContent=rows.map(r=>C.equation(r,self.model)).join('\n');accepted=JSON.stringify(rows);self.status('');}catch(e){self.status(e.message);}};
    if(!creating)preview.textContent=C.equation(value,self.model);
    return function(){var rows=read();if(mode.value!=='manual'&&JSON.stringify(rows)!==accepted){preview.textContent=rows.map(C.equation).join('\n');accepted=JSON.stringify(rows);throw new Error('Review the generated equations and click Create condition again.');}return mode.value==='manual'?rows[0]:{helperRows:rows,defaults:{mode:mode.value,masterNode:number(master),endNode:number(end),mask:checks.map(function(c){return Number(c.checked);})}};};
  };
  P.conditionEditors.releases=function(form,value){
    var self=this,member=field(form,'elementId','Frame',value.elementId),end=field(form,'end','End',value.end||1,[[1,'1 · first node'],[2,'2 · second node']]),components=M.editorAdapter.info(this.model.elementType).releases,checks=mask(form,components,[]);
    function refresh(){var rows=self.model.releases.filter(function(r){return r.elementId===number(member)&&r.end===number(end);});checks.forEach(function(input,i){input.checked=rows.some(function(r){return r.component===components[i];});});self.updateConditionFocus();}
    member.onchange=refresh;end.onchange=refresh;refresh();
    el('p','Moments in the local axes at the frame end. An empty mask removes all releases at the selected end.',form);
    return function(){return {elementId:number(member),end:number(end),releaseMask:checks.flatMap(function(input,i){return input.checked?[components[i]]:[];})};};
  };
}(window));
