(function (global) {
  'use strict';
  var M = global.MKEFPre;
  function el(tag, text, parent) { var e = document.createElement(tag); if (text !== undefined) e.textContent = text; if (parent) parent.appendChild(e); return e; }
  M.SpatialPanel = function (select, changed) {
    this.select = select; this.pages = {node:0, element:0}; this.model = null; this.scene = null;
    this.changed = changed; this.tool = 'select';
    this.memberDefaults = {}; this.memberDraft = {}; this.pendingMember = null;
    this.root = document.getElementById('spatial-tables');
    this.host = document.getElementById('spatial-viewport');
    this.form = el('div',undefined,document.getElementById('selection-details').parentElement); this.form.id = 'spatial-edit-form'; this.form.hidden = true;
    this.analysisForm = el('div',undefined,document.getElementById('analysis-type').parentElement); this.analysisForm.id='spatial-analysis-form';this.analysisForm.hidden=true;
    var self = this;
    this.deleteButton=el('button','Delete',this.form.parentElement);this.deleteButton.id='selection-delete';this.deleteButton.type='button';
    var original=el('button','Download original',this.root);original.id='spatial-original';original.type='button';
    this.deleteDialog = el('dialog',undefined,document.body); this.deleteDialog.id = 'spatial-delete-dialog';
    this.deleteButton.addEventListener('click',function () { self.confirmDelete(); });
    document.querySelectorAll('[data-spatial-tool]').forEach(function (button) { button.addEventListener('click',function () { self.setTool(button.dataset.spatialTool); }); });
    ['spatial-snap-grid','spatial-snap-nodes'].forEach(function (id) { document.getElementById(id).addEventListener('change',function () { if (self.scene) self.scene.setSnap({grid:document.getElementById('spatial-snap-grid').checked,nodes:document.getElementById('spatial-snap-nodes').checked}); }); });
    document.getElementById('spatial-original').addEventListener('click',function () {
      if (!self.model || !self.model.sourceText) return;
      var url = URL.createObjectURL(new Blob([self.model.sourceText],{type:'text/plain;charset=utf-8'})), link = el('a');
      link.href = url; link.download = self.model.name+'-original.txt'; link.click(); setTimeout(function () { URL.revokeObjectURL(url); },0);
    });
    document.querySelectorAll('[data-view]').forEach(function (button) {
      button.addEventListener('click',function () {
        if (!self.scene) return;
        var view = button.dataset.view;
        if (view === 'fit') self.scene.fitModel(); else if (view === 'origin') self.scene.resetOrigin(); else self.scene.setView(view);
      });
    });
    ['spatial-plane','spatial-offset','spatial-step','spatial-grid'].forEach(function (id) {
      document.getElementById(id).addEventListener('change',function () {
        if (!self.scene) return;
        var offset = document.getElementById('spatial-offset'), step = document.getElementById('spatial-step');
        var a = offset.value.trim() ? Number(offset.value) : NaN, b = step.value.trim() ? Number(step.value) : NaN;
        offset.setCustomValidity(Number.isFinite(a) ? '' : 'Enter a finite offset.'); step.setCustomValidity(Number.isFinite(b) && b > 0 ? '' : 'Enter a positive finite step.');
        if (!offset.reportValidity() || !step.reportValidity()) return;
        self.scene.setWorkPlane({plane:document.getElementById('spatial-plane').value,offset:a,step:b,visible:document.getElementById('spatial-grid').checked});
      });
    });
  };
  M.SpatialPanel.prototype.leave = function () {
    if (this.scene) this.scene.dispose(); this.scene = null; this.model = null;
    this.host.replaceChildren(); this.unavailable = false;
    this.form.hidden = true; this.history = null; this.tool = 'select';
    this.analysisForm.hidden=true;
    this.deleteDialog.close();
  };
  M.SpatialPanel.prototype.status = function (text) {
    if (!text) this.operationError = null;
    document.getElementById('spatial-edit-status').textContent = text || '';
    if (this.catalogMessage) this.catalogMessage.textContent = text || '';
  };
  // Pointer hints must not replace a rejected operation's diagnostic.
  M.SpatialPanel.prototype.previewStatus = function (text) { if (!this.operationError) this.status(text); };
  M.SpatialPanel.prototype.changeAnalysis = function(type) {
    const old=this.model.analysis,confirmed=!(old.type==='static'&&type==='transient'&&this.model.loads.some(r=>r.type===10))||global.confirm('Type 10 loads will act on the first step only in transient analysis. Continue?');
    if(!confirmed){document.getElementById('analysis-type').value=old.type;return;}
    this.analysisInvalid=false;
    if(!this.execute(h=>h.setAnalysis({...old,type},confirmed))){document.getElementById('analysis-type').value=old.type;this.analysisError.textContent=this.operationError;}
  };
  M.SpatialPanel.prototype.renderAnalysis = function () {
    const self=this,analysis=this.model.analysis;this.analysisForm.replaceChildren();this.analysisForm.hidden=false;this.analysisInvalid=false;
    const form=el('form',undefined,this.analysisForm),inputs={};form.onsubmit=e=>e.preventDefault();
    this.analysisError=el('p','',this.analysisForm);this.analysisError.id='analysis-error';this.analysisError.setAttribute('role','alert');
    if(analysis.type!=='transient')return;
    for(const [key,title] of [['timeStep','Time step'],['duration','Duration'],['monitorNode','Monitor node']]){const input=el('input',undefined,el('label',title,form));input.type='number';input.step=key==='monitorNode'?'1':'any';input.required=true;input.value=analysis[key];input.id='spatial-analysis-'+key;inputs[key]=input;}
    const select=el('select',undefined,el('label','Monitor component',form));select.id='spatial-analysis-monitorDOF';M.editorAdapter.info(this.model.elementType).dofLabels.forEach((name,i)=>{el('option',name,select).value=i+1;});select.value=analysis.monitorDOF;inputs.monitorDOF=select;
    form.onchange=function(){const next={...analysis};for(const key in inputs)next[key]=inputs[key].value.trim()?Number(inputs[key].value):NaN;
      const errors=M.editorAdapter.validate({...self.model,analysis:next},false);self.analysisInvalid=!!errors.length;
      if(errors.length){self.analysisError.textContent=errors[0];self.changed();return;}
      self.execute(h=>h.setAnalysis(next,true));
    };
  };
  M.SpatialPanel.prototype.setTool = function (tool) {
    this.assignSequential = false;
    this.tool = tool; if (this.scene) this.scene.setTool(tool); this.status('');
    this.select(tool === 'node' || tool === 'member' ? null : this.selection);
  };
  M.SpatialPanel.prototype.cancelPending = function () {
    this.pendingMember = null; this.memberDraft.node1 = ''; this.memberDraft.node2 = '';
    if (this.tool === 'member') ['node1','node2'].forEach(function (key) { var input = document.getElementById('spatial-member-'+key); if (input) input.value = ''; });
  };
  M.SpatialPanel.prototype.readMemberDraft = function () {
    var draft = {}; this.form.querySelectorAll('[data-member-field]').forEach(function (input) { draft[input.dataset.memberField] = input.value; });
    var auto = document.getElementById('spatial-auto-reference'); draft.auto = !auto || auto.checked; return draft;
  };
  M.SpatialPanel.prototype.memberValue = function (draft) {
    function number(key) { return draft[key] === undefined || String(draft[key]).trim() === '' ? NaN : Number(draft[key]); }
    var e = {code:this.model.elementType,node1:number('node1'),node2:number('node2')};
    M.editorAdapter.info(this.model.elementType).properties.forEach(function (key) { e[key] = number(key); });
    if (this.model.elementType === 1436) { e.referenceY = draft.auto ? M.spatialEdit.autoReference(this.model,e.node1,e.node2) : ['refX','refY','refZ'].map(number); e.angleDegrees = number('angleDegrees'); }
    return e;
  };
  M.SpatialPanel.prototype.applyMember = function (draft,index) {
    var self = this;
    return this.execute(function (h) {
      var value = self.memberValue(draft);
      h.member(index,value);
      if (index === null) {
        var defaults = {}; M.editorAdapter.info(self.model.elementType).properties.forEach(function (key) { defaults[key] = value[key]; });
        defaults.sectionId = value.sectionId || 0;
        self.memberDefaults[self.model.elementType] = defaults; self.memberDraft = {}; self.pendingMember = null; if (self.scene) self.scene.setPendingMember(null);
      }
    });
  };
  M.SpatialPanel.prototype.handlePick = function (next) {
    if (this.tool !== 'member' || !next || next.kind !== 'node') { this.select(next); return; }
    this.memberDraft = this.readMemberDraft();
    if (this.pendingMember === null) {
      this.pendingMember = next.index; this.memberDraft.node1 = next.index+1; this.memberDraft.node2 = '';
      if (this.scene) this.scene.setPendingMember(next.index); this.select(next); this.status('Select the second existing node.');
    } else {
      this.memberDraft.node2 = next.index+1;
      if (!this.applyMember(this.memberDraft,null)) { var message = document.getElementById('spatial-edit-status').textContent; this.select(next); this.status(message); }
    }
  };
  M.SpatialPanel.prototype.execute = function (operation) {
    try { operation(this.history); this.status(''); this.changed(); return true; }
    catch (error) { this.operationError = error.message; this.status(error.message); return false; }
  };
  M.SpatialPanel.prototype.travel = function (back) {
    if (this.scene) this.scene.cancelGesture();
    if (this.history.travel(back)) { this.status(''); this.changed(); }
  };
  M.SpatialPanel.prototype.confirmDelete = function () {
    if (!this.selection) return;
    if (global.MKEFSpatialConditions.keys.includes(this.selection.kind)) { var record=this.selection; this.execute(function(h){if(record.kind==='releases'){var r=h.model.releases[record.index];h.setReleaseMask(r.elementId,r.end,[]);}else h.condition(record.kind,record.index,null);}); return; }
    var self = this, selected = this.selection, change = M.spatialEdit.deletion(this.model,selected), dialog = this.deleteDialog;
    dialog.replaceChildren(); el('h3','Delete geometry and dependent records',dialog);
    var report = Object.keys(change.report).filter(function (key) { return change.report[key].length; }).map(function (key) { return key+' ('+change.report[key].length+'): '+change.report[key].join(', '); }).join('\n');
    el('pre',report,dialog); el('p','All listed records will be deleted. Undo restores them together with their references.',dialog);
    var form = el('form',undefined,dialog), monitor;
    if (change.needsMonitor) {
      var label = el('label','Replacement monitor node (existing ID before deletion)',form); monitor = el('input',undefined,label); monitor.id = 'spatial-monitor-replacement'; monitor.type = 'number'; monitor.step = '1'; monitor.min = '1'; monitor.max = String(this.model.nodes.length); monitor.required = true;
    }
    var message = el('p','',form); message.setAttribute('role','alert');
    var cancel = el('button','Cancel',form); cancel.type = 'button'; cancel.onclick = function () { dialog.close(); };
    var apply = el('button','Delete and renumber',form); apply.type = 'submit'; apply.id = 'spatial-delete-confirm';
    form.addEventListener('submit',function (event) {
      event.preventDefault();
      try { self.history.remove(selected,monitor ? Number(monitor.value) : undefined); dialog.close(); if (self.scene) self.scene.cancelGesture(); self.changed(); }
      catch (error) { message.textContent = error.message; }
    });
    dialog.showModal();
  };
  M.SpatialPanel.prototype.renderForm = function () {
    var self = this, selected = this.selection, creating = this.tool === 'node'; this.form.replaceChildren(); this.form.hidden = false;
    if(this.scene)this.scene.setConditionFocus(null);
    if(this.renderConditionForm && this.renderConditionForm())return;
    if (this.tool === 'member' || !creating && selected && selected.kind === 'element') { this.renderMemberForm(); return; }
    if (!creating && (!selected || selected.kind !== 'node')) return;
    var node = creating ? {x:0,y:0,z:0} : this.model.nodes[selected.index];
    if (creating) { var plane = this.scene ? this.scene.workPlane : {plane:document.getElementById('spatial-plane').value,offset:Number(document.getElementById('spatial-offset').value)}; node[{XY:'z',XZ:'y',YZ:'x'}[plane.plane]] = plane.offset; }
    el('h3',(creating ? 'New node · exact ' : 'Edit ')+(this.model.dimension===2?'XY':'XYZ'),this.form);
    var form = el('form',undefined,this.form), inputs = {};
    (this.model.dimension===2?['x','y']:['x','y','z']).forEach(function (key) {
      var label = el('label',key.toUpperCase(),form), input = el('input',undefined,label); input.type = 'number'; input.step = 'any'; input.required = true; input.id = 'spatial-node-'+key; input.value = node[key]; inputs[key] = input;
    });
    var apply = el('button',creating ? 'Create node' : 'Apply XYZ',form); apply.type = 'submit'; apply.id = 'spatial-node-apply';
    form.addEventListener('submit',function (event) {
      event.preventDefault(); var point = {z:0}; Object.keys(inputs).forEach(function (key) { point[key] = inputs[key].value.trim() ? Number(inputs[key].value) : NaN; });
      self.execute(function (h) { if (creating) h.addNode(point); else h.moveNode(selected.index,point); });
    });
    if (!creating) el('p',this.model.dimension===2?'Changing coordinates removes distributed loads on connected elements. Undo restores the operation.':'Conditions are preserved; MPC coefficients are not recalculated.',this.form).className = 'control-note';
  };
  M.SpatialPanel.prototype.renderMemberForm = function () {
    var self = this, creating = this.tool === 'member', index = creating ? null : this.selection.index;
    var value = creating ? Object.assign({node1:'',node2:'',refX:0,refY:1,refZ:0,angleDegrees:0,auto:true},this.memberDefaults[this.model.elementType] || {},this.memberDraft) : Object.assign({},this.model.elements[index],{auto:false});
    if (!creating && value.referenceY) { value.refX = value.referenceY[0]; value.refY = value.referenceY[1]; value.refZ = value.referenceY[2]; }
    el('h3',creating ? 'New element · required properties' : 'Edit element '+(index+1),this.form);
    var form = el('form',undefined,this.form), fields = ['node1','node2'].concat(M.editorAdapter.info(this.model.elementType).properties), auto, inputs = {};
    if (this.model.elementType === 1436) fields = fields.concat(['refX','refY','refZ','angleDegrees']);
    var names = {node1:'Node 1',node2:'Node 2',area:'A',youngsModulus:'E',density:'ρ',shearModulus:'G',torsionConstant:'J',refX:'Reference Y · X',refY:'Reference Y · Y',refZ:'Reference Y · Z',angleDegrees:'Angle, °'};
    fields.forEach(function (key) {
      var label = el('label',names[key] || key,form), input = el('input',undefined,label); input.id = 'spatial-member-'+key; input.dataset.memberField = key; input.type = 'number'; input.step = key === 'node1' || key === 'node2' ? '1' : 'any'; input.required = true; input.value = value[key] === undefined ? '' : value[key]; inputs[key] = input;
      if (key.indexOf('ref') === 0) input.addEventListener('input',function () { if (auto) auto.checked = false; });
    });
    function updateReference() {
      if (!auto || !auto.checked) return;
      var ref = M.spatialEdit.autoReference(self.model,Number(inputs.node1.value),Number(inputs.node2.value)); ['refX','refY','refZ'].forEach(function (key,i) { inputs[key].value = ref[i]; });
    }
    if (creating && this.model.elementType === 1436) {
      var label = el('label','Suggest reference Y automatically',form); auto = el('input',undefined,label); auto.type = 'checkbox'; auto.id = 'spatial-auto-reference'; auto.checked = value.auto; auto.addEventListener('change',updateReference);
      inputs.node1.addEventListener('change',updateReference); inputs.node2.addEventListener('change',updateReference); updateReference();
    }
    var apply = el('button',creating ? 'Create element' : 'Apply properties',form); apply.type = 'submit'; apply.id = 'spatial-member-apply';
    form.addEventListener('submit',function (event) { event.preventDefault(); var draft = self.readMemberDraft(); if (!creating) draft.auto = false; self.memberDraft = creating ? draft : self.memberDraft; self.applyMember(draft,index); });
    el('p',this.model.dimension===2?'Changing numerical properties or endpoints removes this element’s distributed loads. Catalog assignments preserve loads.':'Conditions are preserved; MPC coefficients are not recalculated. J is independent of Iy and Iz.',this.form).className='control-note';
  };
  M.SpatialPanel.prototype.render = function (model, selection) {
    var self = this, changed = this.model !== model;
    ['spatial-plane','spatial-offset'].forEach(id=>document.getElementById(id).disabled=model.dimension===2);document.querySelector('[data-view=iso]').disabled=model.dimension===2;
    if(changed&&model.dimension===2){document.getElementById('spatial-plane').value='XY';document.getElementById('spatial-offset').value=0;}
    if (changed) { this.deleteDialog.close(); this.status(''); this.pages = {node:0, element:0}; this.history = new M.spatialEdit.History(model); this.tool = 'select'; this.memberDraft = {}; this.pendingMember = null; this.renderedRevision = -1; } this.model = model;
    var changedSelection = selection !== this.selection; this.selection = selection;
    this.history.selection = selection;
    if(changed||this.renderedRevision!==this.history.revision) this.renderAnalysis(model.analysis.type);
    this.deleteButton.disabled = !selection;this.deleteButton.hidden=!selection;this.deleteButton.setAttribute('aria-label',selection?'Delete '+selection.kind+' '+(selection.index+1):'Delete');
    document.getElementById('spatial-original').disabled = !model.sourceText;
    document.querySelectorAll('[data-spatial-tool]').forEach(function (button) { button.classList.toggle('active',button.dataset.spatialTool === self.tool); button.setAttribute('aria-pressed',String(button.dataset.spatialTool === self.tool)); });
    var originalButton=document.getElementById('spatial-original');this.root.replaceChildren();
    el('h2','Scene objects',this.root);
    [['node','Nodes',model.nodes],['element','Elements',model.elements]].forEach(function (group) {
      var kind = group[0], items = group[2];
      if (changedSelection && selection && selection.kind === kind) self.pages[kind] = Math.floor(selection.index/200);
      var section = el('section',undefined,self.root); section.className = 'spatial-list';
      el('h3',group[1] + ' (' + items.length + ')',section);
      var nav = el('div',undefined,section), total = Math.max(1,Math.ceil(items.length/200));
      var prev = el('button','Previous',nav), status = el('span',' Page ' + (self.pages[kind]+1) + ' / ' + total + ' ',nav), next = el('button','Next',nav);
      status.setAttribute('role','status'); prev.disabled = self.pages[kind] === 0; next.disabled = self.pages[kind]+1 >= total;
      prev.onclick = function () { self.pages[kind]--; self.render(model,self.selection); }; next.onclick = function () { self.pages[kind]++; self.render(model,self.selection); };
      var scroll = el('div',undefined,section); scroll.className = 'table-scroll';
      var table = el('table',undefined,scroll), head = el('tr',undefined,el('thead',undefined,table));
      var keys = kind === 'node' ? (model.dimension===2?['x','y']:['x','y','z']) : ['node1','node2','sectionId'];
      ['ID'].concat(keys).forEach(function (key) { el('th',key,head); });
      var body = el('tbody',undefined,table); body.id = 'spatial-' + kind + '-rows';
      items.slice(self.pages[kind]*200,(self.pages[kind]+1)*200).forEach(function (item,j) {
        var index = self.pages[kind]*200+j, row = el('tr',undefined,body); row.tabIndex = 0;
        if (selection && selection.kind === kind && selection.index === index) row.className = 'selected-row';
        el('td',String(index+1),row); keys.forEach(function (key) { el('td',key==='sectionId'?(item[key]||'Manual'):String(Number(item[key].toPrecision(8))),row); });
        row.onclick = function () { self.tool='select';if(self.scene)self.scene.setTool('select');self.handlePick({kind:kind,index:index}); };
        row.onkeydown = function (e) { if (e.key === 'Enter') row.onclick(); };
      });
    });
    var summary = el('section',undefined,this.root); summary.id = 'spatial-summary';
    el('h3','Analysis conditions',summary);summary.appendChild(originalButton);
    el('p','Supports: ' + model.supports.length + ' · Nodal loads: ' + model.loads.length + ' · Member loads: ' + model.elementLoads.length + ' · MPC: ' + model.mpcs.length + ' · Releases: ' + model.releases.length,summary);
    var details = el('details',undefined,summary); el('summary','Original TXT (preserved unchanged)',details);
    var source = el('textarea',undefined,details); source.id = 'spatial-source'; source.readOnly = true; source.value = model.sourceText;
    var current = el('details',undefined,summary); el('summary','Current TXT',current);
    var currentText = el('textarea',undefined,current); currentText.id = 'spatial-current'; currentText.readOnly = true;
    try { currentText.value = this.history.text(); } catch (error) { currentText.value = error.message; }
    document.getElementById('node-coordinate-editor').hidden = true;
    var selected = selection && (selection.kind === 'node' ? model.nodes : selection.kind === 'element' ? model.elements : model[selection.kind]||[])[selection.index];
    document.getElementById('selection-details').textContent = selected ? selection.kind + ' ' + (selection.index+1) : 'Nothing selected';
    if (global.MKEFSpatial.SpatialScene && !this.scene && !this.unavailable) {
      try { this.scene = new global.MKEFSpatial.SpatialScene(this.host,{select:function (next) { self.handlePick(next); },cancel:function () { self.cancelPending(); },createNode:function (point) { self.execute(function (h) { h.addNode(point); }); },moveNode:function (index,point) { self.execute(function (h) { h.moveNode(index,point); }); },validateNode:function (index,point) { return M.spatialEdit.pointError(self.model,point,index); },preview:function (text) { self.previewStatus(text); },change:function (view) {
        document.getElementById('spatial-grid-status').textContent = view.workPlane.plane + ' · ' + ({XY:'Z',XZ:'Y',YZ:'X'}[view.workPlane.plane]) + '=' + view.workPlane.offset + ' · displayed spacing ' + Number(view.displayStep.toPrecision(6));
        document.getElementById('spatial-plane').value = view.workPlane.plane;
        document.getElementById('spatial-offset').value = view.workPlane.offset;
        document.getElementById('spatial-step').value = view.workPlane.step;
        document.getElementById('spatial-grid').checked = view.workPlane.visible;
      }},{dimension:model.dimension}); }
      catch (error) { this.unavailable = true; this.host.replaceChildren(); el('p','WebGL 2 unavailable. Geometry display requires WebGL 2. The source and summary remain available. ' + error.message,this.host); }
    }
    if (this.scene) {
      if (changed) { this.scene.setTool('select'); this.scene.setSnap({grid:document.getElementById('spatial-snap-grid').checked,nodes:document.getElementById('spatial-snap-nodes').checked}); }
      if (changed || this.renderedRevision !== this.history.revision) this.scene.setModel({nodes:model.nodes,elements:model.elements},{preserveView:!changed});
      this.scene.setSelection(selection);
      this.scene.callbacks.choose=function(hits){self.chooseCondition(hits);};
    }
    if(this.renderConditions)this.renderConditions(changed || this.renderedRevision !== this.history.revision);
    if (changed || changedSelection || this.renderedRevision !== this.history.revision || this.renderedTool !== this.tool) this.renderForm();
    this.renderedRevision = this.history.revision; this.renderedTool = this.tool;
  };
}(window));
