(function (global) {
  'use strict';
  var M = global.MKEFPre, F = M.editorAdapter || M.spatialFormat, keys = ['nodes','elements','supports','loads','elementLoads','mpcs','releases','sections'];
  function same(a,b) {
    if (a === b) return true;
    if (!a || !b || typeof a !== 'object' || typeof b !== 'object') return false;
    var names = Object.keys(a); return names.length === Object.keys(b).length && names.every(function (key) { return Object.prototype.hasOwnProperty.call(b,key) && same(a[key],b[key]); });
  }
  function clone(a) { return JSON.parse(JSON.stringify(a)); }
  function pointError(m,n,index) {
    if (![n.x,n.y,n.z].every(Number.isFinite)) return 'XYZ must be finite.';
    if (index !== undefined && same(m.nodes[index],n)) return '';
    if (m.nodes.some(function (p,i) { return i !== index && p.x === n.x && p.y === n.y && p.z === n.z; })) return 'Another node occupies this point; automatic merging is disabled.';
    if (index !== undefined) {
      var candidate = Object.assign({},m,{nodes:m.nodes.slice()}); candidate.nodes[index] = n;
      for (var i = 0; i < m.elements.length; i++) {
        var e = m.elements[i]; if (e.node1 !== index+1 && e.node2 !== index+1) continue;
        var error = F.memberError(candidate,e); if (error) return 'Member '+(i+1)+': '+error;
      }
    }
    return '';
  }
  function autoReference(m,a,b) {
    var p = m.nodes[a-1], q = m.nodes[b-1]; if (!p || !q) return [0,1,0];
    var d = [q.x-p.x,q.y-p.y,q.z-p.z], order = [1,2,0];
    order.sort(function (i,j) { return Math.abs(d[i])-Math.abs(d[j]); }); var ref = [0,0,0]; ref[order[0]] = 1; return ref;
  }
  function patch(key,index,before,after) { return {key:key,index:index,before:before,after:after}; }
  function candidate(m,patches,back) {
    var result = Object.assign({},m);
    keys.forEach(function (key) {
      var list = patches.filter(function (p) { return p.key === key; }); if (!list.length) return;
      var byIndex = new Map(), shift = 0;
      list.sort(function (a,b) { return a.index-b.index; }).forEach(function (p) {
        var index = p.index+shift; shift += (p.before === null ? 1 : 0)-(p.after === null ? 1 : 0);
        if (!back) index = p.index;
        var entry = {before:back ? p.after : p.before,after:back ? p.before : p.after};
        if (!byIndex.has(index)) byIndex.set(index,[]); byIndex.get(index).push(entry);
      });
      var rows = [];
      for (var i = 0; i <= m[key].length; i++) {
        var entries = byIndex.get(i) || [], replaced = false;
        entries.forEach(function (p) { if (p.after !== null) rows.push(p.after); if (p.before !== null) replaced = true; });
        if (i < m[key].length && !replaced) rows.push(m[key][i]);
      }
      result[key] = rows;
    });
    ['analysis','name'].forEach(function(key){var change=patches.find(function(p){return p.key===key;});if(change)result[key]=back?change.before:change.after;});
    if(M.editorAdapter)result.elements=M.editorAdapter.resolve(result);else if(global.MKEFSpatialSections)result.elements=global.MKEFSpatialSections.resolve(result);
    return result;
  }
  // Only changed records are retained. Candidate arrays are temporary and share unchanged records.
  function History(model) { this.model = model; this.undo = []; this.redo = []; this.revision = 0; this.saved = 0; this.next = 1; this.selection = null; }
  History.prototype.apply = function (patches,selection) {
    patches = patches.filter(function (p) { return !same(p.before,p.after); }); if (!patches.length) return false;
    var next = candidate(this.model,patches,false); var errors=F.validate(next,false);if(errors.length)throw Error(errors.join("\n"));
    const sections=global.MKEFSpatialSections;
    if(sections&&next.elementType===1436)for(const p of patches){
      if(!p.after)continue;
      if(p.key==='sections'&&!sections.allowedForElement(1436,p.after.type)&&(!p.before||p.before.type!==p.after.type))throw Error(sections.frameRestriction);
      if(p.key==='elements'&&p.after.sectionId&&(!p.before||p.before.sectionId!==p.after.sectionId)){
        const s=next.sections.find(s=>s.id===p.after.sectionId);
        if(s&&!sections.allowedForElement(1436,s.type))throw Error(sections.frameRestriction);
      }
    }
    var command = {patches:patches,beforeSelection:this.selection,afterSelection:selection,beforeRevision:this.revision,afterRevision:this.next++};
    keys.concat('analysis','name').forEach(function (key) { this.model[key] = next[key]; },this);
    this.undo.push(command); if (this.undo.length > 100) this.undo.shift(); this.redo = [];
    this.revision = command.afterRevision; this.selection = selection; return true;
  };
  History.prototype.travel = function (back) {
    var from = back ? this.undo : this.redo, to = back ? this.redo : this.undo; if (!from.length) return false;
    var command = from.pop(), next = candidate(this.model,command.patches,back);
    keys.concat('analysis','name').forEach(function (key) { this.model[key] = next[key]; },this);
    this.revision = back ? command.beforeRevision : command.afterRevision; this.selection = back ? command.beforeSelection : command.afterSelection; to.push(command); return true;
  };
  History.prototype.text = function () { return this.revision === 0 && this.model.sourceText ? this.model.sourceText : F.serialize(this.model); };
  History.prototype.dirty = function () { return this.revision !== this.saved; };
  History.prototype.section=function(index,value){
    if(index!==null&&!this.model.sections[index])throw Error('Section no longer exists.');
    if(value===null){var id=this.model.sections[index].id,used=this.model.elements.flatMap(function(e,i){return e.sectionId===id?[i+1]:[];});if(used.length)throw Error('Section is used by members: '+used.join(', '));}
    return this.apply([patch('sections',index===null?this.model.sections.length:index,index===null?null:this.model.sections[index],value===null?null:clone(value))],this.selection);
  };
  History.prototype.assignSection=function(index,id,angle){
    var old=this.model.elements[index];if(!old)throw Error('Choose an existing member.');
    const s=this.model.sections.find(s=>s.id===id),P=global.MKEFSpatialSections;
    if(s&&P&&!P.allowedForElement(this.model.elementType,s.type))throw Error(P.frameRestriction);
    var e=Object.assign({},old,{sectionId:id});if(!Number.isFinite(angle))throw Error('Angle must be finite.');
    if(this.model.elementType===1436)e.angleDegrees=angle;else e.sectionAngle=angle;
    return this.apply([patch('elements',index,old,e)],{kind:'element',index:index});
  };
  History.prototype.condition = function (key,index,value) {
    if (global.MKEFSpatialConditions.keys.indexOf(key)<0) throw new Error('Unknown condition kind.');
    var at=index===null ? this.model[key].length : index,old=index===null ? null : this.model[key][index];
    if (index!==null&&!old) throw new Error('Condition no longer exists.');
    if(value&&key==='loads'){
      if(!['static','transient'].includes(this.model.analysis.type))throw new Error('Choose a static or transient analysis before assigning loads.');
      if(this.model.analysis.type==='transient'&&value.type===10&&(!old||old.type!==10))throw new Error('Transient type 10 is preserved from import/static models; create type 11, 12 or 13.');
    }
    if(value&&key==='elementLoads'&&this.model.analysis.type!=='static')throw new Error('Distributed loads require static analysis.');
    return this.apply([patch(key,at,old,value===null ? null : clone(value))],value===null ? null : {kind:key,index:at});
  };
  History.prototype.rename = function(name){return this.apply([patch('name',0,this.model.name,name)],this.selection);};
  History.prototype.setAnalysis = function (analysis,confirmed) {
    if (this.model.analysis.type==='static'&&analysis.type==='transient'&&this.model.loads.some(function (r) { return r.type===10; })&&!confirmed) throw new Error('Confirm: type 10 becomes a one-step load in transient analysis.');
    return this.apply([patch('analysis',0,this.model.analysis,clone(analysis))],this.selection);
  };
  History.prototype.addMPCs=function(rows){
    var first=this.model.mpcs.length;return this.apply(rows.map(function(row){return patch('mpcs',first,null,clone(row));}),{kind:'mpcs',index:first});
  };
  History.prototype.setReleaseMask=function(elementId,end,components){
    var m=this.model;
    if(![113,1436].includes(this.model.elementType)||!Number.isInteger(elementId)||!this.model.elements[elementId-1]||![1,2].includes(end)||components.some(function(c){return !(m.dimension===2?['Mz']:['Mx','My','Mz']).includes(c);})||new Set(components).size!==components.length)throw new Error('Choose a frame, end and valid release mask.');
    var patches=[],present=new Set(),kept=[],count=this.model.releases.length;
    this.model.releases.forEach(function(r,i){var target=r.elementId===elementId&&r.end===end;if(target)present.add(r.component);if(target&&!components.includes(r.component))patches.push(patch('releases',i,r,null));else kept.push(r);});
    components.forEach(function(c){if(!present.has(c)){var r={elementId:elementId,end:end,component:c};patches.push(patch('releases',count,null,r));kept.push(r);}});
    var selected=kept.findIndex(function(r){return r.elementId===elementId&&r.end===end;});
    return this.apply(patches,selected<0?{kind:'element',index:elementId-1}:{kind:'releases',index:selected});
  };
  History.prototype.addNode = function (point) {
    var found = this.model.nodes.findIndex(function (n) { return same(n,point); });
    if (found >= 0) { this.selection = {kind:'node',index:found}; return false; }
    var error = pointError(this.model,point); if (error) throw new Error(error);
    var index = this.model.nodes.length; return this.apply([patch('nodes',index,null,clone(point))],{kind:'node',index:index});
  };
  History.prototype.moveNode = function (index,point) {
    if (!this.model.nodes[index]) throw new Error('Select an existing node.');
    var error = pointError(this.model,point,index); if (error) throw new Error(error);
    var patches=[patch('nodes',index,this.model.nodes[index],clone(point))];
    if(this.model.dimension===2&&!same(this.model.nodes[index],point))this.model.elementLoads.forEach(function(r,i){var e=this.model.elements[r.elementId-1];if(e.node1===index+1||e.node2===index+1)patches.push(patch('elementLoads',i,r,null));},this);
    return this.apply(patches,{kind:'node',index:index});
  };
  History.prototype.member = function (index,element) {
    var m = this.model;
    if (element.sectionId && global.MKEFSpatialSections) element = (M.editorAdapter?M.editorAdapter.resolve:global.MKEFSpatialSections.resolve)(Object.assign({},m,{elements:[element]}))[0];
    var error = F.memberError(m,element); if (error) throw new Error(error);
    var old = index === null ? null : m.elements[index]; if (index !== null && !old) throw new Error('Select an existing member.');
    if (!old || old.node1 !== element.node1 || old.node2 !== element.node2) {
      if (m.elements.some(function (e,i) { return i !== index && ((e.node1 === element.node1 && e.node2 === element.node2) || (e.node1 === element.node2 && e.node2 === element.node1)); })) throw new Error('A member already connects these nodes.');
    }
    var at = index === null ? m.elements.length : index;
    var patches=[patch('elements',at,old,clone(element))];
    if(m.dimension===2&&old&&['node1','node2'].concat(M.editorAdapter.info(m.elementType).properties).some(function(k){return old[k]!==element[k];}))m.elementLoads.forEach(function(r,i){if(r.elementId===at+1)patches.push(patch('elementLoads',i,r,null));});
    return this.apply(patches,{kind:'element',index:at});
  };
  function deletion(m,selection,replacement) {
    if (!selection) throw new Error('Select an object to delete.');
    var nodeId = selection.kind === 'node' ? selection.index+1 : 0, ids = new Set(), patches = [], report = {};
    if (nodeId && !m.nodes[nodeId-1] || !nodeId && !m.elements[selection.index]) throw new Error('Object no longer exists.');
    m.elements.forEach(function (e,i) { if (nodeId ? e.node1 === nodeId || e.node2 === nodeId : i === selection.index) ids.add(i+1); });
    var map = new Map(), kept = 0; m.elements.forEach(function (e,i) { if (!ids.has(i+1)) map.set(i+1,++kept); });
    function node(n) { return n > nodeId && nodeId ? n-1 : n; }
    keys.forEach(function (key) {
      report[key] = [];
      m[key].forEach(function (item,i) {
        var remove = key === 'nodes' ? i+1 === nodeId : key === 'elements' ? ids.has(i+1) : key === 'supports' || key === 'loads' ? item.node === nodeId : key === 'mpcs' ? item.depNode === nodeId || item.masters.some(function (a) { return a.node === nodeId; }) : ids.has(item.elementId);
        var after = item;
        if (remove) { after = null; report[key].push(i+1); }
        else if (key === 'elements') after = Object.assign({},item,{node1:node(item.node1),node2:node(item.node2)});
        else if (key === 'supports' || key === 'loads') after = Object.assign({},item,{node:node(item.node)});
        else if (key === 'mpcs') after = Object.assign({},item,{depNode:node(item.depNode),masters:item.masters.map(function (a) { return Object.assign({},a,{node:node(a.node)}); })});
        else if (key === 'releases' || key === 'elementLoads') after = Object.assign({},item,{elementId:map.get(item.elementId)});
        if (!same(item,after)) patches.push(patch(key,i,item,after));
      });
    });
    var needsMonitor = nodeId && m.analysis.type === 'transient' && m.analysis.monitorNode === nodeId;
    if (nodeId && m.analysis.type === 'transient') {
      var a = Object.assign({},m.analysis), target = needsMonitor ? replacement : a.monitorNode;
      if (needsMonitor && (!Number.isInteger(target) || target === nodeId || !m.nodes[target-1])) target = null;
      a.monitorNode = target === null ? null : node(target); patches.push(patch('analysis',0,m.analysis,a));
    }
    return {patches:patches,report:report,needsMonitor:!!needsMonitor,valid:!needsMonitor || Number.isInteger(replacement) && replacement !== nodeId && !!m.nodes[replacement-1]};
  }
  History.prototype.remove = function (selection,replacement) {
    var change = deletion(this.model,selection,replacement); if (!change.valid) throw new Error('Choose a remaining monitor node before deleting this node.');
    return this.apply(change.patches,null);
  };
  M.spatialEdit = {History:History,pointError:pointError,autoReference:autoReference,deletion:deletion};
}(window));
