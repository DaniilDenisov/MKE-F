(function (global) {
  'use strict';
  var M = global.MKEFPre;
  function typeInfo(code) { return global.MKEFSpatial.registry.types.find(function (t) { return t.code === code; }); }
  function isSpatial(code) { var t = typeInfo(code); return !!t && t.dimension === 3; }
  function detects(text) {
    return text.split(/\r?\n/).some(function (line) {
      var match = /^elems_(\d+)$/.exec(line.trim()); return match && isSpatial(Number(match[1]));
    });
  }
  function empty(code) {
    var m = M.caseFormat.newModel('', code); m.dimension = 3; m.sourceText = ''; m.blocks = {}; return m;
  }
  function parse(text) {
    if (new TextEncoder().encode(text).length > M.config.maxFileBytes) throw new Error('File exceeds the 10 MiB limit.');
    var lines = text.replace(/^\uFEFF/, '').split(/\r?\n/), cursor = 0, m = empty(0), seen = {}, info, assignments=[];
    function fail(line, message) { throw new Error('Line ' + line + ': ' + message); }
    function next() {
      while (cursor < lines.length) { var s = lines[cursor++].trim(); if (s && s[0] !== '#') return { text:s, line:cursor }; }
      fail(cursor + 1, 'Unexpected end of file.');
    }
    function nums(r, n) {
      var p = r.text.split(',');
      if (p.length !== n) fail(r.line, 'Expected ' + n + ' fields.');
      return p.map(function (s) { var v = s.trim() === '' ? NaN : Number(s); if (!Number.isFinite(v)) fail(r.line, 'Expected finite numbers.'); return v; });
    }
    function id(value, limit, r, what) { if (!Number.isInteger(value) || value < 1 || value > limit) fail(r.line, 'Invalid ' + what + ' reference.'); }
    function node(value, r) { id(value, m.nodes.length, r, 'node'); }
    function dof(value, r) { id(value, info.dofLabels.length, r, 'DOF'); }
    while (cursor < lines.length) {
      if (!lines[cursor].trim() || lines[cursor].trim()[0] === '#') { cursor++; continue; }
      var marker = next(), key = marker.text, r, v;
      if (key === 'analysis') {
        if (seen.analysis) fail(marker.line, 'Duplicate analysis.'); seen.analysis = true;
        r = next(); var p = r.text.split(',').map(function (s) { return s.trim(); });
        if (['static','modal','transient'].indexOf(p[0]) < 0 || p.length !== (p[0] === 'transient' ? 5 : 1)) fail(r.line, 'Invalid analysis.');
        m.analysis.type = p[0];
        if (p[0] === 'transient') {
          v = nums({text:p.slice(1).join(','), line:r.line},4);
          if (v[0] <= 0 || v[1] <= 0 || !Number.isInteger(v[2]) || v[2] < 1 || !Number.isInteger(v[3]) || v[3] < 1) fail(r.line,'Invalid transient settings.');
          m.analysis = {type:p[0],timeStep:v[0],duration:v[1],monitorNode:v[2],monitorDOF:v[3]};
          m.analysisLine = r.line;
        }
        continue;
      }
      var elem = /^elems_(\d+)$/.exec(key);
      var allowed = ['sections','element_sections','nodes','bcfix_3d','bcforce_stat_3d','bcforce_harm_3d','bcforce_pulse_3d','bcforce_step_3d','eload_uniform_3d','eload_linear_3d','mpc','releases'];
      if (!elem && allowed.indexOf(key) < 0) fail(marker.line, 'Unsupported block "' + key + '".');
      if (elem) {
        var code = Number(elem[1]);
        if (!isSpatial(code) || (info && code !== info.code)) fail(marker.line,'Mixed or unsupported element families.');
        if (!seen.nodes) fail(marker.line,'Nodes must precede elements.');
        info = typeInfo(code); m.elementType = code;
      } else if (!['nodes','sections','element_sections'].includes(key) && !info) fail(marker.line,'Elements must precede conditions.');
      if (key === 'nodes' && seen.nodes) fail(marker.line,'Duplicate nodes.');
      seen[key] = true;
      r = next(); var count = nums(r,1)[0];
      if (!Number.isInteger(count) || count < 0 || count > 100000) fail(r.line,'Invalid block count (limit 100000).');
      if (elem && m.elements.length + count > M.config.maxElements) fail(r.line,'Element limit exceeded.');
      m.blocks[key] = (m.blocks[key] || 0) + count;
      for (var i = 0; i < count; i++) {
        r = next();
        if(key==='sections'){try{var section=JSON.parse(r.text);global.MKEFSpatialSections.properties(section);m.sections.push(section);}catch(error){fail(r.line,error.message);}}
        else if(key==='element_sections'){v=nums(r,3);assignments.push({values:v,line:r.line});}
        else if (key === 'nodes') { v = nums(r,3); m.nodes.push({x:v[0], y:v[1], z:v[2]}); }
        else if (elem) {
          v = nums(r, info.family === 'frame' ? 14 : 6);
          if (v[0] !== info.code) fail(r.line,'Element code differs from block.'); node(v[1],r); node(v[2],r);
          var a = m.nodes[v[1]-1], b = m.nodes[v[2]-1], d = [b.x-a.x,b.y-a.y,b.z-a.z], length = Math.hypot.apply(null,d);
          if (!Number.isFinite(length) || length <= 0) fail(r.line,'Invalid element length.');
          var item = {node1:v[1],node2:v[2],code:info.code};
          info.properties.forEach(function (name,j) { if (v[3+j] <= 0) fail(r.line,'Properties must be positive.'); item[name] = v[3+j]; });
          if (info.family === 'frame') {
            var ref = v.slice(10,13), norm = Math.hypot.apply(null,ref), x = d.map(function (n) { return n/length; });
            var unit = ref.map(function (n) { return n/norm; }), dot = unit.reduce(function (s,n,j) { return s+n*x[j]; },0);
            if (!Number.isFinite(norm) || norm === 0 || Math.hypot.apply(null,unit.map(function (n,j) { return n-dot*x[j]; })) < 1e-10) fail(r.line,'Degenerate reference Y vector.');
            item.referenceY = ref; item.angleDegrees = v[13];
          }
          m.elements.push(item);
        } else if (key === 'bcfix_3d') {
          v = nums(r,7); node(v[0],r);
          if (v.slice(1).some(function (n) { return n !== 0 && n !== 1; }) || (info.family === 'truss' && v.slice(4).some(Boolean))) fail(r.line,'Invalid support mask.');
          m.supports.push({node:v[0],mask:v.slice(1)});
        } else if (key.indexOf('bcforce_') === 0) {
          var type = {bcforce_stat_3d:10,bcforce_harm_3d:11,bcforce_pulse_3d:12,bcforce_step_3d:13}[key];
          v = nums(r,type === 11 ? 9 : 8); node(v[1],r);
          if (v[0] !== type || (type === 11 && v[8] <= 0) || (info.family === 'truss' && v.slice(5,8).some(Boolean))) fail(r.line,'Invalid load type, frequency or truss moment.');
          m.loads.push({type:v[0],node:v[1],components:v.slice(2,8),frequency:type === 11 ? v[8] : null});
        } else if (key.indexOf('eload_') === 0) {
          var linear = key === 'eload_linear_3d'; v = nums(r,linear ? 9 : 6); id(v[1],m.elements.length,r,'element');
          if (info.family !== 'frame' || v[0] !== (linear ? 21 : 20) || [1,2].indexOf(v[2]) < 0) fail(r.line,'Invalid member load.');
          m.elementLoads.push({type:v[0],elementId:v[1],coordinateSystem:v[2],components:v.slice(3)});
        } else if (key === 'releases') {
          var fields = r.text.split(',').map(function (s) { return s.trim(); });
          if (fields.length !== 3 || info.releases.indexOf(fields[2]) < 0) fail(r.line,'Invalid release component.');
          id(Number(fields[0]),m.elements.length,r,'element'); id(Number(fields[1]),2,r,'end'); m.releases.push({elementId:Number(fields[0]),end:Number(fields[1]),component:fields[2]});
        } else if (key === 'mpc') {
          var n = Number(r.text.split(',')[3]);
          if (!Number.isInteger(n) || n < 1) fail(r.line,'Invalid MPC master count.');
          v = nums(r,4+3*n); node(v[0],r); dof(v[1],r);
          if (v[2] !== 0) fail(r.line,'Only homogeneous MPC supported.');
          for (var j = 0; j < n; j++) { node(v[4+3*j],r); dof(v[5+3*j],r); }
          var constraint = {depNode:v[0],depDOF:v[1],rhs:v[2],masters:[]};
          for (var t = 0; t < n; t++) constraint.masters.push({node:v[4+3*t],dof:v[5+3*t],coefficient:v[6+3*t]});
          m.mpcs.push(constraint);
        }
      }
    }
    if (!seen.nodes || !info) fail(cursor,'Missing nodes or spatial element block.');
    if (m.analysis.type === 'transient') { r = {line:m.analysisLine}; node(m.analysis.monitorNode,r); dof(m.analysis.monitorDOF,r); }
    if (m.elementLoads.length && m.analysis.type && m.analysis.type !== 'static') fail(cursor,'Distributed loads are static only.');
    if ((m.analysis.type === 'static' && m.loads.some(function (v) { return v.type !== 10; })) || (m.analysis.type === 'modal' && m.loads.length)) fail(cursor,'Nodal loads do not match the analysis type.');
    var assigned=new Set();assignments.forEach(function(row){var v=row.values,id=v[0],sid=v[1];if(!Number.isInteger(id)||!m.elements[id-1]||assigned.has(id)||!Number.isInteger(sid)||sid<0||(sid===0&&m.elementType!==1133))fail(row.line,'Invalid or duplicate section assignment.');assigned.add(id);var e=m.elements[id-1];e.sectionId=sid;if(m.elementType===1436)e.angleDegrees=v[2];else e.sectionAngle=v[2];});
    if(global.MKEFSpatialSections){try{m.elements=global.MKEFSpatialSections.resolve(m);}catch(error){fail(cursor,error.message);}}
    var problems=validate(m,false); if (problems.length) fail(cursor,problems.slice(0,12).join('\n'));
    m.sourceText = text; return m;
  }
  function memberError(m,e) {
    var info = typeInfo(m.elementType), a = m.nodes[e.node1-1], b = m.nodes[e.node2-1];
    if (!Number.isInteger(e.node1) || !Number.isInteger(e.node2) || !a || !b || e.node1 === e.node2) return 'Choose two distinct existing nodes.';
    if (info.properties.some(function (key) { return !Number.isFinite(e[key]) || e[key] <= 0; })) return 'All material and section properties must be positive and finite.';
    var d = [b.x-a.x,b.y-a.y,b.z-a.z], length = Math.hypot.apply(null,d);
    if (!Number.isFinite(length) || length <= 0) return 'Member length must be positive and finite.';
    if (info.family === 'frame') {
      var ref = e.referenceY, norm = ref && Math.hypot.apply(null,ref);
      if (!ref || ref.length !== 3 || !ref.every(Number.isFinite) || !Number.isFinite(norm) || norm === 0 || !Number.isFinite(e.angleDegrees)) return 'Enter a finite reference Y vector and angle.';
      var x = d.map(function (v) { return v/length; }), unit = ref.map(function (v) { return v/norm; });
      var dot = unit.reduce(function (s,v,i) { return s+v*x[i]; },0);
      if (Math.hypot.apply(null,unit.map(function (v,i) { return v-dot*x[i]; })) < 1e-10) return 'Reference Y is parallel or nearly parallel to the member; change the reference vector.';
    }
    return '';
  }
  function validate(m, forCalculation) {
    var errors = [];
    if (!isSpatial(m.elementType)) return ['Choose a spatial element family.'];
    if(global.MKEFSpatialSections)try{var resolved=global.MKEFSpatialSections.resolve(m);resolved.forEach(function(e,i){['area','Iy','Iz','torsionConstant'].forEach(function(k){if(e[k]!==m.elements[i][k])throw Error('Member '+(i+1)+': catalog property mismatch.');});});}catch(error){errors.push(error.message);}
    if (m.nodes.length > M.config.maxNodes || m.elements.length > M.config.maxElements) errors.push('Limit: 100000 nodes and members.');
    m.nodes.forEach(function (n,i) { if (![n.x,n.y,n.z].every(Number.isFinite)) errors.push('Node '+(i+1)+': XYZ must be finite.'); });
    m.elements.forEach(function (e,i) { var error = memberError(m,e); if (error) errors.push('Member '+(i+1)+': '+error); });
    if (!errors.length) errors=errors.concat(global.MKEFSpatialConditions.validate(m).errors);
    if (forCalculation && (!m.nodes.length || !m.elements.length || !m.analysis.type)) errors.push('Geometry draft: nodes, members and an analysis configuration are required for calculation.');
    if(forCalculation&&global.MKEFSpatialSections){var sectionError=global.MKEFSpatialSections.calculationError(m);if(sectionError)errors.push(sectionError);}
    return errors;
  }
  function serialize(m) {
    var errors = validate(m,false); if (errors.length) throw new Error(errors.slice(0,20).join('\n')+(errors.length>20?'\n… '+errors.length+' diagnostics.':''));
    var out = [], info = typeInfo(m.elementType);
    function block(name,rows,always) { if (rows.length || always) { out.push(name,String(rows.length)); rows.forEach(function (r) { out.push(r.join(',')); }); } }
    if (m.analysis.type) {
      out.push('analysis',m.analysis.type === 'transient' ? ['transient',m.analysis.timeStep,m.analysis.duration,m.analysis.monitorNode,m.analysis.monitorDOF].join(',') : m.analysis.type);
    }
    block('nodes',m.nodes.map(function (n) { return [n.x,n.y,n.z]; }),true);
    block('elems_'+m.elementType,m.elements.map(function (e) { return [m.elementType,e.node1,e.node2].concat(info.properties.map(function (k) { return e[k]; }),info.family === 'frame' ? e.referenceY.concat(e.angleDegrees) : []); }),true);
    if(m.sections.length){out.push('sections',String(m.sections.length));m.sections.forEach(function(s){out.push(JSON.stringify(s));});}
    // sectionId=0 is a spatial truss orientation-only record after detaching.
    block('element_sections',m.elements.flatMap(function(e,i){return e.sectionId||m.elementType===1133&&e.sectionAngle?[[i+1,e.sectionId||0,m.elementType===1436?e.angleDegrees:e.sectionAngle||0]]:[];}));
    block('bcfix_3d',m.supports.map(function (r) { return [r.node].concat(r.mask); }));
    block('mpc',m.mpcs.map(function (r) { var row = [r.depNode,r.depDOF,r.rhs,r.masters.length]; r.masters.forEach(function (a) { row.push(a.node,a.dof,a.coefficient); }); return row; }));
    block('releases',m.releases.map(function (r) { return [r.elementId,r.end,r.component]; }));
    // Keep record order, using repeated blocks for noncontiguous types.
    [['loads',{10:'bcforce_stat_3d',11:'bcforce_harm_3d',12:'bcforce_pulse_3d',13:'bcforce_step_3d'}],['elementLoads',{20:'eload_uniform_3d',21:'eload_linear_3d'}]].forEach(function (group) {
      var rows = m[group[0]];
      for (var i = 0; i < rows.length;) {
        var type = rows[i].type, end = i+1; while (end < rows.length && rows[end].type === type) end++;
        block(group[1][type],rows.slice(i,end).map(function (r) { return group[0] === 'loads' ? [r.type,r.node].concat(r.components,r.type === 11 ? [r.frequency] : []) : [r.type,r.elementId,r.coordinateSystem].concat(r.components); })); i = end;
      }
    });
    var text = out.join('\n')+'\n';
    if (new TextEncoder().encode(text).length > M.config.maxFileBytes) throw new Error('Edited TXT exceeds the 10 MiB limit.');
    return text;
  }
  M.spatialFormat = {parse:parse, detects:detects, empty:empty, isSpatial:isSpatial, typeInfo:typeInfo, memberError:memberError,validate:validate,serialize:serialize};
}(window));
