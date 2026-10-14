(function (M) {
  'use strict';
  var S=window.MKEFSections;
  function node(tag,text) {var n=document.createElement(tag);if(text!==undefined)n.textContent=text;return n;}
  function nextId(model) {return Math.max.apply(null,(model.sections || []).map(function (s) {return s.id;}).concat([0]))+1;}
  function defaults(type) {return type==='circle' ? {d:.1} : type==='rectangle' ? {b:.1,h:.2} : {b:.1,h:.2,tw:.008,tf:.012};}
  M.renderSections=function (model,commit,selection) {
    var host=document.getElementById('section-catalog'); host.replaceChildren();
    var add=node('button','Add section');add.id='add-section';add.type='button';host.appendChild(add);
    add.onclick=function () {commit(function () {model.sections=model.sections || [];var id=nextId(model);model.sections.push({id:id,name:'Section '+id,type:'rectangle',dimensions:defaults('rectangle')});});};
    (model.sections || []).forEach(function (s) {
      var card=node('article');card.className='section-card';card.dataset.sectionId=s.id;host.appendChild(card);
      var controls=node('div');card.appendChild(controls);
      function input(label,value,type,change) {var l=node('label',label),i=node('input');i.type=type;i.value=value;if(type==='number')i.step='any';i.setAttribute('aria-label',label+' · section '+s.id);l.appendChild(i);controls.appendChild(l);i.onchange=function () {commit(function () {change(type==='number' ? (i.value.trim()==='' ? NaN : Number(i.value)) : i.value);S.sync(model);});};return i;}
      input('Name',s.name,'text',function (v) {s.name=v;});
      var label=node('label','Shape'),select=node('select');Object.keys(S.types).forEach(function (type) {var o=node('option',S.types[type]);o.value=type;select.appendChild(o);});select.value=s.type;select.setAttribute('aria-label','Shape · section '+s.id);label.appendChild(select);controls.appendChild(label);
      select.onchange=function () {commit(function () {s.type=select.value;s.dimensions=defaults(s.type);S.sync(model);});};
      S.fields(s.type).forEach(function (key) {input(key,s.dimensions[key],'number',function (v) {s.dimensions[key]=v;});});
      var p=S.properties(s),info=node('p','A='+p.area.toPrecision(6)+' · Iy='+p.Iy.toPrecision(6)+' · Iz='+p.Iz.toPrecision(6)+' · centroid (y,z)=('+p.centroidY.toPrecision(5)+', '+p.centroidZ.toPrecision(5)+')');info.className='control-note';controls.appendChild(info);
      var assigned=[];model.elements.forEach(function (e,i) {if(e.sectionId===s.id)assigned.push('E'+(i+1));});
      var assignments=node('p',assigned.length ? 'Assigned to: '+assigned.join(', ') : 'Not assigned to any element.');assignments.className='control-note section-assignments';controls.appendChild(assignments);
      var assign=node('button','Assign to selected element'),copy=node('button','Copy'),del=node('button','Delete');assign.type=copy.type=del.type='button';
      var selected=selection && selection.kind==='element' && model.elements[selection.index];
      assign.disabled=!selected;assign.title=selected ? 'Assign to E'+(selection.index+1) : 'Select an element on the diagram or in the Elements table first.';
      assign.onclick=function () {if(!selected)return;commit(function () {selected.sectionId=s.id;selected.sectionAngle=selected.sectionAngle || 0;S.sync(model);});};
      controls.appendChild(assign);controls.appendChild(copy);controls.appendChild(del);
      copy.onclick=function () {commit(function () {var c=JSON.parse(JSON.stringify(s));c.id=nextId(model);c.name+=' copy';model.sections.push(c);});};
      del.onclick=function () {commit(function () {var used=[];model.elements.forEach(function (e,i) {if(e.sectionId===s.id)used.push(i+1);});if(used.length)throw new Error('Section is assigned to elements '+used.join(', ')+'. Remove those assignments first.');model.sections=model.sections.filter(function (a) {return a.id!==s.id;});});};
      var svg=document.createElementNS('http://www.w3.org/2000/svg','svg');svg.id='section-preview-'+s.id;svg.classList.add('section-preview');card.appendChild(svg);window.MKEFSectionDrawing.preview(svg,s,0);
    });
  };
}(window.MKEFPre));
