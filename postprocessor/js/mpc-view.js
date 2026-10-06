(function (M) {
  'use strict';
  var V=M.mpcView={};
  function clear(el) { while (el.firstChild) el.removeChild(el.firstChild); }
  function format(v) { return Number(v.toPrecision(7)).toString(); }
  V.mount=function (dataset,renderer,panel) {
    clear(panel); panel.hidden=dataset.raw.version !== 3;
    if (panel.hidden) return null;
    var a=dataset.raw.analysis, model=dataset.raw.model, state={index:0};
    var title=document.createElement('h2'); title.textContent='MPC constraints and forces'; panel.appendChild(title);
    var note=document.createElement('p'); note.className='control-note';
    note.textContent='r = support reactions + Cᵀλ; dependent coefficient is +1. '+(a.type==='modal'?'Values use the stored eigenvector amplitude, independent of the visual deformation scale.':''); panel.appendChild(note);
    function toggle(id,label,layer) {
      var host=document.createElement('label'), input=document.createElement('input'); input.type='checkbox'; input.checked=true; input.id=id;
      host.appendChild(input); host.appendChild(document.createTextNode(' '+label)); panel.appendChild(host);
      input.addEventListener('change',function () { layer.style.display=input.checked?'':'none'; });
      return input;
    }
    toggle('show-mpc-links','MPC links',renderer.layers.mpcs);
    toggle('show-spc-forces','Support reactions',renderer.layers['spc-forces']);
    toggle('show-mpc-forces','MPC forces',renderer.layers['mpc-forces']);
    var select=document.createElement('select'); select.id='mpc-equation'; select.setAttribute('aria-label','MPC equation');
    model.mpcs.forEach(function (m,i) { var option=document.createElement('option'); option.value=String(i); option.textContent='MPC '+m.id+': '+window.MKEFMPC.equation(m); select.appendChild(option); }); panel.appendChild(select);
    var value=document.createElement('p'); value.id='mpc-multiplier'; panel.appendChild(value);
    var table=document.createElement('table'), head=document.createElement('tr');
    ['Node','DOF','Support reaction','MPC force'].forEach(function (label) { var th=document.createElement('th'); th.textContent=label; head.appendChild(th); }); table.appendChild(head);
    var body=document.createElement('tbody'); table.appendChild(body); var scroll=document.createElement('div'); scroll.style.overflow='auto'; scroll.style.maxHeight='260px'; scroll.appendChild(table); panel.appendChild(scroll);
    var chart=M.svgElement('svg',{id:'mpc-history-chart',role:'img','aria-label':'MPC multiplier history'}); panel.appendChild(chart);
    model.mpcs.forEach(function (m,i) {
      var p=dataset.nodesById.get(m.depNode), group=M.svgElement('g',{'data-mpc-id':m.id,stroke:'#7c3aed','stroke-width':renderer.size*.002,fill:'none'});
      var tooltip=M.svgElement('title'); tooltip.textContent='MPC '+m.id+': '+window.MKEFMPC.equation(m); group.appendChild(tooltip);
      m.masters.forEach(function (term) { var q=dataset.nodesById.get(term.node); group.appendChild(M.svgElement('line',{x1:p.x,y1:-p.y,x2:q.x,y2:-q.y,'stroke-dasharray':renderer.size*.025+' '+renderer.size*.012})); });
      group.appendChild(M.svgElement('circle',{cx:p.x,cy:-p.y,r:renderer.size*.02})); renderer.layers.mpcs.appendChild(group);
    });
    function vector(name,index) {
      if (!a[name]) return null;
      if (a.type==='static') return a[name];
      if (a.type==='modal') return a[name].map(function (row) { return row[index]; });
      var values=new Array(model.nodes.length*model.dofPerNode).fill(null);
      a.globalDOFIds.forEach(function (id,row) { values[id-1]=a[name][row][index]; }); return values;
    }
    function drawForces(layer,values,color,label) {
      clear(layer); if (!values) return;
      var max=Math.max.apply(null,values.map(function (v) { return Math.abs(v||0); }).concat([1e-30])), size=renderer.size*.12;
      model.nodes.forEach(function (node,row) {
        var ids=model.dofMap[row], fx=values[ids[0]-1]||0, fy=values[ids[1]-1]||0, moment=ids.length>2?(values[ids[2]-1]||0):0;
        var magnitude=Math.hypot(fx,fy), path='';
        if (magnitude>max*1e-10) {
          var dx=fx/magnitude*size, dy=-fy/magnitude*size, x=node.x,y=-node.y;
          path='M '+(x-dx)+' '+(y-dy)+' L '+x+' '+y+' M '+(x-dx*.25-dy*.12)+' '+(y-dy*.25+dx*.12)+' L '+x+' '+y+' L '+(x-dx*.25+dy*.12)+' '+(y-dy*.25-dx*.12);
        }
        if (Math.abs(moment)>max*1e-10) {
          var radius=size*.45, sx=node.x-radius*.7, sy=-node.y-radius*.7, direction=moment>0?-1:1;
          path+=' M '+(node.x+radius)+' '+(-node.y)+' A '+radius+' '+radius+' 0 1 '+(moment>0?0:1)+' '+sx+' '+sy+' l '+(radius*.3)+' '+(direction*radius*.05)+' M '+sx+' '+sy+' l '+(direction*radius*.05)+' '+(radius*.3);
        }
        if (path) { var symbol=M.svgElement('path',{d:path,stroke:color,fill:'none','stroke-width':2,'vector-effect':'non-scaling-stroke','data-node-id':node.id}); var tip=M.svgElement('title'); tip.textContent=label+' node '+node.id+': Fx='+format(fx)+', Fy='+format(fy)+', Mz='+format(moment); symbol.appendChild(tip); layer.appendChild(symbol); }
      });
    }
    state.update=function (index) {
      state.index=index||0;
      var equation=Number(select.value), m=model.mpcs[equation], unit=m.depDOF===3?dataset.raw.metadata.units.moment:dataset.raw.metadata.units.force;
      var lambda=a.mpcMultipliers ? (a.type==='static'?a.mpcMultipliers[equation]:a.mpcMultipliers[equation][state.index]) : null;
      value.textContent='MPC '+m.id+' · λ = '+(lambda===null?'not exported':format(lambda)+(unit?' '+unit:''));
      var supports=vector('supportReactions',state.index), forces=vector('mpcForces',state.index);
      drawForces(renderer.layers['spc-forces'],supports,'#c56a00','Support reaction');
      drawForces(renderer.layers['mpc-forces'],forces,'#7c3aed','MPC force');
      clear(body);
      var relevant=new Set(); model.mpcs.forEach(function (q) { relevant.add(q.depNode+':'+q.depDOF); q.masters.forEach(function (term) { relevant.add(term.node+':'+term.dof); }); });
      model.supports.forEach(function (s) { window.MKEFSupports.dofs(s,model.dofPerNode).forEach(function (d) { relevant.add(s.nodeId+':'+d); }); });
      relevant.forEach(function (key) {
        var pair=key.split(':').map(Number), row=dataset.nodeIndexById.get(pair[0]), id=model.dofMap[row][pair[1]-1], tr=document.createElement('tr');
        [pair[0],model.dofLabels[pair[1]-1],supports&&supports[id-1]!==null?format(supports[id-1]):'not exported',forces&&forces[id-1]!==null?format(forces[id-1]):'not exported'].forEach(function (text) { var td=document.createElement('td');td.textContent=String(text);tr.appendChild(td); }); body.appendChild(tr);
      });
      chart.hidden=a.type!=='transient'||!a.mpcMultipliers;
      chart.style.display=chart.hidden?'none':'block';
      if (!chart.hidden) M.charts.render(chart,{x:a.time,y:a.mpcMultipliers[equation],xLabel:'Time',yLabel:'MPC '+m.id+' · λ'+(unit?' ('+unit+')':''),quantity:'mpcMultipliers'},a.time[state.index]);
    };
    select.addEventListener('change',function () { state.update(state.index); });
    state.update(0); return state;
  };
}(window.MKEFPost));
