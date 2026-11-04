(function(M){
  'use strict';const P=M.SpatialPost.prototype,S=M.spatialResults,{el,field,check,fmt,labels}=M.spatialUI;
  const make=P.makeControls,conditions=P.conditions,details=P.showDetails;
  P.makeControls=function(){make.call(this);this.quantity=field(this.settings,'quantity','Diagram','none',[['none','No diagram'],['N','N'],...(this.model.dofPerNode===6?['Vy','Vz','Mx','My','Mz'].map(x=>[x,x]):[]),['stress','N/A']]);this.quantity.onchange=()=>{if(this.prepareDiagrams)this.prepareDiagrams();else this.update();};
    this.diagramGain=field(this.settings,'diagram-gain','Diagram display multiplier',1);this.diagramGain.min='.01';this.diagramGain.onchange=()=>{const v=Number(this.diagramGain.value);this.diagramGain.setCustomValidity(Number.isFinite(v)&&v>0?'':'Use a positive finite multiplier.');if(this.diagramGain.reportValidity())this.update();};
    this.reactionInputs={};for(const [key,title]of [['supportReactions','Support reactions (purple)'],['mpcForces','MPC forces (ochre)']]){this.reactionInputs[key]=check(this.settings,key,title,false);this.reactionInputs[key].onchange=()=>this.update();}
    this.position=field(this.settings,'position','Position x/L',.5);this.position.min=0;this.position.max=1;this.position.onchange=()=>this.update();
    this.legend=el('p','',this.panel);this.legend.id='spatial-result-legend';this.chart= M.svgElement('svg',{id:'spatial-result-chart',role:'img','aria-label':'Chart along member length'});this.panel.insertBefore(this.chart,this.tables);this.chart.style.width='100%';this.chart.style.height='250px';this.chart.hidden=true;
  };
  P.conditions=function(){const data=conditions.call(this),ds=this.ds;for(const kind of ['supportReactions','mpcForces'])if(this.reactionInputs&&this.reactionInputs[kind].checked)this.model.nodes.forEach((node,i)=>{const values=Array.from({length:6},(_,j)=>j<this.model.dofPerNode?S.reaction(ds,kind,i*this.model.dofPerNode+j,this.index,this.phase):0);if(values.some(v=>v!==null&&Math.abs(v)>0))data.loads.push({node:node.id,kind,ownerIndex:i,type:10,components:values.map(v=>v||0)});});return data;};
  P.decorateFrame=function(frame){
    const quantity=this.quantity.value,rows=[],positions=[],owners=[],colors=[],fill=[],fillOwners=[],fillColors=[],guides=[],guideOwners=[],annotations=[];let maximum=0;
    if(quantity!=='none'&&S.full(this.ds))for(const e of this.model.elements){const result=S.recover(this.ds,e,this.vector,this.acceleration);if(!result)continue;const c=S.coefficient(e,result,quantity);if(!c)continue;for(const t of S.extrema(c))maximum=Math.max(maximum,Math.abs(S.value(c,t)));rows.push({e,c});}
    if(this.diagramMaximum!==undefined)maximum=this.diagramMaximum;
    if(this.analysis.type==='modal'&&this.phase)maximum/=Math.abs(this.phase);
    const factor=Number(this.diagramGain.value),gain=maximum?.15*this.size/maximum*(Number.isFinite(factor)&&factor>0?factor:1):0,samples=Math.max(5,Math.min(Number(this.samples.value),Math.floor(1000000/(2*this.model.elements.length))+1));
    const selected=this.selection&&this.selection.kind==='element'?this.selection.index+1:null;
    if(selected)rows.sort((a,b)=>(a.e.id===selected?-1:b.e.id===selected?1:0));
    for(const {e,c}of rows){const a=this.model.nodes[e.nodeIds[0]-1],b=this.model.nodes[e.nodeIds[1]-1],normal=e.localBasis[['Vz','My'].includes(quantity)?2:1];let last=null;
      const ts=new Set(Array.from({length:samples},(_,i)=>i/(samples-1))),extrema=S.extrema(c);if(e.id===selected)extrema.forEach(t=>ts.add(t));
      const color=v=>v>=0?[.75,.15,.15]:[.15,.3,.8];
      for(const t of [...ts].sort((a,b)=>a-b)){const value=S.value(c,t),p=[a.x+(b.x-a.x)*t,a.y+(b.y-a.y)*t,a.z+(b.z-a.z)*t],offset=p.map((v,i)=>v+gain*value*normal[i]);
        if(last){positions.push(...last.point,...offset);owners.push(e.id-1);colors.push(...color(last.value),...color(value));
          if(fill.length+18<=3000000){fill.push(...last.base,...last.point,...offset,...last.base,...offset,...p);fillOwners.push(e.id-1,e.id-1);fillColors.push(...color(last.value),...color(last.value),...color(value),...color(last.value),...color(value),...color(value));}}
        if((e.id===selected||rows.length<500)&&(t===0||t===1)){guides.push(...p,...offset);guideOwners.push(e.id-1);}
        if(e.id===selected&&(t===0||t===1||extrema.includes(t)))annotations.push({position:offset,text:quantity+' '+fmt(Math.abs(value)<maximum*1e-10?0:value)+' '+this.unit(quantity),color:value>=0?'#a32626':'#235cab'});
        last={point:offset,base:p,value};}
    }
    frame.diagrams={positions,owners,colors,visible:quantity!=='none'};this.legend.textContent=quantity==='none'?'':quantity+' · |max| = '+fmt(maximum)+' '+this.unit(quantity)+' · diagram scale ×'+fmt(gain)+' · + red / − blue · diagram in original local axes';
    frame.diagramFill={positions:fill,owners:fillOwners,colors:fillColors,kind:'triangles',opacity:.22,visible:quantity!=='none'};
    frame.diagramGuides={positions:guides,owners:guideOwners,color:'#668099',visible:quantity!=='none'};frame.labels=annotations.slice(0,8);
    if(this.scene){this.scene.setConditions(this.conditions());this.scene.setConditionLayers({supportReactions:this.reactionInputs.supportReactions.checked,mpcForces:this.reactionInputs.mpcForces.checked});}
  };
  P.unit=function(quantity){const units=this.ds.raw.metadata.units;return quantity==='stress'?units.stress||'':quantity.startsWith('M')?units.moment||'':units.force||'';};
  P.showDetails=function(){details.call(this);const s=this.selection;if(!s)return;const lines=[];
    if(s.kind==='node'||s.kind==='supportReactions'||s.kind==='mpcForces'){const i=s.index;for(const field of ['supportReactions','mpcForces']){const values=labels.slice(0,this.model.dofPerNode).map((name,j)=>name+'='+fmt(S.reaction(this.ds,field,i*this.model.dofPerNode+j,this.index,this.phase)));lines.push((field==='supportReactions'?'Support reactions':'MPC forces')+': '+values.join(', '));}}
    if(s.kind==='element'){const e=this.model.elements[s.index],r=S.recover(this.ds,e,this.vector,this.acceleration),x=Number(this.position.value);if(r&&Number.isFinite(x)&&x>=0&&x<=1)for(const quantity of ['N',...(this.model.dofPerNode===6?['Vy','Vz','Mx','My','Mz']:[]),'stress']){const c=S.coefficient(e,r,quantity),ext=S.extrema(c).map(t=>({t,value:S.value(c,t)})),min=ext.reduce((a,b)=>a.value<b.value?a:b),max=ext.reduce((a,b)=>a.value>b.value?a:b);lines.push(quantity+'('+x+')='+fmt(S.value(c,x))+'; min='+fmt(min.value)+' @'+fmt(min.t)+'; max='+fmt(max.value)+' @'+fmt(max.t));}}
    if(s.kind==='mpcs'){const r=this.model.mpcs[s.index],a=this.analysis,values=a.mpcMultipliers;lines.push(labels[r.depDOF-1]+'('+r.depNode+') = '+r.masters.map(t=>t[2]+'·'+labels[t[1]-1]+'('+t[0]+')').join(' + '));if(values){const lambda=a.type==='static'?values[s.index]:values[s.index][this.index]*this.phase;lines.push('λ = '+fmt(lambda));const names=['Fx','Fy','Fz','Mx','My','Mz'];lines.push('Physical forces of this MPC: '+[[r.depNode,r.depDOF,1],...r.masters.map(t=>[t[0],t[1],-t[2]])].map(t=>names[t[1]-1]+'('+t[0]+')='+fmt(t[2]*lambda)).join(', '));}else lines.push('The multiplier and forces for this MPC were not exported.');}
    if(this.analysis.type==='modal')lines.push('Modal quantities: saved amplitude × current phase.');if(this.analysis.type!=='static'&&s.kind==='element')lines.push('Elastic moments at a released end may be balanced by inertial moments.');this.details.textContent+='\n'+lines.join('\n');
  };
  P.afterUpdate=function(){this.drawChart();};
  P.drawChart=function(){
    const s=this.selection,quantity=this.quantity.value;if(!s||s.kind!=='element'||quantity==='none'){this.chart.hidden=true;this.series=null;return;}
    const e=this.model.elements[s.index],r=S.recover(this.ds,e,this.vector,this.acceleration),c=r&&S.coefficient(e,r,quantity);if(!c){this.chart.hidden=true;this.series=null;return;}const x=[...new Set([...Array.from({length:201},(_,i)=>i/200),...S.extrema(c)])].sort((a,b)=>a-b);this.series={x,y:x.map(t=>S.value(c,t)),xLabel:'x/L',yLabel:quantity+' '+this.unit(quantity)};this.chart.hidden=false;M.charts.render(this.chart,this.series);
  };
}(window.MKEFPost));
