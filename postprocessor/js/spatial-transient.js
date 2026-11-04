(function(M){
  'use strict';const P=M.SpatialPost.prototype,S=M.spatialResults,{el,field,fmt,labels}=M.spatialUI;
  const make=P.makeControls,prepare=P.prepare,draw=P.drawChart,after=P.afterUpdate,play=P.play,conditions=P.conditions,select=P.select;
  P.makeControls=function(){make.call(this);if(this.analysis.type!=='transient')return;
    this.timeInput=field(this.settings,'time','Saved frame (zero-based)',0);this.timeInput.type='range';this.timeInput.min=0;this.timeInput.max=this.analysis.time.length-1;this.timeInput.step=1;
    this.frameNumber=field(this.settings,'frame','Frame number',0);this.frameNumber.min=0;this.frameNumber.max=this.analysis.time.length-1;this.frameNumber.step=1;
    this.timeInput.oninput=()=>this.setFrame(Number(this.timeInput.value));this.frameNumber.onchange=()=>this.setFrame(Number(this.frameNumber.value));
    const nav=el('div',undefined,this.settings);for(const delta of [-1,1]){const b=el('button',delta<0?'Previous frame':'Next frame',nav);b.onclick=()=>this.setFrame(Math.max(0,Math.min(this.analysis.time.length-1,this.index+delta*Number(this.stride.value))));}
    this.speed=field(this.settings,'speed','Playback speed',1,[[.25,'0.25×'],[.5,'0.5×'],[1,'1×'],[2,'2×'],[4,'4×']]);this.stride=field(this.settings,'stride','Displayed frame stride',1);this.stride.min=1;this.stride.step=1;
    this.playButton=el('button','Play',this.settings);this.playButton.id='spatial-result-play';this.playButton.onclick=()=>this.play();this.timeNote=el('p','',this.settings);
    this.prepareButton=el('button','Prepare scales',this.settings);this.prepareButton.onclick=()=>this.prepare();this.cancelPreparation=el('button','Cancel preparation',this.settings);this.cancelPreparation.onclick=()=>{this.token++;this.prepared=false;this.preparing=false;this.cancelPreparation.hidden=true;this.status.textContent='Preparation canceled. Data is preserved.';};this.cancelPreparation.hidden=true;
    this.historyKind=field(this.settings,'history','Chart','length',[['length','Along member length'],['displacements','Displacement'],['velocities','Velocity'],['accelerations','Acceleration'],['supportReactions','Support reaction'],['mpcForces','MPC force'],['mpcMultiplier','MPC multiplier'],['elementForce','Section force'],['spectrum','Displacement spectrum']]);
    this.historyNode=field(this.settings,'history-node','Chart node',1);this.historyComponent=field(this.settings,'history-component','Physical component',0,labels.slice(0,this.model.dofPerNode).map((s,i)=>[i,s]));this.historyElement=field(this.settings,'history-element','Chart member',1);this.historyMPC=field(this.settings,'history-mpc','MPC number',1);
    for(const input of [this.historyKind,this.historyNode,this.historyComponent,this.historyElement,this.historyMPC])input.onchange=()=>this.drawChart();
  };
  P.prepare=function(){if(this.analysis.type!=='transient')return prepare.call(this);this.prepareHistory().catch(e=>{if(!this.disposed){this.status.textContent=e.message;this.preparing=false;this.cancelPreparation.hidden=true;}});};
  P.prepareDiagrams=function(){if(this.analysis.type==='transient')this.prepare();else this.update();};
  P.prepareHistory=async function(){
    const token=++this.token;this.stop();this.prepared=false;this.preparing=true;this.cancelPreparation.hidden=false;this.playButton.disabled=true;
    let translation=0,rotation=0,force=0,stressMaximum=0,operations=0;const quantity=this.quantity.value,full=S.full(this.ds),count=this.analysis.displacements?this.analysis.time.length:0;
    for(let t=0;t<count;t++){const vector=S.vector(this.ds,'displacements',t);for(const e of this.model.elements){const u=S.local(e,vector);if(u&&!u.some(v=>v===null)){const b=S.bounds(e,u);translation=Math.max(translation,b.translation);rotation=Math.max(rotation,b.rotation);const r=S.recover(this.ds,e,vector);if(quantity!=='none'){const c=S.coefficient(e,r,quantity);for(const xi of S.extrema(c))force=Math.max(force,Math.abs(S.value(c,xi)));}if(S.sectionStress){const stress=S.sectionStress(e,r);if(stress)stressMaximum=Math.max(stressMaximum,Math.abs(stress.minimum.value),Math.abs(stress.maximum.value));}}
      if(++operations%500===0){this.status.textContent='Preparing scales: '+Math.floor(100*(t+1)/count)+'%';await new Promise(resolve=>setTimeout(resolve,0));if(token!==this.token)return;}
    }}
    if(token!==this.token)return;this.maximum=translation;this.maxRotation=rotation;this.diagramMaximum=force;this.stressMaximum=stressMaximum;this.prepared=true;this.preparing=false;this.cancelPreparation.hidden=true;this.playButton.disabled=!full||!this.scene;
    this.update();if(this.scene&&!this.initialFit){this.scene.fitModel();this.initialFit=true;}
  };
  P.setFrame=function(i){if(!Number.isInteger(i)||i<0||i>=this.analysis.time.length){this.status.textContent='Invalid frame number.';return;}this.stop();this.index=i;this.update();};
  P.afterUpdate=function(){after.call(this);if(this.timeNote){this.timeInput.value=this.index;this.frameNumber.value=this.index;this.timeNote.textContent='t = '+fmt(this.analysis.time[this.index])+' '+(this.ds.raw.metadata.units.time||'time units')+' · frame '+this.index+'. Entire history: 5 s playback at 1×. Export: '+this.analysis.sampling.exportedSampleCount+' / '+this.analysis.sampling.originalSampleCount+' samples, stride '+this.analysis.sampling.timeStride+'.';}};
  P.play=function(){if(this.analysis.type!=='transient')return play.call(this);if(this.animation!==undefined){this.stop();return;}if(!this.prepared||!S.full(this.ds)||!this.scene||this.scene.lost)return;
    const stride=Number(this.stride.value);if(!Number.isInteger(stride)||stride<1){this.status.textContent='Frame stride must be a positive integer.';return;}
    if(this.index===this.analysis.time.length-1)this.index=0;const started=performance.now(),first=this.index;this.playButton.textContent='Pause';
    const tick=now=>{const f=M.transientView.playbackFrame(this.analysis.time,first,now-started,Number(this.speed.value),stride,5000);if(this.index!==f.index){this.index=f.index;this.update();}if(f.complete)this.stop();else this.animation=requestAnimationFrame(tick);};this.animation=requestAnimationFrame(tick);
  };
  P.conditions=function(){const data=conditions.call(this);if(this.analysis.type!=='transient')return data;const a=this.analysis,t=a.time[this.index],dt=a.timeStep||(a.time[a.time.length-1]-a.time[0])/(a.sampling.originalSampleCount-1);
    data.loads.forEach(r=>{if(r.kind)return;const factor=r.type===11?Math.sin(2*Math.PI*r.frequency*t):r.type===13?1:Math.abs(t-dt)<=Number.EPSILON*Math.max(Number.MIN_VALUE,Math.abs(t),Math.abs(dt))*8?1:0;r.components=r.components.map(v=>v*factor);});return data;
  };
  P.select=function(s){if(this.historyNode&&s){if(s.kind==='node')this.historyNode.value=s.index+1;if(s.kind==='element')this.historyElement.value=s.index+1;if(s.kind==='mpcs')this.historyMPC.value=s.index+1;}select.call(this,s);};
  P.drawChart=function(){
    if(!this.historyKind||this.historyKind.value==='length')return draw.call(this);
    const kind=this.historyKind.value,a=this.analysis,d=this.model.dofPerNode,node=Number(this.historyNode.value),component=Number(this.historyComponent.value),row=(node-1)*d+component,eid=Number(this.historyElement.value),mid=Number(this.historyMPC.value)-1,xi=Number(this.position.value),quantity=this.quantity.value;
    const key=JSON.stringify([kind,node,component,eid,mid,xi,quantity]);let series=this.historyCache&&this.historyCache.key===key?this.historyCache.series:null;
    if(!series){let y=null,x=a.time,caption=kind+' · '+labels[component]+'('+node+')';
      if(['displacements','velocities','accelerations'].includes(kind)&&Number.isInteger(node)&&node>=1&&node<=this.model.nodes.length){const values=a[kind],terms=this.ds.physicalRows[row];if(values&&terms.length&&terms.every(([id])=>this.ds.systemRows.has(id)))y=a.time.map((_,j)=>terms.reduce((sum,[id,c])=>sum+c*values[this.ds.systemRows.get(id)][j],0));}
      else if(['supportReactions','mpcForces'].includes(kind)&&Number.isInteger(node)&&node>=1&&node<=this.model.nodes.length){y=a.time.map((_,j)=>S.reaction(this.ds,kind,row,j));if(y.some(v=>v===null))y=null;}
      else if(kind==='spectrum'&&this.ds.physicalRows[row]){const spectrum=S.spectrum(this.ds,row);if(spectrum){x=spectrum.x;y=spectrum.y;}}
      else if(kind==='mpcMultiplier'&&a.mpcMultipliers&&a.mpcMultipliers[mid]){y=a.mpcMultipliers[mid];caption='λ · MPC '+(mid+1);}
      else if(kind==='elementForce'&&this.model.elements[eid-1]&&quantity!=='none'&&Number.isFinite(xi)&&xi>=0&&xi<=1&&a.displacements){const e=this.model.elements[eid-1];if(e.globalDOFs.every(id=>this.ds.systemRows.has(id-1))){y=a.time.map((_,j)=>{const vector={};e.globalDOFs.forEach(id=>{vector[id-1]=a.displacements[this.ds.systemRows.get(id-1)][j];});return S.value(S.coefficient(e,S.recover(this.ds,e,vector),quantity),xi);});caption=quantity+' · E'+eid+' · x/L='+xi;}}
      if(!y||y.some(v=>!Number.isFinite(v))){this.chart.hidden=true;this.series=null;this.legend.textContent='Chart unavailable: select an existing object and an exported component. For member forces, also select a diagram.';return;}
      series={x,y,xLabel:kind==='spectrum'?'Frequency, Hz':'Analysis time',yLabel:caption};this.historyCache={key,series};
    }
    this.series=series;this.chart.hidden=false;M.charts.render(this.chart,series,kind==='spectrum'?NaN:a.time[this.index],kind==='spectrum'?undefined:i=>this.setFrame(i));
  };
}(window.MKEFPost));
