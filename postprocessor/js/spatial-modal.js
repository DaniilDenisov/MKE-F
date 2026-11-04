(function(M){
  'use strict';const P=M.SpatialPost.prototype,{el,field,fmt}=M.spatialUI,make=P.makeControls,after=P.afterUpdate;
  P.makeControls=function(){make.call(this);if(this.analysis.type!=='modal')return;
    this.modeInput=field(this.settings,'mode','Mode number',1);this.modeInput.min=1;this.modeInput.max=this.analysis.frequenciesHz.length;this.modeInput.step=1;
    this.modeInput.onchange=()=>this.changeMode(Number(this.modeInput.value)-1);
    this.phaseInput=field(this.settings,'phase','Phase (multiplier)',1);this.phaseInput.type='range';this.phaseInput.min=-1;this.phaseInput.max=1;this.phaseInput.step=.01;
    this.phaseInput.oninput=()=>{this.stop();this.phase=Number(this.phaseInput.value);this.update();};
    this.speed=field(this.settings,'speed','Playback speed',1,[[.25,'0.25×'],[.5,'0.5×'],[1,'1×'],[2,'2×'],[4,'4×']]);
    this.playButton=el('button','Play',this.settings);this.playButton.id='spatial-result-play';this.playButton.onclick=()=>this.play();this.modeNote=el('p','',this.settings);
    const details=el('details',undefined,this.settings);el('summary','Frequency list',details);const list=el('div',undefined,details);let page=0;
    const draw=()=>{list.replaceChildren();const prev=el('button','Previous',list),next=el('button','Next',list);prev.disabled=page===0;next.disabled=(page+1)*200>=this.analysis.frequenciesHz.length;prev.onclick=()=>{page--;draw();};next.onclick=()=>{page++;draw();};const body=el('tbody',undefined,el('table',undefined,list));this.analysis.frequenciesHz.slice(page*200,(page+1)*200).forEach((f,j)=>{const i=page*200+j,tr=el('tr',undefined,body);el('td',i+1,tr);el('td',fmt(f)+' Hz',tr);tr.tabIndex=0;tr.onclick=()=>this.changeMode(i);tr.onkeydown=e=>{if(e.key==='Enter')tr.click();};});};draw();
  };
  P.changeMode=function(index){if(!Number.isInteger(index)||index<0||index>=this.analysis.frequenciesHz.length){this.status.textContent='Invalid mode number.';return;}this.stop();this.index=index;this.phase=1;this.phaseInput.value=1;this.modeInput.value=index+1;this.prepare();};
  P.afterUpdate=function(){after.call(this);if(this.modeNote)this.modeNote.textContent='f = '+fmt(this.analysis.frequenciesHz[this.index])+' Hz · phase '+fmt(this.phase)+' · 2 s playback cycle at 1×. Forces and reactions refer to the saved mode amplitude.';};
  P.stop=function(){if(this.animation!==undefined){cancelAnimationFrame(this.animation);this.animation=undefined;}if(this.playButton)this.playButton.textContent='Play';};
  P.play=function(){if(this.animation!==undefined){this.stop();return;}if(!this.scene||this.scene.lost)return;this.playButton.textContent='Pause';const started=performance.now(),angle=Math.acos(Math.max(-1,Math.min(1,this.phase)));
    const tick=now=>{this.phase=Math.cos(angle+(now-started)*2*Math.PI/2000*Number(this.speed.value));this.phaseInput.value=this.phase;this.update();this.animation=requestAnimationFrame(tick);};this.animation=requestAnimationFrame(tick);
  };
}(window.MKEFPost));
