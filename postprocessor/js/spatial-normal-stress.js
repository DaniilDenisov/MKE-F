(function(M,V){
  'use strict';
  const P=M.SpatialPost.prototype,make=P.makeControls,decorate=P.decorateFrame,after=P.afterUpdate,dispose=P.dispose,samplePositions=P.memberSamplePositions;
  P.makeControls=function(){
    make.call(this);if(this.analysis.type!=='static')return;
    this.normalStress=new V.Session(this.ds,this.settings,this.panel);
    this.normalStressColor=M.spatialUI.field(this.settings,'color','Deformed member color','none',[['none','Plain'],['normalStress','Maximum absolute normal stress']]);
    this.settings.insertBefore(this.normalStressColor.parentElement,this.quantity.parentElement);this.settings.insertBefore(this.normalStress.controls,this.quantity.parentElement);
    this.normalStressColor.onchange=()=>{this.normalStress.setActive(this.normalStressColor.value==='normalStress');this.update();};
  };
  P.decorateFrame=function(frame){
    decorate?.call(this,frame);if(!this.normalStress?.active)return;
    const colors=[];frame.deformed.owners.forEach((owner,i)=>{const e=this.model.elements[owner],r=this.normalStress.results.get(e.id);
      for(const t of frame.deformed.parameters.slice(2*i,2*i+2)){const v=V.at(r,t);colors.push(...(v===null?[135/255,151/255,167/255]:V.rgb(v,this.normalStress.maximum)));}});
    frame.deformed.colors=colors;
  };
  P.memberSamplePositions=function(e,n){
    const ts=samplePositions?samplePositions.call(this,e,n):Array.from({length:n},(_,i)=>i/(n-1)),r=this.normalStress?.active&&this.normalStress.results.get(e.id);
    return r?.status==='ok'?[...new Set(ts.concat(r.positions))].sort((a,b)=>a-b):ts;
  };
  P.afterUpdate=function(){after?.call(this);if(this.normalStress)this.normalStress.select(this.selection?.kind==='element'?this.model.elements[this.selection.index].id:null);};
  P.dispose=function(){this.normalStress?.dispose();dispose.call(this);};
}(window.MKEFPost,window.MKEFPost.normalStressResults));
