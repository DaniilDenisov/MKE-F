(function(M,V){
  'use strict';
  const el=(tag,text,parent)=>{const item=document.createElement(tag);if(text!==undefined)item.textContent=text;if(parent)parent.append(item);return item;};
  const fmt=v=>Number.isFinite(v)?Number(v.toPrecision(7)).toString():'—';
  const unit=(ds,key)=>ds.raw.metadata?.units?.[key]||key+' units not specified';
  V.stops=[{offset:0,color:'#440154'},{offset:.25,color:'#3b528b'},{offset:.5,color:'#21918c'},{offset:.75,color:'#5ec962'},{offset:1,color:'#fde725'}];
  V.rgb=function(value,maximum){
    const t=maximum?Math.max(0,Math.min(1,value/maximum)):0,i=Math.min(3,Math.floor(t*4)),u=t*4-i;
    const rgb=c=>[1,3,5].map(p=>parseInt(c.slice(p,p+2),16)/255),a=rgb(V.stops[i].color),b=rgb(V.stops[i+1].color);
    return a.map((v,j)=>v+(b[j]-v)*u);
  };
  V.color=(value,maximum)=>'#'+V.rgb(value,maximum).map(c=>Math.round(c*255).toString(16).padStart(2,'0')).join('');
  V.csv=function(session){
    const rows=[['element','xi','sigmaMin','sigmaMax','sigmaAbsMax','sigmaX','eta','y','z','stressUnit','lengthUnit','status','reason']];
    for(const e of session.ds.raw.model.elements){const r=session.results.get(e.id),units=[unit(session.ds,'stress'),unit(session.ds,'length')];
      if(r.status!=='ok')rows.push([e.id,...Array(8).fill(''),...units,r.status,r.reason]);
      else for(const p of r.samples)rows.push([e.id,p.xi,p.sigmaMin,p.sigmaMax,p.sigmaAbsMax,p.sigmaX,p.eta??'',p.y??'',p.z??'',...units,r.status,'']);
    }
    const escape=v=>'"'+String(v).replace(/"/g,'""')+'"';return rows.map(r=>r.map(escape).join(',')).join('\r\n');
  };
  class Session {
    constructor(ds,settingsHost,panelHost){
      this.ds=ds;this.results=new Map();this.maximum=0;this.active=false;this.started=false;this.selected=null;
      this.controls=el('section',undefined,settingsHost);this.controls.className='normal-stress-settings';this.controls.hidden=true;
      el('h3','Maximum absolute normal stress',this.controls);
      el('p','Normal stresses from axial force and bending only. Transverse shear and torsion are excluded. Pure torsion can give zero normal stress; zero does not mean a stress-free or verified safe member.',this.controls);
      this.exportButton=el('button','Export normal stress CSV',this.controls);this.exportButton.id='normal-stress-export-csv';this.exportButton.onclick=()=>{
        const url=URL.createObjectURL(new Blob(['\uFEFF'+V.csv(this)],{type:'text/csv;charset=utf-8'})),a=el('a');a.href=url;a.download='MKE-F-normal-stress.csv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
      };
      this.panel=el('section',undefined,panelHost);this.panel.id='normal-stress-panel';this.panel.hidden=true;panelHost.prepend(this.panel);
      el('h2','Maximum absolute normal stress',this.panel);this.legend=el('div',undefined,this.panel);this.legend.id='normal-stress-legend';
      this.summary=el('p','',this.panel);this.summary.id='normal-stress-status';this.exclusions=el('details',undefined,this.panel);el('summary','Unavailable members',this.exclusions);this.exclusionList=el('pre','',this.exclusions);
      this.details=el('pre','',this.panel);this.details.id='normal-stress-details';
      this.chart=M.svgElement('svg',{id:'normal-stress-chart',role:'img','aria-label':'Maximum absolute normal stress along the member',tabindex:0});this.panel.append(this.chart);
      this.sample=el('pre','',this.panel);this.sample.id='normal-stress-sample';
    }
    setActive(value){
      if(this.active===value)return;this.active=value;this.controls.hidden=this.panel.hidden=!value;
      if(value&&!this.started){Object.assign(this,V.collect(this.ds));this.started=true;}
      if(value)this.render();
    }
    select(id){if(this.selected===id)return;this.selected=id;if(this.active)this.renderSelection();}
    render(){
      this.legend.replaceChildren();const total=this.ds.raw.model.elements.length,excluded=[...this.results].filter(([,r])=>r.status!=='ok');
      this.summary.textContent='Supported: '+(this.results.size-excluded.length)+' / '+total+' · unavailable: '+excluded.length+'. Scale uses the whole model; grey means unavailable. Member values at shared nodes are not averaged.';
      this.exclusions.hidden=!excluded.length;this.exclusionList.textContent=excluded.map(([id,r])=>'E'+id+' · '+r.reason).join('\n');
      if(this.results.size>excluded.length){el('p','max |σx| · 0 … '+fmt(this.maximum)+' · '+unit(this.ds,'stress'),this.legend);const bar=el('div',undefined,this.legend);bar.className='legend-bar';bar.style.background='linear-gradient(90deg,'+V.stops.map(s=>s.color+' '+s.offset*100+'%').join(',')+')';}
      this.renderSelection();
    }
    renderSelection(){
      const r=this.results.get(this.selected);this.chart.hidden=!r||r.status!=='ok';this.sample.textContent='';
      if(!r){this.details.textContent='Select a member to inspect stresses and the maximum location.';M.charts.render(this.chart,null);return;}
      if(r.status!=='ok'){this.details.textContent='E'+this.selected+' · Normal stress unavailable: '+r.reason;M.charts.render(this.chart,null);return;}
      this.details.textContent='E'+this.selected+' · extrema along member\nσmin='+fmt(r.sigmaMin)+'; σmax='+fmt(r.sigmaMax)+'\nmax |σx|='+fmt(r.sigmaAbsMax)+' · '+unit(this.ds,'stress')+'\n'+this.pointText(r.peak);
      M.charts.render(this.chart,{x:r.samples.map(p=>p.xi),y:r.samples.map(p=>p.sigmaAbsMax),xLabel:'x/L',yLabel:'max |σx| · '+unit(this.ds,'stress'),compact:true},r.peak.xi,index=>{
        const p=r.samples[index];this.sample.textContent='Section at x/L='+fmt(p.xi)+'\nσmin='+fmt(p.sigmaMin)+'; σmax='+fmt(p.sigmaMax)+'; max |σx|='+fmt(p.sigmaAbsMax)+' · '+unit(this.ds,'stress')+'\n'+this.pointText(p);
      });
    }
    pointText(p){return 'Governing '+(p.sigmaX<0?'compression':p.sigmaX>0?'tension':'zero normal stress')+': σx='+fmt(p.sigmaX)+' at x/L='+fmt(p.xi)+'\n'+(p.eta!==null?'η='+fmt(p.eta)+' · '+unit(this.ds,'length'):p.y!==null?'y='+fmt(p.y)+', z='+fmt(p.z)+' · '+unit(this.ds,'length'):'Uniform axial stress; no unique critical point in the section.');}
    dispose(){V.clear(this.ds);this.results.clear();this.controls.remove();this.panel.remove();}
  }
  V.Session=Session;
  V.planar=function(ds){
    const session=ds.normalStress,colors=new Map(),empty={positions:[],owners:[],visible:false};
    return {colors,diagrams:empty,diagramFill:empty,maximum:session.maximum,legend:{kind:'text',text:'max |σx| · normal stresses from axial force and bending only. See the normal stress scale and member details.',summary:'Maximum absolute normal stress'}};
  };
  const originalResult=M.planarPresentation.result;
  M.planarPresentation.result=function(ds,mode,count){return mode==='normalStress'?V.planar(ds):originalResult(ds,mode,count);};
  const nodeDetails=M.staticResults.nodeResultDetails;
  M.staticResults.nodeResultDetails=function(ds,id,mode){if(mode!=='normalStress')return nodeDetails(ds,id,mode);return ds.raw.model.elements.filter(e=>e.nodeIds.includes(id)).map(e=>{const v=V.at(ds.normalStress.results.get(e.id),e.nodeIds[0]===id?0:1);return 'E'+e.id+' · max |σx| = '+(v===null?'unavailable':fmt(v))+' · '+unit(ds,'stress');});};
}(window.MKEFPost,window.MKEFPost.normalStressResults));
