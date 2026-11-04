(function(M){
  'use strict';
  const P=M.SpatialPost.prototype,{el,fmt}=M.spatialUI,make=P.makeControls;
  function download(blob,name){const url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
  P.makeControls=function(){make.call(this);const box=el('div',undefined,this.panel);
    for(const [kind,title]of [['svg','Chart SVG'],['csv','Chart CSV']]){const b=el('button',title,box);b.id='spatial-export-'+kind;b.onclick=()=>{if(!this.series){this.status.textContent='Select an available chart first.';return;}this.exportChart(kind);};}
  };
  P.caption=function(){return [this.ds.raw.metadata.title||'MKE-F 3D',this.analysis.type==='modal'?'Mode '+(this.index+1)+' · f='+fmt(this.analysis.frequenciesHz[this.index])+' Hz · phase '+fmt(this.phase):this.analysis.type==='transient'?'t='+fmt(this.analysis.time[this.index])+' · frame '+this.index:'Static analysis',this.scaleNote.textContent,this.legend.textContent].filter(Boolean);};
  P.exportChart=function(kind){
    if(kind==='csv'){const quote=s=>'"'+String(s).replace(/"/g,'""')+'"',s=this.series;download(new Blob(['\uFEFF'+[quote(s.xLabel)+','+quote(s.yLabel),...s.x.map((v,i)=>String(v)+','+String(s.y[i]))].join('\r\n')],{type:'text/csv;charset=utf-8'}),'MKE-F-3D-chart.csv');return;}
    const svg=this.chart.cloneNode(true);svg.removeAttribute('hidden');svg.setAttribute('xmlns',M.svgNS);svg.setAttribute('width','800');svg.setAttribute('height','240');svg.setAttribute('class','chart-export');svg.removeAttribute('style');
    const style=M.svgElement('style');style.textContent=M.exporting.svgStyles;svg.prepend(style);const metadata=M.svgElement('metadata');metadata.textContent=this.caption().join('\n');svg.prepend(metadata);
    download(new Blob([new XMLSerializer().serializeToString(svg)],{type:'image/svg+xml;charset=utf-8'}),'MKE-F-3D-chart.svg');
  };
  P.exportPNG=async function(resolution=2){
    if(!this.scene||this.scene.lost)throw new Error('Scene PNG export requires WebGL 2.');
    this.stop();const s=this.scene,canvas=document.createElement('canvas'),lines=this.caption(),width=Math.max(s.width,800),sceneHeight=Math.round(s.pixels*width/s.width),height=sceneHeight+28+22*lines.length;
    canvas.width=width*resolution;canvas.height=height*resolution;const ctx=canvas.getContext('2d');ctx.scale(resolution,resolution);ctx.fillStyle='#fff';ctx.fillRect(0,0,width,height);
    // Capture synchronously after rendering; no permanent preserveDrawingBuffer or second WebGL context.
    s.depth();s.updateGrid();s.conditions.update();s.sections.update();s.renderer.render(s.scene,s.camera);s.labels();ctx.drawImage(s.canvas,0,0,width,sceneHeight);ctx.drawImage(s.overlay,0,0,width,sceneHeight);ctx.fillStyle='#263847';ctx.font='14px sans-serif';
    lines.forEach((line,i)=>ctx.fillText(line,16,sceneHeight+24+i*22,width-32));
    const blob=await new Promise(resolve=>canvas.toBlob(resolve,'image/png'));if(!blob)throw new Error('Could not create PNG.');download(blob,'MKE-F-3D.png');
  };
}(window.MKEFPost));
