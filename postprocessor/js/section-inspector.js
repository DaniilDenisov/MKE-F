(function(M,S,P){
  'use strict';
  const api='/api/v1/sections/',supported=new Set(['rectangle','circle','square_tube','rectangular_tube','circular_tube']);
  const fmt=v=>Number(v.toPrecision(7)).toString(),create=(tag,text,parent)=>{const e=document.createElement(tag);if(text!==undefined)e.textContent=text;if(parent)parent.append(e);return e;};
  const save=(blob,name)=>{const url=URL.createObjectURL(blob),a=create('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);};
  const release=id=>{if(id)fetch(api+'jobs/'+id,{method:'DELETE',keepalive:true}).catch(()=>{});};
  class Inspector {
    constructor(){
      this.xi=.5;this.cache=new Map();this.cacheBytes=0;this.active=null;this.pending=null;
      this.panel=document.getElementById('section-result-panel');this.panel.replaceChildren();this.panel.className='section-inspector';
      document.getElementById('result-sidebar').prepend(this.panel);
      const heading=create('div',undefined,this.panel);heading.className='section-inspector-heading';
      this.title=create('h2','Section stresses',heading);this.expand=create('button','Expand',heading);this.expand.id='section-expand';this.expand.onclick=()=>this.open();
      this.close=create('button','Close expanded view',heading);this.close.hidden=true;this.close.onclick=()=>this.dialog.close();
      this.caption=create('p','',this.panel);this.caption.id='section-caption';
      const controls=create('div',undefined,this.panel);controls.className='section-coordinate-controls';
      const input=(parent,title,id,type='number')=>{const label=create('label',title,parent),e=create('input',undefined,label);e.id=id;e.type=type;e.step='any';return e;};
      this.range=input(controls,'Position x/L','section-position','range');this.range.min=0;this.range.max=1;this.range.step=.001;
      this.relative=input(controls,'x/L','section-position-number');this.relative.min=0;this.relative.max=1;
      this.distance=input(controls,'x','section-distance');this.distance.min=0;
      for(const field of [this.range,this.relative,this.distance])field.oninput=()=>{
        const value=field.value.trim()===''?NaN:Number(field.value),xi=field===this.distance?value/this.data.length:value;
        field.setCustomValidity(Number.isFinite(xi)&&xi>=0&&xi<=1?'':'Use a coordinate within the member.');
        if(field.checkValidity())this.setPosition(xi);
      };
      this.nu=input(this.panel,'Poisson’s ratio ν (2D frames)','section-nu');this.nu.value=.3;
      this.nu.oninput=()=>{const source=document.getElementById('member-nu');source.value=this.nu.value;source.dispatchEvent(new Event('input'));};
      document.getElementById('member-nu').addEventListener('input',()=>{this.nu.value=document.getElementById('member-nu').value;if(this.element)this.refresh();});
      this.maps=create('div',undefined,this.panel);this.maps.className='section-maps';
      const normal=create('article',undefined,this.maps);create('h3','Normal stress σx',normal);
      this.svg=M.svgElement('svg',{id:'section-result-svg',xmlns:'http://www.w3.org/2000/svg',role:'img','aria-label':'Normal stress across the section'});normal.append(this.svg);
      this.normalInfo=create('p','',normal);this.normalInfo.id='section-result-info';
      const vm=create('article',undefined,this.maps);create('h3','Von Mises stress σVM',vm);
      this.canvas=create('canvas',undefined,vm);this.canvas.id='section-vm-field';this.canvas.setAttribute('aria-label','Von Mises stress across the section');
      this.vmInfo=create('p','',vm);this.vmInfo.id='section-vm-info';this.vmInfo.setAttribute('role','status');
      this.retry=create('button','Retry section calculation',vm);this.retry.id='section-retry';this.retry.onclick=()=>this.schedule(true);
      this.wire=input(vm,'Show mesh','section-mesh','checkbox');this.wire.onchange=()=>this.draw();
      this.details=create('details',undefined,this.panel);create('summary','Calculation details',this.details);this.report=create('pre','',this.details);
      this.point=create('details',undefined,this.panel);create('summary','Normal stress at a material point',this.point);
      this.pointY=input(this.point,'Material point · Y','section-point-y');this.pointZ=input(this.point,'Material point · Z','section-point-z');this.pointY.value=this.pointZ.value=0;
      this.pointButton=create('button','Show point stress',this.point);this.pointButton.id='section-point';this.pointButton.onclick=()=>this.showPoint();
      this.pointInfo=create('p','',this.point);
      this.svg.onclick=event=>{if(!this.layout)return;const point=this.svg.createSVGPoint();point.x=event.clientX;point.y=event.clientY;const q=point.matrixTransform(this.svg.getScreenCTM().inverse()),f=this.layout.frame;
        this.pointY.value=fmt((f.cy-q.y)/f.scale);this.pointZ.value=fmt((q.x-f.cx)/f.scale);this.point.open=true;this.showPoint();};
      const actions=create('div',undefined,this.panel);actions.className='section-export-actions';
      this.exportSVG=create('button','SVG · normal stress',actions);this.exportSVG.id='export-section-svg';this.exportSVG.onclick=()=>this.saveSVG();
      this.exportCSV=create('button','CSV · normal stress',actions);this.exportCSV.id='section-export-csv';this.exportCSV.onclick=()=>this.saveCSV();
      this.exportPNG=create('button','PNG · both maps',actions);this.exportPNG.id='section-export-png';this.exportPNG.onclick=()=>this.savePNG().catch(e=>{this.vmInfo.textContent=e.message;});
      this.dialog=create('dialog',undefined,document.body);this.dialog.className='section-expanded';this.dialog.setAttribute('aria-label','Section stresses');
      this.dialog.addEventListener('close',()=>{document.getElementById('result-sidebar').prepend(this.panel);this.close.hidden=true;this.expand.hidden=false;this.draw();if(!this.panel.hidden)this.expand.focus();});
      this.observer=new ResizeObserver(()=>{if(!this.panel.hidden)this.draw();});this.observer.observe(this.maps);
      window.addEventListener('pagehide',()=>this.cancel());
    }
    load(ds){
      this.hide();this.ds=ds;this.xi=.5;this.cache.clear();this.cacheBytes=0;this.nu.value=document.getElementById('member-nu').value;
    }
    select(e,viewer){
      if(!e||this.ds?.raw.analysis.type!=='static'){this.hide();return;}
      this.viewer=viewer||null;
      if(this.element===e){if(viewer&&Number(viewer.position.value)!==this.xi)this.setPosition(Number(viewer.position.value),false);return;}
      this.element=e;this.panel.hidden=false;
      const spatial=this.ds.raw.model.dimension===3;
      this.svg.id=spatial?'spatial-section-svg':'section-result-svg';this.pointY.id=spatial?'spatial-result-point-y':'section-point-y';this.pointZ.id=spatial?'spatial-result-point-z':'section-point-z';
      this.pointButton.id=spatial?'spatial-section-point':'section-point';this.exportSVG.id=spatial?'spatial-section-export-svg':'export-section-svg';this.exportCSV.id=spatial?'spatial-section-export-csv':'section-export-csv';
      this.point.id=spatial?'spatial-section-cut':'section-point-controls';
      if(viewer)viewer.position.value=this.xi;
      this.refresh();
    }
    refresh(){
      this.title.textContent=`Section stresses · E${this.element.id}`;this.caption.textContent='';
      try{
        this.data=M.sectionMap.adapt(this.ds,this.element,Number(this.nu.value));
        for(const field of [this.range,this.relative,this.distance])field.disabled=false;
        this.caption.textContent=(this.data.section?.name||'No section geometry')+` · orientation ${this.data.angle}° · x from node ${this.element.nodeIds[0]}`;
        this.nu.parentElement.hidden=this.ds.raw.model.dimension===3||this.data.truss;
        this.point.hidden=!this.data.profile;this.exportSVG.disabled=this.exportCSV.disabled=!this.data.profile;
        this.sync();this.schedule();
      }catch(error){this.data=null;this.entry=null;this.pending=null;this.currentKey=null;clearTimeout(this.timer);this.svg.replaceChildren();this.svg.hidden=this.canvas.hidden=true;this.point.hidden=true;this.normalInfo.textContent=error.message;this.vmInfo.textContent='Section stresses unavailable.';this.retry.hidden=true;this.exportSVG.disabled=this.exportCSV.disabled=this.exportPNG.disabled=true;for(const field of [this.range,this.relative,this.distance])field.disabled=true;}
    }
    sync(){
      for(const e of [this.range,this.relative,this.distance])e.setCustomValidity('');
      this.range.value=this.relative.value=String(this.xi);this.distance.value=String(this.xi*this.data.length);this.distance.max=this.data.length;
      this.distance.parentElement.firstChild.textContent=`x (${this.data.lengthUnit})`;
    }
    setPosition(xi,notify=true){
      if(!this.data||!Number.isFinite(xi)||xi<0||xi>1)return;
      if(this.xi===xi){this.sync();return;}this.xi=xi;this.sync();this.schedule();
      if(this.viewer&&notify){this.viewer.position.value=xi;this.viewer.update();}
    }
    request(){
      const d=this.data,forces=d.coeff.map(c=>S.evaluate(c,this.xi));
      if(d.truss||forces.slice(1,4).every(v=>v===0))return null;
      if(!d.profile)throw Error('Assign a section in the preprocessor and rerun to view stress maps.');
      if(!supported.has(d.section.type))throw Error('Von Mises with shear or torsion is not supported for this profile. Normal stress remains available.');
      if(!Number.isFinite(d.nu)||!(d.nu>-1&&d.nu<.5)||this.ds.raw.model.dimension===2&&this.nu.value.trim()==='')throw Error('Enter a Poisson’s ratio between −1 and 0.5.');
      if(location.protocol==='file:')throw Error('Open the postprocessor through Docker Compose to calculate von Mises with shear or torsion.');
      const dimensions=Object.fromEntries(S.fields(d.section.type).filter(k=>d.section.dimensions[k]!==undefined).map(k=>[k,d.section.dimensions[k]]));
      const body={section:{type:d.section.type,dimensions},nu:d.nu,forces};
      if(d.element.properties.torsionConstant)body.J_model=d.element.properties.torsionConstant;
      return body;
    }
    schedule(immediate=false){
      clearTimeout(this.timer);this.pending=null;this.entry=null;this.local=false;this.currentKey=null;this.retry.hidden=true;this.status='';
      try{
        const body=this.request();
        if(!body){this.local=true;this.status=this.data.profile?'Axial and bending stress · σVM = |σx|':'No section geometry. Axial stress values are available.';}
        else{
          this.currentKey=JSON.stringify(body);
          if(this.cache.has(this.currentKey)){this.entry=this.cache.get(this.currentKey);this.cache.delete(this.currentKey);this.cache.set(this.currentKey,this.entry);}
          else {this.status='Calculating von Mises for the selected section…';const task={key:this.currentKey,body,ds:this.ds};this.timer=setTimeout(()=>{this.pending=task;this.pump();},immediate?0:350);}
        }
      }catch(error){this.status=error.message;this.retry.hidden=false;}
      this.draw();
    }
    async json(url,options){
      const response=await fetch(url,options);
      if(!response.ok)throw Error(response.status===429?'Section service queue is full. Retry when ready.':response.status===404?'Section field expired. Retry the calculation.':`Section service unavailable (HTTP ${response.status}).`);
      return response;
    }
    async pump(){
      if(this.active||!this.pending)return;
      const task=this.pending;this.pending=null;this.active=task;task.controller=new AbortController();
      try{
        // Do not abort the POST: a late identifier must still be released.
        const submitted=await (await this.json(api+'diagnostic',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(task.body)})).json();task.id=submitted.id;
        if(this.active!==task){release(task.id);return;}
        let job;
        while(this.active===task){
          job=await (await this.json(api+'jobs/'+task.id,{signal:task.controller.signal})).json();
          if(!['queued','running'].includes(job.state))break;
          await new Promise(resolve=>setTimeout(resolve,300));
        }
        if(this.active!==task)return;
        if(job.state==='failed')throw Error(job.result?.reason||'Section calculation failed.');
        const result=job.result;
        if(result?.display?.status!=='ok')throw Error(`Section field unavailable: ${result?.reason||result?.status||job.state}.`);
        const buffer=await (await this.json(api+'jobs/'+task.id+'/field',{signal:task.controller.signal})).arrayBuffer();
        if(this.active!==task)return;
        const entry={result,field:new window.SectionField(buffer),bytes:buffer.byteLength};
        while(this.cache.size&&this.cacheBytes+entry.bytes>48*1024**2){const key=this.cache.keys().next().value;this.cacheBytes-=this.cache.get(key).bytes;this.cache.delete(key);}
        if(entry.bytes<=48*1024**2){if(this.cache.has(task.key))this.cacheBytes-=this.cache.get(task.key).bytes;this.cache.set(task.key,entry);this.cacheBytes+=entry.bytes;}
        if(this.ds===task.ds&&this.currentKey===task.key){this.entry=entry;this.retry.hidden=true;this.draw();}
      }catch(error){if(this.active===task&&this.ds===task.ds&&this.currentKey===task.key){this.status=error instanceof TypeError?'Cannot reach the section service. Normal stress remains available.':error.message;this.retry.hidden=false;this.draw();}}
      finally{
        release(task.id);
        if(this.active===task){this.active=null;if(this.pending){if(this.cache.has(this.pending.key)&&this.pending.key===this.currentKey){this.entry=this.cache.get(this.pending.key);this.pending=null;this.draw();}else this.pump();}}
      }
    }
    draw(){
      if(!this.data||this.panel.hidden)return;
      const d=this.data,at=d.normal.at(this.xi);
      this.normalInfo.textContent=`σmin=${fmt(at.sigmaMin)} · σmax=${fmt(at.sigmaMax)}\n${d.unit} · x/L=${fmt(this.xi)} · x=${fmt(this.xi*d.length)} ${d.lengthUnit}`;
      this.svg.hidden=!d.profile;this.canvas.hidden=!d.profile||(!this.local&&!this.entry);this.wire.parentElement.hidden=!this.entry;
      this.exportPNG.disabled=this.canvas.hidden;
      if(d.profile){
        this.layout=M.sectionMap.normal(this.svg,d,this.xi,d.normal.sigmaAbsMax,this.svg.getBoundingClientRect().width);
        if(!this.canvas.hidden)M.sectionMap.vm(this.canvas,d,this.xi,this.layout,this.entry,this.wire.checked);
      }else this.layout=null;
      if(this.entry){
        const r=this.entry.result,peak=r.maximum?.mesh_peak||r.maximum,singular=!!r.maximum?.singularPeak,confirmed=r.maximumStatus==='ok';
        this.vmInfo.textContent=(Number.isFinite(peak?.value)?`${singular?'Mesh peak · sharp inner corner':confirmed?'Maximum σVM':'Current mesh peak'} = ${fmt(peak.value)} ${d.unit}. `:'')+
          (singular?'No finite global maximum is reported. ':confirmed?'':'Maximum convergence is not confirmed. ')+(r.fieldStatus==='ok'?'':'Field convergence is not confirmed. ');
        this.report.textContent=JSON.stringify(r,null,2);this.details.hidden=false;
      }else{
        this.vmInfo.textContent=this.local?`Maximum σVM = ${fmt(at.sigmaAbsMax)} ${d.unit}. ${this.status}`:this.status;
        this.details.hidden=true;
      }
      this.showPoint();
    }
    showPoint(){
      if(!this.data?.profile)return;
      const y=this.pointY.value.trim()?Number(this.pointY.value):NaN,z=this.pointZ.value.trim()?Number(this.pointZ.value):NaN;
      const inside=P.contains(this.data.profile,y,z);this.svg.querySelector('.section-material-point')?.remove();
      if(inside&&this.layout){const [x,yPixel]=this.layout.frame.project(y,z);this.svg.append(M.svgElement('circle',{class:'section-material-point',cx:x,cy:yPixel,r:3,fill:'#111',stroke:'white'}));}
      this.pointInfo.textContent=inside?`σx(${fmt(y)}, ${fmt(z)}) = ${fmt(this.data.field.point(this.xi,y,z).value)} ${this.data.unit}`:'Point outside the material or inside a hole: stress is undefined.';
    }
    open(){this.dialog.append(this.panel);this.expand.hidden=true;this.close.hidden=false;this.dialog.showModal();this.draw();this.close.focus();}
    cancel(){clearTimeout(this.timer);this.pending=null;const task=this.active;this.active=null;if(task){task.controller.abort();release(task.id);}}
    hide(){
      this.cancel();this.element=null;this.data=null;this.currentKey=null;this.entry=null;this.panel.hidden=true;if(this.dialog.open)this.dialog.close();
      // The modal/transient viewer owns the spatial IDs when the static panel is hidden.
      this.svg.id='section-result-svg';this.pointY.id='section-point-y';this.pointZ.id='section-point-z';this.point.id='section-point-controls';
      this.pointButton.id='section-point';this.exportSVG.id='export-section-svg';this.exportCSV.id='section-export-csv';
    }
    svgText(){const svg=this.svg.cloneNode(true),meta=M.svgElement('metadata');meta.textContent=this.title.textContent+' · '+this.caption.textContent+' · '+this.normalInfo.textContent;svg.prepend(meta);return new XMLSerializer().serializeToString(svg);}
    saveSVG(){if(this.data?.profile)save(new Blob([this.svgText()],{type:'image/svg+xml'}),`section-E${this.element.id}.svg`);}
    saveCSV(){if(!this.data?.profile)return;const rows=this.data.profile.loops.flatMap((loop,boundary)=>loop.map(v=>({boundary:boundary+1,...this.data.field.point(this.xi,v.y,v.z)}))),at=this.data.field.at(this.xi);rows.push({boundary:0,...at.minimum},{boundary:0,...at.maximum});
      save(new Blob(['boundary,x/L,y,z,sigmaX\r\n'+rows.map(v=>[v.boundary,v.xi,v.y,v.z,v.value].join(',')).join('\r\n')],{type:'text/csv'}),`section-E${this.element.id}.csv`);}
    async savePNG(){
      if(this.canvas.hidden||!this.layout)return;
      // Capture all state before image decoding; an arriving result cannot mix sections.
      const {w,h}=this.layout,canvas=create('canvas');canvas.width=2*w+36;canvas.height=h+170;const ctx=canvas.getContext('2d');
      const title=this.title.textContent,caption=this.caption.textContent,normal=this.normalInfo.textContent,vm=this.vmInfo.textContent,name=`section-E${this.element.id}.png`;
      ctx.fillStyle='white';ctx.fillRect(0,0,canvas.width,canvas.height);ctx.drawImage(this.canvas,w+24,110);
      const url=URL.createObjectURL(new Blob([this.svgText()],{type:'image/svg+xml'}));
      try{const img=new Image();await new Promise((resolve,reject)=>{img.onload=resolve;img.onerror=()=>reject(Error('Could not export section maps.'));img.src=url;});ctx.drawImage(img,12,110);}finally{URL.revokeObjectURL(url);}
      ctx.fillStyle='#17313e';ctx.font='bold 17px sans-serif';ctx.fillText(title,12,25);ctx.font='13px sans-serif';ctx.fillText(caption,12,48);ctx.fillText(normal.replace('\n',' · '),12,70);
      ctx.font='bold 15px sans-serif';ctx.fillText('Normal stress σx',12,100);ctx.fillText('Von Mises stress σVM',w+24,100);
      ctx.font='12px sans-serif';let line='',y=h+135;for(const word of vm.split(' ')){if(ctx.measureText(line+word).width>canvas.width-24){ctx.fillText(line,12,y);y+=16;line='';}line+=word+' ';}ctx.fillText(line,12,y);
      const blob=await new Promise(resolve=>canvas.toBlob(resolve));if(blob)save(blob,name);
    }
  }
  M.sectionInspector={load(ds){if(!this.view)this.view=new Inspector();this.view.load(ds);},select(e,viewer){this.view?.select(e,viewer);},hide(){this.view?.hide();}};
  M.sectionResults.show=function(ds,selection){if(ds)M.sectionInspector.select(selection?.kind==='element'?ds.elementsById.get(selection.id):null);};
})(window.MKEFPost,window.MKEFSections,window.MKEFSpatialSections);
