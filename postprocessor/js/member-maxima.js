(function(M,sections){
  'use strict';
  const el=id=>document.getElementById('member-'+id);
  const api='/api/v1/sections/', finite=Number.isFinite;
  const supported=new Set(['rectangle','circle','rectangular_tube','square_tube','circular_tube']);
  let dataset=null,active=null,rows=[];

  function bodyFor(ds,e,nu2d){
    const result=ds.resultsByElementId.get(e.id);
    if(!result)throw Error('Missing member forces');
    const spatial=ds.raw.model.dimension===3,truss=[112,1133].includes(e.type);
    let axial,coeff;
    if(spatial){
      coeff=result.sectionForceCoefficients;
      axial=coeff?.[0]??[-result.localEndForces?.[0]];
    }else{
      const q=M.staticResults.elementIntensities(ds,e),L=q.length,f=result.localEndForces;
      axial=[-f[0],-q.x*L,-q.dx*L/2];
      if(!truss){
        const v=[-f[1],-q.y*L,-q.dy*L/2],m=[-f[2],f[1]*L,q.y*L*L/2,q.dy*L*L/6];
        const {c,s}=sections.rotation(e.sectionAngle??0);
        coeff=[axial,v.map(x=>c*x),v.map(x=>-s*x),[0],m.map(x=>s*x),m.map(x=>c*x)];
      }
    }
    if(truss){
      if(!axial.every(finite)||!(e.properties.area>0))throw Error('Missing axial force or area');
      return {axial,area:e.properties.area};
    }
    const profile=ds.raw.model.sections?.find(p=>p.id===e.sectionId);
    if(!profile)throw Error('Section geometry required');
    if(!supported.has(profile.type))throw Error('Unsupported profile: '+profile.type);
    if(!coeff||coeff.length!==6||!coeff.every(row=>Array.isArray(row)&&row.length&&row.every(finite)))throw Error('Missing force coefficients');
    const nu=spatial?e.properties.youngsModulus/(2*e.properties.shearModulus)-1:nu2d;
    if(!(nu>-1&&nu<.5))throw Error('Poisson’s ratio must be between −1 and 0.5');
    const dimensions=Object.fromEntries(sections.fields(profile.type).filter(k=>profile.dimensions[k]!==undefined).map(k=>[k,profile.dimensions[k]]));
    return {section:{type:profile.type,dimensions},nu,coefficients:coeff,J_model:e.properties.torsionConstant,include_mesh_peak:true};
  }

  function release(id){
    if(id)fetch(api+'jobs/'+id,{method:'DELETE',keepalive:true}).catch(()=>{});
  }
  function cancel(){
    if(!active)return;
    active.controller.abort();release(active.id);active=null;
    el('progress').textContent='Calculation cancelled.';
    for(const row of rows)if(row.state.textContent==='Calculating')row.state.textContent='Cancelled';
    el('cancel').disabled=true;el('run').disabled=false;
  }
  function ratio(row){
    const strength=el('yield').valueAsNumber;
    row.ratio.textContent=finite(row.value)&&strength>0&&finite(strength)?(row.value/strength).toPrecision(4):'—';
  }
  function newRow(element){
    const tr=el('rows').insertRow();
    tr.insertCell().textContent='E'+element.id;
    const row={element,stress:tr.insertCell(),xi:tr.insertCell(),ratio:tr.insertCell(),state:tr.insertCell()};
    for(const cell of [row.stress,row.xi,row.ratio])cell.textContent='—';
    row.state.textContent='Not calculated';
    return row;
  }
  async function json(response){
    if(!response.ok)throw Error(response.status===429?'Section service queue is full. Retry later.':'Section service unavailable (HTTP '+response.status+').');
    return response.json();
  }
  async function run(){
    if(active||!dataset)return;
    const task={controller:new AbortController(),id:null}, ds=dataset,nu=el('nu').valueAsNumber;
    active=task;el('run').disabled=true;el('cancel').disabled=false;
    el('rows').replaceChildren();rows=ds.raw.model.elements.map(newRow);
    try{
      for(const [index,row] of rows.entries()){
        el('progress').textContent='Member '+(index+1)+' of '+rows.length;
        let body;
        try{body=bodyFor(ds,row.element,nu);}catch(error){row.state.textContent=error.message;continue;}
        row.state.textContent='Calculating';
        const response=await fetch(api+'member',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});
        if(response.status===422){
          if(active!==task)return;
          row.state.textContent='Invalid section or force data';continue;
        }
        const submitted=await json(response);
        if(active!==task){release(submitted.id);return;}
        task.id=submitted.id;
        let job;
        do{
          job=await json(await fetch(api+'jobs/'+task.id,{signal:task.controller.signal}));
          if(active!==task)return;
          if(job.state!=='running'&&job.state!=='queued')break;
          await new Promise(resolve=>setTimeout(resolve,300));
          if(active!==task)return;
        }while(true);
        task.id=null;
        const r=job.result;
        if(job.state==='completed'&&r.status==='ok'&&finite(r.maximum?.value)){
          row.value=r.maximum.value;row.stress.textContent=row.value.toPrecision(6);
          row.xi.textContent=r.maximum.xi.toPrecision(4);row.state.textContent='OK';ratio(row);
        }else if(r?.maximum?.singularPeak && finite(r.maximum.mesh_peak?.value)){
          const peak=r.maximum.mesh_peak;
          row.meshPeak=true;
          row.stress.textContent=peak.value.toPrecision(6)+' (mesh)';
          row.xi.textContent=peak.xi.toPrecision(4);
          row.state.textContent='Sharp corner · mesh level '+r.history.at(-1).level;
          if(r.reason)row.state.textContent+=' · '+r.reason;
        }else{
          row.state.textContent=r?.singularPeak?'Unavailable · sharp inner corner':(r?.reason||r?.status||job.state);
        }
        const dims=body.section?.dimensions;
        if(dims?.inner_radius)row.state.textContent+=' · rᵢ='+dims.inner_radius+', rₒ='+Number((dims.inner_radius+dims.t).toPrecision(6));
      }
      const diagnosticCount=rows.filter(row=>row.meshPeak).length;
      el('progress').textContent=rows.filter(row=>finite(row.value)).length+' of '+rows.length+' confirmed member maxima.'+
        (diagnosticCount?' Mesh peaks shown for '+diagnosticCount+'. Stress at a sharp inner corner depends on mesh refinement; a yield ratio is unavailable.':'');
    }catch(error){
      if(active===task){
        release(task.id);
        const reason=error instanceof TypeError?'Could not reach the section service.':error.message;
        el('progress').textContent=reason;
        const unfinished=rows.find(row=>row.state.textContent==='Calculating');
        if(unfinished)unfinished.state.textContent=reason;
      }
    }finally{
      if(active===task){active=null;el('cancel').disabled=true;el('run').disabled=false;}
    }
  }

  M.memberMaxima={load(ds){
    cancel();dataset=ds;rows=[];
    const panel=el('maxima');panel.hidden=ds.raw.analysis.type!=='static';panel.open=false;
    el('rows').replaceChildren();el('yield').value='';el('nu').value='0.3';
    const units=ds.raw.metadata.units||{};
    el('stress-unit').textContent=units.stress||(units.force&&units.length?units.force+'/'+units.length+'²':'model force / length²');
    el('stress-heading').textContent='Peak σVM ('+el('stress-unit').textContent+')';
    el('nu-label').hidden=ds.raw.model.dimension===3||ds.raw.model.elements.every(e=>e.type===112);
    el('run').disabled=location.protocol==='file:';
    el('progress').textContent=el('run').disabled?'Open the postprocessor through Docker Compose to calculate stresses.':'Stress and reference yield strength use the same model units.';
    el('run').onclick=run;el('cancel').onclick=cancel;el('yield').oninput=()=>rows.forEach(ratio);
    el('nu').oninput=()=>{cancel();rows=[];el('rows').replaceChildren();el('progress').textContent='Calculate again for the updated Poisson ratio.';};
  }};
  window.addEventListener('pagehide',cancel);
})(window.MKEFPost,window.MKEFSections);
