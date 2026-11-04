(()=>{
  'use strict';
  const el=id=>document.getElementById(id), base='/api/v1/sections/';
  const labels={ok:'Расчёт завершён.',nonconverged:'Сходимость не подтверждена.',
    resource_limit:'Достигнут лимит ресурсов.',cancelled:'Расчёт отменён.',unavailable:'Глобальный максимум недоступен.'};
  let current=null,field=null,peak=null;
  const fmt=v=>Number(v.toPrecision(7)).toString();
  function render(){
    if(!field)return;
    const [lo,hi]=field.draw(el('field'),Number(el('component').value),el('mesh').checked,peak);
    el('low').textContent=fmt(lo);el('high').textContent=fmt(hi);
  }
  function release(id){
    if(id)fetch(base+'jobs/'+id,{method:'DELETE',keepalive:true}).catch(()=>{});
  }
  function stop(){
    if(current){current.controller.abort();release(current.id);current=null;}
    el('cancel').disabled=true;
  }
  async function checked(url,options){
    const response=await fetch(url,options);
    if(!response.ok){
      if(response.status===429)throw Error('Очередь заполнена. Повторите расчёт позже.');
      if(response.status===422)throw Error('Проверьте тип профиля, размеры, ν и шесть усилий.');
      if(response.status===404)throw Error('Данные задания больше недоступны. Повторите расчёт.');
      throw Error('Сервис сечений недоступен (HTTP '+response.status+').');
    }
    return response;
  }
  el('component').onchange=el('mesh').onchange=render;
  el('cancel').onclick=()=>{stop();el('progress').textContent=labels.cancelled;};
  window.addEventListener('pagehide',stop);
  if(location.protocol==='file:'){
    el('run').disabled=true;
    el('progress').textContent='Откройте страницу через Docker Compose: http://localhost:8080/reference/section-fem-lab.html';
    return;
  }
  el('inputs').onsubmit=async event=>{
    event.preventDefault();stop();
    const task={id:null,controller:new AbortController()};current=task;
    field=peak=null;el('plot').hidden=el('details').hidden=true;el('critical').textContent='';
    el('cancel').disabled=false;el('progress').textContent='Отправка задания…';
    try{
      const body={section:JSON.parse(el('section').value),nu:Number(el('nu').value),forces:JSON.parse(el('forces').value)};
      // Let a late POST return its id, so a superseded submission can be cancelled.
      const submitted=await (await checked(base+'diagnostic',{method:'POST',
        headers:{'Content-Type':'application/json'},body:JSON.stringify(body)})).json();
      task.id=submitted.id;
      if(current!==task){release(task.id);return;}
      const url=base+'jobs/'+task.id, options={signal:task.controller.signal};
      let job;
      while(current===task){
        job=await (await checked(url,options)).json();
        if(current!==task)return;
        if(!['queued','running'].includes(job.state))break;
        el('progress').textContent=job.state==='queued'?'Задание в очереди…':'Расчёт поля…';
        await new Promise(resolve=>setTimeout(resolve,250));
      }
      if(current!==task)return;
      const r=job.result;
      if(job.state==='failed')throw Error(r.reason||'Ошибка расчёта.');
      el('progress').textContent=labels[r.status]||r.status;
      el('result').textContent=JSON.stringify(r,null,2);el('details').hidden=false;
      peak=r.maximum?.mesh_peak||r.maximum;
      if(Number.isFinite(peak?.value)){
        const title=r.maximum.singularPeak?'Диагностический сеточный пик · глобальный максимум unavailable':
          r.maximumStatus==='ok'?'Максимум по Мизесу':'Пик текущего уровня сетки';
        el('critical').textContent=title+' · σVM = '+fmt(peak.value)+' · y = '+fmt(peak.y)+' · z = '+fmt(peak.z)+
          ' · σx = '+fmt(peak.sigmaX)+' · τxy = '+fmt(peak.tauXY)+' · τxz = '+fmt(peak.tauXZ);
      }else peak=null;
      if(r.display?.status==='ok'){
        const buffer=await (await checked(url+'/field',options)).arrayBuffer();
        if(current!==task)return;
        field=new SectionField(buffer);el('plot').hidden=false;render();
      }else if(r.display)el('progress').textContent+=' Поле для просмотра: '+(labels[r.display.status]||r.display.status);
    }catch(error){
      if(current===task)el('progress').textContent=error instanceof TypeError?'Нет связи с сервисом сечений.':error.message;
    }finally{
      if(current===task){el('cancel').disabled=true;current=null;}
    }
  };
})();
