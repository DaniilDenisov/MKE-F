(function(global){
 'use strict';
 var M=global.MKEFPre,F=M.editorAdapter,state,spatialPanel,selection=null,dirty=false;
 var elements={},apiReady=false,activeJobId=null,pollTimer=null;
 var activeJobStorageKey='mkef-active-job',submitting=false,activeSubmission=null,completedSubmission=null;
 function byId(id){return document.getElementById(id);}
 function setError(message){elements.error.textContent=message||'';elements.error.hidden=!message;}
 function resultSnapshotNote(){const context=activeSubmission||completedSubmission;if(!context)return;let note=byId('spatial-run-snapshot');if(!note){note=document.createElement('p');note.id='spatial-run-snapshot';elements.runPanel.append(note);}const changed=context.model!==state||context.revision!==spatialPanel.history.revision;note.textContent='Analysis snapshot «'+context.name+'», revision '+context.revision+'.'+(changed?' The current model differs from the submitted snapshot.':'');}
  function storedJobId() {
    try { return global.sessionStorage.getItem(activeJobStorageKey); } catch (_) { return null; }
  }
  function storeJobId(value) {
    try { if (value) global.sessionStorage.setItem(activeJobStorageKey, value); else global.sessionStorage.removeItem(activeJobStorageKey); } catch (_) { /* storage is optional */ }
  }
  function runStatusLabel(status) {
    return { queued: 'Queued', running: 'Running in GNU Octave', succeeded: 'Completed', failed: 'Calculation failed', canceled: 'Canceled', timed_out: 'Timed out' }[status] || status;
  }
  function updateRunPanel(job) {
    elements.runPanel.hidden = false;
    elements.runStatus.textContent = runStatusLabel(job.status);
    elements.runElapsed.textContent = Number(job.elapsedSeconds || 0).toFixed(1) + ' s';
    elements.runLog.textContent = job.logTail || (job.status === 'succeeded' ? 'Calculation completed. The solver produced no log messages.' : 'No solver output yet.');
    elements.cancelRun.disabled = ['queued', 'running'].indexOf(job.status) < 0;
  }
  function clearPollTimer() { if (pollTimer) { global.clearTimeout(pollTimer); pollTimer = null; } }
  function finishActiveJob(job) {
    clearPollTimer(); storeJobId(null); activeJobId = null; render();
    if (job.status === 'failed' || job.status === 'timed_out') setError(job.error && job.error.message ? job.error.message : 'The calculation failed.');
  }

  async function pollActiveJob() {
    if (!activeJobId) return;
    var pollingId=activeJobId;
    try {
      var job = await global.MKEFApi.getJob(pollingId);
      if(activeJobId!==pollingId)return;
      updateRunPanel(job);
      if (job.status === 'succeeded') {
        var completedId = activeJobId;
        {completedSubmission=activeSubmission||{name:"Recovered job",revision:null,model:null};activeSubmission=null;clearPollTimer();storeJobId(null);activeJobId=null;
          var link=byId('spatial-open-result');if(!link){link=document.createElement('a');link.id='spatial-open-result';link.className='button';elements.runPanel.appendChild(link);}link.textContent='Open result';link.href='../postprocessor/?job='+encodeURIComponent(completedId);link.target='_blank';link.rel='noopener';link.hidden=false;render();resultSnapshotNote();return;}

      }
      if (['failed', 'canceled', 'timed_out'].indexOf(job.status) >= 0) { finishActiveJob(job); return; }
      pollTimer = global.setTimeout(pollActiveJob, 1000);
    } catch (error) {
      if(activeJobId!==pollingId)return;
      if (error.status === 404) { clearPollTimer(); storeJobId(null); activeJobId = null; render(); }
      else pollTimer=global.setTimeout(pollActiveJob,1000);
      setError(error.message);
    }
  }

  async function submitCurrentCase() {
    if (!apiReady || activeJobId || submitting || spatialPanel.analysisInvalid) return;
    var caseText;
    try { var validation=F.validate(state,true);if(validation.length)throw new Error(validation.join('\n'));caseText=F.serialize(state); } catch (error) { setError(error.message); return; }
    submitting=true;activeSubmission={dimension:state.dimension||2,model:state,name:state.name,revision:spatialPanel.history.revision};completedSubmission=null;var prior=byId('spatial-open-result');if(prior)prior.hidden=true;
    elements.run.disabled = true; setError('');
    elements.runPanel.hidden = false; elements.runStatus.textContent = 'Submitting'; elements.runElapsed.textContent = '0.0 s'; elements.runLog.textContent = 'No solver output yet.'; elements.cancelRun.disabled = true;
    try {
      var job = await global.MKEFApi.createJob(state.name || 'Case', caseText);
      try{global.sessionStorage.setItem('mkef-job-context',JSON.stringify({id:job.id,dimension:activeSubmission.dimension,name:activeSubmission.name,revision:activeSubmission.revision}));}catch(_){}
      activeJobId = job.id; storeJobId(activeJobId); updateRunPanel(job); pollActiveJob();
    } catch (error) { elements.runPanel.hidden = true; activeSubmission=null;setError(error.message); }
    finally{submitting=false;render();resultSnapshotNote();}
  }

  async function cancelActiveJob() {
    if (!activeJobId) return;
    var cancelingId=activeJobId;
    elements.cancelRun.disabled = true;
    try { var job = await global.MKEFApi.cancelJob(cancelingId);if(activeJobId!==cancelingId)return;updateRunPanel(job); finishActiveJob(job); }
    catch (error) {if(activeJobId!==cancelingId)return;setError(error.message); elements.cancelRun.disabled = false; }
  }

  async function configureSolver() {
    if (!global.MKEFApi || !global.MKEFApi.isHosted()) {
      elements.run.title = 'Run analysis is available when the application is started with Docker Compose.';
      elements.solverMessage.textContent = elements.run.title;
      elements.solverMessage.hidden = false;
      render(); return;
    }
    try {
      var health = await global.MKEFApi.health();
      apiReady = health.status === 'ready';
      elements.run.title = apiReady ? 'Run this case with ' + health.octaveVersion : 'The local solver is unavailable.';
      elements.solverMessage.textContent = apiReady ? '' : elements.run.title;
      elements.solverMessage.hidden = apiReady;
      var candidate = storedJobId();
      if (candidate && global.MKEFApi.isValidJobId(candidate)) { activeJobId = candidate;try{activeSubmission=JSON.parse(global.sessionStorage.getItem('mkef-job-context'));if(activeSubmission&&activeSubmission.id!==candidate)activeSubmission=null;}catch(_){}pollActiveJob(); }
    } catch (error) {
      apiReady = false; elements.run.title = error.message;
      elements.solverMessage.textContent = 'Local solver unavailable: ' + error.message;
      elements.solverMessage.hidden = false;
    }
    render();
  }


 function render(){
  document.body.classList.add('spatial-mode');
  elements.elementType.value=String(state.elementType);elements.analysisType.value=state.analysis.type;elements.caseName.value=state.name;
  elements.elementType.disabled=!!(state.nodes.length||state.elements.length);
  document.querySelectorAll('[data-spatial-tool]').forEach(b=>{b.disabled=b.dataset.spatialTool!=='select'&&!state.elementType;});
  for(const id of ['spatial-controls','spatial-tables'])byId(id).hidden=false;
  byId('spatial-viewport').hidden=false;
  spatialPanel.render(state,selection);
  const errors=F.validate(state,true);elements.run.disabled=!apiReady||!!activeJobId||submitting||!!errors.length||spatialPanel.analysisInvalid;
  elements.download.disabled=!!errors.length||spatialPanel.analysisInvalid;
  elements.undo.disabled=!spatialPanel.history.undo.length;elements.redo.disabled=!spatialPanel.history.redo.length;
  dirty=spatialPanel.history.dirty();elements.caseStatus.textContent=state.name+' · '+state.dimension+'D'+(dirty?' · modified':'');
  elements.validationStatus.textContent=errors.length?'Draft: '+errors[0]:'Model is ready to run';
  elements.solverMessage.textContent=apiReady?'':elements.run.title;elements.solverMessage.hidden=apiReady;resultSnapshotNote();
 }
 function replace(model){spatialPanel.leave();state=model;selection=null;setError('');render();}
 function discard(){return !dirty||global.confirm('Discard the current unsaved changes?');}
 function openFile(file){if(!file||!discard())return;if(file.size>M.config.maxFileBytes){setError('The selected file exceeds the 10 MiB limit.');return;}const reader=new FileReader();reader.onload=()=>{try{const text=new TextDecoder('utf-8',{fatal:true,ignoreBOM:true}).decode(reader.result),model=F.parse(text);model.name=file.name.replace(/\.txt$/i,'');replace(model);}catch(e){setError(e.message);}};reader.onerror=()=>setError('Could not read file.');reader.readAsArrayBuffer(file);}
 document.addEventListener('DOMContentLoaded',()=>{
  const ids={caseStatus:'case-status',fileInput:'file-input',download:'download-case',undo:'undo',redo:'redo',run:'run-analysis',solverMessage:'solver-message',runPanel:'run-panel',runStatus:'run-status',runElapsed:'run-elapsed',runLog:'run-log',cancelRun:'cancel-run',error:'error-panel',notice:'notice-panel',caseName:'case-name',analysisType:'analysis-type',elementType:'element-type',validationStatus:'validation-status'};
  for(const key in ids)elements[key]=byId(ids[key]);
  for(const type of global.MKEFSpatial.registry.types.filter(t=>t.dimension===3)){const o=document.createElement('option');o.value=type.code;o.textContent=type.code+' · 3D '+type.family;elements.elementType.append(o);}
  spatialPanel=new M.SpatialPanel(next=>{selection=next;render();},()=>{selection=spatialPanel.history.selection;render();});
  state=F.empty(0);render();
  const exampleSelect=byId('load-example');
  for(const analysis of ['static','modal','transient'])for(const dimension of [2,3]){
    const group=document.createElement('optgroup');group.label=analysis[0].toUpperCase()+analysis.slice(1)+' · '+dimension+'D';
    for(const example of M.examples.filter(e=>e.analysis===analysis&&e.dimension===dimension)){
      const option=document.createElement('option');option.value=example.id;option.textContent=example.label;group.append(option);
    }
    exampleSelect.append(group);
  }
  exampleSelect.onchange=()=>{
    const example=M.examples.find(e=>e.id===exampleSelect.value);exampleSelect.value='';
    if(!example||!discard())return;
    try{const model=F.parse(example.text);model.name=example.name;replace(model);}catch(e){setError(e.message);}
  };
  byId('new-case').onclick=()=>{if(discard())replace(F.empty(0));};
  elements.fileInput.onchange=()=>{openFile(elements.fileInput.files[0]);elements.fileInput.value='';};
  elements.elementType.onchange=()=>{if(state.nodes.length||state.elements.length){render();return;}const next=F.empty(Number(elements.elementType.value));next.name=state.name;next.analysis={...state.analysis};replace(next);};
  elements.analysisType.onchange=()=>spatialPanel.changeAnalysis(elements.analysisType.value);
  elements.caseName.onchange=()=>spatialPanel.execute(h=>h.rename(elements.caseName.value.trim()||'Case'));
  elements.undo.onclick=()=>spatialPanel.travel(true);elements.redo.onclick=()=>spatialPanel.travel(false);
  byId('reset-view').onclick=()=>spatialPanel.scene&&spatialPanel.scene.fitModel();
  elements.download.onclick=()=>{try{const errors=F.validate(state,true);if(errors.length)throw new Error(errors.join('\n'));const text=spatialPanel.history.text(),url=URL.createObjectURL(new Blob([text],{type:'text/plain;charset=utf-8'})),a=document.createElement('a');a.href=url;a.download=state.name.replace(/[^A-Za-z0-9._ -]/g,'_')+'.txt';a.click();setTimeout(()=>URL.revokeObjectURL(url),0);spatialPanel.history.saved=spatialPanel.history.revision;render();}catch(e){setError(e.message);}};
  elements.run.onclick=submitCurrentCase;elements.cancelRun.onclick=cancelActiveJob;
  document.addEventListener('keydown',event=>{if(event.isComposing||document.querySelector('dialog[open]'))return;const typing=event.target.closest('input,textarea,select,[contenteditable="true"]');if(!typing&&(event.ctrlKey||event.metaKey)&&['z','y'].includes(event.key.toLowerCase())){event.preventDefault();spatialPanel.travel(event.key.toLowerCase()==='z'&&!event.shiftKey);}else if(!typing&&event.key==='Delete'&&!event.ctrlKey&&!event.metaKey&&!event.altKey){event.preventDefault();spatialPanel.confirmDelete();}else if(event.key==='Escape'){event.preventDefault();spatialPanel.setTool('select');selection=null;render();}});
  const drop=byId('drop-zone');drop.addEventListener('dragover',e=>e.preventDefault());drop.addEventListener('drop',e=>{e.preventDefault();openFile(e.dataTransfer.files[0]);});
  global.addEventListener('beforeunload',e=>{if(!dirty)return;e.preventDefault();e.returnValue='';});configureSolver();
 });
}(window));
