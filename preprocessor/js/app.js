(function (global) {
  'use strict';
  var M = global.MKEFPre, F = M.caseFormat;
  var state, renderer, undoStack = [], redoStack = [], dirty = false, selection = null;
  var elements = {};
  var apiReady = false, activeJobId = null, pollTimer = null, navigatingToResult = false;
  var activeJobStorageKey = 'mkef-active-job';

  function byId(id) { return document.getElementById(id); }
  function clone(value) {
    if (Array.isArray(value)) return value.map(clone);
    if (value && typeof value === 'object') {
      var result = {}; Object.keys(value).forEach(function (key) { result[key] = clone(value[key]); }); return result;
    }
    return value;
  }
  function clear(element) { while (element.firstChild) element.removeChild(element.firstChild); }
  function numeric(input) { return input.value.trim() === '' ? NaN : Number(input.value); }
  function option(value, label) { var item = document.createElement('option'); item.value = String(value); item.textContent = label; return item; }
  function cell(row) { var td = document.createElement('td'); row.appendChild(td); return td; }
  function inputCell(row, value, type, change, disabled) {
    var input = document.createElement('input'); input.type = type || 'number'; input.step = type === 'text' ? '' : 'any'; input.value = String(value); input.disabled = !!disabled;
    input.addEventListener('change', function () { change(type === 'text' ? input.value : numeric(input)); }); cell(row).appendChild(input); return input;
  }
  function selectCell(row, value, choices, change) {
    var select = document.createElement('select'); choices.forEach(function (choice) { select.appendChild(option(choice[0], choice[1])); }); select.value = String(value); select.addEventListener('change', function () { change(Number(select.value)); }); cell(row).appendChild(select); return select;
  }
  function idCell(row, value) { var td = cell(row); td.textContent = String(value); }
  function deleteCell(row, action, label) { var button = document.createElement('button'); button.type = 'button'; button.textContent = 'Delete'; button.setAttribute('aria-label', 'Delete ' + label); button.addEventListener('click', action); cell(row).appendChild(button); }
  function snapshot() { return clone(state); }

  function commit(operation) {
    undoStack.push(snapshot()); if (undoStack.length > 100) undoStack.shift(); redoStack = [];
    operation(); dirty = true; render();
  }

  function setNotice(message) { elements.notice.textContent = message || ''; elements.notice.hidden = !message; }
  function setError(message) { elements.error.textContent = message || ''; elements.error.hidden = !message; }

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
    elements.runLog.textContent = job.logTail || 'No solver output yet.';
    elements.cancelRun.disabled = ['queued', 'running'].indexOf(job.status) < 0;
  }
  function clearPollTimer() { if (pollTimer) { global.clearTimeout(pollTimer); pollTimer = null; } }
  function finishActiveJob(job) {
    clearPollTimer(); storeJobId(null); activeJobId = null; render();
    if (job.status === 'failed' || job.status === 'timed_out') setError(job.error && job.error.message ? job.error.message : 'The calculation failed.');
  }

  async function pollActiveJob() {
    if (!activeJobId) return;
    try {
      var job = await global.MKEFApi.getJob(activeJobId);
      updateRunPanel(job);
      if (job.status === 'succeeded') {
        var completedId = activeJobId;
        clearPollTimer(); storeJobId(null); activeJobId = null; navigatingToResult = true;
        global.location.assign('../postprocessor/?job=' + encodeURIComponent(completedId));
        return;
      }
      if (['failed', 'canceled', 'timed_out'].indexOf(job.status) >= 0) { finishActiveJob(job); return; }
      pollTimer = global.setTimeout(pollActiveJob, 1000);
    } catch (error) {
      if (error.status === 404) { clearPollTimer(); storeJobId(null); activeJobId = null; render(); }
      setError(error.message);
    }
  }

  async function submitCurrentCase() {
    if (!apiReady || activeJobId) return;
    var caseText;
    try { caseText = F.serialize(state); } catch (error) { setError(error.message); return; }
    elements.run.disabled = true; setError('');
    elements.runPanel.hidden = false; elements.runStatus.textContent = 'Submitting'; elements.runElapsed.textContent = '0.0 s'; elements.runLog.textContent = 'No solver output yet.'; elements.cancelRun.disabled = true;
    try {
      var job = await global.MKEFApi.createJob(state.name || 'Case', caseText);
      activeJobId = job.id; storeJobId(activeJobId); updateRunPanel(job); pollActiveJob();
    } catch (error) { elements.runPanel.hidden = true; render(); setError(error.message); }
  }

  async function cancelActiveJob() {
    if (!activeJobId) return;
    elements.cancelRun.disabled = true;
    try { var job = await global.MKEFApi.cancelJob(activeJobId); updateRunPanel(job); finishActiveJob(job); }
    catch (error) { setError(error.message); elements.cancelRun.disabled = false; }
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
      if (candidate && global.MKEFApi.isValidJobId(candidate)) { activeJobId = candidate; pollActiveJob(); }
    } catch (error) {
      apiReady = false; elements.run.title = error.message;
      elements.solverMessage.textContent = 'Local solver unavailable: ' + error.message;
      elements.solverMessage.hidden = false;
    }
    render();
  }

  function currentNode() { return selection && selection.kind === 'node' ? selection.index + 1 : 1; }
  function placementPoint(point) {
    var x = point.x, y = point.y;
    if (elements.snapEnabled.checked) {
      var spacing = numeric(elements.gridSpacing);
      if (Number.isFinite(spacing) && spacing > 0) { x = Math.round(x / spacing) * spacing; y = Math.round(y / spacing) * spacing; }
    }
    return { x: Number(x.toPrecision(12)), y: Number(y.toPrecision(12)) };
  }
  function defaultElement(node1, node2) {
    return { node1: node1, node2: node2, area: numeric(elements.defaultArea), youngsModulus: numeric(elements.defaultYoung), density: numeric(elements.defaultDensity), momentOfInertia: state.elementType === 113 ? numeric(elements.defaultInertia) : null };
  }

  function addNode(point) {
    var placed = placementPoint(point);
    commit(function () { state.nodes.push(placed); selection = { kind: 'node', index: state.nodes.length - 1 }; });
  }

  function addElement(firstIndex, secondIndex) {
    if ([112, 113].indexOf(state.elementType) < 0) { setNotice('Choose an element family before creating the element.'); return false; }
    commit(function () { state.elements.push(defaultElement(firstIndex + 1, secondIndex + 1)); selection = { kind: 'element', index: state.elements.length - 1 }; });
    return true;
  }

  function activateCanvasTool(tool, pendingNode) {
    document.querySelectorAll('[data-tool]').forEach(function (button) { button.classList.toggle('active', button.getAttribute('data-tool') === tool); });
    renderer.setTool(tool, pendingNode);
  }

  function selectItem(next) { selection = next; render(); }

  function deleteNode(index) {
    var id = index + 1, dependencies = [];
    if (state.elements.some(function (item) { return item.node1 === id || item.node2 === id; })) dependencies.push('elements');
    if (state.supports.some(function (item) { return item.node === id; })) dependencies.push('supports');
    if (state.loads.some(function (item) { return item.node === id; })) dependencies.push('loads');
    if (state.analysis.type === 'transient' && state.analysis.monitorNode === id) dependencies.push('transient monitor');
    if (dependencies.length) { setNotice('Node ' + id + ' is still referenced by ' + dependencies.join(', ') + '. Remove those references first.'); return; }
    commit(function () {
      state.nodes.splice(index, 1);
      state.elements.forEach(function (item) { if (item.node1 > id) item.node1 -= 1; if (item.node2 > id) item.node2 -= 1; });
      state.supports.forEach(function (item) { if (item.node > id) item.node -= 1; });
      state.loads.forEach(function (item) { if (item.node > id) item.node -= 1; });
      if (state.analysis.monitorNode > id) state.analysis.monitorNode -= 1;
      selection = null;
    });
  }

  function renderNodes() {
    clear(elements.nodesBody);
    state.nodes.forEach(function (node, index) {
      var row = document.createElement('tr'); if (selection && selection.kind === 'node' && selection.index === index) row.className = 'selected-row'; row.addEventListener('click', function (event) { if (['INPUT', 'SELECT', 'BUTTON'].indexOf(event.target.tagName) < 0) selectItem({ kind: 'node', index: index }); }); idCell(row, index + 1);
      inputCell(row, node.x, 'number', function (value) { commit(function () { state.nodes[index].x = value; }); });
      inputCell(row, node.y, 'number', function (value) { commit(function () { state.nodes[index].y = value; }); });
      deleteCell(row, function () { deleteNode(index); }, 'node ' + (index + 1)); elements.nodesBody.appendChild(row);
    });
    elements.nodeCount.textContent = '(' + state.nodes.length + ')';
  }

  function renderElements() {
    clear(elements.elementsHead); clear(elements.elementsBody);
    var header = document.createElement('tr'), headings = ['ID', 'Node 1', 'Node 2', 'A', 'E', 'rho']; if (state.elementType === 113) headings.push('I'); headings.push('');
    headings.forEach(function (heading) { var th = document.createElement('th'); th.textContent = heading; header.appendChild(th); }); elements.elementsHead.appendChild(header);
    state.elements.forEach(function (item, index) {
      var row = document.createElement('tr'); if (selection && selection.kind === 'element' && selection.index === index) row.className = 'selected-row'; row.addEventListener('click', function (event) { if (['INPUT', 'SELECT', 'BUTTON'].indexOf(event.target.tagName) < 0) selectItem({ kind: 'element', index: index }); }); idCell(row, index + 1);
      ['node1', 'node2', 'area', 'youngsModulus', 'density'].forEach(function (field) { inputCell(row, item[field], 'number', function (value) { commit(function () { state.elements[index][field] = value; }); }); });
      if (state.elementType === 113) inputCell(row, item.momentOfInertia, 'number', function (value) { commit(function () { state.elements[index].momentOfInertia = value; }); });
      deleteCell(row, function () { commit(function () { state.elements.splice(index, 1); selection = null; }); }, 'element ' + (index + 1)); elements.elementsBody.appendChild(row);
    });
    elements.elementCount.textContent = '(' + state.elements.length + ')';
  }

  function renderSupports() {
    clear(elements.supportsBody);
    state.supports.forEach(function (support, index) {
      var row = document.createElement('tr'); idCell(row, index + 1);
      inputCell(row, support.node, 'number', function (value) { commit(function () { state.supports[index].node = value; }); });
      var supportDOFs = state.elementType === 112 ? 2 : 3;
      selectCell(row, support.type, [1, 2, 3, 4].map(function (type) { return [type, type + ' · ' + global.MKEFSupportMarkers.label(type, supportDOFs)]; }), function (value) { commit(function () { state.supports[index].type = value; }); });
      deleteCell(row, function () { commit(function () { state.supports.splice(index, 1); }); }, 'support ' + (index + 1)); elements.supportsBody.appendChild(row);
    });
    elements.supportCount.textContent = '(' + state.supports.length + ')';
  }

  function renderLoads() {
    clear(elements.loadsBody);
    state.loads.forEach(function (load, index) {
      var row = document.createElement('tr'); idCell(row, index + 1);
      var choices = state.analysis.type === 'static' ? [[10, '10 · static']] : [[10, '10 · legacy pulse'], [11, '11 · harmonic'], [12, '12 · pulse'], [13, '13 · step']];
      if (!choices.some(function (choice) { return choice[0] === load.type; })) choices.push([load.type, load.type + ' · incompatible']);
      selectCell(row, load.type, choices, function (value) { commit(function () { state.loads[index].type = value; if (value === 11 && !Number.isFinite(state.loads[index].frequency)) state.loads[index].frequency = 1; }); });
      ['node', 'fx', 'fy', 'mz'].forEach(function (field) { inputCell(row, load[field], 'number', function (value) { commit(function () { state.loads[index][field] = value; }); }, field === 'mz' && state.elementType === 112); });
      inputCell(row, load.type === 11 ? load.frequency : '', 'number', function (value) { commit(function () { state.loads[index].frequency = value; }); }, load.type !== 11);
      deleteCell(row, function () { commit(function () { state.loads.splice(index, 1); }); }, 'load ' + (index + 1)); elements.loadsBody.appendChild(row);
    });
    elements.loadCount.textContent = '(' + state.loads.length + ')';
  }

  function renderSelection() {
    elements.nodeCoordinateEditor.hidden = true;
    if (!selection) elements.selectionDetails.textContent = 'Nothing selected';
    else if (selection.kind === 'node' && state.nodes[selection.index]) {
      var node = state.nodes[selection.index]; elements.selectionDetails.textContent = 'Node ' + (selection.index + 1); elements.nodeCoordinateEditor.hidden = false; elements.selectedNodeX.value = node.x; elements.selectedNodeY.value = node.y;
    } else if (selection.kind === 'element' && state.elements[selection.index]) {
      var item = state.elements[selection.index]; elements.selectionDetails.textContent = 'Element ' + (selection.index + 1) + '\nNodes ' + item.node1 + ' — ' + item.node2 + '\nType ' + state.elementType;
    } else { selection = null; elements.selectionDetails.textContent = 'Nothing selected'; }
  }

  function renderAnalysis() {
    elements.analysisType.value = state.analysis.type; elements.elementType.value = String(state.elementType); elements.caseName.value = state.name;
    elements.transientSettings.hidden = state.analysis.type !== 'transient';
    elements.timeStep.value = state.analysis.timeStep; elements.duration.value = state.analysis.duration; elements.monitorNode.value = state.analysis.monitorNode;
    clear(elements.monitorDOF); var labels = state.elementType === 112 ? ['ux', 'uy'] : (state.elementType === 113 ? ['ux', 'uy', 'thetaZ'] : []); labels.forEach(function (label, index) { elements.monitorDOF.appendChild(option(index + 1, (index + 1) + ' · ' + label)); }); elements.monitorDOF.value = String(state.analysis.monitorDOF);
    elements.defaultInertia.hidden = state.elementType !== 113; elements.defaultInertiaLabel.hidden = state.elementType !== 113;
    elements.addLoad.disabled = ['static', 'transient'].indexOf(state.analysis.type) < 0;
  }

  function render() {
    setNotice(''); renderAnalysis(); renderNodes(); renderElements(); renderSupports(); renderLoads(); renderSelection();
    var gridSize = numeric(elements.gridSpacing), validGridSize = Number.isFinite(gridSize) && gridSize > 0;
    elements.gridSpacing.setCustomValidity(validGridSize ? '' : 'Grid size must be positive and finite.');
    renderer.setGrid(elements.gridEnabled.checked, validGridSize ? gridSize : 1);
    renderer.selection = selection; renderer.draw(state);
    var errors = F.validate(state), text = '';
    if (!errors.length) { text = F.serialize(state); elements.validationStatus.textContent = 'Valid ' + state.analysis.type + ' case · ready to download'; elements.validationStatus.className = 'status-valid'; elements.download.disabled = false; }
    else { text = '# Resolve validation errors to generate a case file.\n'; elements.validationStatus.textContent = errors.length + ' validation issue' + (errors.length === 1 ? '' : 's'); elements.validationStatus.className = ''; elements.download.disabled = true; }
    elements.preview.value = text; setError(errors.slice(0, 6).join('\n') + (errors.length > 6 ? '\n…' : ''));
    elements.run.disabled = errors.length > 0 || !apiReady || !!activeJobId;
    elements.undo.disabled = !undoStack.length; elements.redo.disabled = !redoStack.length; elements.caseStatus.textContent = dirty ? state.name + ' · modified' : state.name;
  }

  function confirmDiscard() { return !dirty || global.confirm('Discard the current unsaved changes?'); }

  function loadText(text, filename) {
    try {
      var parsed = F.parse(text); parsed.name = (filename || 'Case').replace(/\.txt$/i, '') || 'Case';
      state = parsed; undoStack = []; redoStack = []; dirty = false; selection = null; setError(''); render(); renderer.fit(); renderer.draw(state);
      if (!state.analysis.type) setNotice('This legacy case has no analysis section. Choose a task before exporting.');
    } catch (error) { setError(error.message); }
  }

  function openFile(file) {
    if (!file || !confirmDiscard()) return;
    if (file.size > M.config.maxFileBytes) { setError('The selected file exceeds the 10 MiB limit.'); return; }
    var reader = new FileReader(); reader.onload = function () { loadText(String(reader.result), file.name); }; reader.onerror = function () { setError('The selected file could not be read.'); }; reader.readAsText(file);
  }

  function bind() {
    elements = {
      caseStatus: byId('case-status'), fileInput: byId('file-input'), download: byId('download-case'), undo: byId('undo'), redo: byId('redo'),
      run: byId('run-analysis'), solverMessage: byId('solver-message'), runPanel: byId('run-panel'), runStatus: byId('run-status'), runElapsed: byId('run-elapsed'), runLog: byId('run-log'), cancelRun: byId('cancel-run'),
      error: byId('error-panel'), notice: byId('notice-panel'), caseName: byId('case-name'), analysisType: byId('analysis-type'), elementType: byId('element-type'),
      transientSettings: byId('transient-settings'), timeStep: byId('time-step'), duration: byId('duration'), monitorNode: byId('monitor-node'), monitorDOF: byId('monitor-dof'),
      gridEnabled: byId('grid-enabled'), snapEnabled: byId('snap-enabled'), gridSpacing: byId('grid-spacing'), defaultArea: byId('default-area'), defaultYoung: byId('default-young'), defaultDensity: byId('default-density'), defaultInertia: byId('default-inertia'), defaultInertiaLabel: byId('default-inertia-label'),
      selectionDetails: byId('selection-details'), nodeCoordinateEditor: byId('node-coordinate-editor'), selectedNodeX: byId('selected-node-x'), selectedNodeY: byId('selected-node-y'), nodesBody: byId('nodes-body'), elementsHead: byId('elements-head'), elementsBody: byId('elements-body'), supportsBody: byId('supports-body'), loadsBody: byId('loads-body'), nodeCount: byId('node-count'), elementCount: byId('element-count'), supportCount: byId('support-count'), loadCount: byId('load-count'), preview: byId('case-preview'), validationStatus: byId('validation-status'), addLoad: byId('add-load')
    };
    renderer = new M.Renderer(byId('viewport'), { addNode: addNode, addElement: addElement, select: selectItem, placementPoint: placementPoint, elementPending: function () { setNotice('Element Add mode: select second node for element'); } });
    state = F.newModel('', 0); render(); renderer.fit(); renderer.draw(state);
    byId('new-case').addEventListener('click', function () { if (!confirmDiscard()) return; state = F.newModel('', 0); undoStack = []; redoStack = []; dirty = false; selection = null; render(); renderer.fit(); });
    elements.fileInput.addEventListener('change', function () { openFile(elements.fileInput.files[0]); elements.fileInput.value = ''; });
    elements.run.addEventListener('click', submitCurrentCase);
    elements.cancelRun.addEventListener('click', cancelActiveJob);
    elements.download.addEventListener('click', function () { try { var text = F.serialize(state), blob = new Blob([text], { type: 'text/plain;charset=utf-8' }), url = URL.createObjectURL(blob), link = document.createElement('a'); link.href = url; link.download = (state.name.replace(/[^A-Za-z0-9._ -]/g, '_') || 'Case') + '.txt'; link.click(); setTimeout(function () { URL.revokeObjectURL(url); }, 0); dirty = false; render(); } catch (error) { setError(error.message); } });
    elements.undo.addEventListener('click', function () { if (!undoStack.length) return; redoStack.push(snapshot()); state = undoStack.pop(); selection = null; dirty = true; render(); });
    elements.redo.addEventListener('click', function () { if (!redoStack.length) return; undoStack.push(snapshot()); state = redoStack.pop(); selection = null; dirty = true; render(); });
    byId('reset-view').addEventListener('click', function () { renderer.fit(); renderer.draw(state); });
    elements.caseName.addEventListener('change', function () { state.name = elements.caseName.value.trim() || 'Case'; dirty = true; render(); });
    elements.analysisType.addEventListener('change', function () { var value = elements.analysisType.value; commit(function () { state.analysis.type = value; }); });
    elements.elementType.addEventListener('change', function () { var value = Number(elements.elementType.value); if (state.elements.length) { global.alert('Delete all elements before changing the element family.'); elements.elementType.value = String(state.elementType); return; } if (value === 112 && (state.loads.some(function (load) { return load.mz !== 0; }) || (state.analysis.type === 'transient' && state.analysis.monitorDOF > 2))) { global.alert('Set all load moments to zero and choose monitor DOF 1 or 2 before changing to a truss.'); elements.elementType.value = String(state.elementType); return; } commit(function () { state.elementType = value; }); });
    [['timeStep', elements.timeStep], ['duration', elements.duration], ['monitorNode', elements.monitorNode]].forEach(function (binding) { binding[1].addEventListener('change', function () { var value = numeric(binding[1]); commit(function () { state.analysis[binding[0]] = value; }); }); });
    elements.monitorDOF.addEventListener('change', function () { commit(function () { state.analysis.monitorDOF = Number(elements.monitorDOF.value); }); });
    elements.gridEnabled.addEventListener('change', function () { renderer.setGrid(elements.gridEnabled.checked, numeric(elements.gridSpacing)); });
    elements.gridSpacing.addEventListener('change', function () { var value = numeric(elements.gridSpacing), valid = Number.isFinite(value) && value > 0; elements.gridSpacing.setCustomValidity(valid ? '' : 'Grid size must be positive and finite.'); if (valid) { renderer.setGrid(elements.gridEnabled.checked, value); renderer.refreshNodePreview(); } });
    elements.snapEnabled.addEventListener('change', function () { renderer.refreshNodePreview(); });
    elements.selectedNodeX.addEventListener('change', function () { if (!selection || selection.kind !== 'node') return; var index = selection.index, value = numeric(elements.selectedNodeX); commit(function () { state.nodes[index].x = value; }); });
    elements.selectedNodeY.addEventListener('change', function () { if (!selection || selection.kind !== 'node') return; var index = selection.index, value = numeric(elements.selectedNodeY); commit(function () { state.nodes[index].y = value; }); });
    document.querySelectorAll('[data-tool]').forEach(function (button) { button.addEventListener('click', function () { activateCanvasTool(button.getAttribute('data-tool')); }); });
    byId('add-node').addEventListener('click', function () { addNode({ x: 0, y: 0 }); });
    byId('add-element').addEventListener('click', function () {
      if (state.nodes.length < 2) { setNotice('Add at least two nodes first.'); return; }
      if (selection && selection.kind === 'node') {
        activateCanvasTool('member', selection.index);
        setNotice('Element Add mode: select second node for element');
      } else {
        activateCanvasTool('member', null);
        setNotice('Element Add mode: select first node for element');
      }
    });
    byId('add-support').addEventListener('click', function () { commit(function () { state.supports.push({ type: 1, node: currentNode() }); }); });
    elements.addLoad.addEventListener('click', function () { if (state.analysis.type === 'modal') return; commit(function () { state.loads.push({ type: state.analysis.type === 'static' ? 10 : 13, node: currentNode(), fx: 0, fy: -1, mz: 0, frequency: null }); }); });
    var drop = byId('drop-zone'); drop.addEventListener('dragover', function (event) { event.preventDefault(); }); drop.addEventListener('drop', function (event) { event.preventDefault(); if (event.dataTransfer.files.length) openFile(event.dataTransfer.files[0]); });
    global.addEventListener('beforeunload', function (event) { if (!dirty || navigatingToResult) return; event.preventDefault(); event.returnValue = ''; });
    configureSolver();
  }

  document.addEventListener('DOMContentLoaded', bind);
}(window));
