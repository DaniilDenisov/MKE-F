(function (M) {
  'use strict';
  var elements = {};
  var dataset = null;
  var renderer = null;
  var displayInfo = null;

  function byId(id) { return document.getElementById(id); }
  function settings() {
    var manual = Number(elements.manualScale.value), samples = Number(elements.sampleCount.value);
    if (!Number.isFinite(manual) || manual <= 0) throw new Error('Manual deformation scale must be a finite positive value.');
    if (!Number.isInteger(samples) || samples < M.config.minFrameSamples || samples > M.config.maxFrameSamples) throw new Error('Frame samples must be an integer from 5 to 201.');
    return {
      showOriginal: elements.showOriginal.checked, showDeformed: elements.showDeformed.checked,
      showNodes: elements.showNodes.checked, showNodeLabels: elements.showNodeLabels.checked,
      showElementLabels: elements.showElementLabels.checked, scaleMode: elements.scaleMode.value,
      manualScale: manual, samples: samples,
      showSupports: elements.showSupports.checked, showLoads: elements.showLoads.checked,
      showReactions: elements.showReactions.checked, trussResult: elements.trussResult.value,
      diagram: elements.diagramResult.value
    };
  }
  function showError(error) {
    elements.error.textContent = error && error.message ? error.message : String(error);
    elements.error.hidden = false;
  }
  function clearMessages() { elements.error.hidden = true; elements.notice.hidden = true; }
  function render() {
    if (!dataset) return;
    clearMessages();
    var analysis = dataset.raw.analysis, displacement;
    if (analysis.type === 'static') displacement = analysis.displacements;
    else {
      displacement = new Array(dataset.raw.model.nodes.length * dataset.raw.model.dofPerNode).fill(0);
      elements.notice.textContent = analysis.type.charAt(0).toUpperCase() + analysis.type.slice(1) + ' data is valid. Its interactive view is planned for a later implementation step.';
      elements.notice.hidden = false;
    }
    try {
      var currentSettings = settings();
      var info = renderer.render(dataset, displacement, currentSettings);
      displayInfo = analysis.type === 'static' ? M.staticResults.render(renderer, dataset, currentSettings, info) : { deformationScale: info.scale, zeroDeformation: info.zero };
      elements.status.textContent = analysis.type + ' · deformation scale ' + formatScale(info.scale) + (info.zero ? ' · zero deformation' : '');
    } catch (error) { showError(error); }
  }
  function formatScale(value) { return Number(value).toPrecision(5).replace(/\.0+$/, ''); }
  function acceptData(parsed) {
    var next = M.validateDataset(parsed);
    dataset = next;
    elements.title.textContent = parsed.metadata.title || 'Untitled dataset';
    elements.dropZone.classList.add('has-data');
    elements.reset.disabled = false;
    var staticAvailable = parsed.analysis.type === 'static';
    elements.exportSvg.disabled = !staticAvailable; elements.exportPng.disabled = !staticAvailable;
    elements.trussResult.disabled = !staticAvailable; elements.diagramResult.disabled = !staticAvailable;
    render();
  }
  function readFile(file) {
    clearMessages();
    if (!file) return;
    if (file.size > M.config.maxFileBytes) { showError(new Error('The selected file exceeds the 100 MiB limit.')); return; }
    var reader = new FileReader();
    reader.onerror = function () { showError(new Error('The selected file could not be read.')); };
    reader.onload = function () {
      try { acceptData(JSON.parse(reader.result)); }
      catch (error) { showError(error instanceof SyntaxError ? new Error('Invalid JSON: ' + error.message) : error); }
    };
    reader.readAsText(file, 'UTF-8');
  }
  function bind() {
    elements.fileInput.addEventListener('change', function () { readFile(this.files[0]); this.value = ''; });
    ['dragenter', 'dragover'].forEach(function (name) { elements.dropZone.addEventListener(name, function (event) { event.preventDefault(); this.classList.add('dragging'); }); });
    ['dragleave', 'drop'].forEach(function (name) { elements.dropZone.addEventListener(name, function (event) { event.preventDefault(); this.classList.remove('dragging'); }); });
    elements.dropZone.addEventListener('drop', function (event) { readFile(event.dataTransfer.files[0]); });
    elements.dropZone.addEventListener('click', function (event) { if (!dataset && event.target === this) elements.fileInput.click(); });
    elements.dropZone.addEventListener('keydown', function (event) { if (!dataset && (event.key === 'Enter' || event.key === ' ')) elements.fileInput.click(); });
    elements.reset.addEventListener('click', function () { renderer.fit(); });
    elements.exportSvg.addEventListener('click', function () { if (dataset) M.exporting.downloadSvg(byId('viewport'), exportContext()); });
    elements.exportPng.addEventListener('click', function () {
      if (!dataset) return;
      M.exporting.downloadPng(byId('viewport'), exportContext(), Number(elements.pngScale.value)).catch(showError);
    });
    elements.scaleMode.addEventListener('change', function () { elements.manualScale.disabled = this.value !== 'manual'; render(); });
    ['showOriginal', 'showDeformed', 'showNodes', 'showNodeLabels', 'showElementLabels', 'showSupports', 'showLoads', 'showReactions', 'trussResult', 'diagramResult', 'manualScale', 'sampleCount'].forEach(function (name) { elements[name].addEventListener('change', render); });
  }
  function exportContext() {
    return { title: dataset.raw.metadata.title, analysisType: dataset.raw.analysis.type,
      scales: displayInfo, diagram: elements.diagramResult.value,
      trussResult: elements.trussResult.value, exportedUtc: new Date().toISOString() };
  }
  document.addEventListener('DOMContentLoaded', function () {
    elements = {
      fileInput: byId('file-input'), reset: byId('reset-view'), exportSvg: byId('export-svg'), exportPng: byId('export-png'), title: byId('dataset-title'),
      dropZone: byId('drop-zone'), error: byId('error-panel'), notice: byId('notice-panel'),
      status: byId('analysis-status'), selection: byId('selection-status'),
      showOriginal: byId('show-original'), showDeformed: byId('show-deformed'), showNodes: byId('show-nodes'),
      showNodeLabels: byId('show-node-labels'), showElementLabels: byId('show-element-labels'),
      showSupports: byId('show-supports'), showLoads: byId('show-loads'), showReactions: byId('show-reactions'),
      trussResult: byId('truss-result'), diagramResult: byId('diagram-result'),
      scaleMode: byId('scale-mode'), manualScale: byId('manual-scale'), sampleCount: byId('sample-count'), pngScale: byId('png-scale')
    };
    renderer = new M.Renderer(byId('viewport'), elements.selection);
    bind();
  });
}(window.MKEFPost));
