(function (M) {
  'use strict';
  var elements = {}, dataset = null, renderer = null, currentVector = null, currentScale = 1, transientScale = 1, playback = null;
  function byId(id) { return document.getElementById(id); }
  function option(value, label) { var item = document.createElement('option'); item.value = value; item.textContent = label; return item; }
  function replaceOptions(select, values) { while (select.firstChild) select.removeChild(select.firstChild); values.forEach(function (item) { select.appendChild(option(item.value, item.label)); }); }
  function number(value) { if (value === 0) return '0'; return Math.abs(value) >= 10000 || Math.abs(value) < .001 ? value.toExponential(5) : Number(value.toPrecision(7)).toString(); }
  function samples() { var value = Number(elements.sampleCount.value); if (!Number.isInteger(value) || value < M.config.minFrameSamples || value > M.config.maxFrameSamples) throw new Error('Frame samples must be an integer from 5 to 201.'); return value; }
  function manualScale() { var value = Number(elements.manualScale.value); if (!Number.isFinite(value) || value <= 0) throw new Error('Manual deformation scale must be a finite positive value.'); return value; }
  function layerSettings() {
    return { showOriginal: elements.showOriginal.checked, showDeformed: elements.showDeformed.checked, showNodes: elements.showNodes.checked,
      showNodeLabels: elements.showNodeLabels.checked, showElementLabels: elements.showElementLabels.checked, showSupports: elements.showSupports.checked,
      showLoads: elements.showLoads.checked, showReactions: elements.showReactions.checked };
  }
  function showError(error) { elements.error.textContent = error && error.message ? error.message : String(error); elements.error.hidden = false; }
  function clearError() { elements.error.hidden = true; }
  function notice(message) { elements.notice.textContent = message || ''; elements.notice.hidden = !message; }
  function setLegend(message) { elements.legend.textContent = message || ''; elements.legend.hidden = !message; }
  function setStatus(message) { elements.status.textContent = message; }
  function setSection(section, visible) { section.hidden = !visible; }

  function updateGeometry() {
    if (!dataset || !currentVector) return;
    try { renderer.updateDeformation(currentVector, currentScale, samples()); updateDetails(renderer.selection); }
    catch (error) { showError(error); }
  }
  function updateScale() {
    if (!dataset || !currentVector) return;
    var type = dataset.raw.analysis.type;
    if (type === 'modal') {
      var amplitude = Number(elements.modalAmplitude.value);
      if (!Number.isFinite(amplitude) || amplitude <= 0) { showError(new Error('Modal display multiplier must be positive.')); return; }
      var modal = M.modalView.display(dataset, Number(elements.modeNumber.value), amplitude, samples());
      currentVector = modal.vector; currentScale = modal.scale;
    } else if (elements.scaleMode.value === 'manual') currentScale = manualScale();
    else if (type === 'transient') currentScale = transientScale;
    else currentScale = M.geometry.automaticScale(dataset, currentVector, samples()).scale;
    updateGeometry();
  }

  function configureStatic() {
    var type = dataset.raw.model.elements[0].type;
    replaceOptions(elements.staticResult, type === 112 ? [
      { value: 'none', label: 'None' }, { value: 'axialForce', label: 'Axial force' }, { value: 'axialStress', label: 'Axial stress' }
    ] : [{ value: 'none', label: 'None' }, { value: 'N', label: 'Axial force N' }, { value: 'V', label: 'Shear force V' }, { value: 'M', label: 'Bending moment M' }]);
    M.staticResults.mount(renderer, dataset);
    currentVector = dataset.raw.analysis.displacements;
    currentScale = M.geometry.automaticScale(dataset, currentVector, samples()).scale;
    setSection(elements.staticControls, true);
    elements.scaleMode.disabled = false; elements.manualScale.disabled = elements.scaleMode.value !== 'manual';
    document.querySelectorAll('.static-layer').forEach(function (item) { item.hidden = false; });
    elements.showLoads.disabled = false; elements.showReactions.disabled = false;
    setLegend(''); notice('');
  }

  function configureModal() {
    replaceOptions(elements.modeNumber, dataset.raw.analysis.frequenciesHz.map(function (frequency, index) { return { value: String(index), label: 'Mode ' + (index + 1) + ' · ' + number(frequency) + ' Hz' }; }));
    elements.modeNumber.value = '0';
    currentVector = M.modalView.displacement(dataset, 0);
    var modal = M.modalView.display(dataset, 0, Number(elements.modalAmplitude.value), samples()); currentScale = modal.scale;
    elements.modeFrequency.textContent = 'f = ' + number(modal.frequencyHz) + ' Hz · ω = ' + number(modal.angularFrequency) + ' rad/s';
    setSection(elements.modalControls, true);
    elements.scaleMode.disabled = true; elements.manualScale.disabled = true;
    document.querySelectorAll('.static-layer').forEach(function (item) { item.hidden = true; });
    elements.showLoads.checked = false; elements.showReactions.checked = false;
    setLegend(modal.zero ? 'This mode has zero translational components; no displaced shape is drawn.' : 'Mode shapes are normalized to 10% of model size before applying the display multiplier.'); notice('');
  }

  function updateHistoryDofOptions(preferredId) {
    if (!dataset || dataset.raw.analysis.type !== 'transient') return;
    var nodeId = Number(elements.historyNode.value), analysis = dataset.raw.analysis;
    var ids = analysis.globalDOFIds.filter(function (id) { return dataset.dofById.get(id).nodeId === nodeId; });
    replaceOptions(elements.historyDof, ids.map(function (id) { var location = dataset.dofById.get(id); return { value: String(id), label: location.label + ' · global DOF ' + id }; }));
    if (preferredId && ids.indexOf(Number(preferredId)) >= 0) elements.historyDof.value = String(preferredId);
  }

  function configureTransient() {
    var analysis = dataset.raw.analysis;
    elements.timeIndex.max = String(analysis.time.length - 1); elements.timeIndex.value = '0';
    var nodeIds = [];
    analysis.globalDOFIds.forEach(function (id) { var nodeId = dataset.dofById.get(id).nodeId; if (nodeIds.indexOf(nodeId) < 0) nodeIds.push(nodeId); });
    replaceOptions(elements.historyNode, nodeIds.map(function (id) { return { value: String(id), label: 'Node ' + id }; }));
    updateHistoryDofOptions();
    var quantities = M.transientView.availableQuantities(dataset);
    replaceOptions(elements.historyQuantity, quantities.map(function (name) { return { value: name, label: name === 'spectrum' ? 'Displacement spectrum' : name }; }));
    transientScale = M.transientView.displayScale(dataset);
    var full = M.transientView.hasFullDisplacements(dataset);
    currentVector = full ? M.transientView.vectorAt(dataset, 'displacements', 0) : null;
    currentScale = elements.scaleMode.value === 'manual' ? manualScale() : transientScale;
    elements.scaleMode.disabled = !full; elements.manualScale.disabled = !full || elements.scaleMode.value !== 'manual';
    elements.showDeformed.disabled = !full;
    if (!full) { elements.showDeformed.checked = false; notice('This reduced-DOF export supports charts only. Structural playback requires displacement histories for every global DOF.'); }
    else notice('');
    setSection(elements.transientControls, true); setSection(elements.chartControls, quantities.length > 0); elements.chartPanel.hidden = quantities.length === 0;
    document.querySelectorAll('.static-layer').forEach(function (item) { item.hidden = true; });
    elements.showLoads.checked = false; elements.showReactions.checked = false;
    setTransientIndex(0, false); setLegend('');
  }

  function resetSections() {
    setSection(elements.staticControls, false); setSection(elements.modalControls, false); setSection(elements.transientControls, false); setSection(elements.chartControls, false);
    elements.chartPanel.hidden = true; M.charts.render(elements.chart, null); elements.showDeformed.disabled = false; elements.scaleMode.disabled = false;
  }

  function acceptData(parsed) {
    stopPlayback(); clearError();
    var next = M.validateDataset(parsed); dataset = next;
    elements.title.textContent = parsed.metadata.title || 'Untitled dataset'; elements.dropZone.classList.add('has-data');
    elements.reset.disabled = false; elements.exportSvg.disabled = false; elements.exportPng.disabled = false;
    resetSections(); renderer.mount(dataset);
    if (parsed.analysis.type === 'static') configureStatic();
    else if (parsed.analysis.type === 'modal') configureModal();
    else configureTransient();
    renderer.setVisibility(layerSettings()); updateGeometry(); renderer.setSelection(null);
    if (parsed.model.elements.length > 2000) notice((elements.notice.hidden ? '' : elements.notice.textContent + '\n') + 'This model exceeds the interactive target of 2,000 elements.');
    updateStatus();
  }

  function readFile(file) {
    clearError(); if (!file) return;
    if (file.size > M.config.maxFileBytes) { showError(new Error('The selected file exceeds the 100 MiB limit.')); return; }
    var reader = new FileReader();
    reader.onerror = function () { showError(new Error('The selected file could not be read.')); };
    reader.onload = function () { try { acceptData(JSON.parse(reader.result)); } catch (error) { showError(error instanceof SyntaxError ? new Error('Invalid JSON: ' + error.message) : error); } };
    reader.readAsText(file, 'UTF-8');
  }

  function applyStaticResult() { if (!dataset || dataset.raw.analysis.type !== 'static') return; setLegend(M.staticResults.applyResult(renderer, dataset, elements.staticResult.value)); renderer.reapplySelection(); updateDetails(renderer.selection); }
  function changeMode() {
    if (!dataset || dataset.raw.analysis.type !== 'modal') return;
    var index = Number(elements.modeNumber.value), modal = M.modalView.display(dataset, index, Number(elements.modalAmplitude.value), samples());
    currentVector = modal.vector; currentScale = modal.scale;
    elements.modeFrequency.textContent = 'f = ' + number(modal.frequencyHz) + ' Hz · ω = ' + number(modal.angularFrequency) + ' rad/s';
    setLegend(modal.zero ? 'This mode has zero translational components; no displaced shape is drawn.' : 'Mode shapes are normalized to 10% of model size before applying the display multiplier.');
    updateGeometry(); updateStatus();
  }

  function setTransientIndex(index, stop) {
    if (!dataset || dataset.raw.analysis.type !== 'transient') return;
    if (stop !== false) stopPlayback();
    var analysis = dataset.raw.analysis, value = Math.max(0, Math.min(analysis.time.length - 1, Number(index)));
    elements.timeIndex.value = String(value);
    if (M.transientView.hasFullDisplacements(dataset)) { currentVector = M.transientView.vectorAt(dataset, 'displacements', value); updateGeometry(); }
    elements.timeValue.textContent = 'sample ' + (value + 1) + ' / ' + analysis.time.length + ' · t = ' + number(analysis.time[value]) + (dataset.raw.metadata.units.time ? ' ' + dataset.raw.metadata.units.time : '');
    updateChart(); updateDetails(renderer.selection); updateStatus();
  }

  function updateChart() {
    if (!dataset || dataset.raw.analysis.type !== 'transient' || !elements.historyDof.value || !elements.historyQuantity.value) return;
    var quantity = elements.historyQuantity.value, series = M.transientView.series(dataset, Number(elements.historyDof.value), quantity);
    var cursor = quantity === 'spectrum' ? NaN : dataset.raw.analysis.time[Number(elements.timeIndex.value)];
    M.charts.render(elements.chart, series, cursor);
  }

  function stopPlayback() {
    if (playback) cancelAnimationFrame(playback.request); playback = null;
    if (elements.playPause) elements.playPause.textContent = 'Play';
  }
  function startPlayback() {
    if (!dataset || dataset.raw.analysis.type !== 'transient' || dataset.raw.analysis.time.length < 2) return;
    if (playback) { stopPlayback(); return; }
    var count = dataset.raw.analysis.time.length, startIndex = Number(elements.timeIndex.value);
    if (startIndex >= count - 1) { startIndex = 0; setTransientIndex(0, false); }
    playback = { startIndex: startIndex, startTime: null, lastDraw: -Infinity, request: 0 };
    elements.playPause.textContent = 'Pause';
    function tick(timestamp) {
      if (!playback) return;
      if (playback.startTime === null) playback.startTime = timestamp;
      var speed = Number(elements.playbackSpeed.value), startProgress = playback.startIndex / (count - 1);
      var progress = Math.min(1, startProgress + (timestamp - playback.startTime) * speed / 5000);
      if (timestamp - playback.lastDraw >= 1000 / 30 || progress === 1) { playback.lastDraw = timestamp; setTransientIndex(Math.min(count - 1, Math.floor(progress * (count - 1))), false); }
      if (progress >= 1) stopPlayback(); else playback.request = requestAnimationFrame(tick);
    }
    playback.request = requestAnimationFrame(tick);
  }

  function valuesAtNode(vector, nodeId) {
    if (!vector) return null;
    var row = dataset.nodeIndexById.get(nodeId), dofs = dataset.raw.model.dofMap[row];
    return dofs.map(function (id) { return vector[id - 1]; });
  }
  function labeledVector(label, values) { return label + ': [' + values.map(number).join(', ') + ']'; }
  function transientNodeValues(field, nodeId, timeIndex) {
    var analysis = dataset.raw.analysis, row = dataset.nodeIndexById.get(nodeId), dofs = dataset.raw.model.dofMap[row];
    return dofs.map(function (id, localIndex) {
      var exportedRow = analysis.globalDOFIds.indexOf(id);
      var value = analysis[field] && exportedRow >= 0 ? number(analysis[field][exportedRow][timeIndex]) : 'not exported';
      return dataset.raw.model.dofLabels[localIndex] + '=' + value;
    });
  }
  function updateDetails(selection) {
    if (!selection || !dataset) { elements.details.textContent = 'Nothing selected'; return; }
    var analysis = dataset.raw.analysis, lines = [];
    if (selection.kind === 'node') {
      var node = dataset.nodesById.get(selection.id), row = dataset.nodeIndexById.get(selection.id), dofs = dataset.raw.model.dofMap[row];
      lines.push('Node ' + node.id, 'Coordinates: [' + number(node.x) + ', ' + number(node.y) + ']', 'Global DOFs: [' + dofs.join(', ') + ']');
      if (analysis.type === 'static') { lines.push(labeledVector('Displacement', valuesAtNode(analysis.displacements, node.id)), labeledVector('Load', valuesAtNode(analysis.loadVector, node.id)), labeledVector('Reaction', valuesAtNode(analysis.reactions, node.id))); }
      else if (analysis.type === 'modal') lines.push(labeledVector('Mode components', valuesAtNode(currentVector, node.id)));
      else {
        var timeIndex = Number(elements.timeIndex.value);
        if (analysis.displacements) lines.push('Displacement: [' + transientNodeValues('displacements', node.id, timeIndex).join(', ') + ']');
        var available = dofs.filter(function (id) { return analysis.globalDOFIds.indexOf(id) >= 0; });
        if (available.length) {
          var current = Number(elements.historyDof.value);
          elements.historyNode.value = String(node.id); updateHistoryDofOptions(available.indexOf(current) >= 0 ? current : available[0]); updateChart();
        }
      }
    } else {
      var element = dataset.elementsById.get(selection.id);
      lines.push('Element ' + element.id, 'Type: ' + element.type, 'Nodes: [' + element.nodeIds.join(', ') + ']');
      if (analysis.type === 'static') { var result = dataset.resultsByElementId.get(element.id); lines.push(labeledVector('Local end forces', result.localEndForces), 'Axial strain: ' + number(result.axialStrain), 'Axial stress: ' + number(result.axialStress), 'Axial force: ' + number(result.axialForce)); }
      else if (currentVector) lines.push(labeledVector('First-node components', valuesAtNode(currentVector, element.nodeIds[0])), labeledVector('Second-node components', valuesAtNode(currentVector, element.nodeIds[1])));
      else lines.push('Structural displacement components were not exported for this reduced-DOF result.');
    }
    elements.details.textContent = lines.join('\n');
  }

  function updateStatus() {
    if (!dataset) return;
    var type = dataset.raw.analysis.type, text = type + ' · deformation scale ' + number(currentScale);
    if (type === 'modal') text += ' · mode ' + (Number(elements.modeNumber.value) + 1);
    if (type === 'transient') text += ' · sample ' + (Number(elements.timeIndex.value) + 1);
    setStatus(text);
  }
  function exportContext() {
    return { title: dataset.raw.metadata.title, analysisType: dataset.raw.analysis.type, deformationScale: currentScale,
      staticResult: elements.staticResult.value, modeIndex: Number(elements.modeNumber.value || 0), timeIndex: Number(elements.timeIndex.value || 0),
      historyDOF: Number(elements.historyDof.value || 0), historyQuantity: elements.historyQuantity.value || '', exportedUtc: new Date().toISOString() };
  }

  function bind() {
    elements.fileInput.addEventListener('change', function () { readFile(this.files[0]); this.value = ''; });
    ['dragenter', 'dragover'].forEach(function (name) { elements.dropZone.addEventListener(name, function (event) { event.preventDefault(); this.classList.add('dragging'); }); });
    ['dragleave', 'drop'].forEach(function (name) { elements.dropZone.addEventListener(name, function (event) { event.preventDefault(); this.classList.remove('dragging'); }); });
    elements.dropZone.addEventListener('drop', function (event) { readFile(event.dataTransfer.files[0]); });
    elements.dropZone.addEventListener('click', function (event) { if (!dataset && event.target === this) elements.fileInput.click(); });
    elements.dropZone.addEventListener('keydown', function (event) { if (!dataset && (event.key === 'Enter' || event.key === ' ')) elements.fileInput.click(); });
    elements.reset.addEventListener('click', function () { renderer.fit(); });
    elements.exportSvg.addEventListener('click', function () { if (dataset) M.exporting.downloadSvg(elements.viewport, exportContext(), elements.chartPanel.hidden ? null : elements.chart); });
    elements.exportPng.addEventListener('click', function () { if (dataset) M.exporting.downloadPng(elements.viewport, exportContext(), Number(elements.pngScale.value), elements.chartPanel.hidden ? null : elements.chart).catch(showError); });
    ['showOriginal', 'showDeformed', 'showNodes', 'showNodeLabels', 'showElementLabels', 'showSupports', 'showLoads', 'showReactions'].forEach(function (name) { elements[name].addEventListener('change', function () { renderer.setVisibility(layerSettings()); }); });
    elements.scaleMode.addEventListener('change', function () { elements.manualScale.disabled = this.value !== 'manual'; updateScale(); updateStatus(); });
    elements.manualScale.addEventListener('change', function () { if (elements.scaleMode.value === 'manual') { updateScale(); updateStatus(); } });
    elements.sampleCount.addEventListener('change', updateScale);
    elements.staticResult.addEventListener('change', applyStaticResult);
    elements.modeNumber.addEventListener('change', changeMode); elements.modalAmplitude.addEventListener('change', changeMode);
    elements.timeIndex.addEventListener('input', function () { setTransientIndex(Number(this.value)); });
    elements.previousFrame.addEventListener('click', function () { setTransientIndex(Number(elements.timeIndex.value) - 1); });
    elements.nextFrame.addEventListener('click', function () { setTransientIndex(Number(elements.timeIndex.value) + 1); });
    elements.playPause.addEventListener('click', startPlayback);
    elements.historyNode.addEventListener('change', function () { updateHistoryDofOptions(); updateChart(); });
    elements.historyDof.addEventListener('change', updateChart); elements.historyQuantity.addEventListener('change', updateChart);
    document.addEventListener('visibilitychange', function () { if (document.hidden) stopPlayback(); });
  }

  document.addEventListener('DOMContentLoaded', function () {
    ['file-input','reset-view','export-svg','export-png','dataset-title','drop-zone','error-panel','notice-panel','analysis-status','selection-details','viewport','result-legend','chart-panel','history-chart','static-controls','modal-controls','transient-controls','chart-controls','static-result','mode-number','modal-amplitude','mode-frequency','time-index','previous-frame','next-frame','play-pause','playback-speed','time-value','history-node','history-dof','history-quantity','show-original','show-deformed','show-nodes','show-node-labels','show-element-labels','show-supports','show-loads','show-reactions','scale-mode','manual-scale','sample-count','png-scale'].forEach(function (id) { var key = id.replace(/-([a-z])/g, function (_, letter) { return letter.toUpperCase(); }); elements[key] = byId(id); });
    elements.title = elements.datasetTitle; elements.reset = elements.resetView; elements.error = elements.errorPanel; elements.notice = elements.noticePanel; elements.status = elements.analysisStatus; elements.details = elements.selectionDetails; elements.legend = elements.resultLegend; elements.chart = elements.historyChart;
    renderer = new M.Renderer(elements.viewport, updateDetails); bind();
  });
}(window.MKEFPost));
