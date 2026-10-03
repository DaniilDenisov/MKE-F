(function (M) {
  'use strict';
  var elements = {}, dataset = null, renderer = null, currentVector = null, currentScale = 1, transientScale = 1, playback = null, currentLegend = null;
  var modalPresentation = null, modalPhaseAngle = 0, modalPhaseFactor = 1;
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
  function setLegend(descriptor) {
    while (elements.legend.firstChild) elements.legend.removeChild(elements.legend.firstChild);
    if (typeof descriptor === 'string') descriptor = descriptor ? { kind: 'text', text: descriptor, summary: descriptor } : { kind: 'none', summary: '' };
    currentLegend = descriptor || { kind: 'none', summary: '' };
    if (currentLegend.kind === 'none') { elements.legend.hidden = true; return; }
    elements.legend.hidden = false;
    if (currentLegend.kind === 'text') { elements.legend.textContent = currentLegend.text; return; }
    var header = document.createElement('div'); header.className = 'legend-header';
    header.textContent = currentLegend.label + (currentLegend.unit ? ' [' + currentLegend.unit + ']' : '');
    var bar = document.createElement('div'); bar.className = 'legend-bar';
    bar.style.background = 'linear-gradient(90deg, ' + currentLegend.stops.map(function (stop) { return stop.color + ' ' + (stop.offset * 100) + '%'; }).join(', ') + ')';
    var ticks = document.createElement('div'); ticks.className = 'legend-ticks' + (currentLegend.ticks.length === 1 ? ' single' : '');
    currentLegend.ticks.forEach(function (value) { var tick = document.createElement('span'); tick.textContent = number(value); ticks.appendChild(tick); });
    elements.legend.appendChild(header); elements.legend.appendChild(bar); elements.legend.appendChild(ticks);
  }
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
      refreshModalPresentation(false);
      return;
    } else if (elements.scaleMode.value === 'manual') currentScale = manualScale();
    else if (type === 'transient') currentScale = transientScale;
    else currentScale = M.geometry.automaticScale(dataset, currentVector, samples()).scale;
    updateGeometry();
  }

  function configureStatic() {
    var type = dataset.raw.model.elements[0].type;
    replaceOptions(elements.staticResult, type === 112 ? [
      { value: 'displacementMagnitude', label: 'Deflection |u|' }, { value: 'none', label: 'None' }, { value: 'axialForce', label: 'Axial force' }, { value: 'axialStress', label: 'Axial stress' }
    ] : [{ value: 'displacementMagnitude', label: 'Deflection |u|' }, { value: 'none', label: 'None' }, { value: 'N', label: 'Axial force N' }, { value: 'V', label: 'Shear force V' }, { value: 'M', label: 'Bending moment M' }]);
    M.staticResults.mount(renderer, dataset);
    currentVector = dataset.raw.analysis.displacements;
    currentScale = M.geometry.automaticScale(dataset, currentVector, samples()).scale;
    setSection(elements.staticControls, true);
    elements.scaleMode.disabled = false; elements.manualScale.disabled = elements.scaleMode.value !== 'manual';
    document.querySelectorAll('.static-layer').forEach(function (item) { item.hidden = false; });
    elements.showLoads.disabled = false; elements.showReactions.disabled = false;
    setLegend(M.staticResults.applyResult(renderer, dataset, elements.staticResult.value, samples())); notice('');
  }

  function configureModal() {
    replaceOptions(elements.modeNumber, dataset.raw.analysis.frequenciesHz.map(function (frequency, index) { return { value: String(index), label: 'Mode ' + (index + 1) + ' · ' + number(frequency) + ' Hz' }; }));
    elements.modeNumber.value = '0';
    if (!Number.isFinite(Number(elements.modalAmplitude.value)) || Number(elements.modalAmplitude.value) <= 0) elements.modalAmplitude.value = '1';
    setSection(elements.modalControls, true);
    elements.scaleMode.disabled = true; elements.manualScale.disabled = true;
    document.querySelectorAll('.static-layer').forEach(function (item) { item.hidden = true; });
    elements.showLoads.checked = false; elements.showReactions.checked = false;
    modalPhaseAngle = 0; modalPhaseFactor = 1;
    refreshModalPresentation(true);
  }

  function modalMultiplier() {
    var value = Number(elements.modalAmplitude.value);
    if (!Number.isFinite(value) || value <= 0) throw new Error('Modal display multiplier must be a finite positive value.');
    return value;
  }

  function modalLegendText() {
    if (!modalPresentation) return '';
    if (modalPresentation.zero) return 'This mode is identically zero; no displaced shape can be drawn.';
    return 'Mode ' + (Number(elements.modeNumber.value) + 1) + ' · f = ' + number(modalPresentation.frequencyHz) +
      ' Hz · normalized arbitrary amplitude · q = ' + number(modalPhaseFactor) + ' · display multiplier ' + number(modalPresentation.displayMultiplier);
  }

  function updateModalControlState() {
    var index = Number(elements.modeNumber.value), count = dataset.raw.analysis.frequenciesHz.length;
    elements.previousMode.disabled = index <= 0;
    elements.nextMode.disabled = index >= count - 1;
    elements.modalPlayPause.disabled = modalPresentation.zero;
    elements.modalPhase.disabled = modalPresentation.zero;
  }

  function applyModalPhase(factor, angle) {
    if (!modalPresentation) return;
    modalPhaseFactor = Math.max(-1, Math.min(1, Number(factor)));
    if (Number.isFinite(angle)) modalPhaseAngle = angle;
    elements.modalPhase.value = String(modalPhaseFactor);
    elements.modalPhaseValue.textContent = 'q = ' + number(modalPhaseFactor) + ' · illustrative phase, not physical time';
    currentVector = modalPresentation.vector;
    currentScale = modalPresentation.zero ? 0 : modalPresentation.peakScale * modalPhaseFactor;
    elements.modeFrequency.textContent = 'f = ' + number(modalPresentation.frequencyHz) + ' Hz · ω = ' +
      number(modalPresentation.angularFrequency) + ' rad/s · peak scale ' + number(modalPresentation.peakScale) +
      ' · effective scale ' + number(currentScale);
    setLegend(modalLegendText()); updateGeometry(); updateStatus();
  }

  function refreshModalPresentation(resetPhase) {
    if (!dataset || dataset.raw.analysis.type !== 'modal') return false;
    var next;
    try { next = M.modalView.display(dataset, Number(elements.modeNumber.value), modalMultiplier(), samples()); }
    catch (error) { showError(error); return false; }
    clearError(); modalPresentation = next; currentVector = next.vector;
    if (resetPhase) { modalPhaseAngle = 0; modalPhaseFactor = 1; }
    updateModalControlState(); applyModalPhase(modalPhaseFactor, modalPhaseAngle);
    return true;
  }

  function selectModalMode(index) {
    if (!dataset || dataset.raw.analysis.type !== 'modal') return;
    var previousIndex = modalPresentation ? modalPresentation.modeIndex : Number(elements.modeNumber.value);
    var wasPlaying = playback && playback.type === 'modal';
    if (wasPlaying) stopPlayback();
    elements.modeNumber.value = String(M.modalView.clampModeIndex(index, dataset.raw.analysis.frequenciesHz.length));
    if (!refreshModalPresentation(true)) {
      elements.modeNumber.value = String(previousIndex); updateModalControlState();
      if (wasPlaying && !modalPresentation.zero) startModalPlayback();
      return;
    }
    if (wasPlaying && !modalPresentation.zero) startModalPlayback();
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
    var quantityLabels = { displacements: 'Displacement', velocities: 'Velocity', accelerations: 'Acceleration', loadHistory: 'Applied load', reactions: 'Dynamic residual / reaction', spectrum: 'Displacement spectrum' };
    replaceOptions(elements.historyQuantity, quantities.map(function (name) { return { value: name, label: quantityLabels[name] }; }));
    transientScale = M.transientView.displayScale(dataset);
    var full = M.transientView.hasFullDisplacements(dataset);
    currentVector = full ? M.transientView.vectorAt(dataset, 'displacements', 0) : null;
    currentScale = elements.scaleMode.value === 'manual' ? manualScale() : transientScale;
    elements.scaleMode.disabled = !full; elements.manualScale.disabled = !full || elements.scaleMode.value !== 'manual';
    elements.showDeformed.disabled = !full;
    if (!full) { elements.showDeformed.checked = false; notice('This reduced-DOF export supports charts only. Structural playback requires displacement histories for every global DOF.'); }
    else notice('');
    var sampling = M.transientView.samplingSummary(dataset), fields = M.transientView.fieldStatus(dataset), metadata =
      'Exported ' + sampling.exported + ' of ' + sampling.original + ' samples · export stride ' + sampling.stride + '. ' +
      'Available: ' + (fields.present.length ? fields.present.join(', ') : 'none') + '. Omitted: ' + (fields.omitted.length ? fields.omitted.join(', ') : 'none') + '.';
    if (sampling.decimated) metadata = 'Warning: the exported time history is decimated. ' + metadata;
    elements.transientMetadata.textContent = metadata;
    elements.transientMetadata.classList.toggle('warning-note', sampling.decimated);
    elements.playPause.disabled = analysis.time.length < 2;
    setSection(elements.transientControls, true); setSection(elements.chartControls, quantities.length > 0); elements.chartPanel.hidden = quantities.length === 0;
    document.querySelectorAll('.static-layer').forEach(function (item) { item.hidden = true; });
    elements.showLoads.checked = false; elements.showReactions.checked = false;
    elements.historyValue.textContent = '';
    setTransientIndex(0, false); setLegend('');
  }

  function resetSections() {
    setSection(elements.staticControls, false); setSection(elements.modalControls, false); setSection(elements.transientControls, false); setSection(elements.chartControls, false);
    elements.chartPanel.hidden = true; M.charts.render(elements.chart, null); elements.showDeformed.disabled = false; elements.scaleMode.disabled = false;
    elements.playPause.disabled = false; elements.previousFrame.disabled = false; elements.nextFrame.disabled = false;
    elements.transientMetadata.textContent = ''; elements.historyValue.textContent = '';
    currentVector = null; currentScale = 1; modalPresentation = null; modalPhaseAngle = 0; modalPhaseFactor = 1;
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

  function applyStaticResult() { if (!dataset || dataset.raw.analysis.type !== 'static') return; setLegend(M.staticResults.applyResult(renderer, dataset, elements.staticResult.value, samples())); renderer.reapplySelection(); updateDetails(renderer.selection); }
  function changeMode() {
    if (!dataset || dataset.raw.analysis.type !== 'modal') return;
    selectModalMode(Number(elements.modeNumber.value));
  }

  function setTransientIndex(index, stop) {
    if (!dataset || dataset.raw.analysis.type !== 'transient') return;
    if (stop !== false) stopPlayback();
    var analysis = dataset.raw.analysis, value = M.transientView.clampFrameIndex(index, analysis.time.length);
    elements.timeIndex.value = String(value);
    elements.previousFrame.disabled = value <= 0;
    elements.nextFrame.disabled = value >= analysis.time.length - 1;
    if (M.transientView.hasFullDisplacements(dataset)) { currentVector = M.transientView.vectorAt(dataset, 'displacements', value); updateGeometry(); }
    elements.timeValue.textContent = 'sample ' + (value + 1) + ' / ' + analysis.time.length + ' · t = ' + number(analysis.time[value]) + (dataset.raw.metadata.units.time ? ' ' + dataset.raw.metadata.units.time : '');
    updateChart(); updateDetails(renderer.selection); updateStatus();
  }

  function setHistoryValue(series, index, xValue, yValue) {
    if (!series) { elements.historyValue.textContent = ''; return; }
    if (series.quantity === 'spectrum') elements.historyValue.textContent = 'f = ' + number(xValue) + ' Hz · amplitude = ' + number(yValue) + (series.unit ? ' ' + series.unit : '');
    else elements.historyValue.textContent = 'sample ' + (index + 1) + ' · ' + series.yLabel + ' = ' + number(yValue);
  }

  function updateChart() {
    if (!dataset || dataset.raw.analysis.type !== 'transient' || !elements.historyDof.value || !elements.historyQuantity.value) return;
    var quantity = elements.historyQuantity.value, series = M.transientView.series(dataset, Number(elements.historyDof.value), quantity);
    var index = Number(elements.timeIndex.value), cursor = quantity === 'spectrum' ? NaN : dataset.raw.analysis.time[index];
    if (quantity === 'spectrum') elements.historyValue.textContent = 'Select a spectrum sample to inspect its exact frequency and amplitude.';
    else setHistoryValue(series, index, series.x[index], series.y[index]);
    M.charts.render(elements.chart, series, cursor, function (sampleIndex, xValue, yValue) {
      if (quantity === 'spectrum') setHistoryValue(series, sampleIndex, xValue, yValue);
      else setTransientIndex(sampleIndex);
    });
  }

  function stopPlayback() {
    if (playback) cancelAnimationFrame(playback.request); playback = null;
    if (elements.playPause) elements.playPause.textContent = 'Play';
    if (elements.modalPlayPause) elements.modalPlayPause.textContent = 'Play';
  }

  function startModalPlayback() {
    if (!dataset || dataset.raw.analysis.type !== 'modal' || !modalPresentation || modalPresentation.zero) return;
    if (playback && playback.type === 'modal') { stopPlayback(); return; }
    stopPlayback();
    playback = { type: 'modal', startAngle: modalPhaseAngle, startTime: null, lastDraw: -Infinity, request: 0 };
    elements.modalPlayPause.textContent = 'Pause';
    function tick(timestamp) {
      if (!playback || playback.type !== 'modal') return;
      if (playback.startTime === null) playback.startTime = timestamp;
      var phase = M.modalView.phaseAt(playback.startAngle, timestamp - playback.startTime, Number(elements.modalPlaybackSpeed.value));
      if (timestamp - playback.lastDraw >= 1000 / 30) {
        playback.lastDraw = timestamp; applyModalPhase(phase.factor, phase.angle);
      }
      playback.request = requestAnimationFrame(tick);
    }
    playback.request = requestAnimationFrame(tick);
  }

  function startTransientPlayback() {
    if (!dataset || dataset.raw.analysis.type !== 'transient' || dataset.raw.analysis.time.length < 2) return;
    if (playback && playback.type === 'transient') { stopPlayback(); return; }
    stopPlayback();
    var count = dataset.raw.analysis.time.length, startIndex = Number(elements.timeIndex.value);
    if (startIndex >= count - 1) { startIndex = 0; setTransientIndex(0, false); }
    playback = { type: 'transient', startIndex: startIndex, startTime: null, lastDraw: -Infinity, request: 0 };
    elements.playPause.textContent = 'Pause';
    function tick(timestamp) {
      if (!playback || playback.type !== 'transient') return;
      if (playback.startTime === null) playback.startTime = timestamp;
      var frame = M.transientView.playbackFrame(dataset.raw.analysis.time, playback.startIndex, timestamp - playback.startTime,
        Number(elements.playbackSpeed.value), Number(elements.frameStride.value), 5000);
      if (timestamp - playback.lastDraw >= 1000 / 30 || frame.complete) { playback.lastDraw = timestamp; setTransientIndex(frame.index, false); }
      if (frame.complete) stopPlayback(); else playback.request = requestAnimationFrame(tick);
    }
    playback.request = requestAnimationFrame(tick);
  }

  function valuesAtNode(vector, nodeId) {
    if (!vector) return null;
    var row = dataset.nodeIndexById.get(nodeId), dofs = dataset.raw.model.dofMap[row];
    return dofs.map(function (id) { return vector[id - 1]; });
  }
  function labeledVector(label, values) { return label + ': [' + values.map(number).join(', ') + ']'; }
  function labeledModalVector(label, values) {
    return label + ': [' + values.map(function (value, index) { return dataset.raw.model.dofLabels[index] + '=' + number(value); }).join(', ') + ']';
  }
  function modalNodeDetails(nodeId) {
    var normalized = valuesAtNode(modalPresentation.vector, nodeId);
    return [
      labeledModalVector('Normalized shape (arbitrary)', normalized),
      'Phase factor q = ' + number(modalPhaseFactor),
      labeledModalVector('Visible normalized shape', normalized.map(function (value) { return value * modalPhaseFactor; }))
    ];
  }
  function transientNodeValues(field, nodeId, timeIndex) {
    var analysis = dataset.raw.analysis, row = dataset.nodeIndexById.get(nodeId), dofs = dataset.raw.model.dofMap[row];
    return dofs.map(function (id, localIndex) {
      var exportedRow = analysis.globalDOFIds.indexOf(id);
      var descriptor = M.transientView.quantityDescriptor(dataset, id, field);
      var value = analysis[field] && exportedRow >= 0 ? number(analysis[field][exportedRow][timeIndex]) + (descriptor.unit ? ' ' + descriptor.unit : '') : 'not exported';
      return dataset.raw.model.dofLabels[localIndex] + '=' + value;
    });
  }
  function updateDetails(selection) {
    if (!selection || !dataset) { elements.details.textContent = 'Nothing selected'; return; }
    var analysis = dataset.raw.analysis, lines = [];
    if (selection.kind === 'node') {
      var node = dataset.nodesById.get(selection.id), row = dataset.nodeIndexById.get(selection.id), dofs = dataset.raw.model.dofMap[row];
      if (analysis.type === 'static') lines = M.staticResults.nodeDetails(dataset, node.id).concat(M.staticResults.nodeResultDetails(dataset, node.id, elements.staticResult.value));
      else if (analysis.type === 'modal') lines = modalNodeDetails(node.id);
      else {
        var timeIndex = Number(elements.timeIndex.value);
        var fieldLabels = { displacements: 'Displacement', velocities: 'Velocity', accelerations: 'Acceleration', loadHistory: 'Applied load', reactions: 'Dynamic residual / reaction' };
        lines.push('Node ' + node.id, 'Sample ' + (timeIndex + 1) + ' · t = ' + number(analysis.time[timeIndex]) + (dataset.raw.metadata.units.time ? ' ' + dataset.raw.metadata.units.time : ''));
        ['displacements', 'velocities', 'accelerations', 'loadHistory', 'reactions'].forEach(function (field) {
          if (analysis[field]) lines.push(fieldLabels[field] + ': [' + transientNodeValues(field, node.id, timeIndex).join(', ') + ']');
        });
        if (analysis.reactions) lines.push('Dynamic residual is M·a + K·u − F; restrained rows are the support reactions.');
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
      else if (analysis.type === 'modal') {
        lines = lines.concat(modalNodeDetails(element.nodeIds[0]).map(function (line) { return 'Node ' + element.nodeIds[0] + ' · ' + line; }));
        lines = lines.concat(modalNodeDetails(element.nodeIds[1]).map(function (line) { return 'Node ' + element.nodeIds[1] + ' · ' + line; }));
      }
      else if (currentVector) lines.push(labeledVector('First-node components', valuesAtNode(currentVector, element.nodeIds[0])), labeledVector('Second-node components', valuesAtNode(currentVector, element.nodeIds[1])));
      else lines.push('Structural displacement components were not exported for this reduced-DOF result.');
    }
    elements.details.textContent = lines.join('\n');
  }

  function updateStatus() {
    if (!dataset) return;
    var type = dataset.raw.analysis.type, text;
    if (type === 'modal') text = 'modal · mode ' + (Number(elements.modeNumber.value) + 1) + ' · q ' + number(modalPhaseFactor) +
      ' · peak scale ' + number(modalPresentation.peakScale) + ' · effective scale ' + number(currentScale);
    else text = type + ' · deformation scale ' + number(currentScale);
    if (type === 'transient') text += ' · sample ' + (Number(elements.timeIndex.value) + 1);
    setStatus(text);
  }
  function exportContext() {
    var context = { title: dataset.raw.metadata.title, analysisType: dataset.raw.analysis.type, deformationScale: currentScale,
      staticResult: elements.staticResult.value, modeIndex: Number(elements.modeNumber.value || 0), timeIndex: Number(elements.timeIndex.value || 0),
      historyDOF: Number(elements.historyDof.value || 0), historyQuantity: elements.historyQuantity.value || '', exportedUtc: new Date().toISOString() };
    if (dataset.raw.analysis.type === 'modal') {
      var modeNumber = Number(elements.modeNumber.value) + 1;
      context.modalModeNumber = modeNumber;
      context.modalFrequencyHz = modalPresentation.frequencyHz;
      context.modalAngularFrequencyRadPerSec = modalPresentation.angularFrequency;
      context.modalPhaseRadians = modalPhaseAngle;
      context.modalPhaseFactor = modalPhaseFactor;
      context.modalDisplayMultiplier = modalPresentation.displayMultiplier;
      context.modalPeakScale = modalPresentation.peakScale;
      context.modalEffectiveScale = currentScale;
      context.captionLines = [
        'Mode ' + modeNumber + ' · f = ' + number(modalPresentation.frequencyHz) + ' Hz · ω = ' + number(modalPresentation.angularFrequency) + ' rad/s',
        'Illustrative phase q = ' + number(modalPhaseFactor) + ' · display multiplier ' + number(context.modalDisplayMultiplier),
        'Peak scale ' + number(modalPresentation.peakScale) + ' · effective scale ' + number(currentScale),
        'Mode-shape amplitudes are normalized and arbitrary; they are not physical displacements.'
      ];
    } else if (dataset.raw.analysis.type === 'transient') {
      var analysis = dataset.raw.analysis, timeIndex = Number(elements.timeIndex.value), dofId = Number(elements.historyDof.value || 0);
      var quantity = elements.historyQuantity.value || '', location = dataset.dofById.get(dofId), descriptor = dofId && quantity ? M.transientView.quantityDescriptor(dataset, dofId, quantity) : null;
      var series = dofId && quantity ? M.transientView.series(dataset, dofId, quantity) : null, selectedIndex = quantity === 'spectrum' ? 0 : timeIndex;
      if (quantity === 'spectrum' && elements.chart.__chartState && elements.chart.__chartState.selectedIndex !== null) selectedIndex = elements.chart.__chartState.selectedIndex;
      var sampling = M.transientView.samplingSummary(dataset);
      context.transientTime = analysis.time[timeIndex];
      context.transientSampleNumber = timeIndex + 1;
      context.transientSampleCount = analysis.time.length;
      context.transientDisplayScale = currentScale;
      context.transientFrameStride = Number(elements.frameStride.value);
      context.transientSampling = sampling;
      context.historySampleIndex = series ? selectedIndex : null;
      context.historyX = series ? series.x[selectedIndex] : null;
      context.historyValue = series ? series.y[selectedIndex] : null;
      context.historyUnit = descriptor ? descriptor.unit : '';
      context.captionLines = [
        'Transient sample ' + (timeIndex + 1) + ' / ' + analysis.time.length + ' · t = ' + number(analysis.time[timeIndex]) + (dataset.raw.metadata.units.time ? ' ' + dataset.raw.metadata.units.time : ''),
        'Deformation scale ' + number(currentScale) + ' · playback frame stride ' + context.transientFrameStride,
        location && descriptor ? 'Chart: node ' + location.nodeId + ' · ' + location.label + ' (global DOF ' + dofId + ') · ' + descriptor.label : 'No history field was exported',
        'Sampling: ' + sampling.exported + ' of ' + sampling.original + ' samples exported · export stride ' + sampling.stride + (sampling.decimated ? ' · decimated history' : '')
      ];
      if (series) context.captionLines.splice(3, 0, (quantity === 'spectrum' ? 'Selected spectrum sample: f = ' + number(series.x[selectedIndex]) + ' Hz' : 'Current history value') + ' · ' + number(series.y[selectedIndex]) + (descriptor.unit ? ' ' + descriptor.unit : ''));
      if (quantity === 'reactions') context.captionLines.push('Dynamic residual is M·a + K·u − F; restrained rows are support reactions.');
    }
    return context;
  }

  function bind() {
    elements.fileInput.addEventListener('change', function () { readFile(this.files[0]); this.value = ''; });
    ['dragenter', 'dragover'].forEach(function (name) { elements.dropZone.addEventListener(name, function (event) { event.preventDefault(); this.classList.add('dragging'); }); });
    ['dragleave', 'drop'].forEach(function (name) { elements.dropZone.addEventListener(name, function (event) { event.preventDefault(); this.classList.remove('dragging'); }); });
    elements.dropZone.addEventListener('drop', function (event) { readFile(event.dataTransfer.files[0]); });
    elements.dropZone.addEventListener('click', function (event) { if (!dataset && event.target === this) elements.fileInput.click(); });
    elements.dropZone.addEventListener('keydown', function (event) { if (!dataset && (event.key === 'Enter' || event.key === ' ')) elements.fileInput.click(); });
    elements.reset.addEventListener('click', function () { renderer.fit(); });
    elements.exportSvg.addEventListener('click', function () { if (dataset) M.exporting.downloadSvg(elements.viewport, exportContext(), elements.chartPanel.hidden ? null : elements.chart, currentLegend); });
    elements.exportPng.addEventListener('click', function () { if (dataset) M.exporting.downloadPng(elements.viewport, exportContext(), Number(elements.pngScale.value), elements.chartPanel.hidden ? null : elements.chart, currentLegend).catch(showError); });
    ['showOriginal', 'showDeformed', 'showNodes', 'showNodeLabels', 'showElementLabels', 'showSupports', 'showLoads', 'showReactions'].forEach(function (name) { elements[name].addEventListener('change', function () { renderer.setVisibility(layerSettings()); }); });
    elements.scaleMode.addEventListener('change', function () { elements.manualScale.disabled = this.value !== 'manual'; updateScale(); updateStatus(); });
    elements.manualScale.addEventListener('change', function () { if (elements.scaleMode.value === 'manual') { updateScale(); updateStatus(); } });
    elements.sampleCount.addEventListener('change', function () { updateScale(); if (dataset && dataset.raw.analysis.type === 'static') applyStaticResult(); });
    elements.staticResult.addEventListener('change', applyStaticResult);
    elements.modeNumber.addEventListener('change', changeMode);
    elements.previousMode.addEventListener('click', function () { selectModalMode(Number(elements.modeNumber.value) - 1); });
    elements.nextMode.addEventListener('click', function () { selectModalMode(Number(elements.modeNumber.value) + 1); });
    elements.modalPlayPause.addEventListener('click', startModalPlayback);
    elements.modalPhase.addEventListener('input', function () {
      stopPlayback();
      var factor = Number(this.value);
      applyModalPhase(factor, M.modalView.phaseAngleForFactor(factor));
    });
    elements.modalPlaybackSpeed.addEventListener('change', function () {
      if (playback && playback.type === 'modal') { stopPlayback(); startModalPlayback(); }
    });
    elements.modalAmplitude.addEventListener('change', function () { refreshModalPresentation(false); });
    elements.timeIndex.addEventListener('input', function () { setTransientIndex(Number(this.value)); });
    elements.previousFrame.addEventListener('click', function () { setTransientIndex(M.transientView.stepFrameIndex(Number(elements.timeIndex.value), -1, dataset.raw.analysis.time.length, Number(elements.frameStride.value))); });
    elements.nextFrame.addEventListener('click', function () { setTransientIndex(M.transientView.stepFrameIndex(Number(elements.timeIndex.value), 1, dataset.raw.analysis.time.length, Number(elements.frameStride.value))); });
    elements.playPause.addEventListener('click', startTransientPlayback);
    elements.playbackSpeed.addEventListener('change', function () { if (playback && playback.type === 'transient') { stopPlayback(); startTransientPlayback(); } });
    elements.frameStride.addEventListener('change', function () { if (playback && playback.type === 'transient') { stopPlayback(); startTransientPlayback(); } });
    elements.historyNode.addEventListener('change', function () { updateHistoryDofOptions(); updateChart(); });
    elements.historyDof.addEventListener('change', updateChart); elements.historyQuantity.addEventListener('change', updateChart);
    document.addEventListener('visibilitychange', function () { if (document.hidden) stopPlayback(); });
  }

  document.addEventListener('DOMContentLoaded', function () {
    ['file-input','reset-view','export-svg','export-png','dataset-title','drop-zone','error-panel','notice-panel','analysis-status','selection-details','viewport','result-legend','chart-panel','history-chart','static-controls','modal-controls','transient-controls','chart-controls','static-result','mode-number','previous-mode','next-mode','modal-play-pause','modal-phase','modal-phase-value','modal-playback-speed','modal-amplitude','mode-frequency','time-index','previous-frame','next-frame','play-pause','playback-speed','frame-stride','time-value','transient-metadata','history-node','history-dof','history-quantity','history-value','show-original','show-deformed','show-nodes','show-node-labels','show-element-labels','show-supports','show-loads','show-reactions','scale-mode','manual-scale','sample-count','png-scale'].forEach(function (id) { var key = id.replace(/-([a-z])/g, function (_, letter) { return letter.toUpperCase(); }); elements[key] = byId(id); });
    elements.title = elements.datasetTitle; elements.reset = elements.resetView; elements.error = elements.errorPanel; elements.notice = elements.noticePanel; elements.status = elements.analysisStatus; elements.details = elements.selectionDetails; elements.legend = elements.resultLegend; elements.chart = elements.historyChart;
    renderer = new M.Renderer(elements.viewport, updateDetails); bind();
  });
}(window.MKEFPost));
