(function (M) {
  'use strict';
  var C = M.charts = {};
  function clear(svg) { while (svg.firstChild) svg.removeChild(svg.firstChild); }
  function text(svg, x, y, value, anchor) { var node = M.svgElement('text', { x: x, y: y, 'text-anchor': anchor || 'start' }); node.textContent = value; svg.appendChild(node); }
  function finiteExtent(values) { var minimum = Math.min.apply(null, values), maximum = Math.max.apply(null, values); if (minimum === maximum) { var pad = Math.max(Math.abs(minimum) * .1, 1); minimum -= pad; maximum += pad; } return [minimum, maximum]; }
  function format(value) { return value === 0 ? '0' : Number(value.toPrecision(4)).toString(); }
  function nearestIndex(values, target) {
    var low = 0, high = values.length - 1;
    while (low < high) { var middle = Math.floor((low + high) / 2); if (values[middle] < target) low = middle + 1; else high = middle; }
    if (low > 0 && Math.abs(values[low - 1] - target) <= Math.abs(values[low] - target)) return low - 1;
    return low;
  }
  function showSample(svg, index) {
    var state = svg.__chartState;
    if (!state || index < 0 || index >= state.series.x.length) return;
    var x = state.sx(state.series.x[index]), y = state.sy(state.series.y[index]);
    state.sample.removeAttribute('hidden'); state.tooltip.removeAttribute('hidden');
    state.sample.setAttribute('cx', x); state.sample.setAttribute('cy', y);
    state.tooltipText.textContent = 'x ' + format(state.series.x[index]) + '   y ' + format(state.series.y[index]);
    state.tooltip.setAttribute('transform', 'translate(' + Math.max(8, Math.min(x + 10, 612)) + ' ' + Math.max(8, Math.min(y - 48, 190)) + ')');
    state.sampleTitle.textContent = state.series.xLabel + ': ' + format(state.series.x[index]) + '; ' + state.series.yLabel + ': ' + format(state.series.y[index]);
    state.hoverIndex = index;
  }
  function hideSample(svg) {
    var state = svg.__chartState;
    if (!state) return;
    state.sample.setAttribute('hidden', ''); state.tooltip.setAttribute('hidden', ''); state.hoverIndex = null;
  }
  function eventIndex(svg, event) {
    var state = svg.__chartState, rect = svg.getBoundingClientRect();
    var viewX = (event.clientX - rect.left) / Math.max(rect.width, 1) * 800;
    var value = state.xExtent[0] + (viewX - state.margin.left) / state.frameWidth * (state.xExtent[1] - state.xExtent[0]);
    return nearestIndex(state.series.x, value);
  }
  function bind(svg) {
    if (svg.__chartBound) return;
    svg.__chartBound = true;
    svg.addEventListener('pointermove', function (event) { if (svg.__chartState) showSample(svg, eventIndex(svg, event)); });
    svg.addEventListener('pointerleave', function () { var state = svg.__chartState; if (!state) return; if (state.selectedIndex === null) hideSample(svg); else showSample(svg, state.selectedIndex); });
    svg.addEventListener('click', function (event) { if (!svg.__chartState) return; var index = eventIndex(svg, event); svg.__chartState.selectedIndex = index; showSample(svg, index); if (svg.__chartState.onSelect) svg.__chartState.onSelect(index, svg.__chartState.series.x[index], svg.__chartState.series.y[index]); });
    svg.addEventListener('keydown', function (event) {
      var state = svg.__chartState; if (!state || (event.key !== 'ArrowLeft' && event.key !== 'ArrowRight')) return;
      event.preventDefault();
      var index = state.selectedIndex === null ? 0 : state.selectedIndex + (event.key === 'ArrowLeft' ? -1 : 1);
      state.selectedIndex = Math.max(0, Math.min(state.series.x.length - 1, index)); showSample(svg, state.selectedIndex);
      if (state.onSelect) state.onSelect(state.selectedIndex, state.series.x[state.selectedIndex], state.series.y[state.selectedIndex]);
    });
  }

  C.render = function (svg, inputSeries, cursorX, onSelect) {
    var selectedIndex = svg.__chartState ? svg.__chartState.selectedIndex : null;
    clear(svg);
    svg.__chartState = null;
    if (!inputSeries || !inputSeries.x.length) return;
    var series = M.transientView.downsample(inputSeries, 2000);
    var width = 800, height = 240, margin = { left: 72, right: 20, top: 20, bottom: 42 };
    var frameWidth = width - margin.left - margin.right, frameHeight = height - margin.top - margin.bottom;
    var xExtent = finiteExtent(inputSeries.x), yExtent = finiteExtent(series.y);
    function sx(value) { return margin.left + (value - xExtent[0]) / (xExtent[1] - xExtent[0]) * frameWidth; }
    function sy(value) { return margin.top + frameHeight - (value - yExtent[0]) / (yExtent[1] - yExtent[0]) * frameHeight; }
    svg.setAttribute('viewBox', '0 0 ' + width + ' ' + height);
    svg.appendChild(M.svgElement('rect', { class: 'chart-frame', x: margin.left, y: margin.top, width: frameWidth, height: frameHeight }));
    var d = series.x.map(function (value, index) { return (index ? 'L ' : 'M ') + sx(value).toFixed(3) + ' ' + sy(series.y[index]).toFixed(3); }).join(' ');
    svg.appendChild(M.svgElement('path', { class: 'chart-series', d: d }));
    if (Number.isFinite(cursorX) && cursorX >= xExtent[0] && cursorX <= xExtent[1]) svg.appendChild(M.svgElement('line', { class: 'chart-cursor', x1: sx(cursorX), x2: sx(cursorX), y1: margin.top, y2: margin.top + frameHeight }));
    text(svg, margin.left, height - 10, inputSeries.xLabel, 'start');
    text(svg, 8, margin.top + 12, inputSeries.yLabel, 'start');
    text(svg, margin.left, margin.top + frameHeight + 18, format(xExtent[0]), 'start');
    text(svg, margin.left + frameWidth, margin.top + frameHeight + 18, format(xExtent[1]), 'end');
    text(svg, margin.left - 8, margin.top + 5, format(yExtent[1]), 'end');
    text(svg, margin.left - 8, margin.top + frameHeight, format(yExtent[0]), 'end');
    var hit = M.svgElement('rect', { class: 'chart-hit-area', x: margin.left, y: margin.top, width: frameWidth, height: frameHeight, tabindex: '0', 'aria-label': 'History samples' });
    var sample = M.svgElement('circle', { class: 'chart-sample', r: 4, hidden: '' }), sampleTitle = M.svgElement('title'); sample.appendChild(sampleTitle);
    var tooltip = M.svgElement('g', { class: 'chart-tooltip', hidden: '' });
    tooltip.appendChild(M.svgElement('rect', { width: 180, height: 32, rx: 4 }));
    var tooltipText = M.svgElement('text', { x: 8, y: 21 }); tooltip.appendChild(tooltipText);
    svg.appendChild(hit); svg.appendChild(sample); svg.appendChild(tooltip);
    svg.__chartState = { series: inputSeries, xExtent: xExtent, margin: margin, frameWidth: frameWidth, selectedIndex: selectedIndex,
      sx: sx, sy: sy, sample: sample, sampleTitle: sampleTitle, tooltip: tooltip, tooltipText: tooltipText, hoverIndex: null, onSelect: onSelect || null };
    bind(svg);
    if (selectedIndex !== null) { svg.__chartState.selectedIndex = Math.min(selectedIndex, inputSeries.x.length - 1); showSample(svg, svg.__chartState.selectedIndex); }
  };
}(window.MKEFPost));
