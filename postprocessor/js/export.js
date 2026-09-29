(function (M) {
  'use strict';
  var E = M.exporting = {};
  E.svgStyles = [
    '.viewport-background{fill:#fff}',
    '.original-element{fill:none;stroke:#7b8795;stroke-width:1.5;stroke-dasharray:6 4;vector-effect:non-scaling-stroke}',
    '.deformed-element{fill:none;stroke:#087f9c;stroke-width:2.7;vector-effect:non-scaling-stroke}',
    '.node{fill:#fff;stroke:#263847;stroke-width:1.4;vector-effect:non-scaling-stroke}',
    '.label{fill:#263847;font-family:sans-serif;paint-order:stroke;stroke:#fff;stroke-linejoin:round}',
    '.node-hit{fill:transparent;stroke:none}',
    '.support-symbol{fill:#fff;stroke:#354657;stroke-width:1.8;vector-effect:non-scaling-stroke}',
    '.load-symbol{fill:none;stroke:#b63a3f;stroke-width:2;vector-effect:non-scaling-stroke}',
    '.reaction-symbol{fill:none;vector-effect:non-scaling-stroke}.reaction-force-symbol{stroke:#6d3ca5;stroke-width:2.8;stroke-dasharray:5 2.5}.reaction-moment-symbol{stroke:#c56a00;stroke-width:3;stroke-dasharray:none}',
    '.diagram-fill{fill:#dceff4;stroke:none}.diagram-line{fill:none;stroke:#176b87;stroke-width:2;vector-effect:non-scaling-stroke}',
    'path.selected,circle.selected{stroke:#ef7d18!important;stroke-width:4!important}',
    '.chart-frame,.chart-axis{fill:none;stroke:#9aa8b7;stroke-width:1}.chart-series{fill:none;stroke:#087f9c;stroke-width:2}.chart-cursor{stroke:#b63a3f;stroke-width:1.5}.chart-hit-area{fill:transparent}.chart-sample{fill:#fff;stroke:#b63a3f;stroke-width:2}.chart-tooltip rect{fill:#263847;opacity:.94}.chart-export text{fill:#354657;font:12px sans-serif}.chart-export .chart-tooltip text{fill:#fff}.chart-export [hidden]{display:none}'
  ].join('');

  function legendNumber(value) { if (value === 0) return '0'; return Math.abs(value) >= 10000 || Math.abs(value) < .001 ? value.toExponential(5) : Number(value.toPrecision(7)).toString(); }
  function appendLegend(root, legend, y, width) {
    var height = 92, padding = 32, barY = y + 39, barWidth = width - 2 * padding;
    var group = M.svgElement('g', { class: 'export-legend', 'data-export-legend': '', transform: 'translate(0 0)' });
    group.appendChild(M.svgElement('rect', { x: 0, y: y, width: width, height: height, fill: '#ffffff', stroke: '#d9e0e9' }));
    var title = M.svgElement('text', { x: padding, y: y + 25, fill: '#18212f', 'font-family': 'sans-serif', 'font-size': 14, 'font-weight': 600 });
    title.textContent = legend.label + (legend.unit ? ' [' + legend.unit + ']' : ''); group.appendChild(title);
    var defs = M.svgElement('defs'), gradient = M.svgElement('linearGradient', { id: 'mkef-export-result-legend', x1: '0%', y1: '0%', x2: '100%', y2: '0%' });
    legend.stops.forEach(function (stop) { gradient.appendChild(M.svgElement('stop', { offset: (stop.offset * 100) + '%', 'stop-color': stop.color })); });
    defs.appendChild(gradient); group.appendChild(defs);
    group.appendChild(M.svgElement('rect', { x: padding, y: barY, width: barWidth, height: 14, rx: 7, fill: 'url(#mkef-export-result-legend)', stroke: '#687386', 'stroke-width': 1 }));
    legend.ticks.forEach(function (value, index) {
      var single = legend.ticks.length === 1, x = single ? width / 2 : padding + index * barWidth / (legend.ticks.length - 1);
      var tick = M.svgElement('text', { x: x, y: barY + 34, fill: '#18212f', 'font-family': 'sans-serif', 'font-size': 12, 'text-anchor': single ? 'middle' : index === 0 ? 'start' : index === legend.ticks.length - 1 ? 'end' : 'middle' });
      tick.textContent = legendNumber(value); group.appendChild(tick);
    });
    root.appendChild(group);
    return height;
  }

  function appendCaption(root, lines, y, width) {
    var padding = 32, lineHeight = 20, height = 28 + lines.length * lineHeight;
    var group = M.svgElement('g', { class: 'export-caption', 'data-export-caption': '' });
    group.appendChild(M.svgElement('rect', { x: 0, y: y, width: width, height: height, fill: '#ffffff', stroke: '#d9e0e9' }));
    lines.forEach(function (line, index) {
      var text = M.svgElement('text', { x: padding, y: y + 24 + index * lineHeight, fill: '#18212f', 'font-family': 'sans-serif', 'font-size': index ? 12 : 14, 'font-weight': index ? 400 : 600 });
      text.textContent = line; group.appendChild(text);
    });
    root.appendChild(group);
    return height;
  }

  E.prepareSvg = function (source, context, chart, legend) {
    var rect = source.getBoundingClientRect(), width = Math.max(Math.round(rect.width), 800), sceneHeight = Math.max(Math.round(rect.height), 600);
    var includeChart = chart && chart.childNodes.length && !chart.closest('[hidden]');
    var includeLegend = legend && legend.kind === 'color';
    var captionLines = Array.isArray(context.captionLines) ? context.captionLines.map(String) : [];
    var captionHeight = captionLines.length ? 28 + captionLines.length * 20 : 0;
    var legendHeight = includeLegend ? 92 : 0;
    var chartHeight = includeChart ? Math.max(Math.round(chart.getBoundingClientRect().height), 240) : 0;
    var totalHeight = sceneHeight + captionHeight + legendHeight + chartHeight;
    var root = M.svgElement('svg', { xmlns: M.svgNS, width: width, height: totalHeight, viewBox: '0 0 ' + width + ' ' + totalHeight });
    var style = M.svgElement('style'); style.textContent = E.svgStyles; root.appendChild(style);
    var metadata = M.svgElement('metadata'); metadata.textContent = JSON.stringify(context); root.appendChild(metadata);
    var scene = source.cloneNode(true);
    scene.removeAttribute('hidden');
    scene.setAttribute('x', 0); scene.setAttribute('y', 0); scene.setAttribute('width', width); scene.setAttribute('height', sceneHeight);
    root.appendChild(scene);
    if (captionLines.length) appendCaption(root, captionLines, sceneHeight, width);
    if (includeLegend) appendLegend(root, legend, sceneHeight + captionHeight, width);
    if (includeChart) {
      var chartClone = chart.cloneNode(true);
      chartClone.setAttribute('class', ((chartClone.getAttribute('class') || '') + ' chart-export').trim());
      chartClone.setAttribute('x', 0); chartClone.setAttribute('y', sceneHeight + captionHeight + legendHeight); chartClone.setAttribute('width', width); chartClone.setAttribute('height', chartHeight);
      root.appendChild(chartClone);
    }
    root.querySelectorAll('[tabindex]').forEach(function (element) { element.removeAttribute('tabindex'); });
    return root;
  };
  E.serialize = function (source, context, chart, legend) { return '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(E.prepareSvg(source, context, chart, legend)); };
  E.filename = function (title, extension) { var base = String(title || 'mkef-result').trim().replace(/[^A-Za-z0-9._-]+/g, '-').replace(/^-+|-+$/g, '') || 'mkef-result'; return base + '.' + extension; };
  E.downloadBlob = function (blob, filename) { var url = URL.createObjectURL(blob), link = document.createElement('a'); link.href = url; link.download = filename; document.body.appendChild(link); link.click(); link.remove(); setTimeout(function () { URL.revokeObjectURL(url); }, 0); };
  E.downloadSvg = function (source, context, chart, legend) { E.downloadBlob(new Blob([E.serialize(source, context, chart, legend)], { type: 'image/svg+xml;charset=utf-8' }), E.filename(context.title, 'svg')); };
  E.pngBlob = function (source, context, scale, chart, legend) {
    return new Promise(function (resolve, reject) {
      var prepared = E.prepareSvg(source, context, chart, legend), width = Number(prepared.getAttribute('width')), height = Number(prepared.getAttribute('height'));
      var blob = new Blob([new XMLSerializer().serializeToString(prepared)], { type: 'image/svg+xml;charset=utf-8' }), url = URL.createObjectURL(blob), image = new Image();
      image.onload = function () {
        try {
          var canvas = document.createElement('canvas'); canvas.width = Math.round(width * scale); canvas.height = Math.round(height * scale);
          var drawing = canvas.getContext('2d'); drawing.setTransform(scale, 0, 0, scale, 0, 0); drawing.drawImage(image, 0, 0, width, height);
          canvas.toBlob(function (png) { URL.revokeObjectURL(url); png ? resolve(png) : reject(new Error('Canvas could not encode PNG.')); }, 'image/png');
        } catch (error) { URL.revokeObjectURL(url); reject(error); }
      };
      image.onerror = function () { URL.revokeObjectURL(url); reject(new Error('The prepared SVG could not be rasterized.')); };
      image.src = url;
    });
  };
  E.downloadPng = function (source, context, scale, chart, legend) { return E.pngBlob(source, context, scale, chart, legend).then(function (blob) { E.downloadBlob(blob, E.filename(context.title, 'png')); }); };
}(window.MKEFPost));
