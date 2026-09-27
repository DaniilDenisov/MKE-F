(function (M) {
  'use strict';
  var E = M.exporting = {};
  E.svgStyles = [
    '.viewport-background{fill:#fff}',
    '.original-element{fill:none;stroke:#7b8795;stroke-width:1.5;stroke-dasharray:6 4;vector-effect:non-scaling-stroke}',
    '.deformed-element{fill:none;stroke:#087f9c;stroke-width:2.7;vector-effect:non-scaling-stroke}',
    '.node{fill:#fff;stroke:#263847;stroke-width:1.4;vector-effect:non-scaling-stroke}',
    '.label{fill:#263847;font-family:sans-serif;paint-order:stroke;stroke:#fff;stroke-linejoin:round}',
    '.support-symbol{fill:#fff;stroke:#354657;stroke-width:1.8;vector-effect:non-scaling-stroke}',
    '.load-symbol{fill:none;stroke:#b63a3f;stroke-width:2;vector-effect:non-scaling-stroke}',
    '.reaction-symbol{fill:none;stroke:#6d3ca5;stroke-width:2;stroke-dasharray:4 2;vector-effect:non-scaling-stroke}',
    '.diagram-fill{fill:#dceff4;stroke:none}.diagram-line{fill:none;stroke:#176b87;stroke-width:2;vector-effect:non-scaling-stroke}',
    'path.selected,circle.selected{stroke:#ef7d18!important;stroke-width:4!important}',
    '.chart-frame,.chart-axis{fill:none;stroke:#9aa8b7;stroke-width:1}.chart-series{fill:none;stroke:#087f9c;stroke-width:2}.chart-cursor{stroke:#b63a3f;stroke-width:1.5}.chart-hit-area{fill:transparent}.chart-sample{fill:#fff;stroke:#b63a3f;stroke-width:2}.chart-tooltip rect{fill:#263847;opacity:.94}.chart-export text{fill:#354657;font:12px sans-serif}.chart-export .chart-tooltip text{fill:#fff}.chart-export [hidden]{display:none}'
  ].join('');

  E.prepareSvg = function (source, context, chart) {
    var rect = source.getBoundingClientRect(), width = Math.max(Math.round(rect.width), 800), sceneHeight = Math.max(Math.round(rect.height), 600);
    var includeChart = chart && chart.childNodes.length && !chart.closest('[hidden]');
    var chartHeight = includeChart ? Math.max(Math.round(chart.getBoundingClientRect().height), 240) : 0;
    var root = M.svgElement('svg', { xmlns: M.svgNS, width: width, height: sceneHeight + chartHeight, viewBox: '0 0 ' + width + ' ' + (sceneHeight + chartHeight) });
    var style = M.svgElement('style'); style.textContent = E.svgStyles; root.appendChild(style);
    var metadata = M.svgElement('metadata'); metadata.textContent = JSON.stringify(context); root.appendChild(metadata);
    var scene = source.cloneNode(true);
    scene.removeAttribute('hidden');
    scene.setAttribute('x', 0); scene.setAttribute('y', 0); scene.setAttribute('width', width); scene.setAttribute('height', sceneHeight);
    root.appendChild(scene);
    if (includeChart) {
      var chartClone = chart.cloneNode(true);
      chartClone.setAttribute('class', ((chartClone.getAttribute('class') || '') + ' chart-export').trim());
      chartClone.setAttribute('x', 0); chartClone.setAttribute('y', sceneHeight); chartClone.setAttribute('width', width); chartClone.setAttribute('height', chartHeight);
      root.appendChild(chartClone);
    }
    root.querySelectorAll('[tabindex]').forEach(function (element) { element.removeAttribute('tabindex'); });
    return root;
  };
  E.serialize = function (source, context, chart) { return '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(E.prepareSvg(source, context, chart)); };
  E.filename = function (title, extension) { var base = String(title || 'mkef-result').trim().replace(/[^A-Za-z0-9._-]+/g, '-').replace(/^-+|-+$/g, '') || 'mkef-result'; return base + '.' + extension; };
  E.downloadBlob = function (blob, filename) { var url = URL.createObjectURL(blob), link = document.createElement('a'); link.href = url; link.download = filename; document.body.appendChild(link); link.click(); link.remove(); setTimeout(function () { URL.revokeObjectURL(url); }, 0); };
  E.downloadSvg = function (source, context, chart) { E.downloadBlob(new Blob([E.serialize(source, context, chart)], { type: 'image/svg+xml;charset=utf-8' }), E.filename(context.title, 'svg')); };
  E.pngBlob = function (source, context, scale, chart) {
    return new Promise(function (resolve, reject) {
      var prepared = E.prepareSvg(source, context, chart), width = Number(prepared.getAttribute('width')), height = Number(prepared.getAttribute('height'));
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
  E.downloadPng = function (source, context, scale, chart) { return E.pngBlob(source, context, scale, chart).then(function (blob) { E.downloadBlob(blob, E.filename(context.title, 'png')); }); };
}(window.MKEFPost));
