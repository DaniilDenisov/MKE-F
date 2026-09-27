(function (M) {
  'use strict';
  var E = M.exporting = {};
  E.svgStyles = [
    '.viewport-background{fill:#fff}',
    '.original-element{fill:none;stroke:#7b8795;stroke-width:1.5;stroke-dasharray:6 4;vector-effect:non-scaling-stroke}',
    '.deformed-element{fill:none;stroke:#087f9c;stroke-width:2.7;vector-effect:non-scaling-stroke}',
    '.node{fill:#fff;stroke:#263847;stroke-width:1.4;vector-effect:non-scaling-stroke}',
    '.label,.result-label{fill:#263847;font:12px sans-serif;paint-order:stroke;stroke:#fff;stroke-width:3px}',
    '.support-symbol,.support-ground{fill:#fff;stroke:#354657;stroke-width:1.8;vector-effect:non-scaling-stroke}',
    '.load-symbol{fill:none;stroke:#b63a3f;stroke-width:2;vector-effect:non-scaling-stroke}',
    '.reaction-symbol{fill:none;stroke:#6d3ca5;stroke-width:2;stroke-dasharray:4 2;vector-effect:non-scaling-stroke}',
    '.diagram-fill{fill:#dceff4;stroke:none}.diagram-line{fill:none;stroke:#176b87;stroke-width:2;vector-effect:non-scaling-stroke}',
    '.legend-background{fill:#fff;stroke:#c8d1dc}.legend-text{fill:#263847;font:11px sans-serif}',
    '.selected{stroke:#ef7d18!important;stroke-width:4!important}'
  ].join('');

  E.prepareSvg = function (source, context) {
    var clone = source.cloneNode(true), rect = source.getBoundingClientRect();
    clone.removeAttribute('hidden');
    clone.setAttribute('xmlns', M.svgNS);
    clone.setAttribute('width', Math.max(Math.round(rect.width), 800));
    clone.setAttribute('height', Math.max(Math.round(rect.height), 600));
    var style = M.svgElement('style'); style.textContent = E.svgStyles;
    clone.insertBefore(style, clone.firstChild);
    var metadata = M.svgElement('metadata'); metadata.textContent = JSON.stringify(context);
    clone.insertBefore(metadata, style.nextSibling);
    clone.querySelectorAll('[tabindex]').forEach(function (element) { element.removeAttribute('tabindex'); });
    return clone;
  };
  E.serialize = function (source, context) {
    return '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(E.prepareSvg(source, context));
  };
  E.filename = function (title, extension) {
    var base = String(title || 'mkef-result').trim().replace(/[^A-Za-z0-9._-]+/g, '-').replace(/^-+|-+$/g, '') || 'mkef-result';
    return base + '.' + extension;
  };
  E.downloadBlob = function (blob, filename) {
    var url = URL.createObjectURL(blob), link = document.createElement('a');
    link.href = url; link.download = filename; document.body.appendChild(link); link.click(); link.remove();
    setTimeout(function () { URL.revokeObjectURL(url); }, 0);
  };
  E.downloadSvg = function (source, context) {
    E.downloadBlob(new Blob([E.serialize(source, context)], { type: 'image/svg+xml;charset=utf-8' }), E.filename(context.title, 'svg'));
  };
  E.pngBlob = function (source, context, scale) {
    return new Promise(function (resolve, reject) {
      var prepared = E.prepareSvg(source, context), width = Number(prepared.getAttribute('width')), height = Number(prepared.getAttribute('height'));
      var blob = new Blob([new XMLSerializer().serializeToString(prepared)], { type: 'image/svg+xml;charset=utf-8' });
      var url = URL.createObjectURL(blob), image = new Image();
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
  E.downloadPng = function (source, context, scale) {
    return E.pngBlob(source, context, scale).then(function (blob) { E.downloadBlob(blob, E.filename(context.title, 'png')); });
  };
}(window.MKEFPost));
