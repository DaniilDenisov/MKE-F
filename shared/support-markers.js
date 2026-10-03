(function (global) {
  'use strict';
  var S = global.MKEFSupportMarkers = global.MKEFSupportMarkers || {};

  S.svgStyles = [
    '.support-marker .support-body{fill:#fff;stroke:#41546b;stroke-width:2px;vector-effect:non-scaling-stroke}',
    '.support-marker .support-detail,.support-marker .support-ground,.support-marker .support-hatch{fill:none;stroke:#41546b;stroke-width:2px;vector-effect:non-scaling-stroke}',
    '.support-marker .support-ground{stroke-width:3px}',
    '.support-marker .support-roller{fill:#fff;stroke:#41546b;stroke-width:2px;vector-effect:non-scaling-stroke}',
    '.support-marker.selected .support-body,.support-marker.selected .support-detail,.support-marker.selected .support-ground,.support-marker.selected .support-hatch,.support-marker.selected .support-roller{stroke:#ef7d18!important;stroke-width:4px!important}'
  ].join('');

  var style = document.createElement('style');
  style.setAttribute('data-mkef-support-marker-styles', '');
  style.textContent = S.svgStyles;
  document.head.appendChild(style);

  S.label = function (type, dofPerNode) {
    if (type === 1) return 'fixed';
    if (type === 2) return 'free ux only';
    if (type === 3) return 'free uy only';
    return dofPerNode === 2 ? 'fixed (truss)' : 'free rotation only';
  };

  S.append = function (svgElement, layer, x, y, type, size, attributes, titleText) {
    var markerAttributes = { class: 'support-marker support-type-' + type, 'data-support-type': type };
    Object.keys(attributes || {}).forEach(function (key) {
      markerAttributes[key] = key === 'class' ? markerAttributes.class + ' ' + attributes[key] : attributes[key];
    });
    var marker = svgElement('g', markerAttributes);
    if (titleText) {
      var title = svgElement('title');
      title.textContent = titleText;
      marker.appendChild(title);
    }
    function horizontalGround(baseY, halfWidth) {
      marker.appendChild(svgElement('line', { x1: x - halfWidth, y1: baseY, x2: x + halfWidth, y2: baseY, class: 'support-ground' }));
      [-0.75, -0.25, 0.25, 0.75].forEach(function (offset) {
        marker.appendChild(svgElement('line', { x1: x + halfWidth * offset, y1: baseY, x2: x + halfWidth * offset - size * 0.28, y2: baseY + size * 0.3, class: 'support-hatch' }));
      });
    }
    function verticalGround(baseX, halfHeight) {
      marker.appendChild(svgElement('line', { x1: baseX, y1: y - halfHeight, x2: baseX, y2: y + halfHeight, class: 'support-ground' }));
      [-0.75, -0.25, 0.25, 0.75].forEach(function (offset) {
        marker.appendChild(svgElement('line', { x1: baseX, y1: y + halfHeight * offset, x2: baseX - size * 0.3, y2: y + halfHeight * offset + size * 0.28, class: 'support-hatch' }));
      });
    }
    if (type === 1) {
      marker.appendChild(svgElement('line', { x1: x, y1: y, x2: x, y2: y + size * 0.2, class: 'support-detail' }));
      horizontalGround(y + size * 0.2, size * 0.85);
    } else if (type === 2) {
      marker.appendChild(svgElement('line', { x1: x, y1: y, x2: x, y2: y + size * 0.38, class: 'support-detail' }));
      marker.appendChild(svgElement('circle', { cx: x, cy: y + size * 0.5, r: size * 0.13, class: 'support-roller' }));
      marker.appendChild(svgElement('path', { d: 'M ' + x + ' ' + (y + size * 0.64) + ' l ' + (-size * 0.72) + ' ' + (size * 0.62) + ' l ' + (size * 1.44) + ' 0 Z', class: 'support-body' }));
      horizontalGround(y + size * 1.36, size * 0.82);
    } else if (type === 3) {
      marker.appendChild(svgElement('line', { x1: x, y1: y, x2: x - size * 0.38, y2: y, class: 'support-detail' }));
      marker.appendChild(svgElement('circle', { cx: x - size * 0.5, cy: y, r: size * 0.13, class: 'support-roller' }));
      marker.appendChild(svgElement('path', { d: 'M ' + (x - size * 0.64) + ' ' + y + ' l ' + (-size * 0.62) + ' ' + (-size * 0.72) + ' l 0 ' + (size * 1.44) + ' Z', class: 'support-body' }));
      verticalGround(x - size * 1.36, size * 0.82);
    } else {
      marker.appendChild(svgElement('path', { d: 'M ' + x + ' ' + y + ' l ' + (-size) + ' ' + size + ' l ' + (2 * size) + ' 0 Z', class: 'support-body' }));
      horizontalGround(y + size * 1.1, size * 1.1);
    }
    layer.appendChild(marker);
    return marker;
  };
}(window));
