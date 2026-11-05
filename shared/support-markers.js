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
}(window));
