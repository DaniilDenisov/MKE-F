(function (global) {
  'use strict';
  var M = global.MKEFPre = global.MKEFPre || {};
  M.config = {
    maxFileBytes: 10 * 1024 * 1024,
    maxNodes: 100000,
    maxElements: 100000
  };
  M.svgNS = 'http://www.w3.org/2000/svg';
  M.svgElement = function (name, attributes) {
    var element = document.createElementNS(M.svgNS, name);
    Object.keys(attributes || {}).forEach(function (key) {
      element.setAttribute(key, String(attributes[key]));
    });
    return element;
  };
}(window));
