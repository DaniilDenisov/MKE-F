(function (global) {
  'use strict';

  // Coordinates are in SVG space. Positive Mz is counterclockwise.
  global.MKEFMomentSymbol = {
    path: function (x, y, radius, value) {
      var direction = value > 0 ? -1 : 1, offset = radius * Math.SQRT1_2;
      var start = { x: x + offset, y: y + direction * offset };
      var end = { x: x + offset, y: y - direction * offset };
      // The 90-degree gap stays on the right for either sign. Both endpoints
      // lie on the circle centered at the node, so SVG cannot shift its center.
      var tx = Math.SQRT1_2, ty = direction * Math.SQRT1_2;
      var length = radius * .36, halfWidth = radius * .20;
      var bx = end.x - tx * length, by = end.y - ty * length;
      return 'M ' + start.x + ' ' + start.y +
        ' A ' + radius + ' ' + radius + ' 0 1 ' + (value > 0 ? 0 : 1) + ' ' + end.x + ' ' + end.y +
        ' M ' + (bx - ty * halfWidth) + ' ' + (by + tx * halfWidth) +
        ' L ' + end.x + ' ' + end.y +
        ' L ' + (bx + ty * halfWidth) + ' ' + (by - tx * halfWidth);
    }
  };
}(window));
