(function (global) {
  'use strict';
  // Embedded native TXT keeps the examples available offline as well as when hosted.
  // SI units: m, N, Pa, kg, s. Each geometry is reused across all three analyses.
  const geometries = [
    {
      id: '2d-cantilever', dimension: 2, label: 'Cantilever beam', monitor: '3,2',
      time: '0.0001,0.05', force: '3,0,-1000,0',
      text: `nodes
3
0,0,0
1,0,0
2,0,0
elems_113
2
113,1,2,0.02,2e11,7850,0.00006666666666666668
113,2,3,0.02,2e11,7850,0.00006666666666666668
sections
1
{"id":1,"name":"Rectangle 100 x 200 mm","type":"rectangle","dimensions":{"b":0.1,"h":0.2}}
element_sections
2
1,1,0
2,1,0
bcfix
1
1,1,0,0,0
`
    },
    {
      id: '2d-triangle', dimension: 2, label: 'Triangle truss', monitor: '3,2',
      time: '0.00001,0.005', force: '3,0,-1000,0',
      text: `nodes
3
0,0,0
2,0,0
1,1,0
elems_112
3
112,1,2,0.0004,2e11,7850
112,1,3,0.0004,2e11,7850
112,2,3,0.0004,2e11,7850
bcfix
2
4,1,0,0,0
6,2,0,0,0
`
    },
    {
      id: '3d-cantilever', dimension: 3, label: 'Cantilever beam', monitor: '3,2',
      time: '0.0001,0.05', force: '3,0,-1000,500,100,0,0',
      text: `nodes
3
0,0,0
1,0,0
2,0,0
elems_1436
2
1436,1,2,0.02,2e11,7850,8e10,0.00001666666666666667,0.00006666666666666668,0.00004573633542393663,0,1,0,0
1436,2,3,0.02,2e11,7850,8e10,0.00001666666666666667,0.00006666666666666668,0.00004573633542393663,0,1,0,0
sections
1
{"id":1,"name":"Rectangle 100 x 200 mm","type":"rectangle","dimensions":{"b":0.1,"h":0.2}}
element_sections
2
1,1,0
2,1,0
bcfix_3d
1
1,1,1,1,1,1,1
`
    },
    {
      id: '3d-tripod', dimension: 3, label: 'Tripod truss', monitor: '4,3',
      time: '0.00001,0.005', force: '4,0,0,-1000,0,0,0',
      text: `nodes
4
0,0,0
2,0,0
1,2,0
1,1,2
elems_1133
3
1133,1,4,0.0004,2e11,7850
1133,2,4,0.0004,2e11,7850
1133,3,4,0.0004,2e11,7850
bcfix_3d
3
1,1,1,1,0,0,0
2,1,1,1,0,0,0
3,1,1,1,0,0,0
`
    }
  ];
  global.MKEFPre.examples = ['static', 'modal', 'transient'].flatMap(analysis =>
    geometries.map(geometry => {
      const transient = analysis === 'transient', suffix = geometry.dimension === 3 ? '_3d' : '';
      const settings = transient ? analysis + ',' + geometry.time + ',' + geometry.monitor : analysis;
      const load = analysis === 'modal' ? '' : 'bcforce_' + (transient ? 'step' : 'stat') + suffix + '\n1\n' + (transient ? '13,' : '10,') + geometry.force + '\n';
      return {
        id: geometry.id + '-' + analysis, dimension: geometry.dimension, analysis,
        label: geometry.label + (transient ? ' · step load' : ''),
        name: geometry.dimension + 'D ' + geometry.label + ' ' + analysis,
        text: '# ' + geometry.dimension + 'D ' + geometry.label + ' · ' + analysis + '\n# SI units: m, N, Pa, kg, s\nanalysis\n' + settings + '\n' + geometry.text + load
      };
    })
  );
}(window));
