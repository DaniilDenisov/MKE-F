// Rebuild the committed, offline-compatible browser asset with npm ci && npm run build:spatial.
const fs = require('node:fs');
const path = require('node:path');
const root = path.resolve(__dirname, '..');
const check = process.argv.includes('--check');
const result = require('esbuild').buildSync({
  entryPoints: [path.join(root, 'shared/spatial/index.js')],
  outfile: path.join(root, 'shared/vendor/spatial.js'),
  bundle: true, format: 'iife', globalName: 'MKEFSpatial', target: 'es2020',
  minify: true, legalComments: 'eof', write: !check
});
const licenseSource = path.join(root, 'node_modules/three/LICENSE'), licenseTarget = path.join(root, 'shared/vendor/THREE-LICENSE.txt');
if (check) {
  const normalize = text => text.replace(/\r\n/g,'\n');
  require('node:assert/strict').equal(normalize(result.outputFiles[0].text),normalize(fs.readFileSync(path.join(root,'shared/vendor/spatial.js'),'utf8')),'Rebuild the spatial bundle after changing scene code or the registry.');
  require('node:assert/strict').equal(normalize(fs.readFileSync(licenseSource,'utf8')),normalize(fs.readFileSync(licenseTarget,'utf8')));
  console.log('Committed spatial bundle, registry and license match the pinned dependencies.');
} else fs.copyFileSync(licenseSource,licenseTarget);
