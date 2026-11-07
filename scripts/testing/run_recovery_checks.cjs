const fs=require('node:fs'),assert=require('node:assert/strict'),{spawnSync}=require('node:child_process');
const path=require('node:path'),root=path.resolve(__dirname,'../..'),manifest=require('./release-tests.json');
process.chdir(root);
// Fresh output is supplied by test_recovery_pipeline.ps1. Missing fixtures fail
// before any glob-based test could accidentally pass with zero cases.
for(const name of ['section-parity','spatial-section-parity','spatial-section-dynamic-parity','spatial-condition-parity']){
  const rows=JSON.parse(fs.readFileSync('output/'+name+'.json'));
  assert(rows.length>0,'Empty '+name);
}
assert.equal(JSON.parse(fs.readFileSync('output/spatial-condition-parity.json')).length,21);
for(const name of ['CaseSpatialFrame','CaseSpatialFrameModal','CaseSpatialFrameTransient',
  'CaseSpatialTruss','CaseSpatialTrussModal','CaseSpatialTrussTransient','CaseSpatialTiltedTransient',
  'CaseSpatialSectionReleasedTransient','static-frame','CaseSectionBeam'])assert(fs.statSync('output/'+name+'.json').size>0,name);
for(const name of manifest.recovery){
  const script=manifest.entryPoints[name];
  assert(script && fs.statSync(script).isFile(),'Missing test entry point: '+name);
  const result=spawnSync(process.execPath,[script],{stdio:'inherit'});
  if(result.error)throw result.error;
  if(result.status!==0)process.exit(result.status||1);
}
fs.writeFileSync('output/recovery/PASS.json',JSON.stringify({completedAt:new Date().toISOString(),checks:manifest.recovery.length},null,2)+'\n');
