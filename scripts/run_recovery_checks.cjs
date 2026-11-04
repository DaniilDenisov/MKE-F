const fs=require('node:fs'),assert=require('node:assert/strict'),{spawnSync}=require('node:child_process');
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
for(const name of ['test_result_numbers','test_recovery_operators','test_recovery_pool','test_planar_recovery','test_sections','test_spatial_sections','test_spatial_post',
  'test_normal_stress','test_section_contract','test_section_contract_browser','test_spatial_conditions','test_recovery_browser','test_spatial_post_integration','test_result_resources']){
  const result=spawnSync(process.execPath,['scripts/'+name+'.cjs'],{stdio:'inherit'});
  if(result.error)throw result.error;
  if(result.status!==0)process.exit(result.status||1);
}
fs.writeFileSync('output/recovery/PASS.json',JSON.stringify({completedAt:new Date().toISOString(),checks:14},null,2)+'\n');
