// Shared harness only: production scripts run unchanged in their browser order.
const fs=require('node:fs'),vm=require('node:vm'),path=require('node:path');
const root=path.resolve(__dirname,'../..');
module.exports=function(returnWindow=false){
  const ctx={window:{MKEFSpatial:{registry:require('../../shared/element-types.json')}},TextEncoder};
  vm.createContext(ctx);
  const files=['shared/sections.js','shared/spatial-sections.js','shared/support-catalog.js',
    'shared/support-semantics.js','shared/mpc.js','shared/releases.js',
    ...['namespace','data-model','spatial-data','spatial-results',
      'spatial-stress','geometry','static-view','section-view','normal-stress-results'].map(f=>'postprocessor/js/'+f+'.js')];
  for(const file of files)vm.runInContext(fs.readFileSync(path.join(root,file),'utf8'),ctx,{filename:file});
  return returnWindow?ctx.window:ctx.window.MKEFPost;
};
