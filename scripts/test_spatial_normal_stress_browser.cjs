const {chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE||'playwright'),assert=require('node:assert/strict'),path=require('node:path');
(async()=>{
 const browser=await chromium.launch({headless:true,args:['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
 try{
  const page=await browser.newPage({viewport:{width:1450,height:1050}}),errors=[];page.on('pageerror',e=>errors.push(e.message));
  await page.goto('file://'+path.resolve('postprocessor/index.html'));
  await page.evaluate(()=>{MKEFPost.SpatialPost=new Proxy(MKEFPost.SpatialPost,{construct(T,args){return window.viewer=Reflect.construct(T,args);}});});
  const load=async name=>{await page.evaluate(()=>window.previousViewer=window.viewer);await page.locator('#file-input').setInputFiles(path.resolve('output/'+name+'.json'));await page.waitForFunction(()=>window.viewer&&viewer!==previousViewer);};
  const activate=()=>page.locator('#spatial-result-color').selectOption('normalStress');
  const select=async()=>{await page.locator('#spatial-result-table').selectOption('element');await page.locator('#spatial-result-rows tr').first().click();};
  await load('CaseSpatialFrameSection');await activate();await select();
  assert.match(await page.locator('#normal-stress-details').innerText(),/y=.*z=/);assert(await page.locator('#normal-stress-chart').isVisible());
  await page.evaluate(()=>{const s=viewer.normalStress;window.savedMaximum=s.maximum;window.savedResult=s.results.get(1);window.savedCamera=JSON.stringify(viewer.scene.getViewState());if(!(s.maximum>0))throw Error('Missing maximum');if(viewer.frame.deformed.colors.length!==viewer.frame.deformed.positions.length)throw Error('Color count');});
  await page.locator('#spatial-result-scale-mode').selectOption('manual');await page.locator('#spatial-result-scale').fill('4');await page.locator('#spatial-result-scale').dispatchEvent('change');
  await page.evaluate(()=>{if(viewer.normalStress.maximum!==savedMaximum||viewer.normalStress.results.get(1)!==savedResult)throw Error('Scale recomputed recovery');if(JSON.stringify(viewer.scene.getViewState())!==savedCamera)throw Error('Camera changed');viewer.layers.deformed=false;viewer.update();if(viewer.normalStress.maximum!==savedMaximum)throw Error('Hidden members changed scale');viewer.layers.deformed=true;viewer.update();});
  await page.locator('[data-result-view="fit"]').click();await page.evaluate(()=>{if(viewer.normalStress.results.get(1)!==savedResult)throw Error('Camera recomputed recovery');});
  await page.locator('#normal-stress-chart').press('ArrowRight');assert.match(await page.locator('#normal-stress-sample').innerText(),/σmin=.*σmax=.*max \|σx\|=/);
  await page.locator('#spatial-result-scale-mode').selectOption('auto');await page.locator('[data-result-view="fit"]').click();await page.evaluate(()=>{document.getElementById('result-sidebar').scrollTop=0;viewer.scene.setSelection(null);});await page.screenshot({path:path.resolve('output/normal-stress-spatial.png')});
  for(const name of ['CaseSpatialTrussSection','CaseSpatialTruss',...['circle','circular_tube','rectangle','square_tube','rectangular_tube'].map(s=>'Mises3D_'+s)]){
   await load(name);assert(await page.locator('#normal-stress-panel').isHidden());await activate();await select();
   await page.evaluate(()=>{const s=viewer.normalStress;if([...s.results.values()].some(r=>r.status!=='ok'))throw Error(JSON.stringify([...s.results]));viewer.frame.deformed.owners.forEach((owner,i)=>{const r=s.results.get(viewer.model.elements[owner].id);for(let j=0;j<2;j++){const t=viewer.frame.deformed.parameters[2*i+j],expected=MKEFPost.normalStressResults.rgb(r.at(t).sigmaAbsMax,s.maximum),actual=viewer.frame.deformed.colors.slice(6*i+3*j,6*i+3*j+3);if(actual.some((v,k)=>v!==expected[k]))throw Error('Color disagrees with section');}});});
  }
  // Zero is a valid first scale color even under pure torsion, never unavailable grey.
  await load('Mises3D_circle');await page.evaluate(()=>{viewer.ds.resultsByElementId.get(1).sectionForceCoefficients=Array.from({length:6},(_,i)=>[i===3?1000:0,0,0,0]);const r=viewer.ds.resultsByElementId.get(1);r.sectionStress=MKEFSpatialSections.stress(MKEFPost.spatialResults.sectionProfile(viewer.model.elements[0]),r.sectionForceCoefficients);});await activate();await select();
  assert.match(await page.locator('#normal-stress-details').innerText(),/zero normal stress/);await page.evaluate(()=>{if(viewer.normalStress.maximum!==0)throw Error('Torsion entered normal stress');const expected=MKEFPost.normalStressResults.rgb(0,0);if(!viewer.frame.deformed.colors.every((v,i)=>v===expected[i%3]))throw Error('Zero color');});
  await load('Mises3D_circle');await page.evaluate(()=>{delete viewer.ds.resultsByElementId.get(1).sectionForceCoefficients;});await activate();await select();assert.match(await page.locator('#normal-stress-details').innerText(),/not exported/);assert.match(await page.evaluate(()=>MKEFPost.normalStressResults.csv(viewer.normalStress)),/unavailable/);
  await load('CaseSpatialFrame');await activate();assert.match(await page.locator('#normal-stress-panel details pre').textContent(),/geometry|not exported/);
  await load('CaseSpatialFrameSection');await activate();assert.match(await page.locator('#normal-stress-details').innerText(),/Select a member/);assert.equal(await page.locator('#normal-stress-sample').innerText(),'');
  for(const name of ['CaseSpatialFrameModalSection','CaseSpatialTrussTransient']){await load(name);assert.equal(await page.locator('#spatial-result-color').count(),0);assert.equal(await page.locator('#normal-stress-panel').count(),0);}
  const offline=await browser.newPage({viewport:{width:1450,height:1050}});offline.on('pageerror',e=>errors.push(e.message));
  await offline.addInitScript(()=>{const get=HTMLCanvasElement.prototype.getContext;HTMLCanvasElement.prototype.getContext=function(kind,...args){return kind==='webgl2'?null:get.call(this,kind,...args);};});
  await offline.goto('file://'+path.resolve('postprocessor/index.html'));await offline.locator('#file-input').setInputFiles(path.resolve('output/CaseSpatialFrameSection.json'));await offline.locator('#spatial-result-color').selectOption('normalStress');await offline.locator('#spatial-result-table').selectOption('element');await offline.locator('#spatial-result-rows tr').first().click();assert(await offline.locator('#normal-stress-chart').isVisible());
  assert.deepEqual(errors,[]);console.log('Spatial normal stress: five shapes, manual trusses, exact color values, cache, pure torsion, unavailable data, reset and no WebGL PASS');
 }finally{await browser.close();}
})().catch(e=>{console.error(e);process.exit(1);});
