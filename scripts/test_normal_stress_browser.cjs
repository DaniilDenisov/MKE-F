const {chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE||'playwright'),assert=require('node:assert/strict'),path=require('node:path'),fs=require('node:fs');
(async()=>{
 const browser=await chromium.launch({headless:true,args:['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
 try{
  const page=await browser.newPage({viewport:{width:1450,height:1050}}),errors=[];page.on('pageerror',e=>errors.push(e.message));
  await page.goto('file://'+path.resolve('postprocessor/index.html'));
  await page.evaluate(()=>{const mount=MKEFPost.WebGLRenderer.prototype.mount;MKEFPost.WebGLRenderer.prototype.mount=function(...args){window.renderer=this;window.dataset=args[0];return mount.apply(this,args);};});
  const load=async name=>{await page.evaluate(()=>window.previousDataset=window.dataset);await page.locator('#file-input').setInputFiles(path.resolve('output/'+name+'.json'));await page.waitForFunction(()=>window.dataset&&dataset!==previousDataset&&(dataset.raw.analysis.type!=='static'||!!dataset.normalStress));};
  const activate=()=>page.locator('#static-result').selectOption('normalStress');
  const select=async()=>{await page.locator('#planar-result-objects select').selectOption('element');await page.locator('#planar-result-rows tr').first().click();};
  await load('CaseSectionBeam');await activate();await select();
  assert.match(await page.locator('#normal-stress-details').innerText(),/σmin=.*σmax=.*\nmax \|σx\|=/);
  assert.match(await page.locator('#normal-stress-details').innerText(),/tension.*x\/L=0.5/);
  assert.match(await page.locator('#normal-stress-legend').innerText(),/stress units not specified/);
  assert.match(await page.locator('.normal-stress-settings').innerText(),/Pure torsion can give zero/);
  await page.evaluate(()=>{const s=dataset.normalStress;window.savedResult=s.results.get(1);window.savedMaximum=s.maximum;renderer.flush();if(renderer.frame.deformed.colors.length!==renderer.frame.deformed.positions.length)throw Error('Color buffer size');if(!savedResult.positions.includes(.5))throw Error('Interior peak missing');if(!renderer.frame.original.colors.every((v,i)=>Math.abs(v-[135/255,151/255,167/255][i%3])<1e-12))throw Error('Original should be neutral');});
  await page.locator('#scale-mode').selectOption('manual');await page.locator('#manual-scale').fill('4');await page.locator('#manual-scale').dispatchEvent('change');
  await page.locator('#show-original').uncheck();await page.locator('#normal-stress-chart').press('ArrowRight');
  assert.match(await page.locator('#normal-stress-sample').innerText(),/Section at x\/L=.*\nσmin=.*σmax=.*max \|σx\|=/);
  await page.evaluate(()=>{if(dataset.normalStress.results.get(1)!==savedResult||dataset.normalStress.maximum!==savedMaximum)throw Error('Display state invalidated recovery');const chart=document.getElementById('normal-stress-chart').__chartState.series;chart.y.forEach((v,i)=>{if(v!==savedResult.at(chart.x[i]).sigmaAbsMax)throw Error('Chart disagreement');});});
  const waiting=page.waitForEvent('download');await page.locator('#normal-stress-export-csv').click();const download=await waiting,stream=await download.createReadStream();let csv='';for await(const chunk of stream)csv+=chunk;
  assert.match(csv,/sigmaMin.*sigmaMax.*sigmaAbsMax.*sigmaX.*eta.*y.*z/);assert.doesNotMatch(csv,/nuSource|tauXY|vonMises/);
  await page.evaluate(()=>{const csv=MKEFPost.normalStressResults.csv(dataset.normalStress).split('\r\n').slice(1);csv.forEach(line=>{const v=line.split(',').map(s=>s.slice(1,-1)),r=dataset.normalStress.results.get(Number(v[0])).at(Number(v[1]));if(Number(v[4])!==r.sigmaAbsMax||Number(v[5])!==r.sigmaX)throw Error('CSV disagreement');});});
  await page.screenshot({path:path.resolve('output/normal-stress-planar.png')});
  for(const name of ['SectionCase_rectangular_tube_0','SectionCase_square_tube_0','SectionCase_circular_tube_90',...['circle','circular_tube','rectangle','square_tube','rectangular_tube'].map(s=>'Mises2D_'+s)]){
   await load(name);assert.equal(await page.locator('#normal-stress-panel').count(),1);await activate();await select();
   await page.evaluate(()=>{if([...dataset.normalStress.results.values()].some(r=>r.status!=='ok'))throw Error(JSON.stringify([...dataset.normalStress.results]));renderer.flush();if(!renderer.frame.deformed.colors.every(Number.isFinite))throw Error('Non-finite colors');});
  }
  const manual=JSON.parse(fs.readFileSync('output/CaseSectionBeam.json'));
  delete manual.model.sections;delete manual.model.elementLoads;delete manual.model.nodalLoads;
  manual.model.numberOfDOFs=4;manual.model.recovery.stiffnessFormulation='linear-bar';
  manual.model.dofPerNode=2;manual.model.dofLabels=['ux','uy'];manual.model.dofMap=[[1,2],[3,4]];manual.model.supports=[];
  const member=manual.model.elements[0];member.type=112;delete member.sectionId;member.properties.area=2;
  member.globalDOFs=[1,2,3,4];member.localFromSystem=[0,1,3,4].map(i=>[0,1,3,4].map(j=>member.localFromSystem[i][j]));
  Object.assign(manual.analysis,{displacements:[0,0,.01,0],loadVector:[0,0,100,0],reactions:[-100,0,0,0],equilibriumResidual:[0,0,0],elementResults:[{elementId:1,type:112,localEndForces:[-100,0,100,0],axialStrain:.005,axialStress:50,axialForce:100}]});
  Object.assign(manual.analysis.elementResults[0],{localDisplacements:[0,0,.01,0],equivalentLocalLoadVector:[0,0,0,0],sectionStress:{minimum:{value:50,xi:0},maximum:{value:50,xi:0},positions:[0,1]},sectionForceCoefficients:[[100,0,0,0],[0,0,0,0],[0,0,0,0]]});
  await page.evaluate(()=>window.previousDataset=window.dataset);await page.locator('#file-input').setInputFiles({name:'manual-truss.json',mimeType:'application/json',buffer:Buffer.from(JSON.stringify(manual))});await page.waitForFunction(()=>dataset!==previousDataset&&!!dataset.normalStress);await activate();await select();
  assert.match(await page.locator('#normal-stress-details').innerText(),/max \|σx\|=50/);assert.match(await page.locator('#normal-stress-details').innerText(),/Uniform axial stress/);
  // A mixed model uses only supported members for the model-wide scale.
  await load('SectionCaseLinearMPC');await activate();assert.match(await page.locator('#normal-stress-status').innerText(),/unavailable: [1-9]/);
  await page.evaluate(()=>{const s=dataset.normalStress,max=Math.max(...[...s.results.values()].filter(r=>r.status==='ok').map(r=>r.sigmaAbsMax));if(s.maximum!==max)throw Error('Mixed scale');renderer.flush();dataset.raw.model.elements.forEach((e,owner)=>{if(s.results.get(e.id).status==='ok')return;renderer.frame.deformed.owners.forEach((id,i)=>{if(id===owner&&!renderer.frame.deformed.colors.slice(i*6,i*6+6).every((v,j)=>Math.abs(v-[135/255,151/255,167/255][j%3])<1e-12))throw Error('Unavailable must be grey');});});});
  await load('CaseSectionBeam');assert(await page.locator('#normal-stress-panel').isHidden());await activate();assert.equal(await page.locator('#normal-stress-sample').innerText(),'');assert.match(await page.locator('#normal-stress-details').innerText(),/Select a member/);
  await page.locator('#static-result').selectOption('M');assert(await page.locator('#normal-stress-panel').isHidden());
  await load('CaseReleaseModal');assert.equal(await page.locator('#normal-stress-panel').count(),0);
  const offline=await browser.newPage({viewport:{width:1450,height:1050}});offline.on('pageerror',e=>errors.push(e.message));
  await offline.addInitScript(()=>{const get=HTMLCanvasElement.prototype.getContext;HTMLCanvasElement.prototype.getContext=function(kind,...args){return kind==='webgl2'?null:get.call(this,kind,...args);};});
  await offline.goto('file://'+path.resolve('postprocessor/index.html'));await offline.locator('#file-input').setInputFiles(path.resolve('output/CaseSectionBeam.json'));await offline.locator('#static-result').selectOption('normalStress');await offline.locator('#planar-result-objects select').selectOption('element');await offline.locator('#planar-result-rows tr').first().click();assert(await offline.locator('#normal-stress-chart').isVisible());
  assert.deepEqual(errors,[]);console.log('Planar normal coloring: extrema, neutral original, cache, chart/CSV, thick tubes, manual trusses, mixed availability, file reset and no WebGL PASS');
 }finally{await browser.close();}
})().catch(e=>{console.error(e);process.exit(1);});
