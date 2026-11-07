const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict'),{chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE||'playwright');
(async()=>{
 const browser=await chromium.launch({headless:true,args:['--no-sandbox','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
 try{
  const page=await browser.newPage({viewport:{width:1500,height:1100}}),errors=[];page.on('pageerror',e=>errors.push(e.message));
  for(const name of ['2D-tee-90','2D-rectangular_tube-270','3D-circle--37.5','3D-rectangular_tube--37.5']){
   const file=path.resolve('output/section-contract/'+name+'.json'),raw=JSON.parse(fs.readFileSync(file)),spatial=raw.model.dimension===3;
   await page.goto('file://'+path.resolve('postprocessor/index.html'));
   await page.evaluate(()=>{
    for(const [o,k]of [[MKEFSections,'properties'],[MKEFSpatialSections,'properties'],[MKEFSpatialSections,'stress'],[MKEFSections,'roots']])o[k]=()=>{throw Error('Unexpected browser calculation: '+k);};
   });
   await page.locator('#file-input').setInputFiles(file);
   await page.waitForFunction(()=>!document.querySelector('#download-json').disabled);
   assert.equal(await page.locator('#error-panel').isVisible(),false,name);
   await page.locator(spatial?'#spatial-result-table':'#planar-result-objects select').selectOption('element');
   await page.locator(spatial?'#spatial-result-rows tr':'#planar-result-rows tr').first().click();
   await page.locator(spatial?'#spatial-result-color':'#static-result').selectOption('normalStress');
   assert.match(await page.locator('#normal-stress-details').innerText(),/σmin=.*σmax=/);
   if(spatial){await page.locator('#spatial-result-quantity').selectOption('sigmaMax');assert(await page.locator('#spatial-section-svg path').count());}
   else{assert(await page.locator('#section-result-svg').isVisible());await page.locator('#section-position-number').fill('.31');await page.locator('#section-position-number').dispatchEvent('input');}
   for(const selector of [spatial?'#spatial-section-export-svg':'#export-section-svg','#normal-stress-export-csv']){
    const wait=page.waitForEvent('download');await page.locator(selector).click();const download=await wait,stream=await download.createReadStream();let text='';for await(const part of stream)text+=part;
    assert(text.length>100);assert(!text.includes('NaN'));assert.match(text,selector.includes('csv')?/sigmaMin/:/<svg/);
   }
   if(!name.includes('tee')){
    const body=await page.evaluate(async()=>{
     const fetch=window.fetch;let body;
     window.fetch=async(url,options)=>{if(String(url).endsWith('/member')){body=JSON.parse(options.body);return new Response(JSON.stringify({id:'test'}));}return new Response(JSON.stringify({state:'completed',result:{status:'ok',maximum:{value:1,xi:.5}}}));};
     try{const button=document.querySelector('#member-run');button.disabled=false;await button.onclick();return body;}finally{window.fetch=fetch;}
    });
    assert.equal(body.section.type,raw.model.sections[0].type);
    if(spatial){assert.equal(body.J_model,raw.model.elements[0].sectionProperties.torsionConstant);assert.deepEqual(body.coefficients,raw.analysis.elementResults[0].sectionForceCoefficients);}
    else{const c=raw.analysis.elementResults[0].sectionForceCoefficients;assert.deepEqual(body.coefficients,JSON.parse(JSON.stringify([c[0],[0,0,0,0],c[1],[0],c[2].map(v=>-v),[0,0,0,0]])));}
   }
  }
  assert.deepEqual(errors,[]);console.log('Solver section browser: disabled properties/extrema formulae, 2D/3D colors, cuts, SVG/CSV and Mises request/J PASS');
 }finally{await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
