// Run inside the isolated rebuilt runtime stack via test_v8_stack.ps1.
const fs=require('node:fs'),assert=require('node:assert/strict'),{chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE||'playwright');
const base='http://localhost:8080',pause=ms=>new Promise(resolve=>setTimeout(resolve,ms));
(async()=>{
  const browser=await chromium.launch({headless:true,args:['--no-sandbox','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
  try{
    const page=await browser.newPage({viewport:{width:1500,height:1100}}),errors=[];
    page.on('pageerror',error=>errors.push(error.message));
    let ready=false;
    for(let i=0;i<120;i++){
      try{const r=await page.request.get(base+'/api/v1/health');ready=r.ok()&&(await r.json()).status==='ready';}catch{}
      if(ready)break;await pause(500);
    }
    assert(ready,'Rebuilt solver never became ready');
    assert.equal((await page.request.get(base+'/postprocessor/js/legacy-spatial-operators.js')).status(),404);
    for(const name of ['CasePreprocessorStatic','CaseMPCModal','CaseReleaseTransient','CaseSpatialFrame','CaseSpatialTrussModal','CaseSpatialFrameTransient','CaseSectionBeam','CaseSpatialFrameSection']){
      let caseText=fs.readFileSync('examples/cases/'+(name==='CaseSpatialFrameSection'?'CaseSpatialFrame':name)+'.txt','utf8');
      if(name==='CaseSpatialFrameSection')caseText+='\nsections\n1\n'+JSON.stringify({id:1,name:'Manual J',type:'rectangle',dimensions:{b:.1,h:.2},torsion:{mode:'manual',value:.000123}})+'\nelement_sections\n1\n1,1,-37.5\n';
      const submit=await page.request.post(base+'/api/v1/jobs',{data:{name,caseText}});
      assert.equal(submit.status(),202,await submit.text());const id=(await submit.json()).id;
      await page.goto(base+'/postprocessor/?job='+id);
      await page.waitForFunction(()=>!document.querySelector('#download-json').disabled,{},{timeout:65000});
      assert.equal(await page.locator('#error-panel').isVisible(),false,name);
      const response=await page.request.get(base+'/api/v1/jobs/'+id+'/result');assert.equal(response.status(),200);
      const raw=await response.json();assert.equal(raw.version,8);assert(raw.model.recovery);
      if(name==='CaseSectionBeam'||name==='CaseSpatialFrameSection'){
        assert(raw.model.elements[0].sectionProperties.area>0);
        assert(raw.analysis.elementResults[0].sectionStress.positions.length>=2);
        if(name==='CaseSpatialFrameSection')assert.equal(raw.model.elements[0].sectionProperties.torsionConstant,.000123);
      }
      const spatial=raw.model.dimension===3;
      await page.locator(spatial?'#spatial-result-table':'#planar-result-objects select').selectOption('element');
      await page.locator(spatial?'#spatial-result-rows tr':'#planar-result-rows tr').first().click();
      const waiting=page.waitForEvent('download');await page.locator('#download-json').click();
      const filename='output/'+name+'.json';await(await waiting).saveAs(filename);assert.deepEqual(JSON.parse(fs.readFileSync(filename)),raw);
      await page.goto(base+'/postprocessor/');await page.locator('#file-input').setInputFiles(filename);
      await page.waitForFunction(()=>!document.querySelector('#download-json').disabled);
      assert.equal(await page.locator('#error-panel').isVisible(),false,'Saved result reopen');
      await page.request.delete(base+'/api/v1/jobs/'+id);
      console.log(name+': TXT -> nginx/API -> Octave v8 -> viewer -> saved/reopened JSON PASS');
    }
    for(const name of ['static-frame','modal-frame','transient-history','pinned-portal-uniform']){
      const response=await page.request.get(base+'/examples/results/'+name+'.json');
      assert.equal(response.status(),200);const buffer=await response.body();
      assert.equal(JSON.parse(buffer).version,8);
      await page.goto(base+'/postprocessor/');
      await page.locator('#file-input').setInputFiles({name:name+'.json',mimeType:'application/json',buffer});
      await page.waitForFunction(()=>!document.querySelector('#download-json').disabled);
      assert.equal(await page.locator('#error-panel').isVisible(),false,name);
      console.log(name+': shipped v8 example PASS');
    }
    for(let version=1;version<=7;version++){
      await page.goto(base+'/postprocessor/');
      await page.locator('#file-input').setInputFiles({name:'old.json',mimeType:'application/json',buffer:Buffer.from(JSON.stringify({format:'mkef-postprocessor',version}))});
      await page.waitForFunction(()=>document.querySelector('#error-panel').textContent.includes('Recalculate'));
      assert(await page.locator('#download-json').isDisabled());
    }
    assert.deepEqual(errors,[]);
    console.log('Rebuilt runtime stack: eight API/browser cases, four shipped examples and v1–v7 rejection PASS');
  }finally{await browser.close();}
})().catch(error=>{console.error(error);process.exitCode=1;});
