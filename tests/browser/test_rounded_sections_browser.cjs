// Run against the source-built Compose stack, after test_spatial_sections.m.
const {chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE||'playwright');
const {captureEditor}=require('../helpers/browser-test-helpers.cjs');
const assert=require('node:assert/strict'),fs=require('node:fs');
const base=process.env.MKEF_POST_BASE_URL||'http://localhost:8080';
(async()=>{
  const browser=await chromium.launch({headless:true,args:['--no-sandbox','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
  try{
    const page=await browser.newPage({viewport:{width:1450,height:1050}}),errors=[],requests=[];
    page.setDefaultTimeout(150000);page.on('pageerror',e=>errors.push(e.message));page.on('dialog',d=>d.accept());
    page.on('request',r=>{if(r.method()==='POST'&&r.url().endsWith('/sections/member'))requests.push(r.postDataJSON());});
    for(const [dim,type]of [[2,'square_tube'],[3,'rectangular_tube']]){
      await page.goto(base+'/preprocessor/');await captureEditor(page);
      const name=`Rounded${dim}D_${type}`;
      await page.locator('#file-input').setInputFiles('output/'+name+'.txt');
      await page.locator('#spatial-catalog-toggle').click();
      const edit=()=>page.locator('#spatial-catalog').getByRole('button',{name:'Edit',exact:true}).first().click();
      await edit();assert.equal(Number(await page.locator('#spatial-section-inner_radius').inputValue()),.005);
      const before=await page.locator('#spatial-current').inputValue();
      await page.fill('#spatial-section-inner_radius','.045');await page.click('#spatial-section-apply');
      assert.match(await page.textContent('#spatial-catalog-message'),/Inside radius/);
      assert.equal(await page.locator('#spatial-current').inputValue(),before);
      await page.fill('#spatial-section-inner_radius','.007');
      const preview=page.locator('#spatial-section-preview path');
      assert.equal(await preview.evaluate(p=>p.isPointInFill(new DOMPoint(0,0))),false);
      assert.equal(await preview.evaluate((p,h)=>p.isPointInFill(new DOMPoint(.05,h/2)),dim===2?.1:.2),false);
      await page.click('#spatial-section-apply');
      const saved=await page.locator('#spatial-current').inputValue();assert.match(saved,/"inner_radius":0.007/);
      const expected=await page.evaluate(()=>MKEFSpatialSections.properties(editorPanel.model.sections[0]));
      await page.click('#spatial-catalog-close');await page.click('#undo');
      assert.equal(await page.locator('#spatial-current').inputValue(),before);
      await page.click('#redo');assert.equal(await page.locator('#spatial-current').inputValue(),saved);
      await page.locator('#file-input').setInputFiles({name:name+'.txt',mimeType:'text/plain',buffer:Buffer.from(saved)});
      await page.waitForFunction(()=>editorPanel.model.sections[0].dimensions.inner_radius===.007);
      const beforeRequests=requests.length;
      await page.click('#run-analysis');await page.waitForSelector('#spatial-open-result:visible');
      const href=await page.locator('#spatial-open-result').getAttribute('href');
      const job=new URL(href,base+'/preprocessor/').searchParams.get('job');
      const result=await (await page.request.get(base+'/api/v1/jobs/'+job+'/result')).json();
      assert.equal(result.model.sections[0].dimensions.inner_radius,.007);
      for(const key of dim===3?['area','Iy','Iz','torsionConstant']:['area','momentOfInertia']){
        assert.ok(Math.abs(result.model.elements[0].properties[key]/expected[key]-1)<1e-10,key);
      }
      fs.writeFileSync('output/'+name+'-catalog.json',JSON.stringify(result));
      await page.goto(base+'/postprocessor/?job='+job);await page.waitForSelector('#member-maxima summary');
      assert.equal(requests.length,beforeRequests);assert.equal(await page.locator('#member-rounded').count(),0);
      await page.locator('#member-maxima summary').click();await page.fill('#member-yield','250000000');
      await page.click('#member-run');await page.waitForFunction(()=>!document.querySelector('#member-run').disabled);
      const row=await page.locator('#member-rows tr').first().locator('td').allTextContents();
      assert.match(row[4],/^OK · rᵢ=0.007/);assert.ok(Number(row[1])>0&&Number(row[3])>0,row.join(' | '));
      assert.equal(requests.at(-1).section.dimensions.inner_radius,.007);
      if(dim===3)assert.equal(requests.at(-1).J_model,result.model.elements[0].properties.torsionConstant);
      await page.screenshot({path:'output/'+name+'-catalog.png'});
    }
    assert.deepEqual(errors,[]);console.log('Rounded catalog: validation, preview hole, undo/redo, TXT, structural solve, result JSON and Python maxima PASS');
  }finally{await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
