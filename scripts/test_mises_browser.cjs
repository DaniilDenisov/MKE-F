const {chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE||'playwright');
const assert=require('node:assert/strict'),fs=require('node:fs');
const hosted=process.env.MKEF_POST_BASE_URL;
(async()=>{
  const browser=await chromium.launch({headless:true,args:['--no-sandbox','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
  try{
    const page=await browser.newPage(),errors=[],posts=[];
    page.on('pageerror',e=>errors.push(e.message));
    page.on('request',r=>{if(r.url().endsWith('/sections/member')&&r.method()==='POST')posts.push(r.postDataJSON());});
    page.setDefaultTimeout(150000);
    await page.goto((hosted||'file:///app')+'/postprocessor/index.html');
    const load=async name=>{
      const data=JSON.parse(fs.readFileSync('output/'+name+'.json'));
      data.metadata.title=name;
      await page.locator('#file-input').setInputFiles({name:name+'.json',mimeType:'application/json',buffer:Buffer.from(JSON.stringify(data))});
      await page.waitForFunction(n=>document.querySelector('#dataset-title').textContent===n,name);
      return data;
    };
    const run=async()=>{
      await page.locator('#member-maxima summary').click();
      await page.click('#member-run');
      await page.waitForFunction(()=>!document.querySelector('#member-run').disabled);
      return page.locator('#member-rows tr').first().locator('td').allTextContents();
    };
    await load('CaseSectionBeam');assert.equal(posts.length,0);
    if(hosted){
      const row=await run();assert.equal(row[4],'OK',row.join(' | '));
      assert.ok(Math.abs(Number(row[1])/750000-1)<.01);assert.ok(Math.abs(Number(row[2])-.5)<.01);
      assert.equal(row[3],'—');
      await page.fill('#member-yield','1500000');
      assert.ok(Math.abs(Number(await page.locator('#member-rows td').nth(3).textContent())-.5)<.005);
      assert.equal(posts.length,1);
      await page.fill('#member-yield','0');assert.equal(await page.locator('#member-rows td').nth(3).textContent(),'—');
      await page.fill('#member-nu','.25');assert.equal(await page.locator('#member-rows tr').count(),0);

      const spatial=await load('CaseSpatialFrameSection'),r=await run();
      assert.deepEqual(posts.at(-1).coefficients,spatial.analysis.elementResults[0].sectionForceCoefficients);
      assert.equal(posts.at(-1).nu,.25);assert.equal(r[4],'OK',r.join(' | '));
      const truss=await load('CaseSpatialTrussSection'),t=await run();
      const expected=Math.abs(truss.analysis.elementResults[0].localEndForces[0])/truss.model.elements[0].properties.area;
      assert.equal(t[4],'OK');assert.ok(Math.abs(Number(t[1])-expected)<=1e-5*Math.max(1,expected));
      await load('Mises2D_square_tube');const sharp=await run();
      assert.match(sharp[4],/Sharp corner/);assert.match(sharp[1],/\(mesh\)/);
      assert.ok(parseFloat(sharp[1])>0 && Number(sharp[2])>=0 && Number(sharp[2])<=1);
      await page.fill('#member-yield','250000000');assert.equal(await page.locator('#member-rows td').nth(3).textContent(),'—');
      assert.equal(await page.locator('#member-rounded').count(),0);
      for(const name of ['Rounded2D_square_tube','Rounded3D_rectangular_tube']){
        const rounded=await load(name),row=await run();
        assert.equal(posts.at(-1).section.dimensions.inner_radius,rounded.model.sections[0].dimensions.inner_radius);
        assert.match(row[4],/^OK · rᵢ=/);assert.ok(Number(row[1])>0);
        await page.fill('#member-yield','250000000');assert.ok(Number(await page.locator('#member-rows td').nth(3).textContent())>0);
      }
      // A cancelled submission can return its job id after the next file loads.
      await load('CaseSectionBeam');let release,ready;
      const gate=new Promise(r=>release=r),posted=new Promise(r=>ready=r);
      await page.route('**/sections/member',async route=>{
        const response=await route.fetch();ready((await response.json()).id);await gate;await route.fulfill({response});
      },{times:1});
      await page.locator('#member-maxima summary').click();await page.click('#member-run');
      const id=await posted;await page.click('#member-cancel');await load('CaseSpatialFrameModalSection');
      const removed=page.waitForResponse(r=>r.url().endsWith('/jobs/'+id)&&r.request().method()==='DELETE');
      release();await removed;
    }else{
      await page.locator('#member-maxima summary').click();assert.ok(await page.isDisabled('#member-run'));
      assert.match(await page.textContent('#member-progress'),/Docker Compose/);
      await load('CaseSpatialFrameModalSection');
    }
    assert.equal(await page.isVisible('#member-maxima'),false);
    assert.deepEqual(errors,[]);console.log('Member maxima, dataset lifecycle and yield comparison PASS',hosted||'offline');
  }finally{await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
