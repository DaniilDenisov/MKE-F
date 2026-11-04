// Run in the test Nginx network namespace. All numerical requests go to the real service.
const {chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE||'playwright');
const assert=require('node:assert/strict');
(async()=>{
  const browser=await chromium.launch({headless:true,args:['--no-sandbox']});
  try{
    const page=await browser.newPage(), errors=[],requests=[];
    page.setDefaultTimeout(150000);
    page.on('pageerror',e=>errors.push(e.message));
    page.on('request',r=>{if(r.url().includes('/api/v1/sections/'))requests.push([r.method(),r.url()]);});
    await page.goto('http://localhost:8080/reference/section-fem-lab.html');
    assert.equal(requests.length,0);
    await page.fill('#nu','.3');
    const run=()=>page.click('#run');
    const finished=()=>page.waitForFunction(()=>document.querySelector('#cancel').disabled);
    const mode=process.argv[2];
    if(mode==='down'){
      await run();await finished();
      assert.match(await page.textContent('#progress'),/недоступен/);
      assert.equal(await page.isVisible('#plot'),false);
    }else{
      const tube={type:'square_tube',dimensions:{b:1,t:.15}};
      await page.fill('#section',JSON.stringify(tube));await page.fill('#forces','[0,0,0,1,0,0]');
      await run();await finished();
      assert.equal(await page.isVisible('#plot'),true,await page.textContent('#progress'));
      assert.match(await page.textContent('#critical'),/Диагностический сеточный пик.*unavailable/);
      const fetched=requests.length;
      await page.selectOption('#component','1');
      assert.ok(Number(await page.textContent('#low'))<0 && Number(await page.textContent('#high'))>0);
      const inspect=()=>page.evaluate(()=>{
        const c=document.querySelector('#field'),ctx=c.getContext('2d');
        return {hole:ctx.getImageData(380,240,1,1).data[3],wall:ctx.getImageData(555,240,1,1).data[3],url:c.toDataURL()};
      });
      const plain=await inspect();assert.equal(plain.hole,0);assert.equal(plain.wall,255);
      await page.check('#mesh');assert.notEqual((await inspect()).url,plain.url);
      assert.equal(requests.length,fetched);
      assert.equal(requests.filter(([m,u])=>m==='GET'&&u.endsWith('/field')).length,1);

      // Hold an actual submission response until a newer run has finished.
      let release,posted;
      const gate=new Promise(r=>release=r),ready=new Promise(r=>posted=r);
      await page.route('**/sections/diagnostic',async route=>{
        const response=await route.fetch();posted((await response.json()).id);
        await gate;await route.fulfill({response});
      },{times:1});
      await run();const stale=await ready;
      await page.click('#cancel');assert.match(await page.textContent('#progress'),/отменён/);
      await page.fill('#forces','[1,0,0,0,0,0]');await run();await finished();
      const finite=await page.textContent('#critical');assert.match(finite,/Максимум по Мизесу/);
      const deleted=page.waitForResponse(r=>r.url().endsWith('/jobs/'+stale)&&r.request().method()==='DELETE');
      release();await deleted;
      assert.equal(await page.textContent('#critical'),finite);
      await page.goto('file:///app/reference/section-fem-lab.html');
      assert.equal(await page.isDisabled('#run'),true);assert.match(await page.textContent('#progress'),/Docker Compose/);
    }
    assert.deepEqual(errors,[]);
    console.log('Diagnostic browser PASS',mode);
  }finally{await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
