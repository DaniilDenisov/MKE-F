// First run setup; run_nafems_challenge5('quick') in Octave.
const {chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE || 'playwright');
const assert=require('node:assert/strict'), fs=require('node:fs'), path=require('node:path');
const {pathToFileURL}=require('node:url');
const root=path.resolve(__dirname,'..');
const output=path.join(root,'output','nafems-challenge-5','browser');
fs.mkdirSync(output,{recursive:true});
(async()=>{
  const browser=await chromium.launch({headless:true});
  try {
    const page=await browser.newPage({viewport:{width:1440,height:1000}}), errors=[];
    page.on('pageerror',e=>errors.push(e.message));
    page.on('dialog',d=>d.accept());
    await page.goto(pathToFileURL(path.join(root,'preprocessor/index.html')).href);
    const source=path.join(root,'examples/cases/nafems-challenge-5/frame_n016.txt');
    await page.locator('#file-input').setInputFiles(source);
    await page.waitForFunction(()=>document.querySelectorAll('.end-release').length===4);
    assert.equal(await page.locator('#error-panel').isVisible(),false);
    const snapshot=()=>page.evaluate(()=>{
      const f=window.MKEFPre.caseFormat, m=f.parse(document.querySelector('#case-preview').value);
      return JSON.parse(JSON.stringify(m,(key,value)=>key==='sourceLine'?undefined:value));
    });
    const before=await snapshot();
    assert.equal(before.nodes.length,33); assert.equal(before.elements.length,32);
    assert.equal(before.releases.length,4);
    const downloadPromise=page.waitForEvent('download');
    await page.locator('#download-case').click();
    const download=await downloadPromise, saved=path.join(output,'frame_n016_roundtrip.txt');
    await download.saveAs(saved);
    await page.locator('#file-input').setInputFiles(saved);
    assert.deepEqual(await snapshot(),before);
    await page.screenshot({path:path.join(output,'preprocessor.png'),fullPage:true});
    await page.goto(pathToFileURL(path.join(root,'postprocessor/index.html')).href);
    const json=path.join(root,'output/nafems-challenge-5/quick/frame_n016.json');
    const raw=JSON.parse(fs.readFileSync(json,'utf8'));
    assert.equal(raw.version,4);
    assert.equal(raw.model.dofRegistry.filter(d=>d.kind==='elementEnd').length,4);
    await page.locator('#file-input').setInputFiles(json);
    await page.waitForFunction(()=>document.querySelectorAll('#mode-number option').length>3);
    assert.equal(await page.locator('#error-panel').isVisible(),false);
    assert.equal(await page.locator('[data-layer="original-geometry"] .end-release').count(),4);
    for(const index of [0,1,6,7,11,12]) {
      const value=await page.locator('#mode-number option').nth(index).getAttribute('value');
      await page.locator('#mode-number').selectOption(value);
      assert(await page.evaluate(({raw,index})=>!window.MKEFPost.modalView.display(window.MKEFPost.validateDataset(raw),index,1,21).zero,{raw,index}));
    }
    await page.screenshot({path:path.join(output,'postprocessor.png'),fullPage:true});
    for(const language of ['RU','EN']) {
      await page.goto(pathToFileURL(path.join(root,'reference',language,'06a-nafems-challenge-5.html')).href);
      assert.equal(await page.locator('main.chapter a[href$=".txt"]').count(),20);
      assert.equal(await page.locator('nav.page-nav').count(),1);
      assert.equal(await page.locator('[aria-labelledby="comparison-title"] tbody tr').count(),7);
      assert.match(await page.locator('[aria-labelledby="comparison-title"]').innerText(),/0\.8\.1/);
      assert.match(await page.locator('[aria-labelledby="comparison-title"]').innerText(),/0\.8\.2/);
      assert.match(await page.locator('[aria-labelledby="comparison-title"]').innerText(),/Unknown section marker/);
      assert.match(await page.locator('main').innerText(),/run_nafems_challenge5/);
      await page.screenshot({path:path.join(output,'reference-'+language+'.png'),fullPage:true});
      await page.setViewportSize({width:390,height:1000});
      assert(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1));
      await page.screenshot({path:path.join(output,'reference-'+language+'-mobile.png'),fullPage:true});
      await page.setViewportSize({width:1440,height:1000});
    }
    for(const plot of ['axial_convergence','bending_convergence','frame_convergence','spectrum','frame_mode_001']) {
      await page.goto(pathToFileURL(path.join(root,'output/nafems-challenge-5/full',plot+'.svg')).href);
      await page.screenshot({path:path.join(output,plot+'.png')});
    }
    assert.deepEqual(errors,[]);
    console.log('PASS NAFEMS saved-case download/reimport, releases, v4 modes and RU/EN instructions');
  } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
