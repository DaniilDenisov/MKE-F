// Generate output/CaseLinear*.json and CaseTriangleFrame.json with test_linear_loads.
const {chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE || 'playwright');
const {pathToFileURL}=require('url');
const path=require('path'), fs=require('fs'), assert=require('assert/strict');
const root=path.resolve(__dirname,'..');
const close=(actual,expected)=>assert(Math.abs(actual-expected)<1e-7*Math.max(1,Math.abs(expected)),`${actual} != ${expected}`);
(async()=>{
  const browser=await chromium.launch({headless:true});
  try {
    const page=await browser.newPage({viewport:{width:1440,height:1000}}), errors=[];
    page.on('pageerror',e=>errors.push(e.message));
    page.on('dialog',d=>d.accept());
    await page.goto(pathToFileURL(path.join(root,'preprocessor/index.html')).href);
    await page.locator('#file-input').setInputFiles(path.join(root,'examples/cases/CaseTriangleFrame.txt'));
    await page.waitForFunction(()=>document.querySelectorAll('#element-loads-body input').length===5);
    await page.locator('#element-loads-body').evaluate(e=>e.closest('details').open=true);
    assert.equal(await page.locator('#error-panel').isVisible(),false);
    assert.equal(await page.locator('.element-load-envelope').count(),1);
    const start=page.getByRole('spinbutton',{name:'Load 1 qy1 at node 1',exact:true});
    await start.fill('-300'); await start.dispatchEvent('change');
    assert.match(await page.locator('#case-preview').inputValue(),/21,1,2,0,-300,0,-1000/);
    await page.locator('#undo').click(); assert.equal(await start.inputValue(),'0');
    await page.locator('#redo').click(); assert.equal(await start.inputValue(),'-300');
    await page.locator('#element-loads-body select').first().selectOption('20');
    assert.match(await page.locator('#case-preview').inputValue(),/20,1,2,0,-650/);
    assert.equal(await page.locator('#element-loads-body input').count(),3);
    await page.locator('#undo').click(); assert.equal(await start.inputValue(),'-300');
    await page.locator('#element-loads-body select').nth(1).selectOption('1');
    assert.match(await page.locator('#case-preview').inputValue(),/21,1,1,0,-300,0,-1000/);
    await page.screenshot({path:path.join(root,'output/linear-preprocessor.png'),fullPage:true});
    await page.setViewportSize({width:390,height:844});
    assert(await start.isVisible());
    await page.screenshot({path:path.join(root,'output/linear-preprocessor-mobile.png'),fullPage:true});
    await page.setViewportSize({width:1440,height:1000});
    // Sign-changing endpoint values must produce finite graphics in both directions.
    await page.getByRole('spinbutton',{name:'Load 1 qy2 at node 2',exact:true}).fill('300');
    await page.getByRole('spinbutton',{name:'Load 1 qy2 at node 2',exact:true}).dispatchEvent('change');
    const directions=await page.locator('.element-load-symbol line').evaluateAll(lines=>lines.map(l=>Number(l.getAttribute('y2'))-Number(l.getAttribute('y1'))));
    assert(directions.some(y=>y>0) && directions.some(y=>y<0) && directions.every(Number.isFinite));
    console.log('linear preprocessor: PASS (endpoints, type/axes, undo/redo, sign change, mobile)');

    for (const name of ['CaseLinearFrame','CaseTriangleFrame','CaseLinearRelease','CaseLinearMPC']) {
      const filename=path.join(root,'output',name+'.json'), raw=JSON.parse(fs.readFileSync(filename,'utf8'));
      await page.goto(pathToFileURL(path.join(root,'postprocessor/index.html')).href);
      await page.locator('#file-input').setInputFiles(filename);
      await page.waitForFunction(()=>document.querySelectorAll('.element-load-envelope').length>0);
      assert.equal(await page.locator('#error-panel').isVisible(),false);
      await page.locator('#static-result').selectOption('M');
      const checks=await page.evaluate(raw=>{
        const M=window.MKEFPost,d=M.validateDataset(raw),S=M.staticResults;
        return d.raw.model.elements.map(e=>{
          const f=d.resultsByElementId.get(e.id).localEndForces;
          return {end:[S.elementDiagram(d,e,'N',1),S.elementDiagram(d,e,'V',1),S.elementDiagram(d,e,'M',1)],expected:[f[3],f[4],f[5]],maxM:Math.max(...S.diagramSamples(d,e,'M',6).map(s=>Math.abs(s.value)))};
        });
      },raw);
      checks.forEach(c=>c.end.forEach((v,i)=>close(v,c.expected[i])));
      assert(await page.locator('.element-load-label').count()>0);
      const both=await page.locator('.element-load-envelope').count();
      await page.locator('#show-original').uncheck();
      assert.equal(await page.locator('.element-load-envelope').count(),both/2);
      await page.locator('#show-original').check();
      await page.screenshot({path:path.join(root,'output',name+'-browser.png'),fullPage:true});
      console.log(name+': PASS (v5 import, end forces, diagrams, original/deformed loads)');
    }
    assert.deepEqual(errors,[]);
  } finally { await browser.close(); }
})().catch(error=>{console.error(error);process.exitCode=1;});
