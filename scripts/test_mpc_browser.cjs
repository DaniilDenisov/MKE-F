const {captureEditor,capturePost,pickNode}=require('./browser-test-helpers.cjs');
// Generate fixtures first: setup; export_mpc_examples (GNU Octave).
const {chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE || 'playwright');
const {pathToFileURL}=require('url');
const path=require('path'), fs=require('fs'), assert=require('assert/strict');
const root=path.resolve(__dirname,'..');
(async()=>{
  const browser=await chromium.launch({headless:true,args:['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
  try {
    const page=await browser.newPage({viewport:{width:1500,height:1000}}), errors=[];
    page.on('pageerror',e=>errors.push(e.message));
    await page.goto(pathToFileURL(path.join(root,'preprocessor/index.html')).href);await captureEditor(page);
    await page.locator('#file-input').setInputFiles(path.join(root,'examples/cases/CaseMPCTruss.txt'));
    await page.waitForSelector('#condition-mpcs-rows tr');
    await page.locator('#condition-mpcs-rows tr').click();await page.locator('#selection-delete').click();
    await page.locator('[data-condition-tool=mpcs]').click();await page.locator('#condition-mode').selectOption('axis');
    await page.locator('#condition-node').fill('3');await page.locator('#condition-masterNode').fill('1');await page.locator('#condition-endNode').fill('2');await page.locator('#condition-preview').click();await page.locator('#condition-apply').click();
    assert.equal(await page.locator('#condition-mpcs-rows tr').count(),1);
    assert.match(await page.locator('#spatial-current').inputValue(),/3,2,0,2,1,2,0.5,2,2,0.5/);
    await page.locator('#spatial-node-rows tr').nth(2).click();await page.locator('#selection-delete').click();await page.locator('#spatial-delete-confirm').click();
    assert.equal(await page.locator('#condition-mpcs-rows tr').count(),0);
    await page.locator('#undo').click();assert.equal(await page.locator('#condition-mpcs-rows tr').count(),1);await page.locator('#redo').click();assert.equal(await page.locator('#condition-mpcs-rows tr').count(),0);await page.locator('#undo').click();
    await page.screenshot({path:path.join(root,'output/mpc-preprocessor.png'),fullPage:true});
    await page.goto(pathToFileURL(path.join(root,'postprocessor/index.html')).href);await capturePost(page);
    for (const name of ['CaseMPCFrame','CaseMPCModal','CaseMPCTransient','CaseMPCUniform','CaseMPCTruss','CaseMPCPartial']) {
      const json=JSON.parse(fs.readFileSync(path.join(root,'output',name+'.json'),'utf8'));
      json.metadata.title=name;
      await page.locator('#file-input').setInputFiles({name:name+'.json',mimeType:'application/json',buffer:Buffer.from(JSON.stringify(json))});
      await page.waitForFunction(title=>document.getElementById('dataset-title').textContent===title,name);
      await page.waitForFunction(()=>!document.getElementById('mpc-results').hidden);
      assert.equal(await page.locator('#error-panel').isVisible(),false,name+': '+await page.locator('#error-panel').textContent());
      assert.match(await page.locator('#mpc-multiplier').textContent(),/λ = /);
      await page.locator('#show-mpc-forces').uncheck(); await page.locator('#show-mpc-forces').check();
      assert.equal(await page.evaluate(()=>resultRenderer.scene.conditions.data.mpcs.length),1);
      const validation=await page.evaluate(data=>{
        const clone=JSON.parse(JSON.stringify(data)); clone.model.mpcs[0].rhs=1;
        try { window.MKEFPost.validateDataset(clone); return false; } catch(e) { return e.name==='DataError'; }
      },json);
      assert(validation,'Malformed MPC must report DataError');
      if (json.analysis.type==='transient') {
        await page.locator('#next-frame').click();
        assert(await page.locator('#mpc-history-chart .chart-series').count()>0);
      }
      if (name==='CaseMPCUniform') {
        await page.locator('#static-result').selectOption('M');
        assert(await page.evaluate(()=>resultRenderer.extra.distributed.positions.length>0));
      }
      if (name==='CaseMPCFrame') await page.screenshot({path:path.join(root,'output/mpc-postprocessor.png'),fullPage:true});
    }
    // Reopening a legacy file must reset the MPC panel.
    await page.locator('#file-input').setInputFiles(path.join(root,'examples/results/static-frame.json'));
    await page.waitForFunction(()=>document.getElementById('mpc-results').hidden);
    assert.equal(await page.locator('#mpc-results').isVisible(),false);
    assert.deepEqual(errors,[]);
    console.log('PASS MPC offline editor, helper, cascade, undo/redo, six Octave v3 exports, forces, histories and legacy reopening');
  } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
