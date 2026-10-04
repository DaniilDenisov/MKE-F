// Run with a locally installed Playwright or MKEF_PLAYWRIGHT_MODULE pointing to it.
const { chromium } = require(process.env.MKEF_PLAYWRIGHT_MODULE || 'playwright');
const fs = require('fs');
const path = require('path');
const http = require('http');
const assert = require('assert/strict');
const root = path.resolve(__dirname, '..');
const server = http.createServer((req, res) => {
  const pathname = decodeURIComponent(new URL(req.url, 'http://localhost').pathname);
  let file = path.resolve(root, '.' + pathname);
  if (!file.startsWith(root + path.sep)) { res.writeHead(403); res.end(); return; }
  if (fs.existsSync(file) && fs.statSync(file).isDirectory()) file = path.join(file, 'index.html');
  if (!fs.existsSync(file)) { res.writeHead(404); res.end(); return; }
  res.setHeader('Content-Type', { '.html':'text/html', '.js':'text/javascript', '.css':'text/css', '.json':'application/json' }[path.extname(file)] || 'text/plain');
  res.end(fs.readFileSync(file));
});
(async () => {
  await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
  const base = 'http://127.0.0.1:' + server.address().port;
  const browser = await chromium.launch({ headless:true });
  try {
    const page = await browser.newPage({ viewport:{width:1440,height:1000} });
    const errors=[]; page.on('pageerror', error => errors.push(error.message));
    for (const suite of ['preprocessor','postprocessor','shared']) {
      await page.goto(base+'/'+suite+'/tests/index.html');
      if (suite==='preprocessor') await page.waitForFunction(() => document.getElementById('test-status').dataset.status !== 'running');
      else await page.waitForFunction(() => ['passed','failed'].includes(document.body.getAttribute('data-test-status')));
      const report=await page.locator('#test-output').innerText();
      console.log(suite+': '+report.split('\n').filter(line=>line.startsWith('PASS')).length+' passed'+(report.includes('FAIL')?'\n'+report:''));
      assert(!report.includes('FAIL'),report);
    }
    await page.goto(base+'/preprocessor/index.html');
    await page.locator('#file-input').setInputFiles(path.join(root,'examples/cases/CaseUniformFrame.txt'));
    await page.waitForFunction(() => document.querySelectorAll('#element-loads-body tr').length===1);
    await page.locator('details').filter({has:page.locator('#nodes-body')}).locator('summary').click().catch(()=>{});
    // Canvas selection gives access to the second coordinate editor.
    await page.locator('#viewport [data-node-id="2"]').dispatchEvent('click');
    await page.locator('#selected-node-x').fill('3'); await page.locator('#selected-node-x').dispatchEvent('change');
    assert.equal(await page.locator('#element-loads-body tr').count(),0);
    await page.locator('#undo').click(); assert.equal(await page.locator('#element-loads-body tr').count(),1);
    await page.locator('#redo').click(); assert.equal(await page.locator('#element-loads-body tr').count(),0);
    await page.locator('#undo').click();
    // Open nodes table and edit the same field through its table binding.
    await page.locator('#nodes-body').evaluate(el=>{el.closest('details').open=true;});
    const coordinate=page.locator('#nodes-body tr').nth(1).locator('input').first();
    await coordinate.fill('4'); await coordinate.dispatchEvent('change');
    assert.equal(await page.locator('#element-loads-body tr').count(),0); await page.locator('#undo').click();
    await page.locator('#nodes-body tr').first().getByRole('button',{name:'Delete node 1',exact:true}).click();
    assert.equal(await page.locator('#elements-body tr').count(),0);
    assert.equal(await page.locator('#supports-body tr').count(),0);
    assert.equal(await page.locator('#nodes-body tr').count(),1);
    await page.locator('#undo').click();
    assert.equal(await page.locator('#element-loads-body tr').count(),1);
    await page.locator('#analysis-type').selectOption('modal');
    assert(await page.locator('#download-case').isDisabled()); assert.equal(await page.locator('#element-loads-body tr').count(),1);
    await page.locator('#analysis-type').selectOption('static');
    assert(!(await page.locator('#download-case').isDisabled()));
    await page.locator('#viewport [data-element-id="1"].model-element').dispatchEvent('click');
    await page.locator('#add-element-load').click(); assert.equal(await page.locator('#element-loads-body tr').count(),2);
    await page.locator('#undo').click();
    fs.mkdirSync(path.join(root,'output'),{recursive:true});
    await page.screenshot({path:path.join(root,'output/uniform-preprocessor.png'),fullPage:true});
    console.log('preprocessor integration: PASS (both editors, cascade, undo/redo, analysis switch, creation)');
    const resultPath=path.join(root,'output/uniform-result.json');
    if (!fs.existsSync(resultPath)) throw new Error('Export output/uniform-result.json using CaseUniformFrame.txt before running this suite.');
    await page.goto(base+'/postprocessor/index.html');
    await page.locator('#file-input').setInputFiles(resultPath);
    await page.waitForFunction(()=>document.querySelectorAll('.element-load-symbol').length>0);
    await page.locator('#static-result').selectOption('M');
    assert.equal(await page.locator('.load-symbol[data-node-id]').count(),0);
    const count=await page.locator('.element-load-symbol').count(); assert.equal(count,14);
    await page.locator('#show-original').uncheck(); assert.equal(await page.locator('.element-load-symbol').count(),7);
    await page.locator('#show-original').check();
    await page.screenshot({path:path.join(root,'output/uniform-postprocessor.png'),fullPage:true});
    console.log('postprocessor integration: PASS (Octave JSON v2, curved moment diagram, original/deformed loads)');
    assert.deepEqual(errors,[]);
  } finally { await browser.close(); }
})().catch(error=>{console.error(error);process.exitCode=1;}).finally(()=>server.close());
