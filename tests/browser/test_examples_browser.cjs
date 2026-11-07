const {chromium} = require(process.env.MKEF_PLAYWRIGHT_MODULE || 'playwright');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const http = require('node:http');
const root = path.resolve(__dirname, '../..');
const output = path.join(root, 'output', 'load-examples');

(async () => {
  fs.mkdirSync(output, {recursive: true});
  const server = http.createServer((req, res) => {
    const file = path.resolve(root, '.' + new URL(req.url, 'http://localhost').pathname);
    if (!file.startsWith(root + path.sep) || !fs.existsSync(file) || !fs.statSync(file).isFile()) {
      res.writeHead(404).end(); return;
    }
    res.setHeader('Content-Type', {'.html':'text/html', '.js':'text/javascript', '.css':'text/css', '.png':'image/png'}[path.extname(file)] || 'text/plain');
    fs.createReadStream(file).pipe(res);
  });
  await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
  const browser = await chromium.launch({headless:true, args:['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
  try {
    const page = await browser.newPage({viewport:{width:1600,height:1000}}), errors = [], requests = [];
    page.on('pageerror', e => errors.push(e.message));
    page.on('request', r => { if (/^https?:/.test(r.url())) requests.push(r.url()); });
    await page.goto('file://' + path.join(root, 'preprocessor/index.html'));
    const select = page.getByRole('combobox', {name:'Load example', exact:true});
    assert.equal(await select.locator('optgroup').count(), 6);
    for (const group of await select.locator('optgroup').all()) assert.equal(await group.locator('option').count(), 2);
    const examples = await page.evaluate(() => MKEFPre.examples.map(e => ({id:e.id, name:e.name, dimension:e.dimension, analysis:e.analysis})));
    assert.equal(examples.length, 12);
    assert.equal(new Set(examples.map(e => e.id)).size, 12);
    for (const example of examples) {
      await select.selectOption(example.id);
      assert.equal(await select.inputValue(), '');
      assert.equal(await page.locator('#case-name').inputValue(), example.name);
      assert.equal(await page.locator('#analysis-type').inputValue(), example.analysis);
      assert.equal(await page.locator('#validation-status').innerText(), 'Model is ready to run');
      assert(await page.locator('#undo').isDisabled());
      assert(await page.locator('#redo').isDisabled());
      assert.equal(await page.locator('#selection-details').innerText(), 'Nothing selected');
      const parsed = await page.evaluate(() => {
        const F=MKEFPre.editorAdapter, model=F.parse(document.getElementById('spatial-current').value);
        return {dimension:model.dimension, errors:F.validate(model,true), loads:model.loads.length, supports:model.supports.length};
      });
      assert.equal(parsed.dimension, example.dimension);
      assert.deepEqual(parsed.errors, []);
      assert(parsed.supports > 0);
      assert.equal(parsed.loads, example.analysis === 'modal' ? 0 : 1);
      const downloading = page.waitForEvent('download');
      await page.locator('#download-case').click();
      const download = await downloading;
      assert.equal(download.suggestedFilename(), example.name + '.txt');
      await download.saveAs(path.join(output, example.id + '.txt'));
    }
    assert.deepEqual(requests, [], 'Built-in examples must work offline');

    // Replacing an edited model must preserve it on cancellation, then reset history on acceptance.
    await page.locator('#case-name').fill('My unsaved model');
    await page.locator('#case-name').dispatchEvent('change');
    const before = await page.locator('#spatial-current').inputValue();
    page.once('dialog', d => d.dismiss());
    await select.selectOption('2d-cantilever-static');
    assert.equal(await page.locator('#case-name').inputValue(), 'My unsaved model');
    assert.equal(await page.locator('#spatial-current').inputValue(), before);
    assert.equal(await select.inputValue(), '');
    assert(!(await page.locator('#undo').isDisabled()));
    page.once('dialog', d => d.accept());
    await select.selectOption('2d-cantilever-static');
    assert.equal(await page.locator('#case-name').inputValue(), '2D Cantilever beam static');
    assert(await page.locator('#undo').isDisabled());
    // A fresh selection restores the original sample even after editing and saving it.
    await page.locator('#case-name').fill('Edited sample');
    await page.locator('#case-name').dispatchEvent('change');
    const saved = page.waitForEvent('download'); await page.locator('#download-case').click(); await saved;
    await select.selectOption('2d-cantilever-static');
    assert.equal(await page.locator('#case-name').inputValue(), '2D Cantilever beam static');
    const reference = page.getByRole('link', {name:'Examples', exact:true});
    assert.match(await reference.getAttribute('href'), /10-input-format.html#all-examples$/);
    for (const width of [1250,1600]) {
      await page.setViewportSize({width,height:1000});
      const box = await select.boundingBox(); assert(box.x >= 0 && box.x + box.width <= width);
      assert(await reference.isVisible());
    }
    await page.screenshot({path:path.join(output, 'selector-2d.png')});
    await select.selectOption('3d-tripod-transient');
    await page.screenshot({path:path.join(output, 'selector-3d.png')});

    // Hosted UI: each selection can be submitted immediately, with the correct native TXT.
    const job = {id:'11111111-1111-4111-8111-111111111111',status:'succeeded'};
    let submitted;
    await page.route('**/api/v1/**', async route => {
      const request = route.request();
      if (request.url().endsWith('/health')) return route.fulfill({json:{status:'ready',octaveVersion:'test'}});
      if (request.method() === 'POST') submitted = request.postDataJSON();
      await route.fulfill({json:job});
    });
    await page.goto('http://127.0.0.1:' + server.address().port + '/preprocessor/index.html');
    for (const example of examples) {
      await select.selectOption(example.id);
      await page.waitForFunction(() => !document.getElementById('run-analysis').disabled);
      const native = await page.evaluate(() => {
        const F = MKEFPre.editorAdapter;
        return F.serialize(F.parse(document.getElementById('spatial-current').value));
      });
      await page.locator('#run-analysis').click();
      await page.waitForFunction(() => document.getElementById('run-status').textContent === 'Completed');
      assert.equal(submitted.name, example.name);
      assert.equal(submitted.caseText, native);
    }
    assert.deepEqual(errors, []);
    console.log('PASS: 12 offline examples, grouped selector, download/round-trip, unsaved-change guard, fresh reload, layout and hosted submission.');
  } finally {
    await browser.close();
    await new Promise(resolve => server.close(resolve));
  }
})().catch(error => { console.error(error); process.exitCode=1; });
