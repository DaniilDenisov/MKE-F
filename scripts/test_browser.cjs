// Run with a locally installed Playwright or MKEF_PLAYWRIGHT_MODULE pointing to it.
const { chromium } = require(process.env.MKEF_PLAYWRIGHT_MODULE || 'playwright');
const fs = require('fs');
const path = require('path');
const http = require('http');
const assert = require('assert/strict');
const root = path.resolve(__dirname, '..');
async function checkConstraintAvailability(page, available) {
  for (const id of ['add-support', 'add-mpc', 'mpc-axis-create']) {
    assert.equal(await page.locator('#' + id).isDisabled(), !available, id);
    if (!available) await page.locator('#' + id).dispatchEvent('click');
  }
  if (!available) {
    assert.equal(await page.locator('#supports-body tr').count(), 0, 'Created a support without elements');
    assert.equal(await page.locator('#mpcs-body > tr').count(), 0, 'Created an MPC without elements');
  }
}
async function checkCanvasTools(page) {
  const hint = page.locator('#canvas-tool-hint');
  await checkConstraintAvailability(page, false);
  assert(await page.locator('#add-load').isDisabled());
  await page.locator('#add-load').dispatchEvent('click');
  assert.equal(await page.locator('#loads-body tr').count(), 0);
  assert(await page.locator('#undo').isDisabled(), 'Blocked constraints changed undo history');
  assert.equal(await page.locator('.settings [data-tool]').count(), 0);
  assert.equal(await page.locator('.table-actions [data-tool]').count(), 5);
  assert.equal(await page.locator('#add-node').count(), 0);
  await page.locator('#grid-spacing').fill('0.5'); await page.locator('#grid-spacing').dispatchEvent('change');
  await page.locator('#place-node').click();
  assert.equal(await page.locator('#nodes-body tr').count(), 0, 'Activating placement inserted a node');
  assert.match(await hint.innerText(), /Click the canvas/);
  for (const point of [{x:-.5,y:0}, {x:.5,y:0}, {x:.5,y:.5}]) {
    const screen = await page.locator('#viewport').evaluate((svg, p) => {
      const screen = new DOMPoint(p.x, -p.y).matrixTransform(svg.getScreenCTM());
      return {x:screen.x,y:screen.y};
    }, point);
    await page.mouse.click(screen.x, screen.y);
  }
  assert.equal(await page.locator('#nodes-body tr').count(), 3);
  assert.deepEqual(await page.locator('#nodes-body tr').evaluateAll(rows => rows.map(row => [...row.querySelectorAll('input')].map(input => +input.value))), [[-.5,0],[.5,0],[.5,.5]]);
  await checkConstraintAvailability(page, false);
  await page.locator('#add-element').click();
  assert.match(await hint.innerText(), /select first node/);
  assert.equal(await page.locator('#selection-details').innerText(), 'Nothing selected');
  assert.equal(await page.locator('.model-node.selected, .model-node.pending, .selected-row').count(), 0);
  await page.locator('.model-node[data-node-id="3"]').dispatchEvent('click');
  assert.match(await hint.innerText(), /select second node/);
  await page.keyboard.press('Escape');
  assert.equal(await page.locator('.model-node.pending').count(), 0);
  assert.match(await hint.innerText(), /select first node/);
  assert.equal(await page.locator('#add-element').getAttribute('aria-pressed'), 'true');
  await page.locator('.model-node[data-node-id="3"]').dispatchEvent('click');
  await page.locator('.model-node[data-node-id="2"]').dispatchEvent('click');
  assert.equal(await page.locator('#elements-body tr').count(), 0);
  assert.match(await hint.innerText(), /select second node/, 'Rejected creation lost the first endpoint');
  await page.locator('#element-type').selectOption('113');
  await page.locator('#analysis-type').selectOption('static');
  await checkConstraintAvailability(page, false);
  await page.locator('.model-node[data-node-id="2"]').dispatchEvent('click');
  assert.equal(await page.locator('#elements-body tr').count(), 1);
  await checkConstraintAvailability(page, true);
  await page.locator('#add-support').click();
  assert.equal(await page.locator('#supports-body tr').count(), 0, 'Entering support mode created a support');
  assert.equal(await page.locator('#selection-details').innerText(), 'Nothing selected');
  assert.equal(await page.locator('#add-support').getAttribute('aria-pressed'), 'true');
  assert.match(await hint.innerText(), /Add support: select a node/);
  assert(await page.locator('#new-member-defaults').isHidden());
  await page.locator('.model-node[data-node-id="3"]').dispatchEvent('click');
  assert.equal(await page.locator('#supports-body tr').count(), 1);
  assert.equal(await page.locator('#supports-body tr input').first().inputValue(), '3');
  assert.match(await page.locator('#supports-body').textContent(), /ux, uy, thetaZ/);
  await page.keyboard.press('Escape');
  assert.equal(await page.locator('#selection-details').innerText(), 'Nothing selected');
  assert.equal(await page.locator('#supports-body tr').count(), 1);
  assert.equal(await page.locator('#add-support').getAttribute('aria-pressed'), 'true');
  await page.locator('#place-node').click();
  await page.locator('.model-node[data-node-id="2"]').dispatchEvent('click');
  assert.equal(await page.locator('#supports-body tr').count(), 1, 'Support placement remained active after tool change');
  await page.locator('#add-support').click();
  assert.equal(await page.locator('#selection-details').innerText(), 'Nothing selected');
  assert.equal(await page.locator('#supports-body tr').count(), 1, 'Re-entering support mode created a support');
  await page.locator('#undo').click();
  await page.locator('#undo').click();
  await checkConstraintAvailability(page, false);
  assert.equal(await page.locator('[data-tool="select"]').getAttribute('aria-pressed'), 'true');
  await page.locator('#redo').click();
  await checkConstraintAvailability(page, true);
  await page.locator('#add-element').click();
  assert.match(await hint.innerText(), /select first node/);
  assert.equal(await page.locator('#add-element').getAttribute('aria-pressed'), 'true');
  await page.locator('.model-node[data-node-id="1"]').dispatchEvent('click');
  assert.match(await hint.innerText(), /select second node/);
  await page.locator('.model-node[data-node-id="2"]').dispatchEvent('click');
  assert.equal(await page.locator('#elements-body tr').count(), 2);
  assert.match(await hint.innerText(), /select first node/);
  await page.locator('#undo').click(); assert.equal(await page.locator('#elements-body tr').count(), 1);
  assert.match(await hint.innerText(), /select first node/);
  await page.locator('#redo').click(); assert.equal(await page.locator('#elements-body tr').count(), 2);
  assert.match(await hint.innerText(), /select first node/);
  await page.locator('.model-node[data-node-id="1"]').dispatchEvent('click');
  await page.locator('#place-node').click();
  assert.equal(await page.locator('.model-node.pending').count(), 0);
  assert.match(await hint.innerText(), /Click the canvas/);
  await page.locator('[data-tool="select"]').click();
  assert.match(await hint.innerText(), /Select a node or element/);
  assert.equal(await page.locator('#place-node').getAttribute('aria-pressed'), 'false');
  assert.equal(await page.locator('#add-element').getAttribute('aria-pressed'), 'false');
  for (const target of ['.model-node[data-node-id="1"]', '.model-element[data-element-id="1"]']) {
    await page.locator(target).dispatchEvent('click');
    assert.notEqual(await page.locator('#selection-details').innerText(), 'Nothing selected');
    await page.keyboard.press('Escape');
    assert.equal(await page.locator('#selection-details').innerText(), 'Nothing selected');
    assert.equal(await page.locator('#viewport .selected, .selected-row').count(), 0);
    assert(await page.locator('#node-coordinate-editor').isHidden());
  }
  await page.locator('.model-element[data-element-id="1"]').dispatchEvent('click');
  await page.locator('#place-node').click();
  assert.equal(await page.locator('#selection-details').innerText(), 'Nothing selected');
  assert.equal(await page.locator('#viewport .selected, .selected-row').count(), 0);
  await checkNodalLoadTool(page);
  fs.mkdirSync(path.join(root,'output'), {recursive:true});
  await page.locator('#add-element').click();
  await page.screenshot({path:path.join(root,'output/preprocessor-canvas-tools.png'),fullPage:true});
  page.once('dialog', dialog => dialog.accept());
  await page.locator('#new-case').click();
  await checkConstraintAvailability(page, false);
  assert.equal(await page.locator('[data-tool="select"]').getAttribute('aria-pressed'), 'true');
  console.log('preprocessor tools: PASS (canvas placement, repeated members, rejected creation, persistent prompts, undo/redo, mode changes)');
}
async function checkNodalLoadTool(page) {
  for (const analysis of ['static', 'transient']) {
    await page.locator('#analysis-type').selectOption(analysis);
    await page.locator('[data-tool="select"]').click();
    await page.locator('.model-node[data-node-id="1"]').dispatchEvent('click');
    await page.locator('#add-load').click();
    assert.equal(await page.locator('#loads-body tr').count(), 0, 'Entering load mode created a load');
    assert.equal(await page.locator('#selection-details').innerText(), 'Nothing selected');
    assert.equal(await page.locator('#add-load').getAttribute('aria-pressed'), 'true');
    assert.match(await page.locator('#canvas-tool-hint').innerText(), /Add nodal load: select a node/);
    assert(await page.locator('#new-member-defaults').isHidden());
    await page.locator('.model-node[data-node-id="3"]').dispatchEvent('click');
    assert.equal(await page.locator('#loads-body tr').count(), 1);
    assert.equal(await page.locator('#loads-body tr input').first().inputValue(), '3');
    assert.equal(await page.locator('#loads-body tr select').inputValue(), analysis === 'static' ? '10' : '13');
    await page.keyboard.press('Escape');
    assert.equal(await page.locator('#selection-details').innerText(), 'Nothing selected');
    assert.equal(await page.locator('#loads-body tr').count(), 1);
    assert.equal(await page.locator('#add-load').getAttribute('aria-pressed'), 'true');
    await page.locator('#undo').click();
    assert.equal(await page.locator('#loads-body tr').count(), 0);
    await page.locator('#redo').click();
    assert.equal(await page.locator('#loads-body tr input').first().inputValue(), '3');
    await page.locator('[data-tool="select"]').click();
    await page.locator('.model-node[data-node-id="2"]').dispatchEvent('click');
    assert.equal(await page.locator('#loads-body tr').count(), 1, 'Load placement remained active after tool change');
    await page.locator('#add-load').click();
    assert.equal(await page.locator('#loads-body tr').count(), 1, 'Re-entering load mode created a load');
    await page.locator('#undo').click();
    await page.locator('#analysis-type').selectOption('modal');
    assert(await page.locator('#add-load').isDisabled());
    assert.equal(await page.locator('[data-tool="select"]').getAttribute('aria-pressed'), 'true');
    await page.locator('#add-load').dispatchEvent('click');
    await page.locator('.model-node[data-node-id="3"]').dispatchEvent('click');
    assert.equal(await page.locator('#loads-body tr').count(), 0, 'Modal analysis allowed load placement');
  }
  await page.locator('#analysis-type').selectOption('static');
  console.log('nodal load tool: PASS (explicit node, static/transient defaults, Escape, undo/redo, mode changes, modal guard)');
}
async function checkTableLayout(page) {
  const panels = page.locator('.tables > details');
  const initialOpen = await panels.evaluateAll(items => items.map(item => item.open));
  for (const viewport of [{width:1680,height:800}, {width:1280,height:720}, {width:390,height:844}]) {
    await page.setViewportSize(viewport);
    for (const mode of ['closed', 'elements', 'all']) {
      await panels.evaluateAll((items, mode) => items.forEach((item, i) => {
        item.open = mode === 'all' || (mode === 'elements' && i === 1);
      }), mode);
      const layout = await page.evaluate(() => {
        const tables = document.querySelector('.tables');
        return {
          height: tables.clientHeight, content: tables.scrollHeight,
          panels: [...tables.children].map(item => {
            const box = item.getBoundingClientRect();
            return {
              open: item.open, height: box.height,
              header: item.querySelector('summary').getBoundingClientRect().height,
              contentBottom: item.lastElementChild.getBoundingClientRect().bottom - box.top
            };
          })
        };
      });
      for (const panel of layout.panels) {
        assert(panel.height >= panel.header + 1, `Clipped header: ${viewport.width}/${mode}`);
        if (!panel.open) assert(Math.abs(panel.height - panel.header - 2) < 1, 'Closed panel stretched beside open panel');
        else assert(panel.contentBottom <= panel.height, 'Open panel clips its contents');
      }
      if (mode === 'all') {
        assert(layout.content > layout.height, 'Expanded panels must overflow the list region');
        const before = await page.locator('.canvas-wrap').boundingBox();
        const scrolled = await page.locator('.tables').evaluate(el => {
          el.scrollTop = el.scrollHeight;
          return el.scrollTop;
        });
        assert(scrolled > 0, 'List region must scroll');
        assert.deepEqual(await page.locator('.canvas-wrap').boundingBox(), before, 'Scrolling lists moves canvas');
      }
    }
  }
  await panels.evaluateAll((items, states) => items.forEach((item, i) => { item.open = states[i]; }), initialOpen);
  await page.locator('.tables').evaluate(el => { el.scrollTop = 0; });
  await page.setViewportSize({width:1440,height:1000});
  console.log('preprocessor layout: PASS (headers, panel contents, independent scrolling, desktop and mobile)');
}
require('child_process').execFileSync(process.execPath, [path.join(__dirname, 'generate-support-catalog.cjs'), '--check']);
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
    await checkCanvasTools(page);
    await page.locator('#file-input').setInputFiles(path.join(root,'examples/cases/CaseUniformFrame.txt'));
    await page.waitForFunction(() => document.querySelectorAll('#element-loads-body tr').length===1);
    await checkTableLayout(page);
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
    const count=await page.locator('.element-load-symbol').count(); assert.equal(count,22);
    await page.locator('#show-original').uncheck(); assert.equal(await page.locator('.element-load-symbol').count(),11);
    await page.locator('#show-original').check();
    await page.screenshot({path:path.join(root,'output/uniform-postprocessor.png'),fullPage:true});
    console.log('postprocessor integration: PASS (Octave JSON v2, curved moment diagram, original/deformed loads)');
    assert.deepEqual(errors,[]);
  } finally { await browser.close(); }
})().catch(error=>{console.error(error);process.exitCode=1;}).finally(()=>server.close());
