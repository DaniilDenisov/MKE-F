// Run with a locally installed Playwright or MKEF_PLAYWRIGHT_MODULE pointing to it.
const { chromium } = require(process.env.MKEF_PLAYWRIGHT_MODULE || 'playwright');
const fs = require('fs');
const path = require('path');
const http = require('http');
const assert = require('assert/strict');
const root = path.resolve(__dirname, '../..');
require('child_process').execFileSync(process.execPath, [path.join(root, 'scripts/build/generate-support-catalog.cjs'), '--check']);
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
  const browser = await chromium.launch({headless:true,args:['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
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
    // UI integration lives in test_unified_editor_browser, test_unified_post_browser and test_unified_navigation_browser.
    assert.deepEqual(errors,[]);
  } finally { await browser.close(); }
})().catch(error=>{console.error(error);process.exitCode=1;}).finally(()=>server.close());
