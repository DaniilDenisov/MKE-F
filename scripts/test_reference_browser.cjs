// Requires local Playwright (or MKEF_PLAYWRIGHT_MODULE). Always checks file://;
// set MKEF_REFERENCE_BASE_URL to additionally check the Docker/Nginx host.
const { chromium } = require(process.env.MKEF_PLAYWRIGHT_MODULE || 'playwright');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const { pathToFileURL, fileURLToPath } = require('node:url');
const root = path.resolve(__dirname, '..');
const reference = path.join(root, 'reference');
const pages = fs.readdirSync(path.join(reference, 'RU')).filter(name => name.endsWith('.html') && name !== '04c-nafems-challenge-5.html');
const visualizers = pages.filter(name => /^(13|14|15)-/.test(name));
const output = path.join(root, 'output', 'reference-browser');
fs.mkdirSync(output, { recursive: true });
const bases = [pathToFileURL(root + path.sep).href];
if (process.env.MKEF_REFERENCE_BASE_URL) bases.push(process.env.MKEF_REFERENCE_BASE_URL.replace(/\/$/, '') + '/');

function assertEnglish(text, description) {
  assert(!/[\u0400-\u04ff]/u.test(text), `Untranslated text: ${description}`);
}

async function numericSnapshot(page, chapter) {
  return page.evaluate(chapter => {
    if (chapter === '13') return { matrix: displayedMatrix(), result: solution() };
    if (chapter === '14') return { local: localMatrices(), K: assembledMatrix('K', 2), M: assembledMatrix('M', 2), effective: effectiveMatrix(), current: displayContext().matrix };
    const matrix = matrixFrom(allTriplets());
    const reduced = matrix.slice(1).map(row => row.slice(1));
    const load = Array(state.elementCount).fill(0);
    load[load.length - 1] = 30;
    return { matrix, visible: matrixFrom(rawTriplets()), solution: solveLinearSystem(reduced, load) };
  }, chapter);
}

async function exercise(page, chapter, english) {
  const snapshots = [];
  const steps = chapter === '13' ? 7 : chapter === '14' ? 8 : 5;
  const check = async () => {
    if (english) assertEnglish(await page.locator('body').innerText(), `${chapter} dynamic state`);
    snapshots.push(await numericSnapshot(page, chapter));
  };
  for (let step = 0; step < steps; step++) {
    await page.locator(`[data-stage="${step}"]`).click();
    await check();
    if (chapter !== '13') {
      await page.locator('[data-matrix="M"]').click();
      await check();
      await page.locator('[data-matrix="K"]').click();
    }
  }
  const rangeId = chapter === '13' ? 'stiffnessRange1' : chapter === '14' ? 'timeStepRange' : 'elementRange';
  for (const limit of ['min', 'max']) {
    await page.locator('#' + rangeId).evaluate((input, limit) => {
      input.value = input[limit];
      input.dispatchEvent(new Event('input', { bubbles: true }));
    }, limit);
    await check();
  }
  if (chapter === '13') {
    await page.locator('[data-view="sparse"]').click();
    await page.locator('#loadRange').fill('60');
    await check();
  }
  if (chapter === '14') {
    await page.locator('[data-matrix="Keff"]').click();
    await check();
    await page.locator('#matrixTable td button').first().click();
    if (english) assertEnglish(await page.locator('#cellDetail').innerText(), 'matrix cell details');
  }
  const elementButton = chapter === '15' ? 'button[data-element="1"]' : '[data-element-button="1"]';
  await page.locator(elementButton).click();
  await check();
  if (chapter !== '15') {
    await page.locator('#resetButton').click();
    await page.locator('#playButton').click();
    await page.waitForFunction(() => state.stage > 0);
    await page.locator('#playButton').click();
    assert.equal(await page.evaluate(() => state.playing), false);
  }
  await page.locator('#resetButton').click();
  assert.equal(await page.evaluate(() => state.stage), 0);
  return snapshots;
}

(async () => {
  const browser = await chromium.launch({ headless: true });
  try {
    const errors = [];
    const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });
    page.on('pageerror', error => errors.push(error.message));
    page.on('console', message => {
      if (message.type() === 'error' && /Content Security Policy|Refused to execute/.test(message.text())) errors.push(message.text());
    });
    for (const base of bases) {
      const mode = base.startsWith('file:') ? 'offline' : 'nginx';
      for (const name of pages) {
        for (const language of ['RU', 'EN']) {
          const response = await page.goto(base + `reference/${language}/${name}`);
          if (response) assert(response.ok(), `${language}/${name}: ${response.status()}`);
          assert.equal(await page.locator('html').getAttribute('lang'), language.toLowerCase());
          assert.equal(await page.locator('.language-switch a[aria-current="true"]').innerText(), language);
          assert.equal(await page.locator('.language-switch a').count(), 2);
          assert(await page.locator('.language-switch').isVisible());
          if (language === 'EN') assertEnglish(await page.content(), name);
          assert.equal(await page.locator('img').evaluateAll(images => images.filter(img => !img.complete || img.naturalWidth === 0).length), 0, `Broken diagrams: ${name}`);
          if (mode === 'offline') {
            const links = await page.locator('a[href]').evaluateAll(anchors => anchors.map(a => a.href));
            for (const href of new Set(links)) {
              const url = new URL(href);
              if (url.protocol !== 'file:') continue;
              const target = fileURLToPath(url);
              assert(fs.existsSync(target), `Broken link in ${language}/${name}: ${href}`);
              if (url.hash && target.endsWith('.html')) {
                const html = fs.readFileSync(target, 'utf8');
                // Legacy redirect pages forward the fragment to the localized page.
                if (/http-equiv=["']refresh["']/i.test(html)) continue;
                const ids = [...html.matchAll(/\bid=["']([^"']+)["']/g)].map(match => match[1]);
                assert(ids.includes(decodeURIComponent(url.hash.slice(1))), `Broken bookmark in ${language}/${name}: ${href}`);
              }
            }
          }
          const other = language === 'RU' ? 'EN' : 'RU';
          await page.locator(`.language-switch a[lang="${other.toLowerCase()}"]`).focus();
          await page.keyboard.press('Enter');
          await page.waitForURL(`**/reference/${other}/${name}`);
          if (visualizers.includes(name)) assert.equal(await page.evaluate(() => state.stage), 0);
        }
      }
      for (const name of pages.filter(name => name !== 'index.html')) {
        const html = fs.readFileSync(path.join(reference, 'RU', name), 'utf8');
        const anchor = (html.match(/<h2 id="([^"]+)"/) || html.match(/\sid="([^"]+)"/))[1];
        await page.goto(base + `reference/${name}?bookmark=1#${anchor}`);
        await page.waitForURL(`**/reference/RU/${name}?bookmark=1#${anchor}`);
        assert.equal(await page.locator(`[id="${anchor}"]`).count(), 1);
      }
      for (const locale of ['', 'RU/', 'EN/']) {
        const target = locale || 'RU/';
        await page.goto(base + `reference/${locale}04c-nafems-challenge-5.html?bookmark=1#model-title`);
        await page.waitForURL(`**/reference/${target}06a-nafems-challenge-5.html?bookmark=1#model-title`);
        assert.equal(await page.locator('#model-title').count(), 1);
      }
      for (const name of visualizers) {
        const results = [];
        for (const language of ['RU', 'EN']) {
          await page.goto(base + `reference/${language}/${name}`);
          results.push(await exercise(page, name.slice(0, 2), language === 'EN'));
        }
        assert.deepEqual(results[1], results[0], `Numerical parity: ${name}`);
      }
      for (const width of [1440, 390]) {
        await page.setViewportSize({ width, height: 1000 });
        for (const name of ['index.html', 'EN/index.html', ...['capabilities.html', 'project-operations.html', '08-result-recovery.html'].flatMap(name => ['RU/' + name, 'EN/' + name]), 'RU/3d-quick-start.html', 'EN/3d-quick-start.html', 'RU/10-input-format.html', 'EN/10-input-format.html', 'EN/02-element-matrices.html', 'RU/04-boundary-loads.html', 'EN/04-boundary-loads.html', 'RU/04a-mpc.html', 'EN/04a-mpc.html', 'RU/04b-end-releases.html', 'EN/04b-end-releases.html', ...visualizers.map(name => 'EN/' + name)]) {
          await page.goto(base + 'reference/' + name);
          const dimensions = await page.evaluate(() => ({ scroll: document.documentElement.scrollWidth, viewport: innerWidth }));
          assert(dimensions.scroll <= dimensions.viewport + 1, `${mode}/${name} overflows at ${width}: ${JSON.stringify(dimensions)}`);
          await page.screenshot({ path: path.join(output, `${mode}-${width}-${name.replaceAll('/', '-')}.png`) });
        }
      }
      await page.setViewportSize({ width: 1440, height: 1000 });
      for (const app of ['preprocessor', 'postprocessor']) {
        await page.goto(base + app + '/index.html');
        await page.locator('.app-nav a[href="../reference/index.html"]').click();
        await page.waitForURL('**/reference/index.html');
        assert.equal(await page.locator('.language-choice').count(), 2);
      }
      await page.goto(base + 'reference/RU/3d-quick-start.html');
      let downloadedCase = path.join(root,'examples/cases/CaseSpatialQuickStart.txt');
      if(mode==='offline') {
        // Chromium opens file:// TXT links directly; download applies over HTTP.
        await page.locator('a[download]').click();
        await page.waitForURL('**/examples/cases/CaseSpatialQuickStart.txt');
        assert.match(await page.locator('body').innerText(),/elems_1436/);
      } else {
        const downloadEvent = page.waitForEvent('download');
        await page.locator('a[download]').click();
        const download = await downloadEvent;
        downloadedCase = path.join(output, 'CaseSpatialQuickStart.txt');
        await download.saveAs(downloadedCase);
      }
      assert.equal(fs.readFileSync(downloadedCase,'utf8'),fs.readFileSync(path.join(root,'examples/cases/CaseSpatialQuickStart.txt'),'utf8'));
      await page.goto(base + 'preprocessor/index.html');
      await page.locator('#file-input').setInputFiles(downloadedCase);
      await page.waitForSelector('#spatial-node-rows tr');
      assert.equal(await page.locator('#spatial-node-rows tr').count(),2);
      assert.equal(await page.locator('#spatial-element-rows tr').count(),1);
      await page.goto(base + 'postprocessor/index.html');
      await page.locator('#file-input').setInputFiles(path.join(root,'output/CaseSpatialQuickStart.json'));
      await page.locator('#spatial-result-rows tr').nth(1).click();
      assert.match(await page.locator('#spatial-result-details').innerText(),/uy = -0.0002/);
      await page.locator('#spatial-result-table').selectOption('element');
      await page.locator('#spatial-result-rows tr').click();
      await page.locator('#spatial-result-position').fill('0');
      await page.locator('#spatial-result-position').dispatchEvent('change');
      assert.match(await page.locator('#section-result-info').innerText(),/σmin=-3000000/);
      assert.match(await page.locator('#section-result-info').innerText(),/σmax=3000000/);
      console.log(`PASS ${mode}: ${pages.length * 2} pages, ${pages.length - 1} legacy redirects, keyboard switches, visualizer numerical parity, responsive layouts and downloadable 3D example`);
    }
    const noScript = await browser.newContext({ javaScriptEnabled: false });
    const fallback = await noScript.newPage();
    await fallback.goto(bases[0] + 'reference/01-model-dofs.html');
    await fallback.locator('#redirect-target').click();
    await fallback.waitForURL('**/RU/01-model-dofs.html');
    await fallback.locator('.language-switch a[lang="en"]').click();
    await fallback.waitForURL('**/EN/01-model-dofs.html');
    await noScript.close();
    assert.deepEqual(errors, [], 'Browser/CSP errors');
    console.log('PASS no-JavaScript fallback and language switching; no browser/CSP errors');
  } finally {
    await browser.close();
  }
})().catch(error => { console.error(error); process.exitCode = 1; });
