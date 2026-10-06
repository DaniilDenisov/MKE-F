// Requires local Playwright (or MKEF_PLAYWRIGHT_MODULE). Always checks file://;
// set MKEF_REFERENCE_BASE_URL to additionally check the Docker/Nginx host.
const { chromium } = require(process.env.MKEF_PLAYWRIGHT_MODULE || 'playwright');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const { pathToFileURL } = require('node:url');
const root = path.resolve(__dirname, '..');
const reference = path.join(root, 'reference');
const pages = fs.readdirSync(path.join(reference, 'RU')).filter(name => name.endsWith('.html'));
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
        for (const name of ['index.html', 'EN/index.html', 'EN/02-element-matrices.html', 'RU/04-boundary-loads.html', 'EN/04-boundary-loads.html', 'RU/04a-mpc.html', 'EN/04a-mpc.html', ...visualizers.map(name => 'EN/' + name)]) {
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
      console.log(`PASS ${mode}: ${pages.length * 2} pages, ${pages.length - 1} legacy redirects, keyboard switches, visualizer numerical parity and responsive layouts`);
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
