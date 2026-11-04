const {chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE||'playwright');
const assert=require('node:assert/strict'),path=require('node:path');
(async()=>{
 const browser=await chromium.launch({headless:true,args:['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
 try{
  const page=await browser.newPage({viewport:{width:1920,height:1080}}),errors=[];page.on('pageerror',e=>errors.push(e.message));page.on('dialog',d=>d.accept());
  await page.goto('file://'+path.resolve('preprocessor/index.html'));
  await page.evaluate(()=>{const render=MKEFPre.SpatialPanel.prototype.render;MKEFPre.SpatialPanel.prototype.render=function(...a){window.panel=this;return render.apply(this,a);};});
  await page.locator('#element-type').selectOption('1436');
  assert.doesNotMatch(await page.locator('body').innerText(),/[\u0400-\u04ff]/);
  await page.locator('[data-spatial-tool="node"]').click();await page.locator('#spatial-node-apply').click();await page.locator('#spatial-node-x').fill('2');await page.locator('#spatial-node-apply').click();
  await page.locator('[data-spatial-tool="member"]').click();
  assert.equal(Number(await page.locator('#spatial-member-youngsModulus').inputValue()),210e9);assert.equal(Number(await page.locator('#spatial-member-density').inputValue()),7850);assert.equal(Number(await page.locator('#spatial-member-shearModulus').inputValue()),210e9/2.6);
  const original=await page.locator('#spatial-current').inputValue(),before=await page.evaluate(()=>JSON.stringify(panel.scene.getViewState()));
  const canvas=await page.locator('.canvas-wrap').boundingBox(),sidebar=await page.locator('#spatial-tables').boundingBox();assert(sidebar.x>=canvas.x+canvas.width);assert(canvas.height>500);assert(sidebar.y<canvas.y);assert.equal(await page.locator('#spatial-tables #spatial-catalog').count(),0);
  await page.locator('#spatial-catalog-toggle').click();await page.locator('#spatial-catalog-tab-material').click();
  assert.doesNotMatch(await page.locator('#spatial-catalog-tool').innerText(),/[\u0400-\u04ff]/);
  await page.locator('#spatial-material-youngsModulus').fill('0');await page.locator('#spatial-material-apply').click();assert.match(await page.locator('#spatial-catalog-message').innerText(),/positive/);
  await page.locator('#spatial-material-youngsModulus').fill('200000000000');await page.locator('#spatial-material-apply').click();
  assert.equal(await page.locator('#spatial-current').inputValue(),original);assert.equal(Number(await page.locator('#spatial-member-youngsModulus').inputValue()),200e9);
  await page.locator('#spatial-material-steel').click();await page.locator('#spatial-material-apply').click();
  await page.locator('#spatial-catalog-tab-sections').click();await page.locator('#spatial-section-new').click();await page.locator('#spatial-section-name').fill('Стальная труба');await page.locator('#spatial-section-type').selectOption('square_tube');
  await page.screenshot({path:path.resolve('output/spatial-catalog-overlay.png')});
  await page.keyboard.press('Escape');assert.equal(await page.locator('#spatial-catalog-tool').evaluate(e=>e.open),false);assert.equal(await page.evaluate(()=>panel.tool),'member');assert.equal(await page.locator('#spatial-current').inputValue(),original);
  await page.locator('#spatial-catalog-toggle').click();assert.equal(await page.locator('#spatial-section-name').inputValue(),'Стальная труба');await page.locator('#spatial-section-apply').click();
  await page.locator('#spatial-catalog-undo').click();assert.equal(await page.evaluate(()=>panel.model.sections.length),0);await page.locator('#spatial-catalog-redo').click();assert.equal(await page.evaluate(()=>panel.model.sections.length),1);
  await page.locator('#spatial-catalog-close').click();await page.locator('#spatial-member-sectionId').selectOption('1');
  await page.locator('#spatial-member-node1').fill('1');await page.locator('#spatial-member-node2').fill('2');await page.locator('#spatial-member-apply').click();
  assert.equal(await page.evaluate(()=>panel.model.elements[0].youngsModulus),210e9);
  await page.locator('#spatial-catalog-toggle').click();await page.locator('#spatial-catalog-tab-material').click();await page.locator('#spatial-material-youngsModulus').fill('190000000000');await page.locator('#spatial-material-apply').click();await page.locator('#spatial-catalog-close').click();
  assert.equal(await page.evaluate(()=>panel.model.elements[0].youngsModulus),210e9);assert.equal(Number(await page.locator('#spatial-member-youngsModulus').inputValue()),190e9);
  // No view reset on opening, editing settings, closing or creating a member.
  const after=await page.evaluate(()=>JSON.stringify(panel.scene.getViewState()));const a=JSON.parse(before),b=JSON.parse(after);assert.deepEqual(a.target,b.target);assert(Math.abs(a.height-b.height)<1e-12);
  await page.screenshot({path:path.resolve('output/spatial-workspace.png')});
  assert.equal(await page.locator('.screen-size-notice').isVisible(),false);
  await page.setViewportSize({width:1280,height:900});await page.locator('#spatial-catalog-toggle').click();
  const dialog=await page.locator('#spatial-catalog-tool').boundingBox();assert(dialog.x>=0&&dialog.x+dialog.width<=1281);assert(await page.locator('#spatial-catalog-close').isVisible());await page.keyboard.press('Escape');
  await page.setViewportSize({width:390,height:844});assert(await page.locator('.screen-size-notice').isVisible());assert.match(await page.locator('.screen-size-notice').innerText(),/MCP/);
  await page.goto('file://'+path.resolve('postprocessor/index.html'));assert(await page.locator('.screen-size-notice').isVisible());
  assert.deepEqual(errors,[]);console.log('Spatial desktop/tablet workspace, material defaults, undo, camera and narrow-screen notice PASS');
 }finally{await browser.close();}
})().catch(e=>{console.error(e);process.exit(1);});
