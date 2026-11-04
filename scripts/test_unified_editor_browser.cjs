const {chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE||'playwright');
const assert=require('node:assert/strict'),path=require('node:path');
(async()=>{const browser=await chromium.launch({headless:true,args:['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});try{
 const page=await browser.newPage({viewport:{width:1600,height:1000}}),errors=[];page.on('pageerror',e=>{errors.push(e.message);console.error(e.stack);});page.on('dialog',d=>d.accept());
 for(const code of [112,113,1133,1436]){
  await page.goto('file://'+path.resolve('preprocessor/index.html'));
  await page.evaluate(()=>{const render=MKEFPre.SpatialPanel.prototype.render;MKEFPre.SpatialPanel.prototype.render=function(...args){window.panel=this;return render.apply(this,args);};});
  await page.locator('#element-type').selectOption(String(code));
  assert.equal(await page.locator('.settings section:first-child #spatial-catalog-toggle').count(),1);
  await page.locator('[data-spatial-tool=node]').click();await page.locator('#spatial-node-apply').click();
  await page.locator('#spatial-node-x').fill('2');await page.locator('#spatial-node-apply').click();
  assert.equal(await page.locator('#spatial-node-rows tr').count(),2);
  await page.locator('#spatial-catalog-toggle').click();await page.locator('#spatial-section-new').click();await page.locator('#spatial-section-apply').click();await page.locator('#spatial-catalog-close').click();
  await page.locator('[data-spatial-tool=member]').click();await page.locator('#spatial-member-sectionId').selectOption('1');
  await page.locator('#spatial-member-node1').fill('1');await page.locator('#spatial-member-node2').fill('2');await page.locator('#spatial-member-apply').click();
  assert.equal(await page.locator('#spatial-element-rows tr').count(),1,await page.locator('#spatial-edit-status').innerText());
  await page.locator('#analysis-type').selectOption('static');
  const result=await page.evaluate(()=>{const a=MKEFPre.editorAdapter,text=a.serialize(panel.model),m=a.parse(text);return {text,dimension:m.dimension,inertia:m.elements[0].momentOfInertia,errors:a.validate(m,true)};});
  assert.deepEqual(result.errors,[]);assert.equal(result.dimension,code<1000?2:3);if(code===113)assert(result.inertia>0);assert.match(result.text,new RegExp('elems_'+code));
  assert.equal(await page.locator('#spatial-analysis-apply').count(),0);
  await page.locator('#analysis-type').selectOption('transient');assert.equal(await page.evaluate(()=>panel.model.analysis.type),'transient');
  await page.locator('#spatial-analysis-timeStep').fill('-1');await page.locator('#spatial-analysis-timeStep').dispatchEvent('change');assert(await page.locator('#run-analysis').isDisabled());assert((await page.locator('#analysis-error').innerText()).length>0);
  await page.locator('#spatial-analysis-timeStep').fill('0.001');await page.locator('#spatial-analysis-timeStep').dispatchEvent('change');assert.equal(await page.evaluate(()=>panel.model.analysis.timeStep),.001);
  await page.locator('#analysis-type').selectOption('static');
  const oldName=await page.locator('#case-name').inputValue();await page.locator('#case-name').fill('Renamed');await page.locator('#case-name').dispatchEvent('change');await page.locator('#undo').click();assert.equal(await page.locator('#case-name').inputValue(),oldName);await page.locator('#redo').click();assert.equal(await page.locator('#case-name').inputValue(),'Renamed');
  const before=await page.evaluate(()=>JSON.stringify(panel.model.elements));await page.locator('#undo').click();await page.locator('#redo').click();assert.equal(await page.evaluate(()=>JSON.stringify(panel.model.elements)),before);
  const canvas=await page.locator('.canvas-wrap').boundingBox(),sidebar=await page.locator('#spatial-tables').boundingBox();assert(sidebar.x>=canvas.x+canvas.width-1);assert(canvas.height>350);
  const modelTools=await page.locator('.model-tools').boundingBox(),viewTools=await page.locator('.view-tools').boundingBox();assert(viewTools.x>modelTools.x);assert(modelTools.width>viewTools.width);
  assert.equal(await page.locator('[data-spatial-tool=move],#spatial-delete,#spatial-notice').count(),0);
  await page.locator('#spatial-node-rows tr').first().click();await page.locator('#spatial-node-x').focus();await page.keyboard.press('Delete');assert.equal(await page.locator('dialog[open]').count(),0);
  await page.locator('#spatial-node-rows tr').first().click();await page.keyboard.press('Delete');await page.locator('#spatial-delete-confirm').click();
  assert.equal(await page.locator('#spatial-element-rows tr').count(),0);assert.equal(await page.locator('#spatial-node-rows tr').count(),1);
  await page.locator('#undo').click();assert.equal(await page.locator('#spatial-element-rows tr').count(),1);assert.equal(await page.locator('#spatial-node-rows tr').count(),2);
  await page.screenshot({path:path.resolve('output/unified-editor-'+code+'.png')});
 }
 await page.locator('#file-input').setInputFiles(path.resolve('examples/cases/CaseSectionBeam.txt'));await page.waitForSelector('#spatial-element-rows tr');await page.locator('#spatial-element-rows tr').click();await page.locator('#assign-section').click();await page.locator('#spatial-section-id').fill('1');await page.locator('#spatial-section-angle').fill('90');await page.locator('#spatial-section-sequence').click();await page.locator('#spatial-element-rows tr').click();assert.equal(await page.evaluate(()=>panel.model.elements[0].sectionAngle),90);assert.equal(await page.locator('#condition-elementLoads-rows tr').count(),1);await page.keyboard.press('Escape');await page.locator('#undo').click();await page.locator('#spatial-element-rows tr').click();assert.equal(await page.evaluate(()=>panel.model.elements[0].sectionAngle),0);
 await page.locator('#file-input').setInputFiles(path.resolve('examples/cases/CasePreprocessorStatic.txt'));await page.waitForFunction(()=>panel.model.name==='CasePreprocessorStatic');const beforeLaw=await page.locator('#spatial-current').inputValue();await page.locator('#analysis-type').selectOption('transient');assert.equal(await page.locator('#analysis-type').inputValue(),'static');assert.equal(await page.locator('#spatial-current').inputValue(),beforeLaw);assert.match(await page.locator('#spatial-edit-status').innerText(),/not valid for transient/);
 assert.deepEqual(errors,[]);console.log('Shared editor: all four families, catalog, assignment, native TXT, history, deletion and analysis validation PASS');
}finally{await browser.close();}})().catch(e=>{console.error(e);process.exitCode=1;});
