const {chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE||'playwright'),assert=require('node:assert/strict'),path=require('node:path');
(async()=>{const browser=await chromium.launch({headless:true,args:['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});try{
const page=await browser.newPage({viewport:{width:1600,height:1000}});await page.goto('file://'+path.resolve('preprocessor/index.html'));
await page.evaluate(()=>{const render=MKEFPre.SpatialPanel.prototype.render;MKEFPre.SpatialPanel.prototype.render=function(...args){window.panel=this;return render.apply(this,args);};});await page.locator('#element-type').selectOption('1436');
await page.evaluate(()=>{panel.scene.setWorkPlane({plane:'XZ',offset:7,step:.5});});
const before=await page.evaluate(()=>panel.scene.getViewState());
for(const axis of ['X','Y','Z']){await page.locator('.orientation-triad [data-axis='+axis+']').click();const view=await page.evaluate(()=>panel.scene.getViewState());assert.deepEqual(view.workPlane,before.workPlane);assert.deepEqual(view.target,before.target);assert.equal(view.height,before.height);assert(Math.abs(view.direction[{X:0,Y:1,Z:2}[axis]])>.999);}
assert.equal(await page.locator('[data-view=XY],[data-view=XZ],[data-view=YZ]').count(),0);
await page.evaluate(()=>{panel.scene.dispose();window.planar=new MKEFSpatial.SpatialScene(document.getElementById('spatial-viewport'),{}, {dimension:2});planar.setModel({nodes:[{x:0,y:0,z:0},{x:2,y:0,z:0}],elements:[{node1:1,node2:2}]});});
const planar=await page.evaluate(()=>{const old=planar.getViewState();planar.setView('iso');let rejected=false;try{planar.setWorkPlane({plane:'XZ'});}catch(e){rejected=true;}return {old,next:planar.getViewState(),rotate:planar.controls.enableRotate,size:planar.nodes.material.size,rejected};});assert.deepEqual(planar.old,planar.next);assert.equal(planar.rotate,false);assert.equal(planar.size,11);assert(planar.rejected);
assert.equal(await page.locator('.orientation-triad button:disabled').count(),3);
const orientation=await page.evaluate(()=>{planar.capture();return ['X','Y','Z'].map(axis=>{const button=planar.triad.querySelector('[data-axis='+axis+']');return {labels:planar.labelFrame.filter(label=>label.text===axis).length,text:button.textContent,background:getComputedStyle(button).backgroundColor,name:button.getAttribute('aria-label')};});});
for(const axis of orientation){assert.equal(axis.labels,1,'One visible/exportable axis label');assert.equal(axis.text,'','No duplicate HTML lettering');assert.equal(axis.background,'rgba(0, 0, 0, 0)');assert(axis.name.startsWith('View '));}
console.log('Shared navigation: triad preserves work plane/target/scale, planar camera locked, 11px nodes PASS');
}finally{await browser.close();}})().catch(e=>{console.error(e);process.exitCode=1;});
