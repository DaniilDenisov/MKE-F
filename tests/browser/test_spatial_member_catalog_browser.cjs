const {captureEditor,pickNode}=require('../helpers/browser-test-helpers.cjs');
const {chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE||'playwright');
const assert=require('node:assert/strict'),path=require('node:path');
(async()=>{
 const browser=await chromium.launch({headless:true,args:['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
 try{
  for(const noGL of [false,true]){
   const page=await browser.newPage({viewport:{width:1440,height:1050}}),errors=[];
   page.on('pageerror',e=>errors.push(e.message));page.on('dialog',d=>d.accept());
   if(noGL)await page.addInitScript(()=>{const get=HTMLCanvasElement.prototype.getContext;HTMLCanvasElement.prototype.getContext=function(kind,...args){return kind==='webgl2'?null:get.call(this,kind,...args);};});
   await page.goto('file://'+path.resolve('preprocessor/index.html'));
   await page.evaluate(()=>{const render=MKEFPre.SpatialPanel.prototype.render;MKEFPre.SpatialPanel.prototype.render=function(...args){window.panel=this;return render.apply(this,args);};});
   await page.locator('#element-type').selectOption('1436');
   await page.locator('#analysis-type').selectOption('static');
   await page.locator('[data-spatial-tool="node"]').click();
   for(const x of [0,2,4]){await page.locator('#spatial-node-x').fill(String(x));await page.locator('#spatial-node-apply').click();}
   assert.equal(await page.locator('#spatial-catalog-tool').evaluate(e=>e.open),false);
   await page.locator('[data-spatial-tool="member"]').click();
   await page.locator('#spatial-member-youngsModulus').fill('2e11');
   const before=await page.locator('#spatial-current').inputValue(),history=await page.evaluate(()=>panel.history.undo.length);
   await page.locator('#spatial-member-catalog').click();await page.locator('#spatial-section-new').click();
   await page.locator('#spatial-section-name').fill('Прямоугольник');
   // Closing the tool preserves an unfinished profile without touching the model.
   await page.locator('#spatial-catalog-close').click();assert.equal(await page.locator('#spatial-current').inputValue(),before);
   assert.equal(await page.evaluate(()=>panel.history.undo.length),history);
   await page.locator('#spatial-member-catalog').click();assert.equal(await page.locator('#spatial-section-name').inputValue(),'Прямоугольник');
   await page.locator('#spatial-section-apply').click();assert.equal(await page.locator('#spatial-member-youngsModulus').inputValue(),'2e11');
   await page.locator('#spatial-catalog-close').click();
   await page.locator('#spatial-member-sectionId').selectOption('1');
   for(const key of ['area','Iy','Iz','torsionConstant']){assert(await page.locator('#spatial-member-'+key).evaluate(e=>e.readOnly));assert(Number(await page.locator('#spatial-member-'+key).inputValue())>0);}
   await page.locator('#spatial-member-sectionId').selectOption('0');assert.equal(await page.locator('#spatial-member-area').evaluate(e=>e.readOnly),false);
   await page.locator('#spatial-member-sectionId').selectOption('1');
   if(!noGL)await page.locator('[data-view="origin"]').click();
   const camera=()=>page.evaluate(()=>{if(!panel.scene)return null;const v=panel.scene.getViewState();return [...v.target,...v.direction,...v.up,v.height];});
   const view=await camera();
   await page.locator('#spatial-member-density').fill('');await page.locator('#spatial-member-shearModulus').fill('');
   let firstNoGL=true;const pick=async index=>{
    if(noGL){await page.locator('#spatial-member-node'+(firstNoGL?1:2)).fill(String(index+1));if(!firstNoGL)await page.locator('#spatial-member-apply').click();firstNoGL=!firstNoGL;return;}
    await page.evaluate(()=>new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve))));
    const point=await page.evaluate(i=>{const s=panel.scene,p=s.local(panel.model.nodes[i]).project(s.camera),r=s.canvas.getBoundingClientRect();return{x:r.left+(p.x+1)*r.width/2,y:r.top+(1-p.y)*r.height/2};},index);
    await page.mouse.click(point.x,point.y);
   };
   await pick(0);await pick(1);
   const error=await page.locator('#spatial-edit-status').innerText();if(noGL)assert(await page.locator('#spatial-member-density').evaluate(e=>!e.validity.valid));else assert.match(error,/positive and finite/);
   assert.equal(await page.locator('#spatial-element-rows tr').count(),0);
   if(!noGL){const b=await page.locator('#spatial-viewport').boundingBox();await page.mouse.move(b.x+20,b.y+20,{steps:3});}
   if(!noGL)assert.equal(await page.locator('#spatial-edit-status').innerText(),error);
   await page.locator('#spatial-member-density').fill('7850');await page.locator('#spatial-member-shearModulus').fill('8e10');
   const count=await page.evaluate(()=>panel.history.undo.length);
   await page.locator('#spatial-member-apply').click();assert.equal(await page.locator('#spatial-element-rows tr').count(),1);
   assert.equal(await page.evaluate(()=>panel.history.undo.length),count+1);
   assert.equal(await page.locator('#spatial-edit-status').innerText(),'');
   assert.equal(await page.locator('#spatial-member-sectionId').inputValue(),'1');
   assert.equal(await page.locator('#spatial-catalog-tool').evaluate(e=>e.open),false);
   const after=await camera();if(view)view.forEach((v,i)=>assert(Math.abs(v-after[i])<1e-12));
   const saved=await page.locator('#spatial-current').inputValue();assert.match(saved,/element_sections/);
   await page.locator('#undo').click();assert.equal(await page.locator('#spatial-element-rows tr').count(),0);
   await page.locator('#redo').click();assert.equal(await page.locator('#spatial-current').inputValue(),saved);
   await pick(1);await pick(2);assert.equal(await page.locator('#spatial-element-rows tr').count(),2,await page.locator('#spatial-edit-status').innerText());
   assert.deepEqual(await page.evaluate(()=>panel.model.elements.map(e=>e.sectionId)),[1,1]);
   const parsed=await page.evaluate(()=>{const m=MKEFPre.spatialFormat.parse(panel.history.text());return m.elements.map(e=>({sectionId:e.sectionId,area:e.area}));});assert.equal(parsed[1].sectionId,1);assert.equal(parsed[0].area,parsed[1].area);
   assert.deepEqual(errors,[]);await page.screenshot({path:path.resolve('output/member-catalog'+(noGL?'-no-webgl':'')+'.png')});
   if(noGL){
    // Import a draft explicitly: exporting an element-free model is prohibited.
    const text='analysis\nstatic\nnodes\n2\n0,0,0\n1,2,3\nelems_1133\n0\nsections\n205\n'+Array.from({length:205},(_,i)=>JSON.stringify({id:i+1,name:'Profile '+(i+1),type:'rectangle',dimensions:{b:.1,h:.2}})).join('\n')+'\n';
    await page.locator('#file-input').setInputFiles({name:'catalog-draft.txt',mimeType:'text/plain',buffer:Buffer.from(text)});
    await page.locator('[data-spatial-tool="member"]').click();
    assert.equal(await page.locator('#spatial-member-sectionId').inputValue(),'0');
    assert.equal(await page.locator('#spatial-member-sectionId option').count(),201);
    await page.locator('#spatial-member-section-search').fill('Profile 205');
    await page.locator('#spatial-member-sectionId').selectOption('205');
    for(const [key,value]of Object.entries({node1:'1',node2:'2',youngsModulus:'2e11',density:'7850'}))await page.locator('#spatial-member-'+key).fill(value);
    // A removed draft choice must not silently become a numeric-only member.
    await page.evaluate(()=>panel.execute(h=>h.section(204,null)));
    await page.locator('#spatial-member-apply').click();
    assert.match(await page.locator('#spatial-edit-status').innerText(),/Unknown section 205/);
    assert.equal(await page.locator('#spatial-element-rows tr').count(),0);
    await page.locator('#undo').click();await page.locator('#spatial-member-apply').click();
    assert.equal(await page.locator('#spatial-element-rows tr').count(),1);
    assert.equal(await page.evaluate(()=>panel.model.elements[0].sectionId),205);
    assert.deepEqual(errors,[]);console.log('Truss creation, bounded catalog search, new-model reset and removed draft section PASS');
   }
   await page.close();console.log('Catalog member creation, persistent errors, collapsed tool, camera and undo/redo '+(noGL?'without WebGL':'with scene clicks')+' PASS');
  }
 }finally{await browser.close();}
})().catch(e=>{console.error(e);process.exit(1);});
