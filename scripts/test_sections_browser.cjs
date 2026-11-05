const {captureEditor,capturePost}=require('./browser-test-helpers.cjs');
const {chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE || 'playwright');
const fs=require('fs'),path=require('path'),http=require('http'),assert=require('assert/strict');
const root=path.resolve(__dirname,'..');
const server=http.createServer((req,res)=>{
  res.setHeader('Content-Security-Policy',"default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; img-src 'self' data: blob:; connect-src 'self'; object-src 'none'; base-uri 'none'; frame-ancestors 'none'");
  const file=path.resolve(root,'.'+decodeURIComponent(req.url.split('?')[0]));
  if(!file.startsWith(root+path.sep)){res.writeHead(403).end();return;}
  fs.readFile(file,(err,data)=>{if(err){res.writeHead(404).end();return;}res.setHeader('Content-Type',file.endsWith('.html')?'text/html':file.endsWith('.js')?'text/javascript':file.endsWith('.css')?'text/css':'application/octet-stream');res.end(data);});
});
(async()=>{
  await new Promise(resolve=>server.listen(0,'127.0.0.1',resolve));
  const base='http://127.0.0.1:'+server.address().port,browser=await chromium.launch({headless:true,args:['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
  const page=await browser.newPage({viewport:{width:1440,height:1000}}),errors=[];page.on('pageerror',e=>errors.push(e.message));
  page.on('dialog',dialog=>dialog.accept());
  try {
    await page.goto(base+'/preprocessor/index.html');
    await captureEditor(page);
    const importCase=async()=>{await page.locator('#file-input').setInputFiles(path.join(root,'examples/cases/CaseSectionBeam.txt'));await page.waitForSelector('#spatial-element-rows tr');await page.locator('#spatial-element-rows tr').first().click();};
    const openCatalog=async()=>{await page.locator('#spatial-catalog-toggle').click();await page.locator('#spatial-catalog-tab-sections').click();};
    const editSection=()=>page.locator('#spatial-catalog').getByRole('button',{name:'Edit',exact:true}).first().click();
    await importCase();const original=await page.evaluate(()=>({...editorPanel.model.elements[0]}));
    assert(await page.locator('#spatial-member-area').evaluate(e=>e.readOnly));await openCatalog();
    await page.locator('#spatial-section-id').fill('0');await page.locator('#spatial-section-assign').click();await page.locator('#spatial-catalog-close').click();
    assert.equal(Number(await page.locator('#spatial-member-area').inputValue()),original.area);assert(!(await page.locator('#spatial-member-area').evaluate(e=>e.readOnly)));
    await page.locator('#spatial-member-area').fill('0.03');await page.locator('#spatial-member-momentOfInertia').fill('0.00002');await page.locator('#spatial-member-apply').click();
    await openCatalog();await editSection();await page.locator('#spatial-section-h').fill('0.4');await page.locator('#spatial-section-apply').click();
    assert.equal(await page.evaluate(()=>editorPanel.model.elements[0].area),.03);await page.locator('#spatial-catalog').getByRole('button',{name:'Delete',exact:true}).click();await page.locator('#spatial-catalog-close').click();
    assert.equal(await page.evaluate(()=>editorPanel.model.sections.length),0);assert(!(await page.locator('#download-case').isDisabled()));
    for(const type of ['square_tube','rectangular_tube','circular_tube']){
      await importCase();await openCatalog();await editSection();await page.locator('#spatial-section-type').selectOption(type);await page.locator('#spatial-section-t').fill('0.01');await page.locator('#spatial-section-t').dispatchEvent('input');await page.locator('#spatial-section-apply').click();
      const expected=type==='square_tube'?.0036:type==='rectangular_tube'?.0056:Math.PI*.0009;
      assert(Math.abs(await page.evaluate(()=>editorPanel.model.elements[0].area)-expected)<1e-14);assert.equal(await page.evaluate(()=>editorPanel.model.elementLoads.length),1);
      await editSection();const shape=page.locator('#spatial-section-preview path');assert.equal(await shape.getAttribute('fill-rule'),'evenodd');assert.equal(await shape.evaluate(e=>e.isPointInFill(new DOMPoint(0,0))),false);
      const saved=await page.locator('#spatial-current').inputValue();await page.locator('#spatial-section-t').fill('0.05');await page.locator('#spatial-section-apply').click();assert.match(await page.locator('#spatial-catalog-message').innerText(),/thickness leaves no clear space/);assert.equal(await page.locator('#spatial-current').inputValue(),saved);
      await page.locator('#spatial-catalog-close').click();await page.locator('#undo').click();await page.locator('#redo').click();assert.equal(await page.locator('#spatial-current').inputValue(),saved);
      await page.locator('#file-input').setInputFiles({name:type+'.txt',mimeType:'text/plain',buffer:Buffer.from(saved)});await page.waitForFunction(type=>editorPanel.model.sections[0].type===type,type);
    }
    await importCase();await openCatalog();await editSection();const hostile='<img src=x onerror="window.sectionInjected=true">';await page.locator('#spatial-section-name').fill(hostile);await page.locator('#spatial-section-apply').click();assert.equal(await page.locator('#spatial-catalog img').count(),0);assert.equal(await page.evaluate(()=>window.sectionInjected),undefined);
    await page.locator('#spatial-catalog').getByRole('button',{name:'Delete',exact:true}).click();assert.match(await page.locator('#spatial-catalog-message').innerText(),/used by members/);
    await editSection();await page.locator('#spatial-section-name').fill('Сечение, «Т»');await page.locator('#spatial-section-type').selectOption('tee');await page.locator('#spatial-section-h').fill('0.3');await page.locator('#spatial-section-apply').click();
    await page.locator('#spatial-section-angle').fill('90');await page.locator('#spatial-section-assign').click();await page.locator('#spatial-catalog').getByRole('button',{name:'Copy',exact:true}).click();await page.locator('#spatial-section-apply').click();await page.locator('#spatial-catalog-close').click();
    const saved=await page.locator('#spatial-current').inputValue();assert.match(saved,/Сечение/);assert.equal(await page.evaluate(()=>editorPanel.model.sections.length),2);assert.equal(await page.evaluate(()=>editorPanel.model.elementLoads.length),1);fs.writeFileSync(path.join(root,'output/SectionEditorRoundtrip.txt'),saved);
    await page.locator('#file-input').setInputFiles({name:'saved.txt',mimeType:'text/plain',buffer:Buffer.from(saved)});await page.waitForFunction(()=>editorPanel.model.name==='saved');assert.equal(await page.evaluate(()=>editorPanel.model.elements[0].sectionAngle),90);assert.equal(await page.evaluate(()=>editorPanel.model.sections.length),2);
    await page.screenshot({path:path.join(root,'output/sections-preprocessor.png')});
    await page.goto(base+'/postprocessor/index.html');await capturePost(page);
    await page.locator('#file-input').setInputFiles(path.join(root,'output/CaseSectionBeam.json'));
    assert(await page.locator('#error-panel').isHidden());
    await page.locator('#static-result').selectOption('sigmaMax');
    assert(await page.evaluate(()=>{resultRenderer.flush();return resultRenderer.frame.diagrams.positions.length>0;}));
    await page.locator('#planar-result-objects select').selectOption('element');await page.locator('#planar-result-rows tr').first().click();
    assert(await page.locator('#section-result-panel').isVisible());
    assert(await page.evaluate(()=>resultRenderer.frame.diagramFill.positions.length>0));
    assert.match(await page.locator('#selection-details').innerText(),/7.5.*e\+5|750000/);
    assert.match(await page.locator('#section-result-info').innerText(),/7.5.*e\+5|750000/);
    await page.locator('#section-position-number').fill('0');await page.locator('#section-position-number').dispatchEvent('input');
    const endValues=(await page.locator('#section-result-info').innerText()).match(/σmin=([^ ]+) · σmax=([^\n]+)/);assert(Math.abs(Number(endValues[1]))<1e-6 && Math.abs(Number(endValues[2]))<1e-6);
    await page.locator('#section-position-number').fill('0.5');await page.locator('#section-position-number').dispatchEvent('input');
    for(const width of [1440,640]) {
      await page.setViewportSize({width,height:1000});
      await page.waitForFunction(()=>{
        const svg=document.querySelector('#section-result-svg'),bounds=svg.getBoundingClientRect();
        return Math.abs(Number(svg.getAttribute('width'))-bounds.width)<1;
      });
      const labels=await page.locator('#section-result-svg .section-stress-labels text').evaluateAll(nodes=>nodes.map(e=>{
        const box=e.getBoundingClientRect(),svg=e.closest('svg').getBoundingClientRect();
        return {size:parseFloat(getComputedStyle(e).fontSize),inside:box.left>=svg.left && box.right<=svg.right && box.top>=svg.top && box.bottom<=svg.bottom};
      }));
      assert.equal(labels.length,5);assert(labels.every(label=>label.size>=14 && label.inside),'Stress labels must be readable and fit the preview');
      await page.locator('#section-result-svg').screenshot({path:path.join(root,'output/section-readable-'+width+'.png')});
    }
    await page.setViewportSize({width:1440,height:1000});
    const download=page.waitForEvent('download');await page.locator('#export-section-svg').click();const svg=await download;await svg.saveAs(path.join(root,'output/section-export.svg'));assert(fs.readFileSync(path.join(root,'output/section-export.svg'),'utf8').includes('linearGradient'));
    const modelDownload=page.waitForEvent('download');await page.locator('#export-svg').click();await (await modelDownload).saveAs(path.join(root,'output/sections-model.svg'));
    const pngDownload=page.waitForEvent('download');await page.locator('#export-png').click();await (await pngDownload).saveAs(path.join(root,'output/sections-model.png'));
    assert(fs.statSync(path.join(root,'output/sections-model.png')).size>1000);
    assert((await page.locator('#viewport').boundingBox()).height>=240,'Section inspector collapsed the model canvas');
    await page.locator('.result-inspector').evaluate(e=>e.scrollTop=0);
    await page.screenshot({path:path.join(root,'output/sections-postprocessor.png'),fullPage:true});
    // Native mouse focus must not create a black SVG ring. dispatchEvent alone
    // does not focus the path and would miss this Chromium regression.
    const unassigned=JSON.parse(fs.readFileSync(path.join(root,'output/CaseSectionBeam.json'),'utf8'));
    unassigned.model.elements[0].sectionId=0;
    unassigned.analysis.elementResults[0].sectionStress=[];
    unassigned.metadata.title='Unassigned section';
    await page.locator('#file-input').setInputFiles({name:'unassigned.json',mimeType:'application/json',buffer:Buffer.from(JSON.stringify(unassigned))});
    await page.locator('#show-loads').uncheck();
    await page.waitForFunction(()=>document.querySelector('#dataset-title').textContent==='Unassigned section');
    assert(await page.locator('#error-panel').isHidden());
    await page.locator('#planar-result-objects select').selectOption('element');await page.locator('#planar-result-rows tr').first().click();
    assert.match(await page.locator('#selection-details').innerText(),/Element 1/);
    assert(await page.locator('#section-result-panel').isVisible());
    assert.match(await page.locator('#section-result-info').innerText(),/geometry is required/);
    assert.equal(await page.locator('#section-result-svg').evaluate(e=>e.childNodes.length),0);
    assert(await page.locator('#export-section-svg').isDisabled());
    assert.match(await page.locator('#selection-details').innerText(),/assign a section and rerun/);
    await page.screenshot({path:path.join(root,'output/selection-focus-fixed.png')});
    await page.locator('#planar-result-rows tr').first().press('Enter');assert.match(await page.locator('#selection-details').innerText(),/Element 1/);
    for(const file of fs.readdirSync(path.join(root,'output')).filter(n=>/^SectionCase.*\.json$/.test(n))) {
      const data=JSON.parse(fs.readFileSync(path.join(root,'output',file),'utf8'));data.metadata.title=file;
      await page.locator('#file-input').setInputFiles({name:file,mimeType:'application/json',buffer:Buffer.from(JSON.stringify(data))});
      await page.waitForFunction(title=>document.querySelector('#dataset-title').textContent===title,file);
      assert(await page.locator('#error-panel').isHidden(),file+': '+await page.locator('#error-panel').innerText());
      if(file.includes('_tube_')) {
        await page.locator('#planar-result-objects select').selectOption('element');await page.locator('#planar-result-rows tr').first().click();
        assert(await page.locator('#section-result-panel').isVisible());
        const shape=page.locator('#section-result-svg > path');
        assert.equal(await shape.getAttribute('fill-rule'),'evenodd');
        assert.equal(await shape.evaluate(e=>{const f=MKEFPost.sectionInspector.view.layout.frame;return e.isPointInFill(new DOMPoint(...f.project(0,0)));}),false);
        assert.match(await shape.getAttribute('fill'),/url\(#section-normal-gradient\)/);
        const download=page.waitForEvent('download');await page.locator('#export-section-svg').click();
        const output=path.join(root,'output',file.replace('.json','.svg'));await (await download).saveAs(output);
        assert.match(fs.readFileSync(output,'utf8'),/fill-rule="evenodd"/);
      }
    }
    // Horizontal, vertical, inclined and short members, one triad per orientation.
    const axes=await page.evaluate(()=>{
      const ns='http://www.w3.org/2000/svg',svg=document.createElementNS(ns,'svg');svg.setAttribute('viewBox','-1 -1 8 8');svg.style.cssText='position:fixed;top:0;left:0;width:640px;height:640px;background:white';document.body.appendChild(svg);
      const model={nodes:[{x:0,y:0},{x:2,y:0},{x:0,y:2},{x:0,y:4},{x:3,y:2},{x:5,y:4},{x:4,y:0},{x:4.1,y:0}],elements:[0,90,180,270].map((a,i)=>({node1:2*i+1,node2:2*i+2,sectionId:1,sectionAngle:a}))};
      window.MKEFSectionDrawing.axes(svg,model,true);
      const labels=[...svg.querySelectorAll('text')].map(e=>e.textContent),positions=[...svg.querySelectorAll('path')].map(e=>e.getAttribute('d'));
      svg.remove();return {labels,positions};
    });
    assert.deepEqual(axes.labels,['x','y','z⊙','x','z','y⊙','x','y','z⊗','x','z','y⊗']);assert(axes.positions.every(p=>!p.includes('NaN')));
    // Both browser unit suites and old dataset semantics remain supported.
    for(const app of ['preprocessor','postprocessor']) {
      await page.goto(base+'/'+app+'/tests/index.html');
      if(app==='preprocessor') {await page.waitForFunction(()=>document.querySelector('#test-status').dataset.status);assert.equal(await page.locator('#test-status').getAttribute('data-status'),'passed',await page.locator('#test-output').innerText());}
      else {await page.waitForFunction(()=>document.body.dataset.testStatus!=='running');assert.equal(await page.locator('body').getAttribute('data-test-status'),'passed',await page.locator('#test-output').innerText());}
    }
    assert.deepEqual(errors,[]);console.log('Sections browser: catalog, undo/redo, roundtrip, axes, stress map, v8 variants, SVG and existing suites PASS');
  } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;}).finally(()=>server.close());
