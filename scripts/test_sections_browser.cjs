const {chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE || 'playwright');
const fs=require('fs'),path=require('path'),http=require('http'),assert=require('assert/strict');
const root=path.resolve(__dirname,'..');
async function checkSectionMode(page) {
  const preview=()=>page.locator('#case-preview').inputValue();
  const tool=page.locator('#assign-section');
  assert.equal(await tool.evaluate(e=>e.previousElementSibling.id),'add-element-load');
  await tool.click();
  assert(await page.locator('#section-tool-settings').isVisible());
  assert(await page.locator('#section-tool-id').isDisabled());
  assert(await page.locator('#section-tool-empty').isVisible());
  await page.locator('#section-tool-catalog').click();
  assert(await page.locator('#sections-editor').evaluate(e=>e.open));
  await page.locator('#add-section').click();
  assert(!(await page.locator('#section-tool-id').isDisabled()));
  await page.locator('#undo').click();
  assert(await page.locator('#section-tool-id').isDisabled());
  await page.keyboard.press('Escape');
  assert.equal(await page.locator('[data-tool="select"]').getAttribute('aria-pressed'),'true');
  // Two unassigned members, with an existing distributed load.
  const input=fs.readFileSync(path.join(root,'examples/cases/CaseSectionBeam.txt'),'utf8').replace(/\r/g,'')
    .replace('nodes\n2\n0,0,0\n2,0,0','nodes\n3\n0,0,0\n2,0,0\n4,0,0')
    .replace('elems_113\n1\n','elems_113\n2\n113,2,3,0.01,200000000000,7800,0.00001\n')
    .replace('element_sections\n1\n1,1,0\n','');
  await page.locator('#file-input').setInputFiles({name:'mode.txt',mimeType:'text/plain',buffer:Buffer.from(input)});
  await page.waitForFunction(()=>document.querySelector('#case-name').value==='mode');
  const initial=await preview();
  await tool.click();
  await page.locator('#section-tool-angle').selectOption('90');
  assert.equal(await preview(),initial,'Activating/configuring the tool must not assign anything');
  await page.locator('.model-node[data-node-id="1"]').dispatchEvent('click');
  assert.equal(await preview(),initial,'Nodes must not receive sections');
  async function clickMember(id) {
    const p=await page.locator('.model-element[data-element-id="'+id+'"]').evaluate(e=>{
      const r=e.getBoundingClientRect();return {x:r.left+r.width*.65,y:r.top+r.height*.65};
    });
    await page.mouse.click(p.x,p.y);
  }
  await clickMember(1);
  assert.match(await preview(),/element_sections\s+1\s+1,1,90/);
  const first=await preview();
  await clickMember(1);
  assert.equal(await preview(),first);
  await page.locator('#undo').click();assert.equal(await preview(),initial,'Repeated assignment must not add a duplicate undo entry');
  await page.locator('#redo').click();assert.equal(await preview(),first);
  await clickMember(2);
  assert.match(await preview(),/element_sections\s+2\s+1,1,90\s+2,1,90/);
  assert.equal(await tool.getAttribute('aria-pressed'),'true');
  await page.locator('#undo').click();assert.equal(await preview(),first);
  await page.locator('#redo').click();
  const assigned=await preview();
  await page.locator('#section-tool-angle').selectOption('180');assert.equal(await preview(),assigned);
  await clickMember(2);assert.match(await preview(),/2,1,180/);
  const model=await page.evaluate(()=>window.MKEFPre.caseFormat.parse(document.querySelector('#case-preview').value));
  assert.equal(model.elementLoads.length,1);assert.equal(model.elementLoads[0].qy,-1000);
  assert(Math.abs(model.elements[0].momentOfInertia-model.elements[1].momentOfInertia/4)<1e-12);
  await page.screenshot({path:path.join(root,'output/section-assignment-mode.png'),fullPage:true});
  await page.locator('.section-card').first().getByRole('button',{name:'Copy',exact:true}).click();
  const copied=await preview();
  await page.locator('#section-tool-id').selectOption('2');
  assert.equal(await preview(),copied,'Choosing another catalog section must wait for a member click');
  await clickMember(1);assert.match(await preview(),/1,2,180/);
  await page.locator('#undo').click();assert.equal(await preview(),copied);
  await page.locator('#undo').click();
  assert.equal(await page.locator('#section-tool-id').inputValue(),'1','Removed catalog choices must fall back to an existing section');
  await page.keyboard.press('Escape');assert(await page.locator('#section-tool-settings').isHidden());
  const finished=await preview();await clickMember(1);assert.equal(await preview(),finished,'Select mode must not apply the section tool');
  // Re-entering starts with the selected member's existing assignment.
  await tool.click();assert.equal(await page.locator('#section-tool-angle').inputValue(),'90');
  await page.keyboard.press('Escape');
  console.log('Section mode: empty catalog, repeated canvas assignment, orientations, unchanged loads, undo/redo, Escape PASS');
}
const server=http.createServer((req,res)=>{
  const file=path.resolve(root,'.'+decodeURIComponent(req.url.split('?')[0]));
  if(!file.startsWith(root+path.sep)){res.writeHead(403).end();return;}
  fs.readFile(file,(err,data)=>{if(err){res.writeHead(404).end();return;}res.setHeader('Content-Type',file.endsWith('.html')?'text/html':file.endsWith('.js')?'text/javascript':file.endsWith('.css')?'text/css':'application/octet-stream');res.end(data);});
});
(async()=>{
  await new Promise(resolve=>server.listen(0,'127.0.0.1',resolve));
  const base='http://127.0.0.1:'+server.address().port,browser=await chromium.launch({headless:true});
  const page=await browser.newPage({viewport:{width:1440,height:1000}}),errors=[];page.on('pageerror',e=>errors.push(e.message));
  page.on('dialog',dialog=>dialog.accept());
  try {
    await page.goto(base+'/preprocessor/index.html');
    await checkSectionMode(page);
    await page.locator('#file-input').setInputFiles(path.join(root,'examples/cases/CaseSectionBeam.txt'));
    await page.waitForFunction(()=>document.querySelector('#case-name').value==='CaseSectionBeam');
    await page.locator('#sections-editor').evaluate(e=>e.open=true);
    await page.locator('#elements-body').evaluate(e=>e.closest('details').open=true);
    assert.equal(await page.locator('.section-card').count(),1);
    assert.equal(await page.locator('.section-axes').count(),1);
    assert.equal(await page.locator('.section-axes text').count(),3);
    const row=page.locator('#elements-body tr').first();
    assert(await row.locator('input').nth(2).isDisabled());assert(await row.locator('input').nth(5).isDisabled());
    assert(await page.getByRole('button',{name:'Assign to selected element',exact:true}).isDisabled());
    await page.locator('.model-element[data-element-id="1"]').dispatchEvent('click');
    await row.locator('select').first().selectOption('0');
    await page.getByRole('button',{name:'Assign to selected element',exact:true}).click();
    assert.equal(await row.locator('select').first().inputValue(),'1');
    assert.match(await page.locator('.section-assignments').innerText(),/Assigned to: E1/);
    assert.equal(await page.locator('#element-loads-body tr').count(),1);
    await page.locator('#undo').click();assert.equal(await row.locator('select').first().inputValue(),'0');
    await page.locator('#redo').click();assert.equal(await row.locator('select').first().inputValue(),'1');
    await row.locator('select').last().selectOption('90');
    assert.equal(await page.locator('#element-loads-body tr').count(),1);
    await page.getByRole('textbox',{name:'Name · section 1',exact:true}).fill('Сечение, «Т»');await page.getByRole('textbox',{name:'Name · section 1',exact:true}).dispatchEvent('change');
    await page.getByLabel('Shape · section 1',{exact:true}).selectOption('tee');
    await page.getByLabel('h · section 1',{exact:true}).fill('0.3');await page.getByLabel('h · section 1',{exact:true}).dispatchEvent('change');
    assert.equal(await page.locator('#element-loads-body tr').count(),1);
    await page.locator('.section-card').getByRole('button',{name:'Delete',exact:true}).click();assert.match(await page.locator('#error-panel').innerText(),/assigned to elements 1/);
    await page.locator('.section-card').getByRole('button',{name:'Copy',exact:true}).click();assert.equal(await page.locator('.section-card').count(),2);
    await page.locator('#undo').click();assert.equal(await page.locator('.section-card').count(),1);
    await page.locator('#redo').click();assert.equal(await page.locator('.section-card').count(),2);
    const saved=await page.locator('#case-preview').inputValue();assert.match(saved,/Сечение/);
    fs.writeFileSync(path.join(root,'output/SectionEditorRoundtrip.txt'),saved);
    // Roundtrip preserves named catalog, orientation, and derived properties.
    await page.locator('#file-input').setInputFiles({name:'saved.txt',mimeType:'text/plain',buffer:Buffer.from(saved)});
    assert.equal(await page.locator('.section-card').count(),2);
    assert.equal(await row.locator('select').last().inputValue(),'90');
    const before=await row.locator('input').nth(2).inputValue();await row.locator('select').first().selectOption('0');
    assert(!(await row.locator('input').nth(2).isDisabled()));assert.equal(await row.locator('input').nth(2).inputValue(),before);
    await page.locator('#undo').click();assert(await row.locator('input').nth(2).isDisabled());
    await page.locator('.tables').evaluate(e=>e.scrollTop=0);await page.locator('#section-catalog').evaluate(e=>e.scrollTop=0);
    await page.screenshot({path:path.join(root,'output/sections-preprocessor.png'),fullPage:true});
    // Triads remain constant in screen units during zoom and point correctly on all directions.
    const axisSize=()=>page.locator('.section-axes text').first().evaluate(e=>e.getBoundingClientRect().height);
    const size=await axisSize();await page.locator('#viewport').hover();await page.mouse.wheel(0,-600);await page.waitForTimeout(150);assert(Math.abs(await axisSize()-size)<1);
    await page.locator('#show-section-axes').uncheck();assert.equal(await page.locator('.section-axes').count(),0);
    await page.setViewportSize({width:640,height:900});await page.screenshot({path:path.join(root,'output/sections-mobile.png'),fullPage:true});
    await page.setViewportSize({width:1440,height:1000});
    await page.goto(base+'/postprocessor/index.html');
    await page.locator('#file-input').setInputFiles(path.join(root,'output/CaseSectionBeam.json'));
    assert(await page.locator('#error-panel').isHidden());
    await page.locator('#static-result').selectOption('sigmaMax');
    assert.equal(await page.locator('.diagram-line').count(),1);
    await page.locator('.original-element[data-element-id="1"]').dispatchEvent('click');
    assert(await page.locator('#section-result-panel').isVisible());
    assert.equal(await page.locator('.diagram-fill').evaluate(e=>getComputedStyle(e).vectorEffect),'non-scaling-stroke');
    assert.match(await page.locator('#selection-details').innerText(),/7.5.*e\+5|750000/);
    assert.match(await page.locator('#section-result-info').innerText(),/7.5.*e\+5|750000/);
    await page.locator('#section-position-number').fill('0');await page.locator('#section-position-number').dispatchEvent('input');
    const endValues=(await page.locator('#section-result-info').innerText()).match(/σmin=([^ ]+) · σmax=([^\n]+)/);assert(Math.abs(Number(endValues[1]))<1e-6 && Math.abs(Number(endValues[2]))<1e-6);
    await page.locator('#section-position-number').fill('0.5');await page.locator('#section-position-number').dispatchEvent('input');
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
    unassigned.model.nodes[1].y=1;
    await page.locator('#file-input').setInputFiles({name:'unassigned.json',mimeType:'application/json',buffer:Buffer.from(JSON.stringify(unassigned))});
    await page.locator('#show-loads').uncheck();
    await page.waitForTimeout(100);
    const target=page.locator('.deformed-element[data-element-id="1"]');
    const clickPoint=await target.evaluate(e=>{const q=e.getPointAtLength(e.getTotalLength()*.5),p=new DOMPoint(q.x,q.y).matrixTransform(e.getScreenCTM());return {x:p.x,y:p.y};});
    await page.mouse.click(clickPoint.x,clickPoint.y);
    assert.match(await page.locator('#selection-details').innerText(),/Element 1/);
    assert.equal(await target.evaluate(e=>e===document.activeElement),true);
    assert.equal(await target.evaluate(e=>getComputedStyle(e).outlineStyle),'none');
    assert(await page.locator('#section-result-panel').isHidden());
    assert.equal(await page.locator('#section-result-svg').evaluate(e=>e.childNodes.length),0);
    assert(await page.locator('#export-section-svg').isDisabled());
    assert.match(await page.locator('#selection-details').innerText(),/assign a section and rerun/);
    await page.screenshot({path:path.join(root,'output/selection-focus-fixed.png')});
    await page.keyboard.press('Tab');await page.keyboard.press('Shift+Tab');
    assert.equal(await target.evaluate(e=>e.matches(':focus-visible')),true);
    assert.equal(await target.evaluate(e=>getComputedStyle(e).outlineStyle),'none');
    assert.equal(await target.evaluate(e=>getComputedStyle(e).vectorEffect),'non-scaling-stroke');
    assert.equal(await target.evaluate(e=>getComputedStyle(e).strokeWidth),'4px');
    for(const file of fs.readdirSync(path.join(root,'output')).filter(n=>/^SectionCase.*\.json$/.test(n))) {
      await page.locator('#file-input').setInputFiles(path.join(root,'output',file));assert(await page.locator('#error-panel').isHidden(),file+': '+await page.locator('#error-panel').innerText());
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
    assert.deepEqual(errors,[]);console.log('Sections browser: catalog, undo/redo, roundtrip, axes, stress map, v6 variants, SVG and existing suites PASS');
  } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;}).finally(()=>server.close());
