// Generate v4 fixtures with setup; addpath('tests'); test_releases in Octave.
const {chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE || 'playwright');
const {pathToFileURL}=require('url');
const path=require('path'), fs=require('fs'), assert=require('assert/strict');
const root=path.resolve(__dirname,'..');
(async()=>{
  const browser=await chromium.launch({headless:true});
  try {
    const page=await browser.newPage({viewport:{width:1440,height:1000}}), errors=[];
    page.on('pageerror',e=>errors.push(e.message));
    page.on('dialog',d=>d.accept());
    await page.goto(pathToFileURL(path.join(root,'preprocessor/index.html')).href);
    await page.locator('#file-input').setInputFiles(path.join(root,'examples/cases/CaseReleaseStatic.txt'));
    await page.waitForFunction(()=>document.querySelectorAll('.end-release').length===2);
    assert.equal(await page.locator('#error-panel').isVisible(),false);
    await page.locator('#elements-body').evaluate(e=>e.closest('details').open=true);
    const release=page.getByLabel('Element 1 release Mz1',{exact:true});
    await release.uncheck(); assert.equal(await page.locator('.end-release').count(),1);
    assert.match(await page.locator('#case-preview').inputValue(),/eload_uniform/);
    await page.locator('#undo').click(); assert.equal(await page.locator('.end-release').count(),2);
    await page.locator('#redo').click(); assert.equal(await page.locator('.end-release').count(),1);
    await release.check();
    const roundtrip=await page.evaluate(()=>{
      const F=window.MKEFPre.caseFormat, text=document.getElementById('case-preview').value;
      const model=F.parse(text); return F.serialize(model)===text;
    }); assert(roundtrip);
    const edit=await page.evaluate(()=>{
      const F=window.MKEFPre.caseFormat,E=window.MKEFPre.modelEdit;
      const m=F.parse(document.getElementById('case-preview').value);
      m.elements.push({...m.elements[0]}); m.releases.push({elementId:2,end:1,component:'Mz'});
      E.updateNode(m,0,'x',-.5); E.updateElement(m,0,'area',.02);
      const kept=m.releases.length===3; E.deleteElement(m,0);
      return kept && m.releases.length===1 && m.releases[0].elementId===1;
    }); assert(edit);
    // Ignore absent rotation restraints, retain the source input and display warning.
    let txt=fs.readFileSync(path.join(root,'examples/cases/CaseReleaseStatic.txt'),'utf8').replace('4,1,0,0,0','1,1,0,0,0');
    await page.locator('#file-input').setInputFiles({name:'ignored.txt',mimeType:'text/plain',buffer:Buffer.from(txt)});
    await page.waitForFunction(()=>document.getElementById('notice-panel').textContent.includes('thetaZ restraint ignored'));
    assert.equal(await page.locator('#error-panel').isVisible(),false);
    assert.match(await page.locator('#case-preview').inputValue(),/1,1,0,0,0/);
    await page.screenshot({path:path.join(root,'output/releases-preprocessor.png'),fullPage:true});
    // Invalid loads are rejected before submission.
    txt+='\nbcforce_stat\n1\n10,1,0,0,1\n';
    const rejected=await page.evaluate(t=>{try{window.MKEFPre.caseFormat.parse(t);return false;}catch(e){return /absent thetaZ/.test(e.message);}},txt);
    assert(rejected);
    await page.goto(pathToFileURL(path.join(root,'postprocessor/index.html')).href);
    for(const name of ['CaseReleaseStatic','CaseReleaseModal','CaseReleaseTransient','CaseReleasePartial','CaseReleaseMPCStatic','CaseReleaseMPCModal','CaseReleaseMPCTransient','CaseReleaseMPCPartial','CaseReleaseMPCSingle','CaseReleaseWarnings']) {
      const file=path.join(root,'output',name+'.json'), raw=JSON.parse(fs.readFileSync(file,'utf8'));
      raw.metadata.title=name;
      await page.locator('#file-input').setInputFiles({name:name+'.json',mimeType:'application/json',buffer:Buffer.from(JSON.stringify(raw))});
      await page.waitForFunction(title=>document.getElementById('dataset-title').textContent===title,raw.metadata.title);
      assert.equal(await page.locator('#error-panel').isVisible(),false,name+': '+await page.locator('#error-panel').textContent());
      assert.equal(await page.locator('[data-layer="original-geometry"] .end-release').count(),2);
      const checks=await page.evaluate(raw=>{
        const M=window.MKEFPost,d=M.validateDataset(raw);
        const bad=JSON.parse(JSON.stringify(raw)); bad.model.elements[0].globalDOFs[2]=bad.model.dofMap[0][0];
        let rejected=false;try{M.validateDataset(bad);}catch(e){rejected=e.name==='DataError';}
        if(raw.analysis.type==='modal') return rejected && !M.modalView.display(d,0,1,21).zero;
        if(raw.analysis.type==='static') return rejected && M.geometry.automaticScale(d,raw.analysis.displacements,21).scale>0;
        return rejected;
      },raw);assert(checks,name);
      await page.locator('.original-element[data-element-id="1"]').focus();
      await page.locator('.original-element[data-element-id="1"]').press('Enter');
      assert.match(await page.locator('#selection-details').textContent(),/End 1 thetaZ/);
      if(name==='CaseReleaseStatic') {
        await page.locator('#static-result').selectOption('M');
        await page.screenshot({path:path.join(root,'output/releases-static.png'),fullPage:true});
      }
      if(raw.analysis.type==='modal') await page.screenshot({path:path.join(root,'output/releases-modal.png'),fullPage:true});
      if(raw.analysis.type==='transient') {
        const internal=raw.model.dofRegistry.find(r=>r.kind==='elementEnd' && raw.analysis.globalDOFIds.includes(r.id));
        const owner='element:'+internal.elementId;
        await page.locator('#history-node').selectOption(owner);
        assert.match(await page.locator('#history-dof option:checked').textContent(),/end .*thetaZ/);
        await page.locator('#next-frame').click();
        assert(await page.locator('#history-chart .chart-series').count()>0);
        assert.equal(await page.locator('#show-deformed').isDisabled(),name.endsWith('Partial'));
      }
      assert.equal(await page.locator('#mpc-results').isVisible(),Boolean(raw.model.mpcs));
      if (raw.model.warnings.length) assert.match(await page.locator('#notice-panel').textContent(),/thetaZ restraint ignored/);
    }
    await page.goto(pathToFileURL(path.join(root,'preprocessor/tests/index.html')).href);
    await page.waitForFunction(()=>document.getElementById('test-status').dataset.status!=='running');
    assert.equal(await page.locator('#test-status').getAttribute('data-status'),'passed',await page.locator('#test-output').textContent());
    await page.goto(pathToFileURL(path.join(root,'postprocessor/tests/index.html')).href);
    await page.waitForFunction(()=>document.body.dataset.testStatus!=='running');
    assert.equal(await page.locator('body').getAttribute('data-test-status'),'passed',await page.locator('#test-output').textContent());
    assert.deepEqual(errors,[]);
    console.log('PASS releases editor, warnings, roundtrip, undo/redo, v4 topology, static/modal shapes, internal histories, partial export and legacy browser suites');
  } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
