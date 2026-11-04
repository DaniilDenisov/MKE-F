const {chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE||'playwright');
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),http=require('node:http');
const root=path.resolve(__dirname,'..'),id='12345678-1234-4234-8234-123456789abc';let submitted,done=false,result,jobStatus='running';
const server=http.createServer((req,res)=>{
  if(req.url.startsWith('/api/')){res.setHeader('Content-Type','application/json');
    if(req.url.endsWith('/health'))return res.end(JSON.stringify({status:'ready',octaveVersion:'test'}));
    if(req.method==='POST'&&req.url.endsWith('/jobs')){let body='';req.on('data',b=>body+=b);return req.on('end',()=>{submitted=JSON.parse(body);res.end(JSON.stringify({id,status:'queued'}));});}
    if(req.url.endsWith('/result'))return res.end(JSON.stringify(result));
    if(req.url.endsWith('/cancel'))jobStatus='canceled';return res.end(JSON.stringify({id,status:done?'succeeded':jobStatus}));
  }
  let file=path.join(root,decodeURIComponent(req.url.split('?')[0]));if(file.endsWith(path.sep))file+='index.html';
  try{res.setHeader('Content-Type',file.endsWith('.js')?'text/javascript':file.endsWith('.css')?'text/css':'text/html');res.setHeader('Content-Security-Policy',"default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; img-src 'self' data: blob:; connect-src 'self'");res.end(fs.readFileSync(file));}catch(_){res.statusCode=404;res.end();}
});
(async()=>{await new Promise(r=>server.listen(0,'127.0.0.1',r));const browser=await chromium.launch({headless:true,args:['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
  try{const context=await browser.newContext({viewport:{width:1450,height:1100}}),page=await context.newPage(),errors=[];context.on('page',p=>p.on('pageerror',e=>errors.push(e.message)));page.on('pageerror',e=>errors.push(e.message));page.on('dialog',d=>d.accept());
    for(const name of ['CaseSpatialFrame','CaseSpatialFrameModal','CaseSpatialFrameTransient']){
      done=false;jobStatus='running';submitted=null;result=JSON.parse(fs.readFileSync(path.join(root,'output',name+'.json')));
      await page.goto('http://127.0.0.1:'+server.address().port+'/preprocessor/');await page.locator('#file-input').setInputFiles(root+'/examples/cases/'+name+'.txt');
      await page.waitForFunction(()=>!document.querySelector('#run-analysis').disabled);await page.locator('#run-analysis').click();await page.waitForFunction(()=>document.querySelector('#run-status').textContent.includes('Running'));
      assert(submitted.caseText.includes('elems_1436'));assert(await page.locator('#run-analysis').isDisabled());
      await page.locator('#spatial-node-rows tr').nth(1).click();await page.locator('#spatial-node-x').fill('3');await page.locator('#spatial-node-apply').click();assert.match(await page.locator('#spatial-run-snapshot').innerText(),/differs/);
      done=true;await page.waitForSelector('#spatial-open-result');assert(!await page.locator('#undo').isDisabled());const popup=context.waitForEvent('page');await page.locator('#spatial-open-result').click();const view=await popup;await view.waitForSelector('#spatial-result-viewport canvas');
      assert(page.url().includes('/preprocessor/'));await page.locator('#undo').click();assert.equal(await page.locator('#spatial-node-rows tr').nth(1).locator('td').nth(1).innerText(),'2');
      if(name==='CaseSpatialFrame'){
        await view.locator('#spatial-result-table').selectOption('element');await view.locator('#spatial-result-rows tr').click();await view.locator('#spatial-result-quantity').selectOption('Mz');
        for(const [selector,ext] of [['#download-json','json'],['#export-png','png'],['#spatial-export-svg','svg'],['#spatial-export-csv','csv']]){
          const waiting=view.waitForEvent('download');await view.locator(selector).click();const download=await waiting,file=root+'/output/spatial-export.'+ext;await download.saveAs(file);const bytes=fs.readFileSync(file);assert(bytes.length>100);
          if(ext==='json')assert.deepEqual(JSON.parse(bytes),result);if(ext==='png')assert.equal(bytes.toString('ascii',1,4),'PNG');if(ext==='csv')assert(bytes.toString().includes('x/L'));if(ext==='svg')assert(bytes.toString().includes('chart-series'));
        }
      }await view.close();
    }
    for(const terminal of ['canceled','failed','timed_out']){
      done=false;jobStatus='running';await page.goto('http://127.0.0.1:'+server.address().port+'/preprocessor/');await page.locator('#file-input').setInputFiles(root+'/examples/cases/CaseSpatialFrame.txt');await page.waitForFunction(()=>!document.querySelector('#run-analysis').disabled);const original=await page.locator('#spatial-current').inputValue();await page.locator('#run-analysis').click();await page.waitForFunction(()=>document.querySelector('#run-status').textContent.includes('Running'));
      if(terminal==='canceled')await page.locator('#cancel-run').click();else jobStatus=terminal;
      await page.waitForFunction(()=>!document.querySelector('#run-analysis').disabled);assert.equal(await page.locator('#spatial-current').inputValue(),original);assert(await page.locator('#undo').isDisabled());assert.equal(await page.locator('#spatial-open-result:visible').count(),0);
    }
    done=false;jobStatus='running';await page.locator('#run-analysis').click();await page.waitForFunction(()=>document.querySelector('#run-status').textContent.includes('Running'));await page.evaluate(()=>sessionStorage.removeItem('mkef-job-context'));await page.reload();await page.waitForFunction(()=>document.querySelector('#run-status').textContent.includes('Running'));done=true;await page.waitForSelector('#spatial-open-result');assert(page.url().includes('/preprocessor/'));
    assert.deepEqual(errors,[]);console.log('Spatial Run: three analyses, immutable submission, editing/undo, cancellation/errors/resume, new tab, PNG/SVG/CSV/JSON passed.');
  }finally{await browser.close();server.close();}
})().catch(e=>{console.error(e);server.close();process.exitCode=1;});
