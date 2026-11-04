const {chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE||'playwright');
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path');
const root=path.resolve(__dirname,'..'),base=root+'/output/recovery/';
(async()=>{
  const browser=await chromium.launch({headless:true,args:['--no-sandbox','--use-gl=angle','--use-angle=swiftshader',
    '--enable-unsafe-swiftshader','--enable-precise-memory-info','--js-flags=--expose-gc']});
  try{
    const page=await browser.newPage({viewport:{width:1450,height:1100}}),errors=[],network=[];
    page.on('pageerror',e=>errors.push(e.message));page.on('request',r=>{if(/^https?:/.test(r.url()))network.push(r.url());});
    await page.goto('file://'+root+'/postprocessor/index.html');
    await page.evaluate(()=>{
      MKEFPost.SpatialPost=new Proxy(MKEFPost.SpatialPost,{construct(T,args){window.viewer=Reflect.construct(T,args);return viewer;}});
    });
    // Exercise the actual file picker and application routing, not only validators.
    const load=async file=>{
      await page.evaluate(()=>window.previousViewer=window.viewer);
      await page.locator('#file-input').setInputFiles(file);
      await page.waitForFunction(()=>window.viewer&&window.viewer!==window.previousViewer);
      assert.equal(await page.locator('#error-panel').isVisible(),false);
    };
    const frozen=JSON.parse(fs.readFileSync(root+'/tests/fixtures/postprocessor-legacy/manifest.json'));
    for(const f of frozen){
      // Reset prevents an earlier file with the same title from satisfying the wait.
      await page.reload();
      await page.locator('#file-input').setInputFiles(root+'/tests/fixtures/postprocessor-legacy/'+f.file);
      await page.waitForFunction(()=>document.querySelector('#error-panel').textContent.includes('Recalculate'));
      assert.equal(await page.locator('#download-json').isDisabled(),true,f.file);
      assert.equal(await page.locator('#spatial-result-viewport canvas').count(),0);
    }
    await page.reload();
    await page.evaluate(()=>{
      MKEFPost.SpatialPost=new Proxy(MKEFPost.SpatialPost,{construct(T,args){window.viewer=Reflect.construct(T,args);return viewer;}});
      MKEFPost.legacySpatialOperators=()=>{throw Error('v8 called legacy matrix formulas');};
    });
    for(const name of ['CaseSpatialFrame','CaseSpatialTruss','CaseSpatialTrussModal','CaseSpatialFrameModal',
      'CaseSpatialTrussTransient','released-transient']){
      await load(base+name+'-v8.json');
      await page.locator('#spatial-result-table').selectOption('element');
      await page.locator('#spatial-result-rows tr').first().click();
      await page.locator('#spatial-result-quantity').selectOption('N');
      assert(await page.locator('#spatial-result-chart').isVisible());
      if(name==='released-transient'){
        await page.waitForFunction(()=>viewer.prepared);
        await page.evaluate(()=>{viewer.setFrame(3);viewer.historyKind.value='elementForce';viewer.quantity.value='Mz';viewer.drawChart();});
        assert(await page.evaluate(()=>viewer.series.y.length===viewer.analysis.time.length&&viewer.series.y.every(Number.isFinite)));
        const download=page.waitForEvent('download');await page.locator('#spatial-export-csv').click();
        await (await download).saveAs(base+'released-force.csv');
        await page.locator('#spatial-result-play').click();
        await page.waitForFunction(()=>viewer.index!==3);
        await page.evaluate(()=>viewer.stop());
      }
    }
    const bad=JSON.parse(fs.readFileSync(base+'CaseSpatialFrameModal-v8.json'));delete bad.model.elements[0].recoveryOperators;
    await page.locator('#file-input').setInputFiles({name:'bad.json',mimeType:'application/json',buffer:Buffer.from(JSON.stringify(bad))});
    await page.waitForFunction(()=>document.querySelector('#error-panel').textContent.includes('missing dynamic recovery operators'));
    assert.equal(await page.evaluate(()=>viewer.analysis.type),'transient');
    await load(base+'released-partial-v8.json');await page.waitForFunction(()=>viewer.prepared);
    assert(await page.locator('#spatial-result-play').isDisabled());
    await load(base+'released-no-displacements-v8.json');await page.waitForFunction(()=>viewer.prepared);
    assert(await page.locator('#spatial-result-play').isDisabled());
    await page.screenshot({path:base+'partial-v8.png'});
    // Deterministic fixture sizes, browser validation time and retained heap.
    // Synthetic repeats isolate per-element format cost; they are not a new solve.
    const measurements=[];
    for(const name of ['CaseSpatialFrame','benchmark-modal','benchmark-transient']){
      for(const version of [8]){
        const original=JSON.parse(fs.readFileSync(base+name+'-v'+version+'.json'));
        for(const count of [1,2000]){
          const stats=await page.evaluate(({original,count})=>{
            if(window.viewer){viewer.dispose();window.viewer=null;}
            window.benchmarkDataset=null;if(window.gc)gc();
            const data=structuredClone(original),element=data.model.elements[0];
            data.model.elements=Array.from({length:count},(_,i)=>({...structuredClone(element),id:i+1}));
            if(data.analysis.elementResults){const r=data.analysis.elementResults[0];data.analysis.elementResults=Array.from({length:count},(_,i)=>({...structuredClone(r),elementId:i+1}));}
            const text=JSON.stringify(data),before=performance.memory?.usedJSHeapSize??null,start=performance.now();
            window.benchmarkDataset=MKEFPost.validateDataset(JSON.parse(text));
            const milliseconds=performance.now()-start;if(window.gc)gc();
            return {jsonBytes:new TextEncoder().encode(text).length,parseAndValidateMs:milliseconds,
              usedHeapBytes:performance.memory?.usedJSHeapSize??null,heapBeforeBytes:before};
          },{original,count});
          measurements.push({name,version,elements:count,...stats});
          assert(stats.jsonBytes<100*1024*1024);
          assert(stats.parseAndValidateMs<30000,'Validation exceeded 30 s diagnostic ceiling');
        }
      }
    }
    fs.writeFileSync(base+'browser-measurements.json',JSON.stringify({browser:browser.version(),measurements},null,2)+'\n');
    assert.deepEqual(errors,[]);assert.deepEqual(network,[]);
    console.log('Browser: rejected v1–v7; v8 static/modal/transient, CSV, animation, malformed/partial results and 6 measurements PASS');
  }finally{await browser.close();}
})().catch(error=>{console.error(error);process.exitCode=1;});
