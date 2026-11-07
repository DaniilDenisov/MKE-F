const fs=require('node:fs'),assert=require('node:assert/strict'),path=require('node:path'),{gzipSync}=require('node:zlib');
const {chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE||'playwright'),context=require('../helpers/recovery-test-context.cjs');
const root=path.resolve(__dirname,'../..'),directory=root+'/output/recovery-pool',read=file=>JSON.parse(fs.readFileSync(file));
const M=context(),S=M.spatialResults;
function compare(inline,table){
  const before=JSON.stringify(table),a=M.validateDataset(inline),b=M.validateDataset(table);
  const count=inline.analysis.type==='modal'?inline.analysis.frequenciesHz.length:inline.analysis.type==='transient'?inline.analysis.time.length:1;
  for(const [i,e] of inline.model.elements.entries()){
    for(const frame of new Set([0,Math.floor(count/2),count-1])){
      const old=S.recover(a,e,S.vector(a,'displacements',frame),S.vector(a,'accelerations',frame));
      const pooled=S.recover(b,table.model.elements[i],S.vector(b,'displacements',frame),S.vector(b,'accelerations',frame));
      assert.deepEqual(pooled,old);
    }
    if(inline.analysis.type!=='static'){
      const op=S.operators(table.model.elements[i],8,table.model);
      assert.strictEqual(op,table.model.recoveryOperators[table.model.elements[i].recoveryOperatorId-1]);
      assert.deepEqual(op.K,e.recoveryOperators.K);assert.deepEqual(op.M,e.recoveryOperators.M);
    }
  }
  assert.equal(JSON.stringify(table),before,'Reader must not expand/mutate table-backed raw JSON');
}
const pairs=read(root+'/output/recovery/manifest.json');assert.equal(pairs.length,12);
for(const row of pairs)compare(read(root+'/output/recovery/'+row.name+'-v8.json'),read(root+'/output/recovery/'+row.name+'-table.json'));
const report=read(directory+'/measurements.json');assert.equal(report.length,6);
for(const row of report){
  const inline=read(directory+'/'+row.name+'-inline.json'),table=read(directory+'/'+row.name+'-table.json');
  assert.equal(table.model.recoveryOperators.length,row.uniqueOperators);
  compare(inline,table);
  const auto=read(directory+'/'+row.name+'-auto.json');M.validateDataset(auto);
  assert(row.variants.auto.bytes<=row.variants.inline.bytes);
}
const source=read(directory+'/solved-grid-modal-table.json');
const invalid=[
  x=>delete x.model.recoveryOperators,
  x=>x.model.recoveryOperators=[],
  x=>x.model.recoveryOperators[0].id=2,
  x=>x.model.recoveryOperators[1].id=1,
  x=>x.model.elements[0].recoveryOperatorId=0,
  x=>x.model.elements[0].recoveryOperatorId=.5,
  x=>x.model.elements[0].recoveryOperatorId=999,
  x=>x.model.elements[0].recoveryOperatorId='1',
  x=>delete x.model.elements[0].recoveryOperatorId,
  x=>x.model.elements[0].recoveryOperators=x.model.recoveryOperators[0],
  x=>x.model.recoveryOperators[0].K.pop(),
  x=>x.model.recoveryOperators[0].M[0][0]=NaN,
  x=>x.model.recoveryOperators[0].K[0][1]=1e20,
  x=>x.model.recoveryOperators.push({...structuredClone(x.model.recoveryOperators[0]),id:3}),
  x=>x.model.recovery.version=1,
  x=>x.model.recovery.version=3
];
for(const mutate of invalid){const bad=structuredClone(source);mutate(bad);assert.throws(()=>M.validateDataset(bad));}
const staticBad=read(root+'/output/recovery/CaseSpatialFrame-v8.json');
staticBad.model.recovery.version=2;staticBad.model.recoveryOperators=source.model.recoveryOperators;
assert.throws(()=>M.validateDataset(staticBad),/static results/);
const inlineBad=read(directory+'/single-modal-inline.json');inlineBad.model.elements[0].recoveryOperatorId=1;
assert.throws(()=>M.validateDataset(inlineBad),/reference requires/);
console.log('Pool recovery: 12 paired cases, solved grids, 6000 synthetic members, shared references, unchanged raw JSON and 18 invalid contracts PASS');

(async()=>{
  const browser=await chromium.launch({headless:true,args:['--no-sandbox','--use-gl=angle','--use-angle=swiftshader',
    '--enable-unsafe-swiftshader','--enable-precise-memory-info','--js-flags=--expose-gc']});
  try{
    const page=await browser.newPage({viewport:{width:1450,height:1100}}),errors=[];
    page.on('pageerror',e=>errors.push(e.message));
    await page.goto('file://'+root+'/postprocessor/index.html');
    await page.evaluate(()=>{
      MKEFPost.SpatialPost=new Proxy(MKEFPost.SpatialPost,{construct(T,args){window.viewer=Reflect.construct(T,args);return viewer;}});
    });
    for(const name of ['solved-grid-modal','solved-grid-transient']){
      await page.evaluate(()=>window.previousViewer=window.viewer);
      await page.locator('#file-input').setInputFiles(directory+'/'+name+'-table.json');
      await page.waitForFunction(()=>window.viewer&&viewer!==window.previousViewer);
      assert.equal(await page.locator('#error-panel').isVisible(),false);
      await page.locator('#spatial-result-table').selectOption('element');await page.locator('#spatial-result-rows tr').first().click();
      await page.locator('#spatial-result-quantity').selectOption('Mz');
      assert(await page.locator('#spatial-result-chart').isVisible());
      if(name.endsWith('transient')){
        await page.waitForFunction(()=>viewer.prepared);await page.evaluate(()=>{viewer.setFrame(3);viewer.historyKind.value='elementForce';viewer.drawChart();});
        assert(await page.evaluate(()=>viewer.series.y.length===viewer.analysis.time.length&&viewer.series.y.every(Number.isFinite)));
      }
    }
    const download=page.waitForEvent('download');await page.locator('#download-json').click();
    await (await download).saveAs(directory+'/roundtrip.json');
    assert.deepEqual(read(directory+'/roundtrip.json'),read(directory+'/solved-grid-transient-table.json'));
    assert.deepEqual(errors,[]);await page.close();
    const measurements=[];
    // Isolated pages and five measured repetitions after warm-up. Text transfer
    // and file I/O are outside timing; GC excludes discarded previous datasets.
    for(const row of report)for(const mode of ['inline','table','auto']){
      const bytes=fs.readFileSync(directory+'/'+row.name+'-'+mode+'.json');
      const p=await browser.newPage();await p.goto('file://'+root+'/postprocessor/index.html');
      await p.evaluate(text=>{window.benchmarkText=text;},bytes.toString('utf8'));
      const cdp=await p.context().newCDPSession(p),times=[],heaps=[];
      for(let i=0;i<6;i++){
        await p.evaluate(()=>{window.benchmarkDataset=null;});
        await cdp.send('HeapProfiler.collectGarbage');
        const before=(await cdp.send('Runtime.getHeapUsage')).usedSize;
        const elapsed=await p.evaluate(()=>{
          const start=performance.now();
          window.benchmarkDataset=MKEFPost.validateDataset(JSON.parse(window.benchmarkText));
          return performance.now()-start;
        });
        await cdp.send('HeapProfiler.collectGarbage');
        const after=(await cdp.send('Runtime.getHeapUsage')).usedSize;
        if(i){times.push(elapsed);heaps.push(after-before);}
      }
      const median=values=>values.slice().sort((a,b)=>a-b)[Math.floor(values.length/2)];
      const result={parseAndValidateMs:median(times),retainedHeapBytes:median(heaps),samplesMs:times};
      assert(result.retainedHeapBytes>0,'CDP must observe the retained parsed dataset');
      assert(result.parseAndValidateMs<30000);await p.close();
      measurements.push({name:row.name,elements:row.elements,uniqueOperators:row.uniqueOperators,
        synthetic:row.synthetic,mode,bytes:bytes.length,gzipBytes:gzipSync(bytes,{level:9}).length,
        exportSeconds:row.variants[mode].exportSeconds,...result});
    }
    fs.writeFileSync(directory+'/browser-measurements.json',JSON.stringify({browser:browser.version(),measurements},null,2)+'\n');
    for(const row of report){const selected=measurements.filter(x=>x.name===row.name),a=selected[0],b=selected[2];
      console.log(`${row.name}: auto JSON ${(100*(1-b.bytes/a.bytes)).toFixed(1)}% smaller; gzip ${(100*(1-b.gzipBytes/a.gzipBytes)).toFixed(1)}%; reader ${a.parseAndValidateMs.toFixed(1)} -> ${b.parseAndValidateMs.toFixed(1)} ms`);}
    console.log('Pooled browser: actual grid charts, JSON download round-trip and 18 measured layouts PASS');
  }finally{await browser.close();}
})().catch(error=>{console.error(error);process.exitCode=1;});
