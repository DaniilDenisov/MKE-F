const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict'),{gzipSync}=require('node:zlib'),{chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE||'playwright');
const read=p=>JSON.parse(fs.readFileSync(p)),base='output/result-resources';
const rows=read(base+'/exports.json').map(r=>({...r,file:base+'/'+r.name+'.json'}));
assert.equal(rows.length,4);
for(const r of read('output/recovery-pool/measurements.json').filter(r=>!r.synthetic&&r.name.startsWith('solved-grid')))
 rows.push({name:r.name,elements:r.elements,exportSeconds:r.variants.auto.exportSeconds,file:'output/recovery-pool/'+r.name+'-auto.json'});
// Explicitly synthetic repeated topology: tests the 2,000-member UI target,
// not a physical 2,000-member solve or its export time.
const repeated=read(base+'/beam-2d-modal.json'),element=repeated.model.elements[0];
repeated.model.elements=Array.from({length:2000},(_,i)=>({...element,id:i+1}));
const repeatedFile=base+'/synthetic-2d-2000.json';fs.writeFileSync(repeatedFile,JSON.stringify(repeated));
rows.push({name:'synthetic-2d-2000',elements:2000,synthetic:true,exportSeconds:null,file:repeatedFile});
(async()=>{
 const browser=await chromium.launch({headless:true,args:['--no-sandbox','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
 try{
  const measurements=[];
  for(const row of rows){
   const page=await browser.newPage(),errors=[];page.on('pageerror',e=>errors.push(e.message));
   await page.goto('file://'+path.resolve('postprocessor/index.html'));
   await page.evaluate(()=>{
    for(const [p,key]of [[MKEFPost.WebGLRenderer.prototype,'flush'],[MKEFPost.SpatialPost.prototype,'update']]){const original=p[key];p[key]=function(...args){window.resourceFrames=(window.resourceFrames||0)+1;return original.apply(this,args);};}
   });
   const file=fs.readFileSync(row.file);assert(file.length<100*1024*1024);
   await page.evaluate(text=>window.resourceText=text,file.toString('utf8'));
   const cdp=await page.context().newCDPSession(page),times=[],heaps=[];
   for(let i=0;i<4;i++){
    await page.evaluate(()=>window.resourceDataset=null);await cdp.send('HeapProfiler.collectGarbage');const before=(await cdp.send('Runtime.getHeapUsage')).usedSize;
    const ms=await page.evaluate(()=>{const t=performance.now();window.resourceDataset=MKEFPost.validateDataset(JSON.parse(resourceText));return performance.now()-t;});
    await cdp.send('HeapProfiler.collectGarbage');const after=(await cdp.send('Runtime.getHeapUsage')).usedSize;
    if(i){times.push(ms);heaps.push(after-before);}
   }
   const spatial=await page.evaluate(()=>resourceDataset.raw.model.dimension===3);
   await page.locator('#file-input').setInputFiles(path.resolve(row.file));await page.waitForFunction(()=>!document.querySelector('#download-json').disabled);
   assert(await page.locator('#error-panel').isHidden());
   const type=await page.evaluate(()=>resourceDataset.raw.analysis.type);let animationFrames=null;
   if(type!=='static'){
    if(spatial)await page.waitForFunction(()=>!document.querySelector('#spatial-result-play').disabled);
    const frames=await page.evaluate(async({spatial,type})=>{
     const id=spatial?'#spatial-result-play':type==='transient'?'#play-pause':'#modal-play-pause';
     // A 21-sample spatial history advances every 250 ms at 1x. Observe
     // a full second so two updates do not depend on a timer boundary.
     const button=document.querySelector(id);window.resourceFrames=0;button.click();await new Promise(resolve=>setTimeout(resolve,1000));
     const frames=window.resourceFrames;if(/Pause/.test(button.textContent))button.click();return frames;
    },{spatial,type});assert(frames>=2,row.name+': animation must update the scene');animationFrames=frames;
   }
   const median=a=>a.slice().sort((x,y)=>x-y)[Math.floor(a.length/2)],result={...row,bytes:file.length,gzipBytes:gzipSync(file).length,parseAndValidateMs:median(times),retainedHeapBytes:median(heaps),animationFrames,animationObservationMs:animationFrames===null?null:1000};
   assert(result.parseAndValidateMs<30000);assert(result.retainedHeapBytes>0);assert.deepEqual(errors,[]);
   measurements.push(result);console.log(`${row.name}: ${(result.bytes/1048576).toFixed(2)} MiB JSON, ${result.parseAndValidateMs.toFixed(1)} ms read, ${(result.retainedHeapBytes/1048576).toFixed(2)} MiB data heap, ${row.exportSeconds?.toFixed(2)??'n/a'} s export`);await page.close();
  }
  fs.writeFileSync(base+'/browser-measurements.json',JSON.stringify({browser:browser.version(),measurements},null,2));
 }finally{await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
