// Fixtures: setup(); addpath('tests/octave'); test_sections(); test_spatial_sections(); test_section_contract(); test_section_inspector();
// Fresh workspace assets, with the real section service proxied through the running web container.
// docker run --network container:mkef-web-1 ... node tests/browser/test_section_inspector_browser.cjs
const {chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE||'playwright');
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),http=require('node:http');
const root=path.resolve(__dirname,'../..'),backend=process.env.MKEF_TEST_BACKEND||'http://127.0.0.1:8080',sectionBackend=process.env.MKEF_SECTION_TEST_BACKEND||backend;
const server=http.createServer((req,res)=>{
  if(req.url.startsWith('/api/')){
    const proxy=http.request(new URL(req.url,req.url.startsWith('/api/v1/sections/')?sectionBackend:backend),{method:req.method,headers:{'content-type':req.headers['content-type']||'application/json'}},up=>{res.writeHead(up.statusCode,up.headers);up.pipe(res);});
    proxy.on('error',()=>res.writeHead(503).end());req.pipe(proxy);return;
  }
  let pathname=decodeURIComponent(req.url.split('?')[0]);if(pathname.endsWith('/'))pathname+='index.html';
  const file=path.resolve(root,'.'+pathname);
  if(!file.startsWith(root+path.sep)){res.writeHead(403).end();return;}
  fs.readFile(file,(error,data)=>{if(error){res.writeHead(404).end();return;}res.setHeader('Content-Type',file.endsWith('.html')?'text/html':file.endsWith('.js')?'text/javascript':file.endsWith('.css')?'text/css':'application/octet-stream');res.end(data);});
});
async function regressions(){
  const {spawn}=require('node:child_process');
  const names=process.env.MKEF_SECTION_REGRESSION_LIST?.split(',')||['test_mises_browser','test_section_diagnostic_browser','test_rounded_sections_browser'];
  for(const script of names)await new Promise((resolve,reject)=>{
    const child=spawn(process.execPath,[path.join(__dirname,script+'.cjs')],{stdio:'inherit',env:{...process.env,MKEF_POST_BASE_URL:`http://127.0.0.1:${server.address().port}`}});
    child.on('error',reject);child.on('exit',code=>code===0?resolve():reject(Error(script+' exited '+code)));
  });
}
(async()=>{
  await new Promise(resolve=>server.listen(0,'127.0.0.1',resolve));
  if(process.env.MKEF_SECTION_REGRESSIONS==='only'){try{await regressions();}finally{server.close();}return;}
  const browser=await chromium.launch({headless:true,args:['--no-sandbox','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
  try{
    const page=await browser.newPage({viewport:{width:1500,height:1050}}),errors=[],posts=[];
    page.setDefaultTimeout(150000);page.on('pageerror',error=>{errors.push(error.message);console.error('Browser:',error.message);});
    page.on('request',r=>{if(r.method()==='POST'&&r.url().endsWith('/sections/diagnostic'))posts.push(r.postDataJSON());});
    const uiBase=process.env.MKEF_SECTION_UI_BASE_URL||`http://127.0.0.1:${server.address().port}`;
    await page.goto(uiBase+'/postprocessor/index.html');
    const load=async(name,select=true)=>{
      await page.locator('#file-input').setInputFiles(path.join(root,'output',name+'.json'));
      assert.equal(await page.locator('#error-panel').isVisible(),false,await page.locator('#error-panel').textContent());
      const spatial=JSON.parse(fs.readFileSync(path.join(root,'output',name+'.json'))).model.dimension===3;
      if(select){await page.locator(spatial?'#spatial-result-table':'#planar-result-objects select').selectOption('element');await page.locator(spatial?'#spatial-result-rows tr':'#planar-result-rows tr').first().click();}console.log('Loaded',name);
    };
    const ready=()=>page.waitForFunction(()=>{const v=MKEFPost.sectionInspector.view;return !v.active&&!v.pending&&!v.canvas.hidden;});
    await load('CaseSpatialTrussSection');await ready();
    assert.equal(posts.length,0,'Truss must not call FEM');
    assert.match(await page.locator('#section-vm-info').innerText(),/σVM = \|σx\|/);
    await page.locator('#section-distance').fill('0.5');assert.equal(Number(await page.locator('#section-position-number').inputValue()),.25);
    assert.equal(Number(await page.locator('#spatial-result-position').inputValue()),.25);
    await page.locator('#spatial-result-position').fill('.7');await page.locator('#spatial-result-position').dispatchEvent('change');assert.equal(Number(await page.locator('#section-position-number').inputValue()),.7);
    await page.locator('#section-expand').click();assert.equal(await page.locator('dialog.section-expanded').evaluate(e=>e.open),true);
    await page.screenshot({path:'output/section-inspector-expanded.png'});
    await page.keyboard.press('Escape');await page.waitForFunction(()=>document.activeElement===document.getElementById('section-expand'));assert.equal(posts.length,0);
    await load('CaseSectionBeam');await ready();
    assert.equal(posts.length,0,'Midspan pure bending uses the analytic map');
    assert.match(await page.locator('#section-result-info').innerText(),/750000/);
    if(process.env.MKEF_SECTION_SMOKE==='1'){
      for(const file of ['postprocessor/js/section-inspector.js','postprocessor/js/section-map.js','shared/section-field.js']){
        const response=await page.request.get(uiBase+'/'+file);assert.equal(response.status(),200);assert.equal(await response.text(),fs.readFileSync(file,'utf8'));
      }
      await page.locator('#section-position-number').fill('.2');await ready();assert.equal(await page.evaluate(()=>!!MKEFPost.sectionInspector.view.entry),true);
      await page.locator('#section-expand').click();await page.screenshot({path:'output/section-inspector-live.png'});
      assert.deepEqual(errors,[]);console.log('Deployed assets, static 2D/3D section panel and real FEM field PASS');return;
    }
    const neutral=await page.evaluate(()=>{
      const v=MKEFPost.sectionInspector.view,{frame}=v.layout,[x,y]=frame.project(0,0),ctx=v.canvas.getContext('2d');
      return Array.from(ctx.getImageData(Math.floor(x),Math.floor(y),1,1).data);
    });
    assert(neutral[0]>200&&neutral[3]===255,'Neutral axis stays at the low end of the Mises scale');
    await page.locator('#section-position-number').fill('.2');await ready();
    assert(posts.length>0);assert.equal(await page.locator('#section-retry').isVisible(),false,await page.locator('#section-vm-info').innerText());
    await page.evaluate(()=>{
      const v=MKEFPost.sectionInspector.view,e=v.entry;
      if(!e)throw Error('Real FEM field missing');
      // The displayed normal stress and the service must agree at every field sample.
      let worst=0,scale=0;
      for(let k=0;k<e.field.triangles.length/6;k++){
        const ids=e.field.triangles.slice(k*6,k*6+3),p=Array.from(ids,i=>[e.field.nodes[2*i],e.field.nodes[2*i+1]]);
        e.field.local.forEach(([u,w],j)=>{const y=p[0][0]+u*(p[1][0]-p[0][0])+w*(p[2][0]-p[0][0]),z=p[0][1]+u*(p[1][1]-p[0][1])+w*(p[2][1]-p[0][1]);
          const normal=v.data.field.point(v.xi,y,z).value,offset=4*(k*15+j),[s,ty,tz,vm]=e.field.values.slice(offset,offset+4);worst=Math.max(worst,Math.abs(normal-s));scale=Math.max(scale,Math.abs(s));
          if(Math.abs(vm-Math.hypot(s,Math.sqrt(3)*Math.hypot(ty,tz)))>1e-10*Math.max(1,vm))throw Error('Mises sample mismatch');});
      }
      if(worst>1e-10*Math.max(1,scale))throw Error('Normal stress adapter mismatch');
    });
    const fetched=posts.length;await page.locator('#section-expand').click();await page.screenshot({path:'output/section-inspector-frame-expanded.png'});await page.keyboard.press('Escape');
    await page.locator('#section-position-number').fill('.5');await ready();await page.locator('#section-position-number').fill('.2');await ready();assert.equal(posts.length,fetched,'Returning to a ready field uses cache');
    const download=page.waitForEvent('download');await page.locator('#section-export-png').click();await (await download).saveAs('output/section-inspector-export.png');
    await page.screenshot({path:'output/section-inspector-desktop.png'});
    assert.equal(posts.length,fetched,'Expand and export must not submit');
    await page.setViewportSize({width:1500,height:1050});
    await load('CaseSpatialFrameSection');await ready();assert.equal(posts.at(-1).nu,.25);
    await page.screenshot({path:'output/section-inspector-spatial.png'});
    // All fresh solver profiles/orientations: no browser geometry property calculation.
    const contracts=fs.readdirSync('output/section-contract').filter(n=>n.endsWith('.json')).map(n=>JSON.parse(fs.readFileSync('output/section-contract/'+n)));
    assert.equal(contracts.length,42);
    await page.evaluate(cases=>{
      const S=MKEFSections,P=MKEFSpatialSections,old=[S.properties,P.properties,P.stress,S.roots];
      S.properties=P.properties=P.stress=S.roots=()=>{throw Error('Unexpected solver-property calculation');};
      try{for(const raw of cases){const ds=MKEFPost.validateDataset(raw);for(const e of raw.model.elements){const d=MKEFPost.sectionMap.adapt(ds,e,.3);if(!d.profile)continue;
        for(const xi of [0,.17,.5,.83,1]){const actual=d.field.at(xi),expected=d.normal.at(xi);for(const [a,b]of [[actual.minimum.value,expected.sigmaMin],[actual.maximum.value,expected.sigmaMax]])if(Math.abs(a-b)>1e-11*Math.max(1,d.normal.sigmaAbsMax))throw Error('Adapter stress mismatch: '+JSON.stringify({dimension:raw.model.dimension,profile:d.section.type,angle:d.angle,xi,a,b}));}
      }}}finally{[S.properties,P.properties,P.stress,S.roots]=old;}
    },contracts);
    console.log('42 solver profiles/orientations: normal stress parity PASS');
    let count=posts.length;
    await load('SectionInspectorCompression');await ready();assert.equal(posts.length,count);
    assert.equal(await page.evaluate(()=>{const v=MKEFPost.sectionInspector.view,a=v.data.normal.at(v.xi);return a.sigmaMax<0&&a.sigmaAbsMax===-a.sigmaMin;}),true);
    await load('SectionInspectorZero');await ready();assert.equal(posts.length,count);assert.match(await page.locator('#section-vm-info').innerText(),/Maximum σVM = 0 /);
    await load('SectionInspectorTorsion');await ready();
    const torsion=await page.evaluate(()=>{const v=MKEFPost.sectionInspector.view;return {sigma:v.data.normal.sigmaAbsMax,vm:v.entry.result.maximum.value,forces:v.request().forces};});
    assert.equal(torsion.sigma,0);assert(Math.abs(torsion.vm/(Math.sqrt(3)*400*.05/(Math.PI*.1**4/32))-1)<.01,'Circular torsion analytical peak');
    assert.deepEqual(torsion.forces.filter((_,i)=>i!==3),[0,0,0,0,0]);
    console.log('Compression, zero load and circular torsion analytical peak PASS');
    await load('SectionCase_square_tube_0');await page.locator('#section-position-number').fill('.2');await ready();
    assert.match(await page.locator('#section-vm-info').innerText(),/Mesh peak.*sharp inner corner/);
    assert.equal(await page.evaluate(()=>{const v=MKEFPost.sectionInspector.view,p=v.layout.frame.project(0,0);return v.canvas.getContext('2d').getImageData(Math.floor(p[0]),Math.floor(p[1]),1,1).data[3];}),0,'Tube hole is empty');
    await load('Rounded3D_rectangular_tube');await ready();assert.equal(posts.at(-1).section.dimensions.inner_radius,.005);
    assert(!/sharp inner corner/.test(await page.locator('#section-vm-info').innerText()));
    console.log('Real sharp/rounded tube fields and empty holes PASS');
    // One active job, one latest pending coordinate. Hold an actual completed response.
    await load('CaseSectionBeam');await ready();let releaseGate,posted;
    const gate=new Promise(r=>releaseGate=r),submitted=new Promise(r=>posted=r);
    await page.route('**/sections/diagnostic',async route=>{const response=await route.fetch();posted((await response.json()).id);await gate;await route.fulfill({response});},{times:1});
    const before=posts.length;await page.locator('#section-position-number').fill('.11');const firstId=await submitted;
    for(const xi of ['.21','.31','.41'])await page.locator('#section-position-number').fill(xi);
    await page.waitForFunction(()=>MKEFPost.sectionInspector.view.pending?.body.forces);
    assert.equal(posts.length,before+1,'Intermediate positions must not submit');assert.equal(await page.locator('#section-vm-field').isVisible(),false,'Stale field must not appear at new coordinate');
    releaseGate();await ready();assert.equal(posts.length,before+2);assert.equal(Number(await page.locator('#section-position-number').inputValue()),.41);
    assert.equal(await page.evaluate(()=>MKEFPost.sectionInspector.view.currentKey===JSON.stringify(MKEFPost.sectionInspector.view.request())),true);
    // Late POST after dataset replacement is deleted, never rendered into the new dataset.
    let releaseLate,latePosted;const lateGate=new Promise(r=>releaseLate=r),lateReady=new Promise(r=>latePosted=r);
    await page.route('**/sections/diagnostic',async route=>{const response=await route.fetch();latePosted((await response.json()).id);await lateGate;await route.fulfill({response});},{times:1});
    await page.locator('#section-position-number').fill('.12');const lateId=await lateReady;await load('CaseSpatialTrussSection');await ready();
    const deleted=page.waitForResponse(r=>r.url().endsWith('/jobs/'+lateId)&&r.request().method()==='DELETE');releaseLate();await deleted;
    assert.match(await page.locator('#section-vm-info').innerText(),/σVM = \|σx\|/);
    console.log('Pending-coordinate coalescing and late-response cleanup PASS',firstId);
    // Service errors do not erase normal stress; retries are explicit.
    await load('CaseSectionBeam');await ready();
    for(const status of [429,503]){
      await page.route('**/sections/diagnostic',route=>route.fulfill({status,body:'{}'}),{times:1});
      const count=posts.length;await page.locator('#section-position-number').fill(status===429?'.13':'.14');
      await page.locator('#section-retry').waitFor({state:'visible'});
      assert.equal(await page.locator('#section-result-svg').isVisible(),true);assert.equal(await page.locator('#section-vm-field').isVisible(),false);
      await page.waitForTimeout(700);assert.equal(posts.length,count+1,'No automatic retry storm');
    }
    await page.locator('#section-retry').click();await ready();
    // Keep a valid field while explicitly marking an unconfirmed numerical result.
    await page.route('**/sections/jobs/*',async route=>{
      const response=await route.fetch();
      if(route.request().method()==='GET'&&!route.request().url().endsWith('/field')){
        const job=await response.json();if(job.state==='completed'&&job.result?.display?.status==='ok'){
          job.result.fieldStatus=job.result.maximumStatus='nonconverged';job.result.status='nonconverged';
          await route.fulfill({response,json:job});return;
        }
      }await route.fulfill({response});
    });
    await page.locator('#section-position-number').fill('.18');await ready();assert.match(await page.locator('#section-vm-info').innerText(),/convergence is not confirmed/);
    await page.unroute('**/sections/jobs/*');
    // The shared 2D nu setting updates both consumers and invalidates only the field request.
    await page.locator('#section-nu').fill('.27');await ready();assert.equal(await page.locator('#member-nu').inputValue(),'.27');assert.equal(posts.at(-1).nu,.27);
    await load('section-contract/2D-tee-90');await page.locator('#section-position-number').fill('.2');
    assert.match(await page.locator('#section-vm-info').innerText(),/not supported for this profile/);assert.equal(await page.locator('#section-result-svg').isVisible(),true);
    // A truss without catalog geometry still has numeric stress, never a fabricated map.
    const truss=JSON.parse(fs.readFileSync('output/CaseSpatialTrussSection.json'));truss.model.elements[0].sectionId=0;
    await page.locator('#file-input').setInputFiles({name:'unassigned-truss.json',mimeType:'application/json',buffer:Buffer.from(JSON.stringify(truss))});
    await page.locator('#spatial-result-table').selectOption('element');await page.locator('#spatial-result-rows tr').first().click();
    assert.match(await page.locator('#section-vm-info').innerText(),/No section geometry/);assert.equal(await page.locator('#section-vm-field').isVisible(),false);
    await load('CaseSpatialFrameModalSection',false);assert.equal(await page.locator('#section-result-panel').isVisible(),false);
    assert.equal(await page.locator('#spatial-section-svg').count(),1,'Dynamic section view must not share IDs with the hidden static view');
    console.log('Unavailable service/profile/geometry and dynamic isolation PASS');
    assert.deepEqual(errors,[]);console.log('Section inspector: local fields, real FEM samples, coordinates, cache, exports and layout PASS');
    if(process.env.MKEF_SECTION_REGRESSIONS==='1')await regressions();
  }finally{await browser.close();server.close();}
})().catch(error=>{console.error(error);process.exitCode=1;server.close();});
