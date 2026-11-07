const {chromium} = require(process.env.MKEF_PLAYWRIGHT_MODULE || 'playwright');
const fs = require('node:fs'), path = require('node:path'), assert = require('node:assert/strict'), http = require('node:http');
const root = path.resolve(__dirname,'../..');
(async () => {
  const browser = await chromium.launch({headless:true,args:['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
  try {
    const page = await browser.newPage({viewport:{width:1400,height:950}}), errors = [];
    page.on('pageerror',e=>errors.push(e.message));
    const external = []; page.on('request',r=> { if (/^https?:/.test(r.url()) && !r.url().startsWith('http://127.0.0.1:')) external.push(r.url()); });
    await page.goto('file://'+root+'/preprocessor/index.html');
    await page.locator('#file-input').setInputFiles(root+'/examples/cases/CaseSpatialFrame.txt');
    await page.waitForFunction(()=>document.querySelector('#case-status').textContent.includes('CaseSpatialFrame'));
    assert.equal(await page.locator('#spatial-viewport canvas').count(),2);
    assert.equal(await page.locator('#spatial-node-rows tr').count(),2);
    assert(await page.locator('#run-analysis').isDisabled());
    await page.locator('#spatial-node-rows tr').nth(1).click();
    assert.match(await page.locator('#selection-details').innerText(),/node 2/);
    assert.equal(await page.locator('#spatial-source').inputValue(),fs.readFileSync(root+'/examples/cases/CaseSpatialFrame.txt','utf8'));
    // A malformed import must not replace the model, selection or original source.
    const source = fs.readFileSync(root+'/examples/cases/CaseSpatialFrame.txt','utf8');
    await page.locator('#file-input').setInputFiles({name:'bad.txt',mimeType:'text/plain',buffer:Buffer.from(source+'\nunknown_block\n0')});
    await page.waitForFunction(()=>document.querySelector('#error-panel').textContent.includes('Unsupported block'));
    assert.equal(await page.locator('#spatial-source').inputValue(),source);
    assert.match(await page.locator('#selection-details').innerText(),/node 2/);
    for (const name of fs.readdirSync(root+'/examples/cases').filter(f=>/^CaseSpatial.*\.txt$/.test(f))) {
      await page.locator('#file-input').setInputFiles(root+'/examples/cases/'+name);
      await page.waitForFunction(name=>document.querySelector('#case-status').textContent.startsWith(name.replace('.txt','')),name);
    }
    // Preserve BOM, CRLF, comments, order and final newline in the actual downloaded bytes.
    const exact = Buffer.from('\ufeff# untouched source\r\n'+source.replace(/\r?\n/g,'\r\n'));
    await page.locator('#file-input').setInputFiles({name:'exact.txt',mimeType:'text/plain',buffer:exact});
    await page.waitForFunction(()=>document.querySelector('#case-status').textContent.startsWith('exact'));
    await page.locator('#spatial-plane').selectOption('XZ'); await page.locator('#spatial-offset').fill('-4'); await page.locator('#spatial-offset').dispatchEvent('change');
    await page.locator('[data-view="plane"]').click(); await page.locator('[data-view="origin"]').click();
    await page.waitForFunction(()=>document.querySelector('#spatial-plane').value === 'XY');
    const downloadEvent = page.waitForEvent('download'); await page.locator('#download-case').click();
    const download = await downloadEvent; assert.deepEqual(fs.readFileSync(await download.path()),exact);
    assert(await page.locator('#undo').isDisabled());
    await page.screenshot({path:root+'/output/spatial-view.png'});
    await page.locator('#file-input').setInputFiles(root+'/examples/cases/CaseSpatialTruss.txt');
    await page.waitForFunction(()=>document.querySelector('#case-status').textContent.includes('CaseSpatialTruss'));
    assert.equal(await page.locator('#spatial-viewport canvas').count(),2);
    await page.locator('#new-case').click(); assert(await page.locator('#spatial-viewport').isVisible());
    assert.equal(await page.locator('#spatial-viewport canvas').count(),2);
    // Test the shared scene directly, independent of app tables or solver state.
    const cameraChecks = await page.evaluate(async () => {
      const host = document.createElement('div'); host.style.cssText = 'position:fixed;inset:0;width:800px;height:600px;z-index:100'; document.body.append(host);
      const scene = new MKEFSpatial.SpatialScene(host), frame = () => new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)));
      const model = {nodes:[{x:1e9,y:2e9,z:3e9},{x:1e9+4,y:2e9+3,z:3e9+2}],elements:[{node1:1,node2:2}]};
      scene.setModel(model); await frame();
      const results = [];
      const check = (ok,name) => { if (!ok) throw new Error(name); results.push(name); };
      const close = (a,b) => a.every((n,i)=>Math.abs(n-b[i])<1e-8);
      const fit = () => scene.corners().every(p=> { p.project(scene.camera); return Math.abs(p.x)<=1.00001 && Math.abs(p.y)<=1.00001 && Math.abs(p.z)<=1; });
      check(fit(),'all bounding-box corners fit in initial isometry');
      for (const plane of ['XY','XZ','YZ']) for (const offset of [-10,10]) {
        const before = scene.getViewState(); scene.setWorkPlane({plane,offset}); await frame(); const after = scene.getViewState();
        check(close(before.target,after.target) && close(before.direction,after.direction) && before.height === after.height,'plane independent '+plane+offset);
      }
      const direction = scene.getViewState().direction; scene.controls.target.x += 100; scene.camera.position.x += 100; scene.fitModel(); await frame();
      check(close(direction,scene.getViewState().direction) && fit(),'fit preserves direction after pan');
      scene.resetOrigin(); await frame(); const origin = scene.getViewState();
      check(close(origin.target,[0,0,0]) && close(origin.direction,[0,0,1]) && origin.workPlane.plane === 'XY' && origin.workPlane.offset === 0 && fit(),'reset origin fits distant model: '+JSON.stringify({origin,fit:fit()}));
      scene.setView('XZ'); scene.fitModel(); await frame(); check(fit(),'XZ fit');
      const beforeResize = scene.getViewState(); host.style.width = '350px'; scene.resize(); await frame();
      check(scene.getViewState().height === beforeResize.height && close(beforeResize.target,scene.getViewState().target),'resize preserves vertical scale and target');
      scene.fitModel(); await frame(); check(fit(),'fit in portrait viewport');
      scene.setWorkPlane({plane:'YZ',offset:-12}); scene.setView('plane'); await frame(); check(scene.getViewState().target[0] === -12,'plane view projects target');
      scene.setModel({nodes:[],elements:[]}); scene.resetOrigin(); await frame(); check(scene.getViewState().height === 20,'empty origin height');
      scene.setModel({nodes:[{x:0,y:0,z:0},{x:0,y:0,z:10}],elements:[{node1:1,node2:2}]}); scene.setView('XY'); scene.fitModel(); await frame();
      check(scene.getViewState().height === 20 && fit(),'end-on member has finite viewport');
      scene.setModel({nodes:[{x:7,y:8,z:9}],elements:[]}); scene.fitModel(); await frame(); check(scene.getViewState().height === 20 && fit(),'single node');
      // Picking uses actual screen positions; an edge midpoint must not select an endpoint.
      const pickModel = {nodes:[{x:0,y:0,z:0},{x:10,y:0,z:0},{x:10,y:10,z:5}],elements:[{node1:1,node2:2},{node1:2,node2:3}]};
      scene.setModel(pickModel); scene.setView('XY'); scene.fitModel(); await frame();
      const screen = n => { const p = scene.local(n).project(scene.camera), b = scene.canvas.getBoundingClientRect(); return [b.left+(p.x+1)*b.width/2,b.top+(1-p.y)*b.height/2]; };
      scene.pick(...screen(pickModel.nodes[0])); check(scene.selection.kind === 'node' && scene.selection.index === 0,'node picking');
      scene.pick(...screen({x:5,y:0,z:0})); check(scene.selection.kind === 'element' && scene.selection.index === 0,'member picking');
      scene.setView('iso'); await frame();
      const depth = pickModel.nodes.map(n=>scene.local(n).project(scene.camera).z); check(new Set(depth).size === 3,'depth differs for spatial nodes');
      const saved = JSON.stringify(scene.getViewState()), extension = scene.renderer.getContext().getExtension('WEBGL_lose_context');
      if (!extension) throw new Error('Missing context loss test extension');
      await new Promise((resolve,reject)=>{ const timer = setTimeout(()=>reject(new Error('context loss timed out')),5000); scene.canvas.addEventListener('webglcontextlost',()=>{clearTimeout(timer);resolve();},{once:true}); extension.loseContext(); });
      check(scene.lost && !scene.message.hidden,'context loss message');
      await new Promise((resolve,reject)=>{ const timer = setTimeout(()=>reject(new Error('context restoration timed out')),5000); scene.canvas.addEventListener('webglcontextrestored',()=>{clearTimeout(timer);resolve();},{once:true}); setTimeout(()=>extension.restoreContext(),100); }); await frame();
      check(!scene.lost && scene.message.hidden && saved === JSON.stringify(scene.getViewState()),'context restoration preserves view');
      const baseline = scene.renderer.info.memory.geometries;
      for (let i = 0; i < 10; i++) { scene.setModel(pickModel); scene.setView(i%2 ? 'XY' : 'YZ'); await frame(); }
      check(scene.renderer.info.memory.geometries <= baseline && scene.scene.children.length <= 5,'bounded geometries across replacements');
      await frame(); check(scene.frame === null,'no idle animation loop');
      window.sceneUnderTest = scene; window.sceneHostUnderTest = host;
      return results;
    });
    console.log(cameraChecks.join('\n'));
    // Real pointer gestures: orbit does not turn into a click selection.
    const canvas = page.locator('canvas[aria-label="3D geometry"]').last(), rect = await canvas.boundingBox();
    await page.evaluate(()=>sceneUnderTest.setSelection({kind:'element',index:1}));
    await page.mouse.move(rect.x+rect.width/2,rect.y+rect.height/2); await page.mouse.down();
    await page.mouse.move(rect.x+rect.width/2+70,rect.y+rect.height/2+30,{steps:5}); await page.mouse.up();
    assert.deepEqual(await page.evaluate(()=>sceneUnderTest.selection),{kind:'element',index:1});
    const dense = await page.evaluate(async () => {
      const scene = sceneUnderTest, frame = () => new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)));
      const nodes = Array.from({length:100000},(_,i)=>({x:1e9+i%100,y:2e9+Math.floor(i/100)%100,z:3e9+Math.floor(i/10000)}));
      const elements = nodes.map((_,i)=>({node1:i+1,node2:(i+1)%nodes.length+1}));
      const start = performance.now(); scene.setModel({nodes,elements}); await frame();
      scene.setSelection({kind:'node',index:99999}); scene.setView('XZ'); scene.fitModel(); await frame();
      const p = scene.local(nodes[0]).project(scene.camera), b = scene.canvas.getBoundingClientRect(); scene.pick(b.x+(p.x+1)*b.width/2,b.y+(1-p.y)*b.height/2); await frame();
      const result = {elapsedMs:Math.round(performance.now()-start),objects:scene.scene.children.length,labels:scene.labelCount,drawCalls:scene.renderer.info.render.calls,picked:!!scene.selection};
      scene.dispose(); result.disposed = scene.frame === null && sceneHostUnderTest.children.length === 0; sceneHostUnderTest.remove(); return result;
    });
    assert(dense.objects <= 5 && dense.labels <= 200 && dense.drawCalls <= 5 && dense.picked && dense.disposed,JSON.stringify(dense));
    console.log('100000 node/member scene: '+JSON.stringify(dense));
    const rows = ['analysis','static','nodes','100000'];
    for (let i=0;i<100000;i++) rows.push([i%100,Math.floor(i/100)%100,Math.floor(i/10000)].join(','));
    rows.push('elems_1133','100000');
    for (let i=0;i<100000;i++) rows.push('1133,'+(i+1)+','+((i+1)%100000+1)+',1,2e11,7800');
    await page.locator('#file-input').setInputFiles({name:'dense.txt',mimeType:'text/plain',buffer:Buffer.from(rows.join('\n'))});
    await page.waitForFunction(()=>document.querySelector('#case-status').textContent.startsWith('dense'));
    assert.equal(await page.locator('#spatial-node-rows tr').count(),200);
    assert.equal(await page.locator('#spatial-element-rows tr').count(),200);
    await page.locator('.spatial-list').first().getByRole('button',{name:'Next',exact:true}).click();
    assert.equal(await page.locator('#spatial-node-rows tr').first().locator('td').first().innerText(),'201');
    // Repeated switching must leave a single scene, including the planar scene.
    for (let i=0;i<5;i++) {
      await page.locator('#new-case').click(); await page.locator('#element-type').selectOption('1436');
      assert.equal(await page.locator('#spatial-viewport canvas').count(),2);
      await page.locator('[data-view="origin"]').click();
    }
    await page.locator('#new-case').click(); assert.equal(await page.locator('#spatial-viewport canvas').count(),2);
    assert.deepEqual(external,[]);
    // Use the production CSP unchanged, including script-src 'self'.
    const csp = fs.readFileSync(root+'/docker/nginx.conf','utf8').match(/Content-Security-Policy\s+"([^"]+)"/)[1];
    const server = http.createServer((req,res)=>{
      const file = path.resolve(root,'.'+new URL(req.url,'http://localhost').pathname);
      if (!file.startsWith(root+path.sep) || !fs.existsSync(file) || !fs.statSync(file).isFile()) { res.writeHead(404); res.end(); return; }
      res.setHeader('Content-Security-Policy',csp); res.setHeader('Content-Type', {'.html':'text/html','.js':'text/javascript','.css':'text/css'}[path.extname(file)] || 'text/plain'); res.end(fs.readFileSync(file));
    });
    await new Promise(resolve=>server.listen(0,'127.0.0.1',resolve));
    try {
      const base = 'http://127.0.0.1:'+server.address().port;
      await page.goto(base+'/preprocessor/index.html');
      await page.evaluate(()=>{ window.cspViolations=[]; document.addEventListener('securitypolicyviolation',e=>cspViolations.push(e.violatedDirective)); });
      await page.locator('#file-input').setInputFiles(root+'/examples/cases/CaseSpatialFrame.txt');
      await page.waitForFunction(()=>document.querySelector('#case-status').textContent.includes('3D'));
      assert.equal(await page.locator('#spatial-viewport canvas').count(),2); assert.deepEqual(await page.evaluate(()=>cspViolations),[]);
      await page.locator('#spatial-node-rows tr').nth(1).click(); await page.locator('#spatial-node-x').fill('3'); await page.locator('#spatial-node-apply').click();
      assert.match(await page.locator('#spatial-current').inputValue(),/3,0,0/); await page.locator('#undo').click(); assert.deepEqual(await page.evaluate(()=>cspViolations),[]);
      await page.locator('[data-condition-tool="loads"]').click();await page.locator('#condition-node').fill('2');await page.locator('#condition-Fy').fill('10');await page.locator('#condition-apply').click();assert.equal(await page.locator('#condition-loads-rows tr').count(),2);await page.locator('#undo').click();assert.deepEqual(await page.evaluate(()=>cspViolations),[]);
      console.log('HTTP with production CSP and no external requests passed.');
    } finally { await new Promise(resolve=>server.close(resolve)); }
    const unavailable = await browser.newPage();
    await unavailable.addInitScript(()=>{ const original = HTMLCanvasElement.prototype.getContext; HTMLCanvasElement.prototype.getContext = function(type,...args) { return type === 'webgl2' ? null : original.call(this,type,...args); }; });
    await unavailable.goto('file://'+root+'/preprocessor/index.html');
    await unavailable.locator('#file-input').setInputFiles(root+'/examples/cases/CaseSpatialFrame.txt');
    await unavailable.waitForFunction(()=>document.querySelector('#spatial-viewport').textContent.includes('WebGL 2 unavailable'));
    assert(!(await unavailable.locator('#download-case').isDisabled()));
    await unavailable.locator('[data-condition-tool="loads"]').click();await unavailable.locator('#condition-node').fill('2');await unavailable.locator('#condition-Fx').fill('400');await unavailable.locator('#condition-apply').click();
    assert.equal(await unavailable.locator('#condition-loads-rows tr').count(),2);await unavailable.locator('#undo').click();assert.equal(await unavailable.locator('#condition-loads-rows tr').count(),1);await unavailable.keyboard.press('Escape');
    await unavailable.locator('#spatial-node-rows tr').first().click();
    await unavailable.locator('#spatial-node-z').fill('1'); await unavailable.locator('#spatial-node-apply').click();
    assert.match(await unavailable.locator('#case-status').innerText(),/modified/); await unavailable.locator('#undo').click(); assert(!(await unavailable.locator('#redo').isDisabled()));
    await unavailable.locator('#file-input').setInputFiles(root+'/examples/cases/CaseUniformFrame.txt');
    await unavailable.waitForFunction(()=>document.querySelector('#case-status').textContent.startsWith('CaseUniformFrame'));
    assert(await unavailable.locator('#spatial-viewport').isVisible()); assert.equal(await unavailable.locator('#spatial-node-rows tr').count(),2); await unavailable.close();
    console.log('Unavailable WebGL 2 preserves import and allows return to 2D.');
    assert.deepEqual(errors,[]);
    console.log('Spatial browser import, scene, selection and 2D return passed.');
  } finally { await browser.close(); }
})().catch(e=>{console.error(e);process.exit(1);});
