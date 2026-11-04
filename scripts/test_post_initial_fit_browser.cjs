const {chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE||'playwright');
const assert=require('node:assert/strict'),path=require('node:path'),fs=require('node:fs');
const root=path.resolve(__dirname,'..');

(async()=>{
  const browser=await chromium.launch({headless:true,args:['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
  try{
    for(const source of ['file','job'])for(const fixture of ['CaseSpatialFrame','CaseSpatialFrameModal','CaseSpatialTiltedTransient']){
      const page=await browser.newPage({viewport:{width:1920,height:940}}),errors=[];
      page.on('pageerror',e=>errors.push(e.message));
      const instrument='MKEFPost.SpatialPost=new Proxy(MKEFPost.SpatialPost,{construct(T,args){window.viewer=Reflect.construct(T,args);return viewer;}});';
      if(source==='job'){
        // Exercise Open Result's hosted loading path, including its status notice.
        await page.route('http://mkef.test/**',route=>{
          const pathname=new URL(route.request().url()).pathname;
          if(pathname.startsWith('/api/v1/jobs/'))return route.fulfill({json:pathname.endsWith('/result')?JSON.parse(fs.readFileSync(path.join(root,'output',fixture+'.json'),'utf8')):{status:'succeeded'}});
          const file=path.join(root,pathname),ext=path.extname(file);
          let body=fs.readFileSync(file);
          if(pathname==='/postprocessor/js/spatial-view.js')body=Buffer.from(body.toString()+'\n'+instrument);
          return route.fulfill({body,contentType:({'.html':'text/html','.js':'text/javascript','.css':'text/css','.svg':'image/svg+xml','.png':'image/png'})[ext]||'application/octet-stream'});
        });
        await page.goto('http://mkef.test/postprocessor/index.html?job=12345678-1234-4234-8234-123456789abc');
        await page.waitForFunction(()=>document.querySelector('#notice-panel').textContent.includes('Results loaded'));
      }else{
        await page.goto('file://'+root+'/postprocessor/index.html');
        await page.evaluate(instrument);
        await page.locator('#file-input').setInputFiles(root+'/output/'+fixture+'.json');
      }
      await page.waitForFunction(()=>window.viewer?.initialFit);
      const settle=()=>page.evaluate(()=>new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve))));
      await settle();
      const state=()=>page.evaluate(()=>{
        const s=viewer.scene,points=s.corners().map(p=>p.project(s.camera));
        return {...s.getViewState(),occupancy:Math.max(
          (Math.max(...points.map(p=>p.x))-Math.min(...points.map(p=>p.x)))/2,
          (Math.max(...points.map(p=>p.y))-Math.min(...points.map(p=>p.y)))/2
        )};
      });
      const initial=await state();
      console.log(source,fixture,'initial height',initial.height,'occupancy',initial.occupancy);
      assert(initial.occupancy>.8&&initial.occupancy<.9,fixture+': visible geometry should fill the initial canvas with padding');
      for(let click=0;click<3;click++){
        await page.locator('[data-result-view="fit"]').click();await settle();
        const fitted=await state();
        assert(Math.abs(fitted.height/initial.height-1)<1e-6,fixture+': Fit all should already be applied on open');
        fitted.target.forEach((v,i)=>assert(Math.abs(v-initial.target[i])<1e-6,fixture+': fit center should be stable'));
      }
      // A user zoom must remain intact during result updates.
      await page.evaluate(()=>{viewer.scene.setHeight(viewer.scene.verticalSize()*1.5);viewer.update();});
      await settle();assert(Math.abs((await state()).height/initial.height-1.5)<1e-6);
      assert.deepEqual(errors,[]);await page.close();
    }
    console.log('Initial spatial result fit and repeated Fit all passed for static, modal and transient results.');
  }finally{await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
