const {chromium} = require(process.env.MKEF_PLAYWRIGHT_MODULE || 'playwright');
const path = require('node:path'), assert = require('node:assert/strict');
const root = path.resolve(__dirname,'../..');
(async()=>{
  const browser = await chromium.launch({headless:true,args:['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
  try {
    const page = await browser.newPage({viewport:{width:1400,height:1000},hasTouch:true}), errors=[]; page.on('pageerror',e=>errors.push(e.message));
    await page.goto('file://'+root+'/preprocessor/index.html');
    const report = await page.evaluate(async()=>{
      const host = document.createElement('div'); host.style.cssText='position:fixed;left:0;top:0;width:900px;height:600px;z-index:100;background:white'; document.body.append(host);
      const scene = new MKEFSpatial.SpatialScene(host), frame = ()=>new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)));
      const check=(condition,message)=>{if(!condition) throw new Error(message);};
      const screen = point=>{const p=scene.local(point).project(scene.camera),rect=scene.canvas.getBoundingClientRect();return {x:rect.x+(p.x+1)*rect.width/2,y:rect.y+(1-p.y)*rect.height/2};};
      const close=(a,b)=>Math.abs(a-b)<1e-7;
      scene.setModel({nodes:[],elements:[]});
      for(const plane of ['XY','XZ','YZ']) for(const offset of [-3,3]) {
        scene.setWorkPlane({plane,offset,step:.25,visible:false}); scene.setView('plane'); await frame();
        const axes={XY:['x','y','z'],XZ:['x','z','y'],YZ:['y','z','x']}[plane],point={};point[axes[0]]=1.26;point[axes[1]]=-2.74;point[axes[2]]=offset;
        const p=screen(point),candidate=scene.pointOnPlane(p.x,p.y);
        check(close(candidate.point[axes[0]],1.25)&&close(candidate.point[axes[1]],-2.75)&&candidate.point[axes[2]]===offset,'plane placement '+plane+offset);
        scene.setHeight(1000); await frame(); const far=screen(point),snapped=scene.pointOnPlane(far.x,far.y);
        check(close(snapped.point[axes[0]],1.25)&&close(snapped.point[axes[1]],-2.75)&&scene.displayStep>.25,'base-step snap remains at distant zoom');
        scene.setHeight(20);
      }
      scene.setModel({nodes:[{x:.3,y:.3,z:1},{x:.3,y:.3,z:0}],elements:[]}); scene.resetOrigin(); scene.setHeight(20); scene.setWorkPlane({step:1}); await frame();
      let p=screen({x:.3,y:.3,z:0}),candidate=scene.pointOnPlane(p.x+2,p.y+2);
      check(candidate.nodeIndex===1&&candidate.point.z===0&&candidate.point.x===.3,'node snap wins over grid and ignores overlapping off-plane node');
      scene.setSnap({nodes:false}); candidate=scene.pointOnPlane(p.x,p.y);check(candidate.point.x===0&&candidate.point.y===0,'node snap can be disabled');
      scene.setSnap({grid:false}); candidate=scene.pointOnPlane(p.x,p.y);check(close(candidate.point.x,.3)&&close(candidate.point.y,.3),'exact plane placement without grid');
      scene.setView('XZ');await frame();check(scene.pointOnPlane(p.x,p.y)===null,'parallel ray rejected');
      const model={nodes:[{x:1e9,y:2e9,z:3e9},{x:1e9+10,y:2e9,z:3e9}],elements:[{node1:1,node2:2}]};
      scene.setModel(model);scene.fitModel();await frame();const before=scene.getViewState();
      scene.setModel({nodes:[model.nodes[0],{x:1e9+30,y:2e9+5,z:3e9+7}],elements:model.elements},{preserveView:true});await frame();const after=scene.getViewState();
      check(before.target.every((v,i)=>close(v,after.target[i]))&&before.direction.every((v,i)=>close(v,after.direction[i]))&&before.height===after.height,'GPU recenter preserves world camera');
      // Keep a real scene for touch and pointer-cancellation tests.
      scene.setModel({nodes:[{x:0,y:0,z:2},{x:5,y:0,z:2}],elements:[{node1:1,node2:2}]});scene.resetOrigin();scene.setHeight(20);scene.setTool('move');
      window.touchEdits=0;scene.callbacks.moveNode=()=>touchEdits++;window.interactionScene=scene;window.interactionHost=host;await frame();
      return {node:screen(scene.model.nodes[0]),view:scene.getViewState()};
    });
    const session=await page.context().newCDPSession(page), point=report.node;
    await session.send('Input.dispatchTouchEvent',{type:'touchStart',touchPoints:[{id:1,x:point.x,y:point.y}]});
    await session.send('Input.dispatchTouchEvent',{type:'touchMove',touchPoints:[{id:1,x:point.x+20,y:point.y+10}]});
    await session.send('Input.dispatchTouchEvent',{type:'touchStart',touchPoints:[{id:1,x:point.x+20,y:point.y+10},{id:2,x:point.x+90,y:point.y+30}]});
    await session.send('Input.dispatchTouchEvent',{type:'touchMove',touchPoints:[{id:1,x:point.x-10,y:point.y},{id:2,x:point.x+120,y:point.y+50}]});
    await session.send('Input.dispatchTouchEvent',{type:'touchEnd',touchPoints:[]});
    assert.equal(await page.evaluate(()=>touchEdits),0);assert(await page.evaluate(()=>interactionScene.controls.enabled));
    assert.notEqual(await page.evaluate(()=>interactionScene.getViewState().height),report.view.height,'two touches must navigate after cancelling node drag');
    const node = await page.evaluate(()=>{
      interactionScene.resetOrigin(); interactionScene.setHeight(20); interactionScene.setTool('move'); interactionScene.depth();
      const p=interactionScene.local(interactionScene.model.nodes[0]).project(interactionScene.camera),r=interactionScene.canvas.getBoundingClientRect();return {x:r.x+(p.x+1)*r.width/2,y:r.y+(1-p.y)*r.height/2};
    });
    await page.mouse.move(node.x,node.y); await page.mouse.down(); await page.mouse.move(node.x+40,node.y+20,{steps:3});
    assert(await page.evaluate(()=>!!interactionScene.gesture));
    await page.evaluate(async()=>{
      const s=interactionScene,ext=s.renderer.getContext().getExtension('WEBGL_lose_context');
      await new Promise((resolve,reject)=>{const t=setTimeout(()=>reject(new Error('loss timeout')),5000);s.canvas.addEventListener('webglcontextlost',()=>{clearTimeout(t);resolve();},{once:true});ext.loseContext();});
      if(s.gesture || !s.controls.enabled || touchEdits) throw new Error('Context loss must cancel the gesture without a model edit');
      await new Promise((resolve,reject)=>{const t=setTimeout(()=>reject(new Error('restore timeout')),5000);s.canvas.addEventListener('webglcontextrestored',()=>{clearTimeout(t);resolve();},{once:true});setTimeout(()=>ext.restoreContext(),100);});
    });
    await page.mouse.up(); assert.equal(await page.evaluate(()=>touchEdits),0);
    await page.evaluate(()=>{interactionScene.dispose();interactionHost.remove();});
    assert.deepEqual(errors,[]); console.log('Spatial interactions: all planes, base grid, snap priority, parallel rays, GPU recenter and touch handoff passed.');
  } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exit(1);});
