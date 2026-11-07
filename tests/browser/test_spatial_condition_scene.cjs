const {chromium}=require(process.env.MKEF_PLAYWRIGHT_MODULE||'playwright');
const assert=require('node:assert/strict'),path=require('node:path');const root=path.resolve(__dirname,'../..');
(async()=>{const browser=await chromium.launch({headless:true,args:['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});try{
  const page=await browser.newPage({viewport:{width:1450,height:1100}}),errors=[];page.on('pageerror',e=>errors.push(e.message));await page.goto('file://'+root+'/preprocessor/index.html');
  const result=await page.evaluate(async()=>{
    const host=document.createElement('div');host.style.cssText='position:fixed;inset:0;width:1000px;height:800px;z-index:100;background:white';document.body.append(host);
    const wait=()=>new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r))),check=(v,message)=>{if(!v)throw Error(message);};let picked=null,overlap=null;
    const s=new MKEFSpatial.SpatialScene(host,{select:value=>picked=value,choose:hits=>overlap=hits}),F=MKEFPre.spatialFormat,C=MKEFSpatialConditions;
    const m=F.empty(1436);m.nodes=[{x:0,y:0,z:0},{x:2,y:3,z:4}];m.elements=[{node1:1,node2:2,referenceY:[0,1,0],angleDegrees:37}];
    s.setModel(m);const basis=C.basis(m,m.elements[0]);
    const data={supports:[{node:1,mask:[1,1,1,1,1,1]}],loads:[{node:2,type:11,frequency:2,components:[100,-200,300,-40,50,-60]},{node:2,type:13,components:[100,0,0,0,0,0]}],elementLoads:[{elementId:1,type:21,coordinateSystem:1,components:[0,200,0,0,-100,0]}],mpcs:[{depNode:2,depDOF:1,masters:[{node:2,dof:2,coefficient:1},{node:1,dof:1,coefficient:1}]}],releases:[{elementId:1,end:2,component:'Mx'}],bases:[basis],inactiveSupports:[],analysisType:'transient'};
    s.setConditions(data);await wait();const before=JSON.stringify(s.getViewState());
    function midpoint(index){const pos=s.conditions.object.geometry.attributes.position.array,p=s.local({x:0,y:0,z:0});p.set((pos[index*6]+pos[index*6+3])/2,(pos[index*6+1]+pos[index*6+4])/2,(pos[index*6+2]+pos[index*6+5])/2).project(s.camera);const rect=s.canvas.getBoundingClientRect();return {x:rect.left+(p.x+1)*rect.width/2,y:rect.top+(1-p.y)*rect.height/2};}
    for(const kind of ['loads','supports','elementLoads','mpcs','releases']){const index=s.conditions.owners.findIndex(o=>o.kind===kind),p=midpoint(index);check(s.conditions.pick(p.x,p.y).some(o=>o.kind===kind),'pick '+kind);}
    const p=midpoint(s.conditions.owners.findIndex(o=>o.kind==='loads'));s.pick(p.x,p.y);check(picked||overlap,'scene forwards condition selection');
    s.setSelection({kind:'releases',index:0});s.setConditionFocus({elementId:1,end:2});await wait();check(!!s.conditions.endMarker,'selected end marker');check(JSON.stringify(s.getViewState())===before,'condition selection preserves camera');
    s.setConditionLayers({loads:false,supports:false,mpcs:false,releases:false,axes:false});s.setConditionFocus(null);await wait();check(!s.conditions.object,'hidden conditions are not pickable');s.setConditionLayers({loads:true,supports:true,mpcs:true,releases:true});
    const ext=s.renderer.getContext().getExtension('WEBGL_lose_context');if(ext){ext.loseContext();await new Promise(r=>setTimeout(r,80));check(s.lost,'context lost');ext.restoreContext();await new Promise(r=>setTimeout(r,180));check(!s.lost,'context restored');check(s.conditions.data===data,'context preserves condition records');}
    // Repeated updates and layer switches retain a bounded number of GPU objects.
    for(let i=0;i<12;i++){s.setConditions({...data});s.setConditionLayers({loads:!!(i%2)});await wait();check(s.scene.children.length<=9,'bounded objects across condition updates');}
    s.setConditionFocus(null);s.setConditionLayers({loads:true,supports:false,mpcs:false,releases:false,axes:false});
    const nodes=Array.from({length:100000},(_,i)=>({x:i%1000,y:Math.floor(i/1000),z:0})),loads=nodes.map((_,i)=>({node:i+1,type:10,components:[0,0,1,0,0,0]}));
    s.setModel({nodes,elements:[]});const start=performance.now();s.setConditions({loads});await wait();check(s.labelCount<=200&&s.conditions.labels.length<=200,'bounded labels');check(s.scene.children.length<=7,'batched 100000 condition records');
    const q=midpoint(0);check(s.conditions.pick(q.x,q.y).length>0,'pick large condition set');const stats={records:loads.length,elapsedMs:Math.round(performance.now()-start),objects:s.scene.children.length,labels:s.labelCount,drawCalls:s.renderer.info.render.calls};
    await wait();check(s.frame===null,'no permanent animation loop');s.dispose();host.remove();return stats;
  });console.log('Spatial condition scene: picking, local axes, context, layers and large model:',result);assert.deepEqual(errors,[]);
}finally{await browser.close();}})().catch(e=>{console.error(e);process.exit(1);});
