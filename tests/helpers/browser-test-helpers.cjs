// Capture the public controllers without adding test-only globals to the applications.
exports.captureEditor=async page=>page.evaluate(()=>{const render=MKEFPre.SpatialPanel.prototype.render;MKEFPre.SpatialPanel.prototype.render=function(...args){window.editorPanel=this;return render.apply(this,args);};});
exports.capturePost=async page=>page.evaluate(()=>{const mount=MKEFPost.WebGLRenderer.prototype.mount;MKEFPost.WebGLRenderer.prototype.mount=function(...args){window.resultRenderer=this;return mount.apply(this,args);};});
exports.pickNode=async(page,index)=>{
 await page.evaluate(()=>new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve))));
 const point=await page.evaluate(index=>{const s=editorPanel.scene,p=s.local(editorPanel.model.nodes[index]).project(s.camera),b=s.canvas.getBoundingClientRect();return {x:b.left+(p.x+1)*b.width/2,y:b.top+(1-p.y)*b.height/2};},index);
 await page.mouse.click(point.x,point.y);
 await page.evaluate(()=>new Promise(resolve=>requestAnimationFrame(resolve)));
};
