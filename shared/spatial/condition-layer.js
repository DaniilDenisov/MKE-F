import {Vector3,Float32BufferAttribute} from 'three';
const colors=['#c33c3c','#26833c','#356dc1'];
const axes=[[1,0,0],[0,1,0],[0,0,1]];
const v=a=>new Vector3(...a);
export class ConditionLayer {
  constructor(scene) { this.view=scene;this.data=null;this.layers={supports:true,loads:true,mpcs:true,releases:true,axes:true};this.labels=[];this.owners=[]; }
  set(data) { this.data=data;this.key=null;this.view.requestRender(); }
  visibility(settings) {Object.assign(this.layers,settings);this.key=null;this.view.requestRender();}
  update() {
    const s=this.view;if(!this.data) return;
    const key=JSON.stringify([s.origin.toArray(),s.camera.position.toArray(),s.camera.quaternion.toArray(),s.verticalSize(),s.width,s.pixels,s.selection,this.layers,this.focus]);if(key===this.key)return;this.key=key;
    s.remove(this.object);s.remove(this.endMarker);this.endMarker=null;this.object=null;this.labels=[];this.owners=[];
    const positions=[],rgb=[],pixel=s.verticalSize()/s.pixels,lanes=new Map();
    const palette=colors.map(c=>parseInt(c.slice(1),16)),tint=number=>[(number>>16&255)/255,(number>>8&255)/255,(number&255)/255];
    const selected=owner=>s.selection&&s.selection.kind===owner.kind&&s.selection.index===owner.index;
    const line=(a,b,color,owner)=>{positions.push(...a.toArray(),...b.toArray());const c=tint(selected(owner)?0xf07816:color);rgb.push(...c,...c);this.owners.push(owner);};
    const local=node=>s.local(node),offscreen=p=>{const q=p.clone().project(s.camera);return Math.abs(q.x)>1.3||Math.abs(q.y)>1.3||Math.abs(q.z)>1;};
    const perpendicular=direction=>{let i=0;if(Math.abs(direction.y)<Math.abs(direction.x))i=1;if(Math.abs(direction.z)<Math.abs(direction.getComponent(i)))i=2;return new Vector3(...axes[i]).cross(direction).normalize();};
    const arrow=(tip,dir,length,color,owner)=>{
      const tail=tip.clone().addScaledVector(dir,-length),wing=perpendicular(dir).multiplyScalar(pixel*4),base=tip.clone().addScaledVector(dir,-Math.min(pixel*9,length*.4));
      line(tail,tip,color,owner);line(tip,base.clone().add(wing),color,owner);line(tip,base.clone().sub(wing),color,owner);
    };
    const arc=(center,dir,sign,bar,color,owner)=>{
      const a=perpendicular(dir),b=dir.clone().cross(a),radius=pixel*13;let last=null;
      for(let i=0;i<=10;i++){const t=sign*i*Math.PI*1.65/10,p=center.clone().addScaledVector(a,Math.cos(t)*radius).addScaledVector(b,Math.sin(t)*radius);if(last)line(last,p,color,owner);last=p;}
      const t=sign*Math.PI*1.65,tangent=a.clone().multiplyScalar(-Math.sin(t)*sign).addScaledVector(b,Math.cos(t)*sign);
      if(bar)line(last.clone().addScaledVector(dir,-pixel*5),last.clone().addScaledVector(dir,pixel*5),color,owner);else arrow(last,tangent,pixel*8,color,owner);
    };
    const baseFor=(node,owner)=>{
      const origin=local(node);if(offscreen(origin))return null;
      const key=owner.kind+':'+[node.x,node.y,node.z].join(','),lane=lanes.get(key)||0;lanes.set(key,lane+1);
      const p=origin.clone().add(new Vector3((lane%5)*12,Math.floor(lane%25/5)*12,0).applyQuaternion(s.camera.quaternion).multiplyScalar(pixel));
      if(lane)line(origin,p,0xa9b3bf,owner);return p;
    };
    const badge=(p,text,owner)=>{if(this.labels.length<199||selected(owner))this.labels.push({point:p,text,owner});};
    const inactive=new Set((this.data.inactiveSupports||[]).map(r=>r.index+':'+r.component));
    if(this.layers.supports)(this.data.supports||[]).forEach((r,index)=>{
      const owner={kind:'supports',index},p=baseFor(s.model.nodes[r.node-1],owner);if(!p)return;
      r.mask.forEach((fixed,j)=>{if(!fixed)return;const dir=v(axes[j%3]),color=inactive.has(index+':'+j)?0x939ba4:palette[j%3];
        if(j<3){const end=p.clone().addScaledVector(dir,-pixel*21),side=perpendicular(dir).multiplyScalar(pixel*6);line(p,end,color,owner);line(end.clone().add(side),end.clone().sub(side),color,owner);}
        else arc(p,dir,1,true,color,owner);
      });badge(p,'S'+(index+1),owner);
    });
    (this.data.loads||[]).forEach((r,index)=>{
      if(!this.layers[r.kind||'loads'])return;
      const owner={kind:r.kind||'loads',index:r.ownerIndex===undefined?index:r.ownerIndex},p=baseFor(s.model.nodes[r.node-1],owner);if(!p)return;
      r.components.forEach((value,j)=>{if(!value)return;const dir=v(axes[j%3]),color=r.kind==='supportReactions'?0x7b3fa1:r.kind==='mpcForces'?0xb66a13:palette[j%3];if(j<3)arrow(p,dir.multiplyScalar(Math.sign(value)),pixel*36,color,owner);else arc(p,dir,Math.sign(value),false,color,owner);});
      if(!r.components.some(Boolean)){line(p.clone().addScalar(-pixel*3),p.clone().addScalar(pixel*3),0x787f88,owner);}
      badge(p,r.kind?(r.kind==='supportReactions'?'Rsupport · N':'RMPC · N')+r.node:'F'+(index+1)+' · '+(r.type===10?'static':r.type===11?'sin '+r.frequency+' Hz':r.type===12?'pulse':'step'),owner);
    });
    if(this.layers.loads)(this.data.elementLoads||[]).forEach((r,index)=>{
      const owner={kind:'elementLoads',index},e=s.model.elements[r.elementId-1],a=local(s.model.nodes[e.node1-1]),b=local(s.model.nodes[e.node2-1]);
      const basis=r.coordinateSystem===1?this.data.bases[r.elementId-1]:axes,ends=r.type===20?r.components.concat(r.components):r.components;
      const maximum=Math.max(...ends.map(Math.abs));if(!maximum)return;
      for(let i=0;i<=6;i++){const t=i/6,p=a.clone().lerp(b,t);if(offscreen(p))continue;
        for(let j=0;j<3;j++){const value=(1-t)*ends[j]+t*ends[j+3];if(!value)continue;arrow(p,v(basis[j]).multiplyScalar(Math.sign(value)),pixel*38*Math.abs(value)/maximum,palette[j],owner);}
        if(i===3)badge(p,'q'+(index+1)+(r.coordinateSystem===1?' local':' global'),owner);
      }
    });
    if(this.layers.mpcs)(this.data.mpcs||[]).forEach((r,index)=>{
      const owner={kind:'mpcs',index},a=local(s.model.nodes[r.depNode-1]),seen=new Set();
      r.masters.forEach(master=>{if(seen.has(master.node))return;seen.add(master.node);const b=local(s.model.nodes[master.node-1]);
        if(master.node===r.depNode){const dx=new Vector3(pixel*10,0,0).applyQuaternion(s.camera.quaternion),dy=new Vector3(0,pixel*10,0).applyQuaternion(s.camera.quaternion);[dx,dy,dx.clone().negate(),dy.clone().negate()].forEach((v,i,vs)=>line(a.clone().add(v),a.clone().add(vs[(i+1)%4]),0x8456ad,owner));}
        else for(let j=0;j<12;j++)line(a.clone().lerp(b,j/12),a.clone().lerp(b,(j+.55)/12),0x8456ad,owner);
      });if(!offscreen(a))badge(a,'MPC '+(index+1),owner);
    });
    if(this.layers.releases)(this.data.releases||[]).forEach((r,index)=>{
      const owner={kind:'releases',index},e=s.model.elements[r.elementId-1],R=this.data.bases[r.elementId-1],j=['Mx','My','Mz'].indexOf(r.component);
      const p=local(s.model.nodes[(r.end===1?e.node1:e.node2)-1]).addScaledVector(v(R[0]),pixel*(r.end===1?14:-14));if(offscreen(p))return;
      arc(p,v(R[j]),1,false,palette[j],owner);badge(p,'E'+r.elementId+'/'+r.end+' '+r.component,owner);
    });
    s.model.elements.forEach((e,index)=>{
      const focus=this.focus&&this.focus.elementId===index+1?this.focus:null;if(!this.layers.axes&&!focus)return;
      const R=this.data.bases[index];if(!R)return;const a=local(s.model.nodes[e.node1-1]),b=local(s.model.nodes[e.node2-1]),p=focus&&focus.end?(focus.end===1?a:b):a.clone().lerp(b,.5);if(offscreen(p))return;
      const owner={kind:'element',index};for(let j=0;j<3;j++)arrow(p.clone().addScaledVector(v(R[j]),pixel*28),v(R[j]),pixel*28,palette[j],owner);
      badge(p,'E'+(index+1)+' · local XYZ',owner);
      if(focus&&focus.end){this.endMarker=s.object(p.toArray(),true,'#e46816',12);this.endMarker.material.depthTest=false;this.endMarker.renderOrder=6;}
    });
    if(positions.length){this.object=s.object(positions,false,'#ffffff');this.object.geometry.setAttribute('color',new Float32BufferAttribute(rgb,3));this.object.material.vertexColors=true;this.object.renderOrder=5;}
  }
  pick(x,y) {
    if(!this.object)return[];const s=this.view,ray=s.pointerRay(x,y);ray.params.Line.threshold=s.verticalSize()/s.pixels*6;
    const hits=ray.intersectObject(this.object),found=new Map();hits.forEach(hit=>{const owner=this.owners[Math.floor(hit.index/2)];if(owner)found.set(owner.kind+':'+owner.index,owner);});return Array.from(found.values());
  }
  drawLabels(ctx,budget) {
    const s=this.view,items=this.labels.slice().sort((a,b)=>Number(!!s.selection&&b.owner.kind===s.selection.kind&&b.owner.index===s.selection.index)-Number(!!s.selection&&a.owner.kind===s.selection.kind&&a.owner.index===s.selection.index));let count=0;
    for(const item of items){if(count>=budget)break;const p=item.point.clone().project(s.camera);if(Math.abs(p.x)>1||Math.abs(p.y)>1||Math.abs(p.z)>1)continue;ctx.fillStyle='#526076';ctx.fillText(item.text,(p.x+1)*s.width/2+7,(1-p.y)*s.pixels/2+16);count++;}return count;
  }
  dispose(){this.view.remove(this.object);this.view.remove(this.endMarker);this.endMarker=null;this.object=null;this.data=null;}
}
