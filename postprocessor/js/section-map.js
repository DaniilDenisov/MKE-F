(function(M,S,P){
  'use strict';
  const V=M.sectionMap={},fmt=v=>Number(v.toPrecision(6)).toString();
  V.frame=function(p,w,h){
    const points=p.loops.flat(),ys=points.map(v=>v.y),zs=points.map(v=>v.z);
    const y0=Math.min(...ys),y1=Math.max(...ys),z0=Math.min(...zs),z1=Math.max(...zs);
    const scale=Math.min((w-70)/(z1-z0),(h-145)/(y1-y0)),cx=w/2-scale*(z0+z1)/2,cy=38+(h-145)/2+scale*(y0+y1)/2;
    const project=(y,z)=>[cx+z*scale,cy-y*scale];
    const path=p.loops.map(loop=>loop.map((v,i)=>(i?'L':'M')+project(v.y,v.z).join(',')).join(' ')+'Z').join(' ');
    return {project,path,scale,cx,cy};
  };
  V.normal=function(svg,data,xi,limit,width){
    const w=Math.max(240,Math.round(width)),h=Math.round(w*.65+155),frame=V.frame(data.profile,w,h),at=data.field.at(xi);
    svg.replaceChildren();svg.setAttribute('width',w);svg.setAttribute('height',h);svg.setAttribute('viewBox',`0 0 ${w} ${h}`);
    const append=(tag,attrs,text,parent=svg)=>{const node=M.svgElement(tag,attrs);if(text!==undefined)node.textContent=text;parent.append(node);return node;};
    append('rect',{width:w,height:h,fill:'white'});
    const a=S.evaluate(data.field.a,xi),b=S.evaluate(data.field.b,xi),c=S.evaluate(data.field.c,xi),norm=b*b+c*c;
    const defs=append('defs',{}),id='section-normal-gradient';let fill='#ffffff';
    if(norm&&limit){
      const start=frame.project((-limit-a)*b/norm,(-limit-a)*c/norm),end=frame.project((limit-a)*b/norm,(limit-a)*c/norm);
      const gradient=append('linearGradient',{id,gradientUnits:'userSpaceOnUse',x1:start[0],y1:start[1],x2:end[0],y2:end[1]},undefined,defs);
      for(const [offset,color]of [[0,'#2166ac'],[.5,'#ffffff'],[1,'#c61b32']])append('stop',{offset,'stop-color':color},undefined,gradient);
      fill=`url(#${id})`;
    }else if(limit){const t=Math.min(1,Math.abs(a)/limit),color=a>=0?[198,27,50]:[33,102,172];fill=`rgb(${color.map(v=>Math.round(255+(v-255)*t)).join(',')})`;}
    append('path',{d:frame.path,'fill-rule':'evenodd',fill,stroke:'#34495e','stroke-width':1});
    for(const [name,color]of [['minimum','#2166ac'],['maximum','#c61b32']]){
      const q=frame.project(at[name].y,at[name].z);append('circle',{cx:q[0],cy:q[1],r:3,fill:color,stroke:'white'});
    }
    const labels=append('g',{class:'section-stress-labels','font-family':'sans-serif','font-size':14,fill:'#17313e'});
    const lines=['y ↑ · z → · viewed from −x',`σmin=${fmt(at.minimum.value)}`,`σmax=${fmt(at.maximum.value)}`,`Scale ±${fmt(limit)} ${data.unit}`, 'Red: tension · blue: compression'];
    lines.forEach((line,i)=>append('text',{x:8,y:i? h-85+(i-1)*22:20},line,labels));
    return {frame,w,h,at};
  };
  V.vm=function(canvas,data,xi,layout,entry,wire){
    const {w,h,frame}=layout;canvas.width=w;canvas.height=h;
    const ctx=canvas.getContext('2d');let high=0,peak=null;
    if(entry?.field){
      for(let i=3;i<entry.field.values.length;i+=4)high=Math.max(high,entry.field.values[i]);
      peak=entry.result.maximum?.mesh_peak||entry.result.maximum;
      entry.field.draw(canvas,3,wire,Number.isFinite(peak?.value)?peak:null,{project:frame.project,limits:[0,high],sequential:true});
    }else{
      const at=data.field.at(xi);high=Math.max(Math.abs(at.minimum.value),Math.abs(at.maximum.value));
      const a=S.evaluate(data.field.a,xi),b=S.evaluate(data.field.b,xi),c=S.evaluate(data.field.c,xi);
      const off=document.createElement('canvas');off.width=w;off.height=h;const other=off.getContext('2d'),image=other.createImageData(w,h);
      // Evaluate |sigma| before colouring every pixel, including the neutral axis.
      for(let y=30;y<h-105;y++)for(let x=0;x<w;x++){
        const value=Math.abs(a+b*(frame.cy-y-.5)/frame.scale+c*(x+.5-frame.cx)/frame.scale),t=high?Math.min(1,value/high):0;
        const color=V.vmColor(t),offset=4*(y*w+x);image.data.set([...color,255],offset);
      }
      other.putImageData(image,0,0);ctx.save();ctx.clip(new Path2D(frame.path),'evenodd');ctx.drawImage(off,0,0);ctx.restore();
    }
    ctx.strokeStyle='#34495e';ctx.lineWidth=1;ctx.stroke(new Path2D(frame.path));
    ctx.fillStyle='#17313e';ctx.font='14px sans-serif';ctx.fillText('y ↑ · z → · viewed from −x',8,20);
    const gradient=ctx.createLinearGradient(12,0,w-12,0);for(const t of [0,.5,1])gradient.addColorStop(t,`rgb(${V.vmColor(t).join(',')})`);
    ctx.fillStyle=gradient;ctx.fillRect(12,h-75,w-24,12);ctx.fillStyle='#17313e';
    ctx.fillText('0',12,h-42);ctx.textAlign='right';ctx.fillText(`${fmt(high)} ${data.unit}`,w-12,h-42);ctx.textAlign='left';ctx.fillText('σVM ≥ 0 · current section scale',8,h-17);
    return high;
  };
  V.vmColor=t=>{const stops=[[241,249,255],[54,144,192],[8,48,107]],k=t<.5?0:1,u=2*t-k;return stops[k].map((v,i)=>Math.round(v+(stops[k+1][i]-v)*u));};
  V.adapt=function(ds,e,nu){
    const spatial=ds.raw.model.dimension===3,truss=[112,1133].includes(e.type),r=ds.resultsByElementId.get(e.id),units=ds.raw.metadata.units||{};
    const section=ds.raw.model.sections?.find(s=>s.id===e.sectionId),normal=M.normalStressResults.get(ds,e);
    if(normal.status!=='ok')throw Error(normal.reason);
    const length=e.length||Math.hypot(...['x','y','z'].map(k=>(ds.nodesById.get(e.nodeIds[1])[k]||0)-(ds.nodesById.get(e.nodeIds[0])[k]||0)));
    const data={element:e,section,normal,length,truss,unit:units.stress||(units.force&&units.length?`${units.force}/${units.length}²`:'model stress units'),lengthUnit:units.length||'model length units',angle:e.sectionAngle||0};
    const n=r.sectionForceCoefficients?.[0]||[-r.localEndForces[0],0,0,0];
    let coeff;
    if(truss)coeff=[n,...Array.from({length:5},()=>[0,0,0,0])];
    else if(spatial)coeff=r.sectionForceCoefficients;
    else {const [axial,v,m]=r.sectionForceCoefficients,{c,s}=S.rotation(e.sectionAngle||0);coeff=[axial,v.map(x=>c*x),v.map(x=>-s*x),[0,0,0,0],m.map(x=>s*x),m.map(x=>c*x)];}
    data.coeff=coeff;data.nu=truss?null:spatial?e.properties.youngsModulus/(2*e.properties.shearModulus)-1:nu;
    if(section){
      data.profile=P.decorate(section,S.resultProperties(section,e));
      data.field=P.stressField(data.profile,coeff);
      // Use accepted axial/bending properties, with the planar solver's convention.
      if(!spatial&&!truss){const p=data.profile,{c,s}=S.rotation(e.sectionAngle||0),m=r.sectionForceCoefficients[2];
        data.field=P.stressField({...p,Iy:p.momentOfInertia,Iz:p.momentOfInertia},[n,coeff[1],coeff[2],coeff[3],m.map(v=>s*v),m.map(v=>c*v)]);
      }
    }
    return data;
  };
})(window.MKEFPost,window.MKEFSections,window.MKEFSpatialSections);
