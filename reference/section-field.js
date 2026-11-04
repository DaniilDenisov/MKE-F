// Display interpolation only. Python supplies the four stresses at every sample.
window.SectionField = class SectionField {
  constructor(buffer) {
    const [version,nodes,elements,divisions] = new Uint32Array(buffer,0,4);
    const samples=(divisions+1)*(divisions+2)/2;
    if(version!==1 || divisions!==4 || !elements ||
       buffer.byteLength!==16+16*nodes+24*elements+32*samples*elements) throw Error('Некорректный формат поля.');
    this.nodes=new Float64Array(buffer,16,2*nodes);
    this.triangles=new Uint32Array(buffer,16+16*nodes,6*elements);
    this.values=new Float64Array(buffer,16+16*nodes+24*elements);
    this.local=[]; this.faces=[];
    const idx=(i,j)=>i*(divisions+1)-i*(i-1)/2+j;
    for(let i=0;i<=divisions;i++)for(let j=0;j<=divisions-i;j++){
      this.local.push([i/divisions,j/divisions]);
      if(i+j<divisions)this.faces.push([idx(i,j),idx(i+1,j),idx(i,j+1)]);
      if(i+j<divisions-1)this.faces.push([idx(i+1,j),idx(i+1,j+1),idx(i,j+1)]);
    }
  }
  draw(canvas,component,wire,peak){
    const ctx=canvas.getContext('2d'), w=canvas.width,h=canvas.height;
    let radiusY=0,radiusZ=0,lo=Infinity,hi=-Infinity;
    for(let i=0;i<this.nodes.length;i+=2){
      radiusY=Math.max(radiusY,Math.abs(this.nodes[i]));radiusZ=Math.max(radiusZ,Math.abs(this.nodes[i+1]));
    }
    for(let i=component;i<this.values.length;i+=4){lo=Math.min(lo,this.values[i]);hi=Math.max(hi,this.values[i]);}
    const scale=Math.min((w-100)/(2*radiusZ),(h-80)/(2*radiusY));
    const project=(y,z)=>[w/2+z*scale,h/2-y*scale];
    const image=ctx.createImageData(w,h), edges=new Path2D();
    const palette=Array.from({length:256},(_,i)=>{
      const stops=[[32,70,121],[241,239,215],[174,64,47]],k=i<128?0:1,t=i/255*2-k;
      return stops[k].map((c,j)=>Math.round(c+(stops[k+1][j]-c)*t));
    });
    for(let e=0;e<this.triangles.length/6;e++){
      const [a,b,c]=Array.from(this.triangles.subarray(6*e,6*e+3),id=>project(this.nodes[2*id],this.nodes[2*id+1]));
      if(wire){edges.moveTo(...a);edges.lineTo(...b);edges.lineTo(...c);edges.closePath();}
      const points=this.local.map(([u,v],k)=>[a[0]+u*(b[0]-a[0])+v*(c[0]-a[0]),
        a[1]+u*(b[1]-a[1])+v*(c[1]-a[1]),this.values[4*(e*this.local.length+k)+component]]);
      for(const ids of this.faces){
        const [p,q,r]=ids.map(k=>points[k]);
        const ux=q[0]-p[0],uy=q[1]-p[1],vx=r[0]-p[0],vy=r[1]-p[1],det=ux*vy-uy*vx;
        const xmax=Math.floor(Math.max(p[0],q[0],r[0])),ymax=Math.floor(Math.max(p[1],q[1],r[1]));
        for(let y=Math.floor(Math.min(p[1],q[1],r[1]));y<=ymax;y++){
          for(let x=Math.floor(Math.min(p[0],q[0],r[0]));x<=xmax;x++){
            const dx=x+.5-p[0],dy=y+.5-p[1],u=(dx*vy-dy*vx)/det,v=(ux*dy-uy*dx)/det;
            if(u<0||v<0||u+v>1)continue;
            const value=p[2]+u*(q[2]-p[2])+v*(r[2]-p[2]);
            const color=palette[hi===lo?128:Math.min(255,Math.max(0,Math.round(255*(value-lo)/(hi-lo))))];
            const offset=4*(y*w+x);
            image.data[offset]=color[0];image.data[offset+1]=color[1];image.data[offset+2]=color[2];image.data[offset+3]=255;
          }
        }
      }
    }
    ctx.putImageData(image,0,0);
    if(wire){ctx.strokeStyle='rgba(15,40,58,.48)';ctx.lineWidth=.65;ctx.stroke(edges);}
    ctx.strokeStyle='#17313e';ctx.fillStyle='#17313e';ctx.lineWidth=1.5;
    ctx.beginPath();ctx.moveTo(18,h-63);ctx.lineTo(18,h-20);ctx.lineTo(63,h-20);ctx.stroke();
    ctx.font='15px sans-serif';ctx.fillText('+y',10,h-68);ctx.fillText('+z',66,h-15);
    if(peak){
      const [x,y]=project(peak.y,peak.z);
      ctx.beginPath();ctx.arc(x,y,6,0,2*Math.PI);ctx.strokeStyle='white';ctx.lineWidth=4;ctx.stroke();
      ctx.strokeStyle='black';ctx.lineWidth=1.5;ctx.stroke();
    }
    return [lo,hi];
  }
};
