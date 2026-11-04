(function(root){
  'use strict';
  var S=root.MKEFSections, P={};
  // Calculation/assignment policy only: legacy geometry and results remain readable.
  P.allowedForElement=(code,type)=>code!==1436||['rectangle','circle','square_tube','rectangular_tube','circular_tube'].includes(type);
  P.frameRestriction='Frames 1436 support only rectangles, solid circles and square, rectangular or circular tubes. Open profiles require shear-center offsets and restrained warping that this model does not include.';
  P.manualNote='Without section geometry, profile applicability is not checked.';
  P.calculationError=function(model){
    if(model.elementType!==1436)return '';
    const catalog=new Map((model.sections||[]).map(s=>[s.id,s])),members=[];
    model.elements.forEach((e,i)=>{const s=catalog.get(e.sectionId);if(s&&!P.allowedForElement(1436,s.type))members.push((i+1)+' ('+s.type+')');});
    return members.length?'Unsupported sections on members '+members.join(', ')+'. '+P.frameRestriction:'';
  };
  P.properties=function(s){
    var p=S.properties(s),d=s.dimensions,t=s.torsion||{mode:'auto'},J;
    if(!t||!['auto','manual'].includes(t.mode))throw Error('Invalid torsion mode.');
    p.approximate=false;p.open=['channel','tee','ibeam'].includes(s.type);
    if(t.mode==='manual'){J=t.value;p.torsionMethod='manual';}
    else if(['circle','circular_tube'].includes(s.type)){J=p.Iy+p.Iz;p.torsionMethod='exact circular';}
    else if(s.type==='rectangle'){
      // Saint-Venant odd series, doi:10.1186/s10033-018-0214-9.
      // Bound the omitted odd terms by the larger all-integer tail.
      var a=Math.max(d.b,d.h),b=Math.min(d.b,d.h),f=192*b/(Math.PI**5*a),sum=0,n=1;
      do{sum+=Math.tanh(n*Math.PI*a/(2*b))/n**5;if(f/(4*n**4)<=1e-12*(1-f*sum))break;n+=2;}while(true);
      J=a*b**3/3*(1-f*sum);p.torsionMethod='Saint-Venant series (relative tail ≤ 1e-12)';
    }else{
      // https://docs.oasys-software.com/structural/gsa/references-theory/other-sections/
      p.approximate=true;
      if(['square_tube','rectangular_tube'].includes(s.type)){
        var h=s.type==='square_tube'?d.b:d.h,rm=p.roundedTube?d.inner_radius+d.t/2:0;
        var enclosed=(d.b-d.t)*(h-d.t)-(4-Math.PI)*rm*rm,perimeter=2*(d.b+h-2*d.t)-(8-2*Math.PI)*rm;
        J=4*enclosed**2*d.t/perimeter;p.torsionMethod='thin-wall closed'+(rm?' · rounded median line':'');
      }else{var nf=s.type==='tee'?1:2;J=(nf*d.b*d.tf**3+(d.h-nf*d.tf)*d.tw**3)/3;p.torsionMethod='thin-wall open';}
    }
    if(!Number.isFinite(J)||J<=0)throw Error('J must be positive and finite.');
    p.torsionConstant=J;p.radius=p.rectangles.length?0:d.d/2;p.innerRadius=s.type==='circular_tube'?d.d/2-d.t:0;
    p.loops=P.loops(p);p.extremeVertices=P.hull(p.vertices);return p;
  };
  P.hull=function(vertices){
    const points=vertices.slice().sort((a,b)=>a.z-b.z||a.y-b.y).filter((v,i,a)=>!i||v.z!==a[i-1].z||v.y!==a[i-1].y),cross=(a,b,c)=>(b.z-a.z)*(c.y-a.y)-(b.y-a.y)*(c.z-a.z);
    function chain(rows){const out=[];for(const v of rows){while(out.length>1&&cross(out[out.length-2],out[out.length-1],v)<=0)out.pop();out.push(v);}out.pop();return out;}
    return chain(points).concat(chain(points.slice().reverse()));
  };
  // Boundary of the rectangle union: shared interior edges are never drawn.
  P.loops=function(p){
    if(p.roundedTube){var d=p.roundedTube;return [S.roundedLoop(d.h,d.b,d.ro),S.roundedLoop(d.h-2*d.t,d.b-2*d.t,d.ri)];}
    if(p.radius)return [p.radius].concat(p.innerRadius?[p.innerRadius]:[]).map(function(r){return Array.from({length:128},function(_,i){var a=2*Math.PI*i/128;return {y:r*Math.sin(a),z:r*Math.cos(a)};});});
    var ys=[...new Set(p.rectangles.flatMap(r=>[r.y,r.y+r.h]))].sort((a,b)=>a-b),zs=[...new Set(p.rectangles.flatMap(r=>[r.z,r.z+r.b]))].sort((a,b)=>a-b),edges=new Map();
    function edge(a,b){var k=a.join(',')+'|'+b.join(','),rev=b.join(',')+'|'+a.join(',');if(edges.has(rev))edges.delete(rev);else edges.set(k,[a,b]);}
    for(var i=0;i<ys.length-1;i++)for(var j=0;j<zs.length-1;j++)if(P.contains(p,(ys[i]+ys[i+1])/2,(zs[j]+zs[j+1])/2)){
      var v=[[i,j],[i,j+1],[i+1,j+1],[i+1,j]];for(var k=0;k<4;k++)edge(v[k],v[(k+1)%4]);
    }
    var next=new Map();edges.forEach(e=>next.set(e[0].join(','),e[1]));var loops=[];
    while(next.size){var start=next.keys().next().value,key=start,loop=[];do{var ij=key.split(',').map(Number);loop.push({y:ys[ij[0]],z:zs[ij[1]]});var end=next.get(key);next.delete(key);key=end.join(',');}while(key!==start);loops.push(loop);}
    return loops;
  };
  P.contains=function(p,y,z){
    if(!Number.isFinite(y)||!Number.isFinite(z))return false;
    if(p.roundedTube){var d=p.roundedTube,eps=1e-12*Math.max(d.h,d.b),distance=Math.hypot(Math.max(0,Math.abs(y)-(d.h/2-d.ro)),Math.max(0,Math.abs(z)-(d.b/2-d.ro)));return distance<=d.ro+eps&&distance>=d.ri-eps;}
    if(p.radius){var r=Math.hypot(y,z),eps=1e-12*p.radius;return r<=p.radius+eps&&r>=p.innerRadius-eps;}
    return p.rectangles.some(r=>y>=r.y&&y<=r.y+r.h&&z>=r.z&&z<=r.z+r.b);
  };
  P.catalog=function(model){
    var sections=model.sections||[],map=new Map();if(sections.length>100000)throw Error('Maximum 100000 sections.');
    sections.forEach(s=>{if(!Number.isInteger(s.id)||s.id<1||map.has(s.id)||typeof s.name!=='string'||!s.name.trim())throw Error('Section IDs must be unique and names nonempty.');map.set(s.id,{section:s,properties:P.properties(s)});});return map;
  };
  // Returns fresh derived records; history retains only the catalog delta.
  P.resolve=function(model){
    var map=P.catalog(model);return model.elements.map(function(e){
      if(!e.sectionId)return e;var item=map.get(e.sectionId);if(!item)throw Error('Unknown section '+e.sectionId+'.');
      var p=item.properties,angle=model.elementType===1436?e.angleDegrees:(e.sectionAngle||0);
      if(!Number.isFinite(angle))throw Error('Section angle must be finite.');
      var values=Object.assign({area:p.area},model.elementType===1436?{Iy:p.Iy,Iz:p.Iz,torsionConstant:p.torsionConstant}:{});
      return Object.keys(values).every(k=>values[k]===e[k])?e:Object.assign({},e,values);
    });
  };
  P.stress=function(p,coeff){
    const a=coeff[0].map(v=>v/p.area),b=coeff[5].map(v=>-v/p.Iz),c=coeff[4].map(v=>v/p.Iy),deriv=p=>p.slice(1).map((v,i)=>v*(i+1));
    const mul=(a,b)=>{const r=Array(a.length+b.length-1).fill(0);a.forEach((v,i)=>b.forEach((w,j)=>r[i+j]+=v*w));return r;},add=(a,b)=>Array.from({length:Math.max(a.length,b.length)},(_,i)=>(a[i]||0)+(b[i]||0));
    const point=(xi,y,z)=>({value:S.evaluate(a,xi)+S.evaluate(b,xi)*y+S.evaluate(c,xi)*z,xi,y,z});
    const vertices=p.extremeVertices||p.vertices;
    const at=xi=>{let points=vertices;if(p.radius||p.roundedTube){const bv=S.evaluate(b,xi),cv=S.evaluate(c,xi),r=Math.hypot(bv,cv),radius=p.radius||p.roundedTube.ro;let y=r?radius*bv/r:radius,z=r?radius*cv/r:0;if(p.roundedTube){y+=Math.sign(bv)*(p.roundedTube.h/2-radius);z+=Math.sign(cv)*(p.roundedTube.b/2-radius);if(!r)y=p.roundedTube.h/2;}points=[{y,z},{y:-y,z:-z}];}
      let minimum={value:Infinity},maximum={value:-Infinity};for(const v of points){const q=point(xi,v.y,v.z);if(q.value<minimum.value)minimum=q;if(q.value>maximum.value)maximum=q;}return {minimum,maximum};};
    let positions=[0,1];
    if(p.radius||p.roundedTube){const D=add(mul(b,b),mul(c,c)),dD=deriv(D),radius=p.radius||p.roundedTube.ro;
      // Squaring generates extraneous candidates: evaluate actual branches at every
      // candidate. Extra interior points cannot exceed the true extrema; roots of
      // D also cover cusps and the identically zero bending case uses a'.
      const centers=p.roundedTube?[-1,1].flatMap(y=>[-1,1].map(z=>({y:y*(p.roundedTube.h/2-radius),z:z*(p.roundedTube.b/2-radius)}))):[{y:0,z:0}];
      for(const v of centers){const da=deriv(a.map((x,j)=>x+b[j]*v.y+c[j]*v.z)),poly=add(mul(mul(da,da),D).map(v=>4*v),mul(dD,dD).map(v=>-radius*radius*v));positions.push(...S.roots(poly),...S.roots(da));}
      positions.push(...S.roots(D),...S.roots(b),...S.roots(c));
    }else for(let i=0;i<vertices.length;i++){const v=vertices[i],next=vertices[(i+1)%vertices.length];positions.push(...S.roots(deriv(a.map((x,j)=>x+b[j]*v.y+c[j]*v.z))));
      // Envelope corners: the supporting vertex changes when an edge has
      // constant stress. Include these cusps for exact sigmaMin/Max plots.
      positions.push(...S.roots(b.map((x,j)=>x*(v.y-next.y)+c[j]*(v.z-next.z))));}
    positions=S.unique(positions);let minimum={value:Infinity},maximum={value:-Infinity};for(const x of positions){const r=at(x);if(r.minimum.value<minimum.value)minimum=r.minimum;if(r.maximum.value>maximum.value)maximum=r.maximum;}
    return {point,at,minimum,maximum,positions,a,b,c};
  };
  P.positions=function(elements,selected,count,xi){
    var ids=[];if(selected!==null&&selected>=0&&elements[selected])ids.push(selected);
    for(var i=0;i<elements.length&&ids.length*count<600;i++)if(i!==selected)ids.push(i);
    return ids.flatMap(function(id){var ts=Array.from({length:count},(_,k)=>(k+.5)/count),nearest=-1;if(id===selected&&Number.isFinite(xi)){nearest=Math.min(count-1,Math.max(0,Math.floor(xi*count)));ts[nearest]=xi;}return ts.map((t,k)=>({index:id,xi:t,selected:id===selected&&k===nearest}));}).slice(0,600);
  };
  P.roll=function(R,angle){var a=(angle%360)*Math.PI/180,c=Math.cos(a),s=Math.sin(a);return [R[0],R[1].map((v,i)=>c*v+s*R[2][i]),R[1].map((v,i)=>-s*v+c*R[2][i])];};
  P.trussBasis=function(a,b){var x=b.map((v,i)=>v-a[i]),L=Math.hypot(...x);x=x.map(v=>v/L);var j=[0,1,2].sort((i,j)=>Math.abs(x[i])-Math.abs(x[j]))[0],y=x.map((v,i)=>(i===j?1:0)-v*x[j]),l=Math.hypot(...y);y=y.map(v=>v/l);return [x,y,[x[1]*y[2]-x[2]*y[1],x[2]*y[0]-x[0]*y[2],x[0]*y[1]-x[1]*y[0]]];};
  P.rotate=function(v,r){var a=Math.hypot(...r);if(!a)return v;var k=r.map(x=>x/a),dot=v.reduce((s,x,i)=>s+x*k[i],0),cross=[k[1]*v[2]-k[2]*v[1],k[2]*v[0]-k[0]*v[2],k[0]*v[1]-k[1]*v[0]];return v.map((x,i)=>x*Math.cos(a)+cross[i]*Math.sin(a)+k[i]*dot*(1-Math.cos(a)));};
  P.transport=function(R,direction){var length=Math.hypot(...direction);if(!length)return R;var x=direction.map(v=>v/length),a=R[0],cross=[a[1]*x[2]-a[2]*x[1],a[2]*x[0]-a[0]*x[2],a[0]*x[1]-a[1]*x[0]],s=Math.hypot(...cross),c=a.reduce((sum,v,i)=>sum+v*x[i],0),axis=s?cross.map(v=>v*Math.atan2(s,c)/s):c<0?R[1].map(v=>v*Math.PI):[0,0,0];return R.map(v=>P.rotate(v,axis));};
  root.MKEFSpatialSections=P;
  if(typeof module!=='undefined')module.exports=P;
})(typeof window==='undefined'?globalThis:window);
