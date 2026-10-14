(function (root) {
  'use strict';
  var S = {};
  S.types = {rectangle:'Rectangle',circle:'Solid circle',channel:'Channel',tee:'T section',ibeam:'I section'};
  S.fields = function (type) { return type === 'circle' ? ['d'] : type === 'rectangle' ? ['b','h'] : ['b','h','tw','tf']; };
  S.rotation = function (angle) {
    var i = [0,90,180,270].indexOf(angle);
    if (i < 0) throw new Error('Section orientation must be 0, 90, 180 or 270 degrees.');
    return {c:[1,0,-1,0][i],s:[0,1,0,-1][i]};
  };
  S.properties = function (section, angle) {
    if (!section || !Object.prototype.hasOwnProperty.call(S.types,section.type)) throw new Error('Unsupported section type.');
    var d=section.dimensions, r=S.rotation(angle === undefined ? 0 : angle);
    if (!d || !S.fields(section.type).every(function (k) { return Number.isFinite(d[k]) && d[k]>0; })) throw new Error('Section dimensions must be positive finite numbers.');
    var rects=[], A=0, cy=0, cz=0, Iy=0, Iz=0;
    function rect(z,y,b,h) { rects.push({z:z,y:y,b:b,h:h}); }
    if (section.type==='circle') { A=Math.PI*d.d*d.d/4; Iy=Iz=Math.PI*Math.pow(d.d,4)/64; }
    else {
      if (section.type==='rectangle') rect(0,0,d.b,d.h);
      else {
        var doubleFlange=section.type!=='tee';
        if (d.tw>=d.b || (doubleFlange ? 2*d.tf : d.tf)>=d.h) throw new Error('Section wall/flange thickness leaves no clear space.');
        rect(0,d.h-d.tf,d.b,d.tf);
        if (doubleFlange) rect(0,0,d.b,d.tf);
        rect(section.type==='channel' ? 0 : (d.b-d.tw)/2,doubleFlange ? d.tf : 0,d.tw,d.h-(doubleFlange ? 2 : 1)*d.tf);
      }
      rects.forEach(function (p) { var a=p.b*p.h; A+=a; cy+=a*(p.y+p.h/2); cz+=a*(p.z+p.b/2); });
      cy/=A; cz/=A;
      rects.forEach(function (p) { p.y-=cy; p.z-=cz; Iy+=p.h*Math.pow(p.b,3)/12+p.b*p.h*Math.pow(p.z+p.b/2,2); Iz+=p.b*Math.pow(p.h,3)/12+p.b*p.h*Math.pow(p.y+p.h/2,2); });
    }
    var vertices=[];
    rects.forEach(function (p) { [[p.z,p.y],[p.z+p.b,p.y],[p.z+p.b,p.y+p.h],[p.z,p.y+p.h]].forEach(function (v) { vertices.push({z:v[0],y:v[1]}); }); });
    var eta=vertices.map(function (p) { return p.y*r.c-p.z*r.s; });
    var result={area:A,centroidY:cy,centroidZ:cz,Iy:Iy,Iz:Iz,momentOfInertia:Iz*r.c*r.c+Iy*r.s*r.s,rectangles:rects,vertices:vertices,
      etaMin:eta.length ? Math.min.apply(null,eta) : -d.d/2,etaMax:eta.length ? Math.max.apply(null,eta) : d.d/2};
    if (![A,Iy,Iz,result.etaMin,result.etaMax].every(Number.isFinite) || A<=0 || Iy<=0 || Iz<=0) throw new Error('Section dimensions overflow or underflow the numerical range.');
    return result;
  };
  S.validate = function (model) {
    var errors=[], ids=new Set();
    (model.sections || []).forEach(function (s) {
      try {
        if (!Number.isInteger(s.id) || s.id<1 || ids.has(s.id)) throw new Error('Section IDs must be unique positive integers.');
        ids.add(s.id);
        if (typeof s.name!=='string' || !s.name.trim()) throw new Error('Section name is required.');
        S.properties(s);
      } catch (e) { errors.push('Section '+s.id+': '+e.message); }
    });
    model.elements.forEach(function (e,i) {
      if (e.sectionId===undefined || e.sectionId===null || e.sectionId===0) return;
      if (!ids.has(e.sectionId)) errors.push('Element '+(e.id || i+1)+': unknown section.');
      try { S.rotation(e.sectionAngle); } catch (err) { errors.push('Element '+(e.id || i+1)+': '+err.message); }
    });
    return errors;
  };
  S.sync = function (model) {
    var errors=S.validate(model); if (errors.length) throw new Error(errors[0]);
    model.elements.forEach(function (e) { if (!e.sectionId) return; var s=model.sections.find(function (s) { return s.id===e.sectionId; }), p=S.properties(s,e.sectionAngle); e.area=p.area; if (model.elementType===113) e.momentOfInertia=p.momentOfInertia; });
  };
  S.evaluate = function (p,x) { return p.reduceRight(function (v,c) { return v*x+c; },0); };
  // Derivative partitioning isolates all real roots in [0,1], including double roots.
  S.roots = function (input) {
    var scale=Math.max.apply(null,input.map(Math.abs)); if (!scale) return [];
    var p=input.map(function (v) { return v/scale; });
    while (p.length>1 && Math.abs(p[p.length-1])<1e-14) p.pop();
    if (p.length<2) return [];
    var cuts=[0].concat(S.roots(p.slice(1).map(function (v,i) { return v*(i+1); })),[1]), result=[];
    cuts.forEach(function (x) { if (Math.abs(S.evaluate(p,x))<1e-12) result.push(x); });
    for (var i=1;i<cuts.length;i++) {
      var a=cuts[i-1],b=cuts[i],fa=S.evaluate(p,a),fb=S.evaluate(p,b);
      if (fa*fb>=0) continue;
      for (var j=0;j<55;j++) { var mid=(a+b)/2,fm=S.evaluate(p,mid); if (fa*fm<=0) b=mid; else {a=mid;fa=fm;} }
      result.push((a+b)/2);
    }
    return S.unique(result);
  };
  S.unique = function (values) { return values.sort(function (a,b) { return a-b; }).filter(function (v,i,a) { return i===0 || v-a[i-1]>1e-10; }); };
  S.stress = function (p,type,L,f,q) {
    q=q || [0,0,0,0];
    var n=[-f[0],-q[0]*L,-(q[2]-q[0])*L/2,0];
    var m=type===112 ? [0,0,0,0] : [-f[2],f[1]*L,q[1]*L*L/2,(q[3]-q[1])*L*L/6];
    var fibers=[p.etaMin,p.etaMax].map(function (eta) { return n.map(function (v,i) { return v/p.area-m[i]*eta/p.momentOfInertia; }); });
    var positions=[0,1].concat(S.roots(m));
    fibers.forEach(function (poly) { positions=positions.concat(S.roots(poly.slice(1).map(function (v,i) { return v*(i+1); }))); });
    positions=S.unique(positions);
    var minimum={value:Infinity}, maximum={value:-Infinity};
    positions.forEach(function (xi) { fibers.forEach(function (poly,i) { var value=S.evaluate(poly,xi), point={value:value,xi:xi,eta:i ? p.etaMax : p.etaMin}; if(value<minimum.value) minimum=point; if(value>maximum.value) maximum=point; }); });
    return {minimum:minimum,maximum:maximum,positions:positions,axial:n,moment:m,fibers:fibers};
  };
  if (typeof module!=='undefined' && module.exports) module.exports=S;
  root.MKEFSections=S;
}(typeof window==='undefined' ? globalThis : window));
