(function (root) {
  'use strict';
  var S = {};
  S.types = {rectangle:'Rectangle',circle:'Solid circle',square_tube:'Square tube',rectangular_tube:'Rectangular tube',circular_tube:'Circular tube',channel:'Channel',tee:'T section',ibeam:'I section'};
  S.fields = function (type) {
    if (type==='square_tube') return ['b','t','inner_radius'];
    if (type==='rectangular_tube') return ['b','h','t','inner_radius'];
    if (type==='circular_tube') return ['d','t'];
    return type === 'circle' ? ['d'] : type === 'rectangle' ? ['b','h'] : ['b','h','tw','tf'];
  };
  S.rotation = function (angle) {
    var i = [0,90,180,270].indexOf(angle);
    if (i < 0) throw new Error('Section orientation must be 0, 90, 180 or 270 degrees.');
    return {c:[1,0,-1,0][i],s:[0,1,0,-1][i]};
  };
  // Geometry primitives in the input profile axes, without area/inertia integration.
  S.geometry = function (section) {
    if (!section || !Object.prototype.hasOwnProperty.call(S.types,section.type)) throw new Error('Unsupported section type.');
    var d=section.dimensions;
    if (!d || !S.fields(section.type).filter(k=>k!=='inner_radius').every(function (k) { return Number.isFinite(d[k]) && d[k]>0; })) throw new Error('Section dimensions must be positive finite numbers.');
    var ri=d.inner_radius===undefined?0:d.inner_radius;
    if (!Number.isFinite(ri) || ri<0 || ri>0 && !['square_tube','rectangular_tube'].includes(section.type)) throw new Error('Inside corner radius must be finite, nonnegative and used only for rectangular tubes.');
    var rects=[];
    function rect(z,y,b,h) { rects.push({z:z,y:y,b:b,h:h}); }
    if (section.type==='circular_tube') {
      if (2*d.t>=d.d) throw new Error('Section wall thickness leaves no clear space.');
    }
    else if (section.type==='circle') { /* Analytic circular boundary. */ }
    else {
      if (section.type==='rectangle') rect(0,0,d.b,d.h);
      else if (section.type==='square_tube' || section.type==='rectangular_tube') {
        var h=section.type==='square_tube' ? d.b : d.h;
        if (2*d.t>=Math.min(d.b,h)) throw new Error('Section wall thickness leaves no clear space.');
        rect(0,0,d.b,d.t); rect(0,h-d.t,d.b,d.t);
        rect(0,d.t,d.t,h-2*d.t); rect(d.b-d.t,d.t,d.t,h-2*d.t);
      }
      else {
        var doubleFlange=section.type!=='tee';
        if (d.tw>=d.b || (doubleFlange ? 2*d.tf : d.tf)>=d.h) throw new Error('Section wall/flange thickness leaves no clear space.');
        rect(0,d.h-d.tf,d.b,d.tf);
        if (doubleFlange) rect(0,0,d.b,d.tf);
        rect(section.type==='channel' ? 0 : (d.b-d.tw)/2,doubleFlange ? d.tf : 0,d.tw,d.h-(doubleFlange ? 2 : 1)*d.tf);
      }
    }
    var roundedTube=null;
    if(ri>0){
      var height=section.type==='square_tube'?d.b:d.h,ro=ri+d.t;
      if(ro>=Math.min(height,d.b)/2*(1-1e-14))throw new Error('Inside radius + wall thickness must be less than half the smaller outside dimension.');
      roundedTube={h:height,b:d.b,ri:ri,ro:ro,t:d.t};
    }
    return {rectangles:rects,roundedTube,radius:rects.length?0:d.d/2,innerRadius:section.type==='circular_tube'?d.d/2-d.t:0};
  };
  // Contours use an explicit centroid, from the editor or the solver snapshot.
  S.contours = function (geometry, centroid, angle) {
    var r=S.rotation(angle===undefined?0:angle),roundedTube=geometry.roundedTube;
    var rects=geometry.rectangles.map(p=>({...p,y:p.y-centroid.centroidY,z:p.z-centroid.centroidZ}));
    var vertices=[];
    rects.forEach(function (p) { [[p.z,p.y],[p.z+p.b,p.y],[p.z+p.b,p.y+p.h],[p.z,p.y+p.h]].forEach(function (v) { vertices.push({z:v[0],y:v[1]}); }); });
    if(roundedTube)vertices=S.roundedLoop(roundedTube.h,roundedTube.b,roundedTube.ro);
    var eta=vertices.map(function (p) { return p.y*r.c-p.z*r.s; });
    return {...geometry,rectangles:rects,vertices,etaMin:eta.length?Math.min(...eta):-geometry.radius,etaMax:eta.length?Math.max(...eta):geometry.radius};
  };
  // Editor preview only. Result readers call resultProperties, never this calculator.
  S.properties=function(section,angle){
    const geometry=S.geometry(section),d=section.dimensions,r=S.rotation(angle===undefined?0:angle);
    let A=0,cy=0,cz=0,Iy=0,Iz=0;
    if(geometry.radius){
      A=section.type==='circle'?Math.PI*d.d*d.d/4:Math.PI*d.t*(d.d-d.t);
      Iy=Iz=A*(d.d*d.d+4*geometry.innerRadius**2)/16;
    }else{
      for(const p of geometry.rectangles){const a=p.b*p.h;A+=a;cy+=a*(p.y+p.h/2);cz+=a*(p.z+p.b/2);}
      cy/=A;cz/=A;
      for(const p of geometry.rectangles){Iy+=p.h*p.b**3/12+p.b*p.h*(p.z+p.b/2-cz)**2;Iz+=p.b*p.h**3/12+p.b*p.h*(p.y+p.h/2-cy)**2;}
      if(geometry.roundedTube){
        function rounded(h,b,r){const cut=(1-Math.PI/4)*r*r,inertia=(span,width)=>width*span**3/12-4*((span/2-r)**2*cut+(span/2-r)*r**3/3+(1/3-Math.PI/16)*r**4);return [h*b-4*cut,inertia(b,h),inertia(h,b)];}
        const t=geometry.roundedTube,outer=rounded(t.h,t.b,t.ro),inner=rounded(t.h-2*t.t,t.b-2*t.t,t.ri);
        [A,Iy,Iz]=outer.map((v,i)=>v-inner[i]);
      }
    }
    const p={area:A,centroidY:cy,centroidZ:cz,Iy,Iz,momentOfInertia:Iz*r.c*r.c+Iy*r.s*r.s},shape=S.contours(geometry,p,angle);
    if(![A,Iy,Iz,shape.etaMin,shape.etaMax].every(Number.isFinite)||A<=0||Iy<=0||Iz<=0)throw Error('Section dimensions overflow or underflow the numerical range.');
    return {...p,...shape};
  };
  // Result consumers use solver values; only drawing geometry is constructed here.
  S.resultProperties=function(section,element){
    const p=element.sectionProperties,spatial=element.type>=1000;
    if(!p||Array.isArray(p))throw Error('Missing solver sectionProperties. Recalculate the original TXT model.');
    for(const key of ['area','centroidY','centroidZ','Iy','Iz','momentOfInertia','etaMin','etaMax'])if(!Number.isFinite(p[key]))throw Error('Invalid section property '+key);
    if(['area','Iy','Iz','momentOfInertia'].some(k=>p[k]<=0)||p.etaMin>=p.etaMax)throw Error('Invalid positive section properties or fiber bounds.');
    const keys=element.type===1436?['area','Iy','Iz','torsionConstant']:element.type===113?['area','momentOfInertia']:['area'];
    for(const k of keys)if(!Number.isFinite(p[k])||p[k]<=0||!Number.isFinite(element.properties[k])||Math.abs(p[k]-element.properties[k])>1e-10*p[k])throw Error('Solver section property mismatch: '+k);
    if(spatial&&(!Number.isFinite(p.torsionConstant)||p.torsionConstant<=0||typeof p.torsionMethod!=='string'||typeof p.approximate!=='boolean'||typeof p.open!=='boolean'))throw Error('Invalid solver torsion properties.');
    const shape=S.contours(S.geometry(section),p,spatial?0:element.sectionAngle),scale=Math.max(Math.abs(shape.etaMin),Math.abs(shape.etaMax));
    if(![shape.etaMin,shape.etaMax].every(Number.isFinite)||Math.abs(shape.etaMin-p.etaMin)>1e-10*scale||Math.abs(shape.etaMax-p.etaMax)>1e-10*scale)throw Error('Section fiber bounds disagree with drawing geometry.');
    return {...p,...shape,etaMin:p.etaMin,etaMax:p.etaMax};
  };
  // Display contours only; areas, inertias and stress extrema use analytic geometry.
  S.roundedLoop=function(h,b,r){
    var points=[];
    for(var corner=0;corner<4;corner++)for(var i=0;i<=32;i++){
      var angle=(corner+i/32)*Math.PI/2;
      points.push({y:(corner<2?1:-1)*(h/2-r)+r*Math.sin(angle),z:(corner===0||corner===3?1:-1)*(b/2-r)+r*Math.cos(angle)});
    }
    return points;
  };
  S.validate = function (model) {
    var errors=[], ids=new Set();
    (model.sections || []).forEach(function (s) {
      try {
        if (!Number.isInteger(s.id) || s.id<1 || ids.has(s.id)) throw new Error('Section IDs must be unique positive integers.');
        ids.add(s.id);
        if (typeof s.name!=='string' || !s.name.trim()) throw new Error('Section name is required.');
        S.geometry(s);
      } catch (e) { errors.push('Section '+s.id+': '+e.message); }
    });
    model.elements.forEach(function (e,i) {
      if (e.sectionId===undefined || e.sectionId===null || e.sectionId===0) return;
      if (!ids.has(e.sectionId)) errors.push('Element '+(e.id || i+1)+': unknown section.');
      try { if(model.dimension===3){if(!Number.isFinite(e.sectionAngle))throw Error('Section angle must be finite.');}else S.rotation(e.sectionAngle); } catch (err) { errors.push('Element '+(e.id || i+1)+': '+err.message); }
    });
    return errors;
  };
  S.sync = function (model) {
    var errors=S.validate(model); if (errors.length) throw new Error(errors[0]);
    model.elements.forEach(function (e) { if (!e.sectionId) return; var s=model.sections.find(function (s) { return s.id===e.sectionId; }), p=S.properties(s,e.sectionAngle); e.area=p.area; if (model.elementType===113) e.momentOfInertia=p.momentOfInertia; });
  };
  S.evaluate = function (p,x) { return p.reduceRight(function (v,c) { return v*x+c; },0); };
  S.validateStress=function(stress,pointValue,contains){
    if(!stress||!Array.isArray(stress.positions)||stress.positions.length<2)throw Error('Missing static section stress summary. Recalculate the original TXT model.');
    const positions=stress.positions;
    if(positions[0]!==0||positions.at(-1)!==1||positions.some((x,i)=>!Number.isFinite(x)||x<0||x>1||i&&x<=positions[i-1]))throw Error('Invalid stress sampling positions.');
    for(const key of ['minimum','maximum']){
      const p=stress[key];
      if(!p||!Number.isFinite(p.value)||!Number.isFinite(p.xi)||p.xi<0||p.xi>1||!positions.some(x=>Math.abs(x-p.xi)<=1e-10)||!contains(p))throw Error('Invalid section stress critical point.');
      const value=pointValue(p);
      if(!Number.isFinite(value)||Math.abs(value-p.value)>1e-8*Math.max(1,Math.abs(value),Math.abs(p.value)))throw Error('Inconsistent section stress critical value.');
    }
    if(stress.minimum.value>stress.maximum.value)throw Error('Inverted section stress extrema.');
  };
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
  // Positions use normalized element coordinates; preserve exact endpoints
  // when merging roots that coincide with them within numerical tolerance.
  S.unique = function (values) {
    return values.map(function (v) { return Math.abs(v)<=1e-10 ? 0 : Math.abs(v-1)<=1e-10 ? 1 : v; })
      .sort(function (a,b) { return a-b; }).filter(function (v,i,a) { return i===0 || v-a[i-1]>1e-10; });
  };
  if (typeof module!=='undefined' && module.exports) module.exports=S;
  root.MKEFSections=S;
}(typeof window==='undefined' ? globalThis : window));
