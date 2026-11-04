// Independent planar stress oracle; never loaded by the application.
module.exports = function (S,p,type,L,f,q) {
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
