(function (root) {
  'use strict';
  var D=root.MKEFSectionDrawing={}, S=root.MKEFSections;
  function el(tag,attrs,text) { var n=document.createElementNS('http://www.w3.org/2000/svg',tag); Object.keys(attrs || {}).forEach(function (k) { n.setAttribute(k,attrs[k]); }); if(text!==undefined) n.textContent=text; return n; }
  D.preview=function (svg,section,angle,stress,xi,unit) {
    svg.replaceChildren();
    var p=S.properties(section,angle), r=S.rotation(angle), vertices=p.vertices;
    function point(z,y) { return {x:z*r.c+y*r.s,y:-(y*r.c-z*r.s)}; }
    var points=vertices.map(function (v) { return point(v.z,v.y); });
    var radius=section.type==='circle' ? section.dimensions.d/2 : 0;
    var xmin=radius ? -radius : Math.min.apply(null,points.map(function (v) {return v.x;})), xmax=radius ? radius : Math.max.apply(null,points.map(function (v) {return v.x;}));
    var ymin=-p.etaMax,ymax=-p.etaMin, size=Math.max(xmax-xmin,ymax-ymin), pad=size*.32;
    svg.setAttribute('viewBox',[xmin-pad,ymin-pad,xmax-xmin+2*pad,ymax-ymin+2*pad+(stress ? size*.23 : 0)].join(' '));
    svg.setAttribute('role','img'); svg.setAttribute('aria-label',section.name+' · '+angle+'°');
    var fill='#dce8f1';
    if(stress) {
      var values=stress.fibers.map(function (poly) { return S.evaluate(poly,xi); }), max=Math.max(Math.abs(stress.minimum.value),Math.abs(stress.maximum.value));
      var defs=el('defs'), grad=el('linearGradient',{id:svg.id+'-stress',gradientUnits:'userSpaceOnUse',x1:0,x2:0,y1:ymin,y2:ymax});
      for(var i=0;i<=16;i++) { var t=i/16, value=values[1]*(1-t)+values[0]*t; var strength=max ? Math.abs(value)/max : 0; grad.appendChild(el('stop',{offset:t,'stop-color':value>=0 ? 'hsl(4 65% '+(96-50*strength)+'%)' : 'hsl(218 65% '+(96-50*strength)+'%)'})); }
      defs.appendChild(grad); svg.appendChild(defs); fill='url(#'+svg.id+'-stress)';
      svg.appendChild(el('text',{x:xmin,y:ymax+pad*.72,'font-size':size*.065,fill:'#23384d'},'σ: '+Math.min.apply(null,values).toPrecision(5)+' … '+Math.max.apply(null,values).toPrecision(5)+' '+(unit || '')));
      var legend=el('linearGradient',{id:svg.id+'-scale',x1:'0%',x2:'100%',y1:'0%',y2:'0%'});
      [[0,'hsl(218 65% 46%)'],[.5,'hsl(0 0% 96%)'],[1,'hsl(4 65% 46%)']].forEach(function (v) {legend.appendChild(el('stop',{offset:v[0],'stop-color':v[1]}));});defs.appendChild(legend);
      var barY=ymax+pad+size*.015;
      var barLeft=xmin-pad*.65,barRight=xmax+pad*.65;
      svg.appendChild(el('rect',{x:barLeft,y:barY,width:barRight-barLeft,height:size*.05,fill:'url(#'+svg.id+'-scale)'}));
      [[barLeft,-max,'start'],[(xmin+xmax)/2,0,'middle'],[barRight,max,'end']].forEach(function (v) {svg.appendChild(el('text',{x:v[0],y:barY+size*.14,'text-anchor':v[2],'font-size':size*.075,fill:'#23384d'},v[1]===0 ? '0' : v[1].toPrecision(3)));});
    }
    var group=el('g',{fill:fill,stroke:'#2b4861','stroke-width':size*.006}); svg.appendChild(group);
    if(radius) group.appendChild(el('circle',{cx:0,cy:0,r:radius}));
    else p.rectangles.forEach(function (v) { var pts=[[v.z,v.y],[v.z+v.b,v.y],[v.z+v.b,v.y+v.h],[v.z,v.y+v.h]].map(function (a) {var q=point(a[0],a[1]);return q.x+','+q.y;}); group.appendChild(el('polygon',{points:pts.join(' ')})); });
    function axis(z,y,label) { var q=point(z,y), length=Math.hypot(q.x,q.y), ux=q.x/length,uy=q.y/length,h=size*.045;
      svg.appendChild(el('path',{d:'M 0 0 L '+q.x+' '+q.y+' m '+(-ux*h-uy*h/2)+' '+(-uy*h+ux*h/2)+' L '+q.x+' '+q.y+' l '+(-ux*h+uy*h/2)+' '+(-uy*h-ux*h/2),fill:'none',stroke:'#172f44','stroke-width':size*.007}));
      svg.appendChild(el('text',{x:q.x+size*.025,y:q.y-size*.025,'font-size':size*.09,fill:'#172f44'},label));
    }
    axis(0,size*.6,'y'); axis(size*.6,0,'z');
    svg.appendChild(el('circle',{cx:0,cy:0,r:size*.012,fill:'#172f44'}));
  };
  D.axes=function (svg,model,enabled) {
    svg.querySelectorAll('.section-axes').forEach(function (n) {n.remove();});
    if(!enabled || !model || !model.elements.some(function (e) {return !!e.sectionId;})) return;
    var matrix=svg.getScreenCTM(), scale=matrix && Math.hypot(matrix.a,matrix.b); if(!(scale>0)) return;
    var nodes=new Map(model.nodes.map(function (n,i) {return [n.id || i+1,n];}));
    var layer=el('g',{'class':'section-axes','pointer-events':'none','font-family':'sans-serif','font-size':12/scale,fill:'#314c63',stroke:'#314c63','stroke-width':1.3/scale}); svg.insertBefore(layer,svg.querySelector('.loads-layer, [data-layer="loads"]'));
    model.elements.forEach(function (e) {
      if(!e.sectionId) return;
      var a=nodes.get(e.node1 || e.nodeIds[0]),b=nodes.get(e.node2 || e.nodeIds[1]); if(!a || !b) return;
      var L=Math.hypot(b.x-a.x,b.y-a.y); if(!L) return;
      var tx=(b.x-a.x)/L,ty=-(b.y-a.y)/L,nx=ty,ny=-tx,r=S.rotation(e.sectionAngle),x=a.x+(b.x-a.x)/3,y=-a.y-(b.y-a.y)/3;
      // On short members put the triad below the member, opposite its ID.
      if(L*scale<170) { var sign=ny>0 || (ny===0 && nx<0) ? 1 : -1; x+=nx*sign*32/scale; y+=ny*sign*32/scale; }
      function arrow(dx,dy,label,color) {
        var len=23/scale,head=5/scale,xx=x+dx*len,yy=y+dy*len;
        layer.appendChild(el('path',{d:'M '+x+' '+y+' L '+xx+' '+yy+' m '+(-dx*head-dy*head/2)+' '+(-dy*head+dx*head/2)+' L '+xx+' '+yy+' l '+(-dx*head+dy*head/2)+' '+(-dy*head-dx*head/2),fill:'none',stroke:color}));
        layer.appendChild(el('text',{x:xx+dx*8/scale,y:yy+dy*8/scale,'dominant-baseline':'middle','text-anchor':'middle',stroke:'none',fill:color},label));
      }
      arrow(tx,ty,'x','#b53d37');
      arrow(nx*(r.c || -r.s),ny*(r.c || -r.s),r.c ? 'y' : 'z','#287b50');
      layer.appendChild(el('text',{x:x-tx*12/scale,y:y-ty*12/scale,'dominant-baseline':'middle','text-anchor':'middle',stroke:'none',fill:'#355fad'},(r.c ? 'z' : 'y')+((r.c || r.s)>0 ? '⊙' : '⊗')));
    });
  };
}(window));
