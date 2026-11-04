import {Vector3} from 'three';
import {ResultLayer} from './result-layer.js';
// At most 600 prepared poses. No catalog, element operators or DOM dependencies.
export class SectionLayer extends ResultLayer {
 constructor(view){super(view,['sectionContours','sectionAxes']);this.data=null;}
 setData(data){this.data=data;this.update();this.view.requestRender();}
 update(){
  if(!this.data)return;const s=this.view,d=this.data,contours=[],owners=[],colors=[],axes=[],axisColors=[],up=new Vector3(),right=new Vector3();
  s.camera.updateMatrixWorld();right.setFromMatrixColumn(s.camera.matrixWorld,0);up.setFromMatrixColumn(s.camera.matrixWorld,1);const pixel=s.verticalSize()/s.pixels;
  const line=(out,a,b)=>out.push(...a,...b);
  for(const pose of d.poses.slice(0,600)){
   const {center,basis}=pose;
   if(d.contours&&pose.loops){const world=v=>center.map((x,i)=>x+d.scale*(v.y*basis[1][i]+v.z*basis[2][i]));
    for(const loop of pose.loops)for(let j=0;j<loop.length;j++){line(contours,world(loop[j]),world(loop[(j+1)%loop.length]));owners.push(pose.element);const color=pose.selected?[.95,.45,.05]:[.15,.42,.52];colors.push(...color,...color);}}
   if(d.axes){const point=center.map((x,i)=>x+pixel*(right.getComponent(i)*22+up.getComponent(i)*16));line(axes,center,point);axisColors.push(.5,.5,.5,.5,.5,.5);
    for(let k=0;k<3;k++){line(axes,point,point.map((x,i)=>x+basis[k][i]*pixel*18));const c=[[.8,.12,.12],[.12,.55,.18],[.15,.3,.85]][k];axisColors.push(...c,...c);}}
  }
  // Screen-sized triads must not feed the previous zoom back into Fit all.
  super.set({sectionContours:{positions:contours,owners,colors},sectionAxes:{positions:axes,owners:[],colors:axisColors,fitBounds:false}},false);
 }
 clear(){super.clear();this.data=null;}
}
