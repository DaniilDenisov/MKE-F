import {Box3,Vector3,Float32BufferAttribute,BufferGeometry,Mesh,MeshBasicMaterial,DoubleSide} from 'three';
// Prepared geometry only. Numerical recovery belongs to the result controller.
export class ResultLayer {
  constructor(view,keys=['original','deformed','diagrams','diagramFill','diagramGuides','sectionAxes','resultNodes','loads','distributed','reactions','spcForces','mpcForces']){this.view=view;this.keys=keys;this.objects={};this.owners={};this.kinds={};this.bounds=new Box3();}
  set(frame,request=true){const s=this.view,point=new Vector3();this.bounds.makeEmpty();
    for(const key of this.keys){
      const data=frame&&frame[key],old=this.objects[key];if(!data){if(old)old.visible=false;continue;}
      let obj=old;const positions=data.positions;if(!obj||obj.geometry.attributes.position.array.length!==positions.length){s.remove(obj);if(data.kind==='triangles'){const geometry=new BufferGeometry();geometry.setAttribute('position',new Float32BufferAttribute(new Float32Array(positions.length),3));obj=new Mesh(geometry,new MeshBasicMaterial({color:data.color||'#9fc4e2',side:DoubleSide,transparent:true,opacity:data.opacity??.25,depthWrite:false}));s.scene.add(obj);}else obj=s.object(new Float32Array(positions.length),key==='resultNodes'||data.kind==='points',data.color||(key==='diagrams'?'#3274ac':key==='sectionAxes'?'#48776c':'#087f9c'),11);this.objects[key]=obj;}
      obj.visible=data.visible!==false;const array=obj.geometry.attributes.position.array;
      for(let i=0;i<positions.length;i+=3){array[i]=positions[i]-s.origin.x;array[i+1]=positions[i+1]-s.origin.y;array[i+2]=positions[i+2]-s.origin.z;if(obj.visible&&data.fitBounds!==false)this.bounds.expandByPoint(point.set(positions[i],positions[i+1],positions[i+2]));}
      obj.geometry.attributes.position.needsUpdate=true;obj.geometry.computeBoundingSphere();this.owners[key]=data.owners;this.kinds[key]=data.ownerKind||'element';
      if(data.colors&&data.colors.length){obj.material.color.set('#ffffff');let attr=obj.geometry.attributes.color;if(!attr||attr.array.length!==data.colors.length){attr=new Float32BufferAttribute(data.colors,3);obj.geometry.setAttribute('color',attr);}else {attr.array.set(data.colors);attr.needsUpdate=true;}if(!obj.material.vertexColors){obj.material.vertexColors=true;obj.material.needsUpdate=true;}}
      else if(obj.material.vertexColors){obj.material.vertexColors=false;obj.material.color.set(data.color||'#087f9c');obj.material.needsUpdate=true;}
    }if(request)s.requestRender();
  }
  pick(x,y){const s=this.view,ray=s.pointerRay(x,y);ray.params.Line.threshold=s.verticalSize()/s.pixels*5;ray.params.Points.threshold=s.verticalSize()/s.pixels*7;
    for(const key of ['resultNodes',...this.keys.filter(k=>!['sectionAxes','resultNodes'].includes(k))]){const obj=this.objects[key];if(!obj||!obj.visible)continue;const hits=ray.intersectObject(obj);if(hits.length){const i=obj.isPoints?hits[0].index:obj.isMesh?hits[0].faceIndex:Math.floor(hits[0].index/2);return {kind:obj.isPoints?'node':this.kinds[key],index:this.owners[key][i]};}}return null;
  }
  clear(){Object.values(this.objects).forEach(o=>this.view.remove(o));this.objects={};this.owners={};this.bounds.makeEmpty();}
}
