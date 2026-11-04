import {
  WebGLRenderer, Scene, Color, OrthographicCamera, Vector2, Vector3, Box3,
  BufferGeometry, Float32BufferAttribute, LineSegments, LineBasicMaterial,
  Points, PointsMaterial, Raycaster, MOUSE, TOUCH
} from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import {ConditionLayer} from './condition-layer.js';
import {ResultLayer} from './result-layer.js';
import {SectionLayer} from './section-layer.js';

// All GPU positions are relative to origin. Public coordinates and view targets remain XYZ.
// No FEM matrices, editing state, network resources or DOM tables belong in this module.
export class SpatialScene {
  constructor(host, callbacks = {}, options = {}) {
    this.host = host; this.callbacks = callbacks; this.disposed = false; this.frame = null;
    this.origin = new Vector3(); this.bounds = new Box3(); this.model = {nodes:[],elements:[]};
    this.selection = null; this.height = 20; this.width = 1; this.pixels = 1;
    this.dimension = options.dimension === 2 ? 2 : 3;
    this.workPlane = {plane:'XY',offset:0,step:1,visible:true};
    this.tool = 'select'; this.snap = {grid:true,nodes:true};
    this.canvas = document.createElement('canvas'); this.canvas.setAttribute('aria-label',this.dimension+'D geometry');
    this.renderer = new WebGLRenderer({canvas:this.canvas,antialias:true,alpha:false});
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1,2));
    this.scene = new Scene(); this.scene.background = new Color('#f4f7fb');
    this.conditions = new ConditionLayer(this);
    this.results = new ResultLayer(this);
    this.sections = new SectionLayer(this);
    this.camera = new OrthographicCamera(-10,10,10,-10,0.01,1000);
    this.camera.up.set(0,0,1); this.camera.position.set(20,20,20);
    this.onChange = () => this.requestRender(); this.makeControls();
    this.overlay = document.createElement('canvas'); this.overlay.className = 'spatial-overlay';
    this.overlay.style.cssText = 'position:absolute;inset:0;pointer-events:none';
    this.message = document.createElement('div'); this.message.className = 'spatial-context-message'; this.message.hidden = true;
    host.append(this.canvas,this.overlay,this.message);
    this.triad=document.createElement('div');this.triad.className='orientation-triad';this.triad.setAttribute('aria-label','Camera orientation');
    this.triad.style.cssText='position:absolute;right:0;top:0;width:112px;height:112px;pointer-events:none';
    // Canvas owns the visible labels (including exports); buttons are accessible
    // hit targets only, so disabled 2D controls cannot duplicate the lettering.
    this.triadButtons=['X','Y','Z'].map((axis,i)=>{const b=document.createElement('button');b.type='button';b.title='View '+['YZ','XZ','XY'][i];b.setAttribute('aria-label',b.title);b.dataset.axis=axis;b.style.cssText='position:absolute;width:28px;height:28px;padding:0;border:0;border-radius:50%;background:transparent;pointer-events:auto';b.onclick=e=>{e.stopPropagation();this.setView(['YZ','XZ','XY'][i]);};b.disabled=this.dimension===2;this.triad.append(b);return b;});host.append(this.triad);
    this.onLost = event => { event.preventDefault(); this.cancelGesture(); this.lost = true; this.message.hidden = false; this.message.textContent = 'WebGL context lost. Waiting for restoration…'; if(this.callbacks.contextLost)this.callbacks.contextLost(); };
    this.onRestored = () => { this.lost = false; this.message.hidden = true; this.requestRender(); };
    this.canvas.addEventListener('webglcontextlost',this.onLost); this.canvas.addEventListener('webglcontextrestored',this.onRestored);
    this.pointers = new Set();
    this.onDown = e => {
      if (this.handoff) return;
      this.pointers.add(e.pointerId);
      if (this.pointers.size === 1) this.pointer = {x:e.clientX,y:e.clientY,id:e.pointerId,drag:false,button:e.button};
      else if (this.pointer) this.pointer.drag = true;
      if (this.pointers.size > 1 && this.gesture) {
        const first = this.gesture.first; this.cancelGesture(); this.handoff = true;
        // Transfer the active first touch to OrbitControls, then let the real second touch through.
        this.canvas.dispatchEvent(new PointerEvent('pointerdown',{...first,bubbles:true})); this.handoff = false; return;
      }
      if (this.tool === 'move' && this.pointers.size === 1 && e.button === 0 && !e.shiftKey && !e.ctrlKey && !e.metaKey) {
        const hit = this.hit(e.clientX,e.clientY);
        if (hit && hit.kind === 'node') {
          const node = this.model.nodes[hit.index], n = this.planeAxes()[2];
          this.gesture = {index:hit.index,offset:[node.x,node.y,node.z][n],candidate:{point:{...node},snap:'Node '+(hit.index+1)},first:{pointerId:e.pointerId,pointerType:e.pointerType,clientX:e.clientX,clientY:e.clientY,button:0,buttons:1}};
          this.controls.enabled = false; this.canvas.setPointerCapture(e.pointerId); e.stopPropagation(); e.preventDefault();
        }
      }
    };
    this.onMove = e => {
      if (this.pointer && Math.hypot(e.clientX-this.pointer.x,e.clientY-this.pointer.y)>4) this.pointer.drag = true;
      if (this.gesture && this.pointer && this.pointer.id === e.pointerId) {
        this.gesture.first.clientX = e.clientX; this.gesture.first.clientY = e.clientY;
        if (this.pointer.drag) { this.gesture.candidate = this.pointOnPlane(e.clientX,e.clientY,this.gesture.offset,this.gesture.index); this.showMoveCandidate(); }
        e.stopPropagation(); return;
      }
      if (this.tool === 'node' && (!this.pointer || !this.pointer.drag)) this.showCandidate(this.pointOnPlane(e.clientX,e.clientY));
      if (this.tool === 'member' && this.pendingMember !== null && this.pendingMember !== undefined && (!this.pointer || !this.pointer.drag)) this.showMemberPreview(e.clientX,e.clientY);
    };
    this.onUp = e => {
      this.pointers.delete(e.pointerId);
      if (this.gesture && this.pointer && this.pointer.id === e.pointerId) {
        const gesture = this.gesture, moved = this.pointer.drag;
        this.cancelGesture();
        if (moved && gesture.candidate) {
          const error = this.callbacks.validateNode ? this.callbacks.validateNode(gesture.index,gesture.candidate.point) : '';
          if (error) { if (this.callbacks.preview) this.callbacks.preview(error); }
          else if (this.callbacks.moveNode) this.callbacks.moveNode(gesture.index,gesture.candidate.point);
        } else if (!moved) this.pick(e.clientX,e.clientY);
        e.stopPropagation(); return;
      }
      if (this.pointer && this.pointer.id === e.pointerId) {
        if (!this.pointer.drag && this.pointer.button === 0) {
          if (this.tool === 'node') {
            const candidate = this.pointOnPlane(e.clientX,e.clientY);
            if (candidate && this.callbacks.createNode) this.callbacks.createNode(candidate.point);
          } else this.pick(e.clientX,e.clientY);
        }
        this.pointer = null;
      }
    };
    this.onCancel = e => { this.pointers.delete(e.pointerId); this.cancelGesture(); };
    this.canvas.addEventListener('pointerdown',this.onDown,true);
    this.canvas.addEventListener('pointermove',this.onMove,true);
    this.canvas.addEventListener('pointerup',this.onUp,true);
    this.canvas.addEventListener('pointercancel',this.onCancel);
    this.observer = new ResizeObserver(() => this.resize()); this.observer.observe(host);
    this.resize();
  }
  object(positions, points, color, size = 11) {
    const geometry = new BufferGeometry(); geometry.setAttribute('position',new Float32BufferAttribute(positions,3));
    const material = points ? new PointsMaterial({color,size,sizeAttenuation:false}) : new LineBasicMaterial({color});
    const object = points ? new Points(geometry,material) : new LineSegments(geometry,material);
    this.scene.add(object); return object;
  }
  makeControls() {
    const target = this.controls ? this.controls.target.clone() : new Vector3();
    if (this.controls) { this.controls.removeEventListener('change',this.onChange); this.controls.dispose(); }
    // OrbitControls caches the up-axis transform at construction. Recreate on up-axis changes.
    this.controls = new OrbitControls(this.camera,this.canvas); this.controls.target.copy(target); this.controls.cursor.copy(target);
    this.controls.enableDamping = false; this.controls.autoRotate = false;
    this.controls.enableRotate = this.dimension !== 2;
    if(this.dimension===2){this.controls.mouseButtons.LEFT=MOUSE.PAN;this.controls.touches.ONE=TOUCH.PAN;}
    this.controls.minZoom = 1e-12; this.controls.maxZoom = 1e12; this.controls.screenSpacePanning = true;
    this.controls.addEventListener('change',this.onChange);
  }
  remove(object) { if (object) { this.scene.remove(object); object.geometry.dispose(); object.material.dispose(); } }
  local(node) { return new Vector3(node.x,node.y,node.z).sub(this.origin); }
  setModel(model, {preserveView = false} = {}) {
    this.results.clear();
    this.sections.clear();
    const previous = this.model, oldOrigin = this.origin.clone();
    this.model = model; this.selection = null;
    this.bounds.makeEmpty();
    for (const n of model.nodes) this.bounds.expandByPoint(new Vector3(n.x,n.y,n.z));
    this.origin.copy(this.bounds.isEmpty() ? new Vector3() : this.bounds.getCenter(new Vector3()));
    this.remove(this.highlight); this.highlight = null; this.clearPreview();
    const shifted = !this.origin.equals(oldOrigin), topology = !this.nodes || model.nodes.length !== previous.nodes.length || model.elements.length !== previous.elements.length;
    if (topology) {
      this.remove(this.nodes); this.remove(this.members);
      this.nodes = this.object(new Float32Array(model.nodes.length*3),true,'#154f79');
      this.members = this.object(new Float32Array(model.elements.length*6),false,'#276e9d');
    }
    const nodes = this.nodes.geometry.attributes.position, members = this.members.geometry.attributes.position;
    model.nodes.forEach((n,i) => { if (shifted || topology || n !== previous.nodes[i]) this.local(n).toArray(nodes.array,i*3); });
    model.elements.forEach((e,i) => {
      if (shifted || topology || e !== previous.elements[i] || model.nodes[e.node1-1] !== previous.nodes[e.node1-1] || model.nodes[e.node2-1] !== previous.nodes[e.node2-1]) {
        this.local(model.nodes[e.node1-1]).toArray(members.array,i*6); this.local(model.nodes[e.node2-1]).toArray(members.array,i*6+3);
      }
    });
    nodes.needsUpdate = true; members.needsUpdate = true; this.nodes.geometry.computeBoundingSphere(); this.members.geometry.computeBoundingSphere();
    if (preserveView) {
      const shift = oldOrigin.sub(this.origin); this.camera.position.add(shift); this.controls.target.add(shift); this.controls.cursor.add(shift); this.depth(); this.requestRender();
    } else { this.controls.target.set(0,0,0); this.setView(this.dimension===2?'XY':'iso'); this.fitModel(); }
  }
  setTool(tool) { this.tool = tool; this.cancelGesture(); }
  setConditions(data) { this.conditions.set(data); }
  setResultFrame(frame) { this.resultLabels=frame&&frame.labels||[];this.results.set(frame); }
  setSectionFrame(frame) { this.sections.setData(frame); }
  visibleBounds() { return this.bounds.clone().union(this.results.bounds).union(this.sections.bounds); }
  setConditionLayers(settings) { this.conditions.visibility(settings); }
  setConditionFocus(focus) { this.conditions.focus=focus;this.requestRender(); }
  setSnap(settings) { this.snap = {...this.snap,...settings}; this.clearPreview(); }
  cancelGesture() {
    const id = this.pointer && this.pointer.id; this.gesture = null; this.pointer = null; this.controls.enabled = true;
    this.pendingMember = null;
    if (id !== null && id !== undefined && this.canvas.hasPointerCapture(id)) this.canvas.releasePointerCapture(id);
    this.clearPreview(); if (this.callbacks.cancel) this.callbacks.cancel();
  }
  setPendingMember(index) { this.pendingMember = index; this.clearPreview(); }
  showMemberPreview(x,y) {
    this.clearPreview(); const node = this.model.nodes[this.pendingMember]; if (!node) return;
    const hit = this.hit(x,y), target = hit && hit.kind === 'node' ? this.model.nodes[hit.index] : this.pointOnPlane(x,y)?.point;
    if (target) { this.previewLines = this.object([...this.local(node).toArray(),...this.local(target).toArray()],false,'#008b67'); this.previewLines.material.depthTest = false; this.previewLines.renderOrder = 19; }
    if (this.callbacks.preview) this.callbacks.preview(hit && hit.kind === 'node' ? 'Second node: '+(hit.index+1) : 'Select the second existing node. This tool does not create nodes.');
  }
  showMoveCandidate() {
    const gesture = this.gesture, candidate = gesture.candidate; this.showCandidate(candidate); if (!candidate) return;
    const error = this.callbacks.validateNode ? this.callbacks.validateNode(gesture.index,candidate.point) : '';
    if (error) { this.preview.material.color.set('#c32332'); if (this.callbacks.preview) this.callbacks.preview(error); }
    const positions = [];
    this.model.elements.forEach(e => {
      if (e.node1 !== gesture.index+1 && e.node2 !== gesture.index+1) return;
      positions.push(...this.local(e.node1 === gesture.index+1 ? candidate.point : this.model.nodes[e.node1-1]).toArray(),...this.local(e.node2 === gesture.index+1 ? candidate.point : this.model.nodes[e.node2-1]).toArray());
    });
    this.previewLines = this.object(positions,false,error ? '#c32332' : '#008b67'); this.previewLines.material.depthTest = false; this.previewLines.renderOrder = 19;
  }
  clearPreview() { this.remove(this.preview); this.preview = null; this.remove(this.previewLines); this.previewLines = null; this.requestRender(); }
  pointerRay(clientX,clientY) {
    const rect = this.canvas.getBoundingClientRect(), ray = new Raycaster();
    ray.setFromCamera(new Vector2((clientX-rect.left)/rect.width*2-1,1-(clientY-rect.top)/rect.height*2),this.camera); return ray;
  }
  pointOnPlane(clientX,clientY,offset = this.workPlane.offset,exclude = -1) {
    const [u,v,n] = this.planeAxes(), ray = this.pointerRay(clientX,clientY).ray, direction = ray.direction.getComponent(n);
    if (Math.abs(direction) < 1e-8) return null;
    const distance = (offset-this.origin.getComponent(n)-ray.origin.getComponent(n))/direction; if (distance < 0) return null;
    const p = ray.at(distance,new Vector3()).add(this.origin), values = p.toArray(), rect = this.canvas.getBoundingClientRect();
    let best = null, bestDistance = 8, bestDepth = Infinity;
    if (this.snap.nodes) this.model.nodes.forEach((node,i) => {
      if (i === exclude) return;
      const xyz = [node.x,node.y,node.z], tolerance = 8*Number.EPSILON*Math.max(1,Math.abs(offset),Math.abs(xyz[n]));
      if (Math.abs(xyz[n]-offset) > tolerance) return;
      const screen = this.local(node).project(this.camera); if (Math.abs(screen.z)>1) return;
      const distance = Math.hypot(rect.left+(screen.x+1)*rect.width/2-clientX,rect.top+(1-screen.y)*rect.height/2-clientY);
      if (distance < bestDistance || distance === bestDistance && screen.z < bestDepth) { best = i; bestDistance = distance; bestDepth = screen.z; }
    });
    if (best !== null) return {point:{...this.model.nodes[best]},snap:'Node '+(best+1),nodeIndex:best};
    if (this.snap.grid) { values[u] = Math.round(values[u]/this.workPlane.step)*this.workPlane.step; values[v] = Math.round(values[v]/this.workPlane.step)*this.workPlane.step; }
    values[n] = offset;
    if (!values.every(Number.isFinite)) return null;
    return {point:{x:values[0],y:values[1],z:values[2]},snap:this.snap.grid ? 'Grid (base step '+this.workPlane.step+')' : 'Working plane'};
  }
  showCandidate(candidate) {
    this.clearPreview();
    if (candidate) { this.preview = this.object(this.local(candidate.point).toArray(),true,'#008b67',10); this.preview.material.depthTest = false; this.preview.renderOrder = 20; }
    if (this.callbacks.preview) this.callbacks.preview(candidate ? candidate.snap+' · XYZ '+[candidate.point.x,candidate.point.y,candidate.point.z].join(', ') : 'Pointer ray is parallel to or behind the working plane. Use “View work plane”.');
  }
  corners() {
    if (this.bounds.isEmpty()) return [];
    const b = this.visibleBounds(), out = [];
    for (const x of [b.min.x,b.max.x]) for (const y of [b.min.y,b.max.y]) for (const z of [b.min.z,b.max.z]) out.push(new Vector3(x,y,z).sub(this.origin));
    return out;
  }
  verticalSize() { return (this.camera.top-this.camera.bottom)/this.camera.zoom; }
  setHeight(height) {
    this.height = height; this.camera.zoom = 1;
    const half = height/2, aspect = this.width/this.pixels;
    this.camera.top = half; this.camera.bottom = -half; this.camera.left = -half*aspect; this.camera.right = half*aspect;
    this.camera.updateProjectionMatrix();
  }
  depth() {
    const visible=this.visibleBounds(),size = visible.isEmpty() ? 0 : visible.getSize(new Vector3()).length();
    const span = Math.max(size,this.controls.target.length(),this.verticalSize(),1e-9);
    const direction = this.camera.position.clone().sub(this.controls.target).normalize();
    if (!direction.lengthSq()) direction.set(1,1,1).normalize();
    this.camera.position.copy(this.controls.target).addScaledVector(direction,span*3);
    this.camera.near = span*1e-5; this.camera.far = span*8;
    this.camera.lookAt(this.controls.target); this.camera.updateProjectionMatrix(); this.camera.updateMatrixWorld();
  }
  fitAtTarget() {
    this.depth();
    const inverse = this.camera.quaternion.clone().invert(); let x = 0, y = 0;
    for (const corner of this.corners()) { corner.sub(this.controls.target).applyQuaternion(inverse); x = Math.max(x,Math.abs(corner.x)); y = Math.max(y,Math.abs(corner.y)); }
    const extent = Math.max(2*y,2*x*this.pixels/this.width);
    const radius = this.bounds.isEmpty() ? 0 : this.bounds.getSize(new Vector3()).length();
    this.setHeight(extent > Math.max(radius*1e-12,1e-15) ? extent*1.2 : 20*this.workPlane.step);
    this.depth(); this.controls.cursor.copy(this.controls.target); this.controls.update(); this.requestRender();
  }
  fitModel() {
    const direction = this.camera.position.clone().sub(this.controls.target).normalize();
    this.controls.target.copy(this.visibleBounds().getCenter(new Vector3()).sub(this.origin)); this.camera.position.copy(this.controls.target).add(direction); this.fitAtTarget();
  }
  setWorkPlane(settings) {
    const next = {...this.workPlane,...settings};
    if(this.dimension===2&&(next.plane!=='XY'||next.offset!==0))throw new Error('Planar models use XY at Z = 0.');
    if (!['XY','XZ','YZ'].includes(next.plane) || !Number.isFinite(next.offset) || !Number.isFinite(next.step) || next.step <= 0) throw new Error('Plane offset must be finite; grid step must be positive and finite.');
    this.cancelGesture(); this.workPlane = next; this.requestRender();
  }
  resetOrigin() {
    this.setWorkPlane({plane:'XY',offset:0});
    this.controls.target.copy(this.origin).negate(); this.setView('XY'); this.fitAtTarget();
  }
  planeAxes() { return {XY:[0,1,2],XZ:[0,2,1],YZ:[1,2,0]}[this.workPlane.plane]; }
  updateGrid() {
    const {step,offset,visible} = this.workPlane, [u,v,n] = this.planeAxes();
    const extent = this.verticalSize()*Math.max(1,this.width/this.pixels)*1.5;
    const power = Math.max(0,Math.ceil(Math.log2(extent)-Math.log2(step)-Math.log2(40)));
    this.displayStep = step*Math.pow(2,power);
    const target = this.controls.target.clone().add(this.origin).toArray(), origin = this.origin.toArray();
    const spacing = this.displayStep;
    const cu = Math.round(target[u]/spacing)*spacing, cv = Math.round(target[v]/spacing)*spacing;
    const radius = Math.min(42,Math.ceil(extent/spacing)+1);
    const key = JSON.stringify([this.workPlane,spacing,cu,cv,radius,origin,extent]);
    if (key === this.gridKey) return; this.gridKey = key;
    this.remove(this.grid); this.grid = null; this.remove(this.axes); this.axes = null;
    if (visible && Number.isFinite(spacing) && spacing > 0) {
      const positions = [];
      const put = (a,b) => { const p = [0,0,0]; p[u] = a-origin[u]; p[v] = b-origin[v]; p[n] = offset-origin[n]; positions.push(...p); };
      for (let i = -radius; i <= radius; i++) {
        put(cu+i*spacing,cv-radius*spacing); put(cu+i*spacing,cv+radius*spacing);
        put(cu-radius*spacing,cv+i*spacing); put(cu+radius*spacing,cv+i*spacing);
      }
      this.grid = this.object(positions,false,'#ced9e4'); this.grid.material.depthWrite = false; this.grid.renderOrder = -1;
    }
    const length = Math.max(this.verticalSize()*.2,step), positions = [], colors = [];
    for (let axis = 0; axis < 3; axis++) {
      const a = this.origin.clone().negate(), b = a.clone(); b.setComponent(axis,b.getComponent(axis)+length);
      positions.push(...a.toArray(),...b.toArray());
      const c = [[.75,.15,.15],[.1,.5,.15],[.15,.3,.8]][axis]; colors.push(...c,...c);
    }
    this.axes = this.object(positions,false,'#ffffff'); this.axes.geometry.setAttribute('color',new Float32BufferAttribute(colors,3)); this.axes.material.vertexColors = true;
  }
  setView(view) {
    if(this.dimension===2&&!['XY','plane'].includes(view))return;
    if (view === 'plane') {
      const n = this.planeAxes()[2], direction = this.camera.position.clone().sub(this.controls.target);
      this.controls.target.setComponent(n,this.workPlane.offset-this.origin.getComponent(n));
      this.camera.position.copy(this.controls.target).add(direction); view = this.workPlane.plane;
    }
    const directions = {XY:[0,0,1],XZ:[0,-1,0],YZ:[1,0,0],iso:[1,1,1]};
    if (!directions[view]) throw new Error('Unknown view.');
    const up = new Vector3(0,view === 'XY' ? 1 : 0,view === 'XY' ? 0 : 1);
    if (!this.camera.up.equals(up)) { this.camera.up.copy(up); this.makeControls(); }
    this.camera.position.copy(this.controls.target).add(new Vector3(...directions[view]).normalize());
    this.depth(); this.controls.cursor.copy(this.controls.target); this.controls.update(); this.requestRender();
  }
  setSelection(selection) {
    this.selection = selection; this.remove(this.highlight); this.highlight = null;
    if (selection && ['node','element'].includes(selection.kind)) {
      const item = (selection.kind === 'node' ? this.model.nodes : this.model.elements)[selection.index];
      if (item) {
        const positions = selection.kind === 'node' ? this.local(item).toArray() : [...this.local(this.model.nodes[item.node1-1]).toArray(),...this.local(this.model.nodes[item.node2-1]).toArray()];
        this.highlight = this.object(positions,selection.kind === 'node','#e46816',11);
        this.highlight.material.depthTest = false; this.highlight.renderOrder = 10;
      }
    }
    this.requestRender();
  }
  hit(clientX,clientY) {
    if (this.lost || !this.nodes) return null;
    const rect = this.canvas.getBoundingClientRect(), ray = new Raycaster();
    ray.params.Points.threshold = this.verticalSize()/this.pixels*9;
    ray.params.Line.threshold = this.verticalSize()/this.pixels*5;
    ray.setFromCamera(new Vector2((clientX-rect.left)/rect.width*2-1,1-(clientY-rect.top)/rect.height*2),this.camera);
    const nodes = this.nodes.visible?ray.intersectObject(this.nodes):[], lines = nodes.length||!this.members.visible ? [] : ray.intersectObject(this.members);
    return nodes.length ? {kind:'node',index:nodes[0].index} : lines.length ? {kind:'element',index:Math.floor(lines[0].index/2)} : null;
  }
  pick(clientX,clientY) {
    if(this.lost)return;
    const geometry=this.sections.pick(clientX,clientY)||this.results.pick(clientX,clientY)||this.hit(clientX,clientY),hits=this.tool==='select'?this.conditions.pick(clientX,clientY):[];
    if(geometry)hits.push(geometry);
    if(hits.length>1&&this.callbacks.choose){this.callbacks.choose(hits);return;}
    const selected = hits[0]||null;
    if (this.callbacks.select) this.callbacks.select(selected); else this.setSelection(selected);
  }
  resize() {
    if (this.disposed) return;
    const height = this.verticalSize(); this.width = Math.max(this.host.clientWidth,1); this.pixels = Math.max(this.host.clientHeight,1);
    this.renderer.setSize(this.width,this.pixels,false); this.overlay.width = this.width; this.overlay.height = this.pixels;
    this.setHeight(height); this.requestRender();
  }
  labels() {
    const context=this.overlay.getContext('2d');this.labelFrame=[];
    const ctx=new Proxy(context,{get:(target,key)=>key==='fillText'?(text,x,y)=>{this.labelFrame.push({text:String(text),x,y,color:target.fillStyle,font:target.font});target.fillText(text,x,y);}:typeof target[key]==='function'?target[key].bind(target):target[key],set:(target,key,value)=>{target[key]=value;return true;}});
    ctx.clearRect(0,0,this.width,this.pixels); ctx.font = '11px sans-serif';
    let count = 0;
    const draw = (kind,index) => {
      const item = (kind === 'node' ? this.model.nodes : this.model.elements)[index]; if (!item) return;
      const p = kind === 'node' ? this.local(item) : this.local(this.model.nodes[item.node1-1]).add(this.local(this.model.nodes[item.node2-1])).multiplyScalar(.5);
      p.project(this.camera); if (Math.abs(p.x)>.98 || Math.abs(p.y)>.98 || Math.abs(p.z)>1) return;
      ctx.fillStyle = this.selection && this.selection.kind === kind && this.selection.index === index ? '#ba4b00' : '#425b70';
      ctx.fillText((kind === 'node' ? 'N' : 'E')+(item.id??(index+1)),(p.x+1)*this.width/2+6,(1-p.y)*this.pixels/2-6); count++;
    };
    if (this.selection && ['node','element'].includes(this.selection.kind)) draw(this.selection.kind,this.selection.index);
    for(const label of this.resultLabels||[]){if(count>=200)break;const p=new Vector3(...label.position).sub(this.origin).project(this.camera);if(Math.abs(p.x)>1||Math.abs(p.y)>1||Math.abs(p.z)>1)continue;const x=(p.x+1)*this.width/2+8,y=(1-p.y)*this.pixels/2-10;ctx.fillStyle='#ffffffe8';ctx.fillRect(x-3,y-12,ctx.measureText(label.text).width+6,16);ctx.fillStyle=label.color||'#263847';ctx.fillText(label.text,x,y);count++;}
    count += this.conditions.drawLabels(ctx,200-count);
    for (const [kind,items] of [['node',this.model.nodes],['element',this.model.elements]]) {
      if(this.labelVisibility&&this.labelVisibility[kind]===false)continue;
      for (let i = 0; i < items.length && count < 200; i++) {
        if (!this.selection || this.selection.kind !== kind || this.selection.index !== i) draw(kind,i);
      }
    }
    this.labelCount = count;
    // Fixed-size orientation indicator is drawn with the camera rotation, independent of pan/zoom.
    const q = this.camera.quaternion.clone().invert(), x = this.width-56, y = 58;
    [['X',[1,0,0],'#bc3030'],['Y',[0,1,0],'#24752e'],['Z',[0,0,1],'#285bbb']].forEach(([name,axis,color],i) => {
      const d = new Vector3(...axis).applyQuaternion(q); ctx.strokeStyle = color; ctx.fillStyle = color;
      ctx.beginPath(); ctx.moveTo(x,y); ctx.lineTo(x+d.x*30,y-d.y*30); ctx.stroke();
      ctx.fillText(name,x+d.x*40-3,y-d.y*40+3);
      const b=this.triadButtons[i];b.style.left=(56+d.x*40-14)+'px';b.style.top=(58-d.y*40-14)+'px';
    });
  }
  requestRender() {
    if (this.disposed || this.frame !== null) return;
    this.frame = requestAnimationFrame(() => {
      this.frame = null; if (this.disposed || this.lost) return;
      this.depth(); if (this.updateGrid) this.updateGrid(); this.conditions.update(); this.sections.update(); this.renderer.render(this.scene,this.camera); this.labels();
      if (this.callbacks.change) this.callbacks.change(this.getViewState());
    });
  }
  getViewState() {
    return {target:this.controls.target.clone().add(this.origin).toArray(),direction:this.camera.position.clone().sub(this.controls.target).normalize().toArray(),up:this.camera.up.toArray(),height:this.verticalSize(),workPlane:{...this.workPlane},displayStep:this.displayStep || this.workPlane.step};
  }
  capture() {
    if(this.lost||this.disposed)throw new Error('Scene export requires an available WebGL context.');
    this.depth();this.updateGrid();this.conditions.update();this.sections.update();this.renderer.render(this.scene,this.camera);this.labels();
  }
  projectPoint(x,y,z) {const p=new Vector3(x,y,z).project(this.camera);return {x:(p.x+1)*this.width/2,y:(1-p.y)*this.pixels/2,z:p.z};}
  dispose() {
    this.disposed = true; if (this.frame !== null) cancelAnimationFrame(this.frame); this.frame = null;
    this.observer.disconnect(); this.controls.removeEventListener('change',this.onChange); this.controls.dispose();
    for (const [event,listener] of [['pointerdown',this.onDown],['pointermove',this.onMove],['pointerup',this.onUp],['pointercancel',this.onCancel],['webglcontextlost',this.onLost],['webglcontextrestored',this.onRestored]]) this.canvas.removeEventListener(event,listener,['pointerdown','pointermove','pointerup'].includes(event));
    for (const object of [...this.scene.children]) this.remove(object);
    this.renderer.dispose(); this.renderer.forceContextLoss(); this.canvas.remove(); this.overlay.remove(); this.message.remove();this.triad.remove();
  }
}
