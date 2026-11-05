(function (global) {
  'use strict';
  const M = global.MKEFPre, planar = M.caseFormat, spatial = M.spatialFormat;
  const copy = value => JSON.parse(JSON.stringify(value));
  const info = code => spatial.typeInfo(code) || {dimension:2,family:'truss',properties:[],dofLabels:[],releases:[],elementLoads:[]};
  // UI records use physical XYZ vectors. MPC DOF numbers remain native to the case.
  function fromNative(model, source) {
    if (model.dimension === 3) return model;
    const m = copy(model); m.dimension = 2; m.sourceText = source || '';
    m.nodes = m.nodes.map(n => ({x:n.x,y:n.y,z:0}));
    m.supports = m.supports.map(r => {const a=global.MKEFSupports.mask(r,m.elementType===112?2:3);return {node:r.node,mask:[+a[0],+a[1],0,0,0,+a[2]]};});
    m.loads = m.loads.map(r => ({node:r.node,type:r.type,components:[r.fx,r.fy,0,0,0,r.mz],frequency:r.frequency}));
    m.elementLoads = m.elementLoads.map(r => ({elementId:r.elementId,type:r.type,coordinateSystem:r.coordinateSystem,components:r.type===21?[r.qx1,r.qy1,0,r.qx2,r.qy2,0]:[r.qx,r.qy,0]}));
    return m;
  }
  function toNative(model) {
    if (model.dimension === 3) return model;
    const m = copy(model); delete m.dimension; delete m.sourceText;
    m.nodes = m.nodes.map(n => ({x:n.x,y:n.y}));
    m.supports = m.supports.map(r => ({node:r.node,fixUx:!!r.mask[0],fixUy:!!r.mask[1],fixThetaZ:!!r.mask[5]}));
    m.loads = m.loads.map(r => ({node:r.node,type:r.type,fx:r.components[0],fy:r.components[1],mz:r.components[5],frequency:r.frequency}));
    m.elementLoads = m.elementLoads.map(r => Object.assign({elementId:r.elementId,type:r.type,coordinateSystem:r.coordinateSystem},r.type===21?{qx1:r.components[0],qy1:r.components[1],qx2:r.components[3],qy2:r.components[4]}:{qx:r.components[0],qy:r.components[1]}));
    return m;
  }
  function resolve(m) {
    if(m.dimension===3)return global.MKEFSpatialSections.resolve(m);
    const native=toNative(m);global.MKEFSections.sync(native);return native.elements;
  }
  function validate(m, calculation) {
    if(m.dimension===3){
      if(!calculation&&!m.nodes.length&&m.analysis.type==='transient'){
        const a=m.analysis,errors=spatial.validate({...m,analysis:{...a,type:''}},false);
        if(!Number.isFinite(a.timeStep)||a.timeStep<=0||!Number.isFinite(a.duration)||a.duration<a.timeStep||!Number.isInteger(a.monitorNode)||a.monitorNode<1||!Number.isInteger(a.monitorDOF)||a.monitorDOF<1||a.monitorDOF>info(m.elementType).dofLabels.length)errors.push('Invalid time settings or monitor component.');
        return errors;
      }
      return spatial.validate(m,calculation);
    }
    const native=toNative(m),errors=planar.validate(native,{requireAnalysis:!!calculation});
    if(m.nodes.some(n=>n.z!==0))errors.push('Planar nodes must have Z = 0.');
    if(m.supports.some(r=>[2,3,4].some(i=>r.mask[i])))errors.push('Out-of-plane support in a planar model.');
    if(m.loads.some(r=>[2,3,4].some(i=>r.components[i])))errors.push('Out-of-plane load in a planar model.');
    if(m.elementLoads.some(r=>r.components[2]||r.components[5]))errors.push('Out-of-plane distributed load.');
    return calculation?errors:errors.filter(e=>!['Add at least one node.','Add at least one element.'].includes(e)&&!(m.elementType===0&&e==='Select truss element 112 or frame element 113.')&&!(m.nodes.length===0&&e==='Transient monitor node is outside the model.'));
  }
  function basis(m,e) {
    if(m.dimension===3)return global.MKEFSpatialConditions.basis(m,e);
    const a=m.nodes[e.node1-1],b=m.nodes[e.node2-1],L=Math.hypot(b.x-a.x,b.y-a.y),c=(b.x-a.x)/L,s=(b.y-a.y)/L;
    return [[c,s,0],[-s,c,0],[0,0,1]];
  }
  const A=M.editorAdapter={
    info,fromNative,toNative,resolve,basis,validate,
    empty:code=>spatial.isSpatial(code)?spatial.empty(code):fromNative(planar.newModel('',code)),
    parse:text=>fromNative(planar.parse(text),text),
    serialize:m=>m.dimension===3?spatial.serialize(m):planar.serialize(toNative(m)),
    memberError(m,e){if(m.dimension===3)return spatial.memberError(m,e);const a=m.nodes[e.node1-1],b=m.nodes[e.node2-1];if(!Number.isInteger(e.node1)||!Number.isInteger(e.node2)||!a||!b||!Math.hypot(b.x-a.x,b.y-a.y))return 'Choose two distinct existing nodes.';return info(m.elementType).properties.some(k=>!Number.isFinite(e[k])||e[k]<=0)?'All material and section properties must be positive and finite.':'';},
    physicalDOFs:m=>m.dimension===2?(m.elementType===112?[0,1]:[0,1,5]):(m.elementType===1133?[0,1,2]:[0,1,2,3,4,5]),
    diagnostics(m){if(m.dimension===3)return global.MKEFSpatialConditions.validate(m);const native=toNative(m);return {errors:validate(m,false),warnings:global.MKEFReleases.warnings(native),bases:m.elements.map(e=>basis(m,e)),inactiveSupports:m.supports.flatMap((r,index)=>r.mask[5]&&global.MKEFReleases.inactive(native,r.node)?[{index,component:5}]:[])};},
    sceneConditions(m,diagnostics){const d=diagnostics||A.diagnostics(m),map=n=>m.dimension===2&&n===3?6:n;return {...m,bases:d.bases,inactiveSupports:d.inactiveSupports,analysisType:m.analysis.type,mpcs:m.mpcs.map(r=>({...r,depDOF:map(r.depDOF),masters:r.masters.map(a=>({...a,dof:map(a.dof)}))}))};}
  };
  // Helpers retain the native planar constraint algebra.
  M.editorConditions={...global.MKEFSpatialConditions,validate:A.diagnostics,basis,
    equation:(r,m)=>m&&m.dimension===2?global.MKEFMPC.equation(r):global.MKEFSpatialConditions.equation(r),
    equalComponents(m,dependent,master,components){if(m.dimension===3)return global.MKEFSpatialConditions.equalComponents(m,dependent,master,components);if(!components.length)throw Error('Choose at least one component.');const rows=components.map(d=>({depNode:dependent,depDOF:d,rhs:0,masters:[{node:master,dof:d,coefficient:1}]})),native=toNative(m);native.mpcs.push(...rows);const errors=global.MKEFMPC.validate(native);if(errors.length)throw Error(errors.join('\n'));return rows;},
    axisHelper:(m,i,j,k)=>m.dimension===3?global.MKEFSpatialConditions.axisHelper(m,i,j,k):[global.MKEFMPC.onAxis(toNative(m),i,j,k)],
    convertLoadAxes(m,row,target){
      const out={...row,components:row.components.slice(),coordinateSystem:target};
      if(row.coordinateSystem===target)return out;
      const e=m.elements[row.elementId-1];if(!e)throw Error('Choose an existing element.');
      out.components=global.MKEFSpatialConditions.transformLoadComponents(row.components,basis(m,e),target);
      return out;
    }
  };
}(window));
