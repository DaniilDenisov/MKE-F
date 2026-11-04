(function(M){
  'use strict';
  const labels=['ux','uy','uz','thetaX','thetaY','thetaZ'];
  function fail(path,message){throw new Error(path+': '+message);}
  function arr(a,p){if(!Array.isArray(a))fail(p,'expected array');return a;}
  function num(x,p){if(!Number.isFinite(x))fail(p,'expected finite number');return x;}
  function integer(x,lo,hi,p){if(!Number.isInteger(x)||x<lo||x>hi)fail(p,'invalid index');return x;}
  function vec(a,n,p){arr(a,p);if(a.length!==n)fail(p,'invalid length');a.forEach(x=>num(x,p));return a;}
  function mat(a,n,k,p){arr(a,p);if(a.length!==n)fail(p,'invalid row count');a.forEach(r=>vec(r,k,p));return a;}
  function ids(a,n,p){arr(a,p);const seen=new Set();a.forEach(x=>{integer(x,1,n,p);if(seen.has(x))fail(p,'duplicate ID');seen.add(x);});return seen;}
  const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
  M.spatialData={validate:function(data){
    if(!data||data.format!=='mkef-postprocessor'||data.version!==8)fail('$','Only result format v8 is supported. Recalculate the original TXT model.');
    const m=data.model,a=data.analysis;if(!m||!a||m.dimension!==3)fail('model','expected spatial model');
    if(!data.metadata||typeof data.metadata.title!=='string'||!data.metadata.units)fail('metadata','missing title/units');
    const nodes=arr(m.nodes,'nodes'),elements=arr(m.elements,'elements'),d=m.dofPerNode,n=m.numberOfDOFs;
    integer(nodes.length,1,100000,'nodes');integer(elements.length,1,100000,'elements');integer(n,1,1800000,'numberOfDOFs');
    if(![3,6].includes(d)||JSON.stringify(m.dofLabels)!==JSON.stringify(labels.slice(0,d)))fail('dofLabels','incorrect components');
    M.validateRecoveryData(data,d);
    nodes.forEach((node,i)=>{if(node.id!==i+1)fail('nodes','IDs must follow physical row order');vec([node.x,node.y,node.z],3,'XYZ');});
    mat(m.dofMap,nodes.length,d,'dofMap');const used=new Set();m.dofMap.forEach(row=>row.forEach((id,j)=>{integer(id,j<3?1:0,n,'dofMap');if(id&&used.has(id))fail('dofMap','duplicate');if(id)used.add(id);}));
    const registry=arr(m.dofRegistry,'dofRegistry');if(registry.length!==n)fail('dofRegistry','wrong count');
    registry.forEach((r,i)=>{if(r.id!==i+1||!['node','elementEnd'].includes(r.kind))fail('dofRegistry','invalid owner');integer(r.nodeId,1,nodes.length,'registry node');if(r.kind==='elementEnd'){integer(r.elementId,1,elements.length,'registry element');integer(r.end,1,2,'registry end');}if(r.direction)vec(r.direction,3,'direction');});
    const G=m.nodeKinematics,rows=Array.from({length:nodes.length*d},()=>[]),pairs=new Set();
    if(!G||G.rows!==rows.length||G.columns!==n)fail('nodeKinematics','incorrect dimensions');
    arr(G.entries,'entries').forEach(r=>{vec(r,3,'entry');integer(r[0],1,rows.length,'physical row');integer(r[1],1,n,'system DOF');const key=r[0]+':'+r[1];if(pairs.has(key))fail('entries','duplicate');pairs.add(key);rows[r[0]-1].push([r[1]-1,r[2]]);});
    // Registry directions own nodal coordinates; private end rotations cannot leak into G.
    const expected=Array.from({length:rows.length},()=>new Map());
    m.dofMap.forEach((row,i)=>row.forEach((id,j)=>{if(!id)return;const r=registry[id-1];if(r.kind!=='node'||r.nodeId!==i+1)fail('dofMap','registry ownership mismatch');if(j<3)expected[i*d+j].set(id-1,1);else {vec(r.direction,3,'rotation direction');r.direction.forEach((v,k)=>{if(v)expected[i*d+3+k].set(id-1,v);});}}));
    rows.forEach((row,i)=>{const actual=new Map(row),keys=new Set([...actual.keys(),...expected[i].keys()]);keys.forEach(k=>{if(Math.abs((actual.get(k)||0)-(expected[i].get(k)||0))>1e-10)fail('nodeKinematics','registry mismatch');});});
    const sectionCatalog=window.MKEFSpatialSections?window.MKEFSpatialSections.resultCatalog(m):new Map();
    elements.forEach((e,i)=>{
      if(e.sectionId){const entry=sectionCatalog.get(e.sectionId);if(!entry)fail('sectionId','unknown section');for(const k of d===6?['area','Iy','Iz','torsionConstant']:['area'])if(Math.abs(e.properties[k]-entry.properties[k])>1e-10*entry.properties[k])fail('section','catalog property mismatch');num(e.sectionAngle,'sectionAngle');}
      if(e.id!==i+1||e.type!==(d===3?1133:1436))fail('elements','inconsistent family or ID');vec(e.nodeIds,2,'nodeIds');e.nodeIds.forEach(id=>integer(id,1,nodes.length,'node'));if(e.nodeIds[0]===e.nodeIds[1])fail('element','same endpoints');
      const p=nodes[e.nodeIds[0]-1],q=nodes[e.nodeIds[1]-1],delta=[q.x-p.x,q.y-p.y,q.z-p.z],L=Math.hypot(...delta);if(!(L>0)||!Number.isFinite(L)||Math.abs(num(e.length,'length')-L)>1e-10*L)fail('length','geometry mismatch');
      const spec=window.MKEFSpatial.registry.types.find(t=>t.code===e.type);spec.properties.forEach(k=>{if(!(num(e.properties&&e.properties[k],k)>0))fail(k,'must be positive');});
      mat(e.localBasis,3,3,'basis');for(let j=0;j<3;j++)for(let k=0;k<3;k++)if(Math.abs(dot(e.localBasis[j],e.localBasis[k])-(j===k?1:0))>1e-10)fail('basis','not orthonormal');
      if(e.localBasis[0].some((v,j)=>Math.abs(v-delta[j]/L)>1e-10))fail('basis','wrong X direction');
      const R=e.localBasis,cross=[R[0][1]*R[1][2]-R[0][2]*R[1][1],R[0][2]*R[1][0]-R[0][0]*R[1][2],R[0][0]*R[1][1]-R[0][1]*R[1][0]];if(dot(cross,R[2])<1-1e-10)fail('basis','not right handed');
      if(d===6){vec(e.referenceY,3,'referenceY');const norm=Math.hypot(...e.referenceY);if(!norm||!Number.isFinite(norm))fail('referenceY','invalid norm');const ref=e.referenceY.map(v=>v/norm),projection=dot(ref,R[0]),y=ref.map((v,j)=>v-projection*R[0][j]),length=Math.hypot(...y);if(length<1e-10)fail('referenceY','parallel to member');const Y=y.map(v=>v/length),X=R[0],Z=[X[1]*Y[2]-X[2]*Y[1],X[2]*Y[0]-X[0]*Y[2],X[0]*Y[1]-X[1]*Y[0]],angle=(num(e.sectionAngle,'sectionAngle')%360)*Math.PI/180;
        if(Y.some((v,j)=>Math.abs(Math.cos(angle)*v+Math.sin(angle)*Z[j]-R[1][j])>1e-10))fail('basis','reference/angle mismatch');}
      ids(e.globalDOFs,n,'globalDOFs');mat(e.localFromSystem,2*d,e.globalDOFs.length,'localFromSystem');
    });
    ['supports','nodalLoads','elementLoads','mpcs','releases','warnings'].forEach(k=>arr(m[k],k));
    if(['supports','nodalLoads','elementLoads','mpcs','releases'].reduce((sum,k)=>sum+m[k].length,0)>100000)fail('conditions','100000 record limit exceeded');
    m.supports.forEach(r=>{integer(r.node,1,nodes.length,'support node');['fixUx','fixUy','fixUz','fixThetaX','fixThetaY','fixThetaZ'].forEach((k,j)=>{if(r[k]!==undefined)integer(r[k],0,j<d?1:0,'support mask');});});
    m.nodalLoads.forEach(r=>{vec(r,9,'nodal load');integer(r[1],1,nodes.length,'load node');if(![10,11,12,13].includes(r[0]))fail('load','invalid law');});
    m.releases.forEach(r=>{integer(r.elementId,1,elements.length,'release element');integer(r.end,1,2,'release end');if(d!==6||!['Mx','My','Mz'].includes(r.component))fail('release','invalid component');});
    const released=new Set(),privateIds=new Map(),connected=Array.from({length:nodes.length},()=>[]);
    m.releases.forEach(r=>{const key=r.elementId+':'+r.end+':'+['Mx','My','Mz'].indexOf(r.component);if(released.has(key))fail('releases','duplicate');released.add(key);});
    registry.forEach(r=>{if(r.kind==='elementEnd'){const e=elements[r.elementId-1],k=labels.indexOf(r.component)-3,key=r.elementId+':'+r.end+':'+k;
      if(!released.has(key)||privateIds.has(key)||r.nodeId!==e.nodeIds[r.end-1]||r.direction.some((v,j)=>Math.abs(v-e.localBasis[k][j])>1e-10))fail('registry','invalid private rotation');privateIds.set(key,r.id-1);}});
    if(privateIds.size!==released.size)fail('registry','missing private rotation');
    elements.forEach(e=>{
      const expectedRows=Array.from({length:2*d},()=>new Map());
      for(let end=0;end<2;end++)for(let j=0;j<d;j++){
        const row=expectedRows[end*d+j],node=e.nodeIds[end]-1,key=e.id+':'+(end+1)+':'+(j-3);
        if(j>=3&&released.has(key))row.set(privateIds.get(key),1);
        else{const axis=e.localBasis[j%3];for(let k=0;k<3;k++)for(const [id,c]of rows[node*d+(j<3?0:3)+k])row.set(id,(row.get(id)||0)+axis[k]*c);
          if(j>=3)connected[node].push(axis);
        }
      }
      expectedRows.forEach((expected,i)=>{const actual=new Map(e.globalDOFs.map((id,j)=>[id-1,e.localFromSystem[i][j]]));for(const id of new Set([...actual.keys(),...expected.keys()]))if(Math.abs((actual.get(id)||0)-(expected.get(id)||0))>1e-10)fail('localFromSystem','kinematics mismatch');});
    });
    if(d===6)m.dofMap.forEach((map,i)=>{const Q=map.slice(3).filter(Boolean).map(id=>registry[id-1].direction);
      for(let j=0;j<Q.length;j++)for(let k=0;k<Q.length;k++)if(Math.abs(dot(Q[j],Q[k])-(j===k?1:0))>1e-10)fail('registry','rotation basis not orthonormal');
      // Compare spans, not SVD basis signs. Both inclusions are required.
      const basis=[];for(const v of connected[i]){let r=v.slice();for(const q of basis){const c=dot(r,q);r=r.map((x,k)=>x-c*q[k]);}const length=Math.hypot(...r);if(length>1e-10)basis.push(r.map(x=>x/length));}
      if(basis.length!==Q.length||basis.some(v=>1-Q.reduce((s,q)=>s+dot(v,q)**2,0)>1e-9))fail('registry','active rotation span mismatch');
    });
    m.mpcs.forEach(r=>{integer(r.depNode,1,nodes.length,'MPC node');integer(r.depDOF,1,d,'MPC DOF');if(r.rhs!==0)fail('MPC','inhomogeneous');arr(r.masters,'masters').forEach(t=>{vec(t,3,'master');integer(t[0],1,nodes.length,'master node');integer(t[1],1,d,'master DOF');});});
    m.elementLoads.forEach(r=>{integer(r.elementId,1,elements.length,'load element');mat(r.ends,3,2,'load ends');if(d!==6||![20,21].includes(r.type)||![1,2].includes(r.coordinateSystem)||a.type!=='static')fail('elementLoads','unsupported load');});
    const kind=a.type;if(!['static','modal','transient'].includes(kind))fail('analysis','unsupported type');
    let columns=1,selected=Array.from({length:n},(_,i)=>i+1);
    if(kind==='modal'){arr(a.frequenciesHz,'frequencies');columns=a.frequenciesHz.length;integer(columns,1,n,'mode count');a.frequenciesHz.forEach(x=>{if(!(num(x,'frequency')>0))fail('frequency','must be positive');});vec(a.eigenvalues,columns,'eigenvalues');vec(a.angularFrequenciesRadPerSec,columns,'omega');mat(a.modeShapes,n,columns,'modeShapes');}
    if(kind==='transient'){arr(a.time,'time');columns=a.time.length;integer(columns,1,10000000,'time count');a.time.forEach((t,i)=>{num(t,'time');if(i&&t<=a.time[i-1])fail('time','not increasing');});ids(a.globalDOFIds,n,'globalDOFIds');selected=a.globalDOFIds;if(!a.sampling||a.sampling.exportedSampleCount!==columns||a.sampling.originalSampleCount<columns)fail('sampling','invalid sample counts');integer(a.sampling.timeStride,1,10000000,'timeStride');integer(a.sampling.originalSampleCount,columns,10000000,'originalSampleCount');if(a.timeStep!==undefined&&!(num(a.timeStep,'timeStep')>0))fail('timeStep','must be positive');}
    for(const field of ['displacements','velocities','accelerations','loadHistory','loadVector','reactions','supportReactions','mpcForces','mpcMultipliers'])if(a[field]!==undefined){const size=field==='mpcMultipliers'?m.mpcs.length:selected.length;kind==='static'?vec(a[field],size,field):mat(a[field],size,columns,field);}
    if(kind==='static'){
      vec(a.displacements,n,'displacements');arr(a.elementResults,'elementResults');if(a.elementResults.length!==elements.length)fail('elementResults','wrong count');a.elementResults.forEach((r,i)=>{if(r.elementId!==i+1||r.type!==elements[i].type)fail('elementResults','ID/type mismatch');['localEndForces','localDisplacements','equivalentLocalLoadVector'].forEach(k=>vec(r[k],2*d,k));if(d===6)mat(r.sectionForceCoefficients,6,4,'section forces');});
    }
    const physical=a.physicalDOFIds||Array.from({length:rows.length},(_,i)=>i+1);ids(physical,rows.length,'physicalDOFIds');if(kind==='transient'&&a.physicalDOFIds){const available=new Set(selected.map(id=>id-1));if(physical.some(id=>rows[id-1].some(([system])=>!available.has(system))))fail('physicalDOFIds','missing required system coordinates');}
    for(const f of ['nodalSupportReactions','nodalMPCForces'])if(a[f]!==undefined){if(kind==='transient'&&!a.physicalDOFIds)fail(f,'missing physicalDOFIds');kind==='static'?vec(a[f],rows.length,f):mat(a[f],kind==='transient'?physical.length:rows.length,columns,f);}
    if((a.displacementAmplitudeSpectrum||a.nodalDisplacementAmplitudeSpectrum)&&!a.spectrumFrequencyHz)fail('spectrum','missing frequency grid');
    if(a.spectrumFrequencyHz){arr(a.spectrumFrequencyHz,'spectrum frequencies').forEach((v,i)=>{num(v,'spectrum frequency');if(v<0||i&&v<=a.spectrumFrequencyHz[i-1])fail('spectrum','invalid frequencies');});if(a.displacementAmplitudeSpectrum)mat(a.displacementAmplitudeSpectrum,selected.length,a.spectrumFrequencyHz.length,'system spectrum');if(a.nodalDisplacementAmplitudeSpectrum){if(!a.physicalDOFIds)fail('spectrum','missing physicalDOFIds');mat(a.nodalDisplacementAmplitudeSpectrum,physical.length,a.spectrumFrequencyHz.length,'physical spectrum');}}
    return {raw:data,physicalRows:rows,systemRows:new Map(selected.map((id,i)=>[id-1,i])),physicalIndex:new Map(physical.map((id,i)=>[id-1,i])),nodesById:new Map(nodes.map(x=>[x.id,x])),elementsById:new Map(elements.map(x=>[x.id,x])),resultsByElementId:new Map((a.elementResults||[]).map(x=>[x.elementId,x]))};
  }};
}(window.MKEFPost));
