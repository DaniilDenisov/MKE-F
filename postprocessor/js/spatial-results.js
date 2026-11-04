(function(M){
  'use strict';
  const S=M.spatialResults={},dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0),mul=(A,v)=>A.map(r=>dot(r,v));
  // Bounded local-operator cache; never assemble or retain a global K/M.
  const operators=new Map();
  function remember(e,value){if(operators.size>=2048)operators.delete(operators.keys().next().value);operators.set(e,value);return value;}
  S.global=(R,v)=>[0,1,2].map(j=>R.reduce((s,r,i)=>s+r[j]*v[i],0));
  S.vector=function(ds,field,index=0,phase=1){
    const a=ds.raw.analysis,n=ds.raw.model.numberOfDOFs,out=Array(n).fill(null);
    if(a.type==='modal'&&field==='accelerations'){const u=S.vector(ds,'displacements',index,phase);return u.map(v=>-a.eigenvalues[index]*v);}
    const values=a[a.type==='modal'&&field==='displacements'?'modeShapes':field];if(!values)return null;
    if(a.type==='static')return values.slice();
    (a.type==='transient'?a.globalDOFIds:Array.from({length:n},(_,i)=>i+1)).forEach((id,row)=>{out[id-1]=values[row][index]*phase;});return out;
  };
  S.physical=function(ds,vector,row){if(!vector)return null;const terms=ds.physicalRows[row];if(!terms.length)return null;let sum=0;for(const [i,v]of terms){if(vector[i]===null||vector[i]===undefined)return null;sum+=v*vector[i];}return sum;};
  S.reaction=function(ds,field,row,index=0,phase=1){const a=ds.raw.analysis,key=field==='supportReactions'?'nodalSupportReactions':'nodalMPCForces',values=a[key];
    if(values){const i=a.type==='transient'?ds.physicalIndex.get(row):row;if(i===undefined)return null;return a.type==='static'?values[i]:values[i][index]*phase;}
    // Truss coordinates are already physical. Frames require exported physical forces.
    if(ds.raw.model.dofPerNode===3){const v=S.vector(ds,field,index,phase);return v&&v[row]!==undefined?v[row]:null;}return null;
  };
  S.local=function(e,vector){if(!vector)return null;return e.localFromSystem.map(row=>{let sum=0;for(let j=0;j<row.length;j++)if(row[j]){const v=vector[e.globalDOFs[j]-1];if(v===null||v===undefined)return null;sum+=row[j]*v;}return sum;});};
  S.operators=function(e){
    if(operators.has(e))return operators.get(e);
    const d=e.type===1133?3:6,n=2*d,K=Array.from({length:n},()=>Array(n).fill(0)),mass=K.map(r=>r.slice()),p=e.properties,L=e.length;
    const block=(A,ids,B,factor=1)=>ids.forEach((i,r)=>ids.forEach((j,c)=>A[i][j]=B[r][c]*factor));
    block(K,[0,d],[[1,-1],[-1,1]],p.area*p.youngsModulus/L);
    if(d===3){for(let j=0;j<3;j++)block(mass,[j,j+3],[[2,1],[1,2]],p.density*p.area*L/6);return remember(e,{K,M:mass});}
    block(K,[3,9],[[1,-1],[-1,1]],p.shearModulus*p.torsionConstant/L);
    const B=[[12,6*L,-12,6*L],[6*L,4*L*L,-6*L,2*L*L],[-12,-6*L,12,-6*L],[6*L,2*L*L,-6*L,4*L*L]],sign=[1,-1,1,-1];
    block(K,[1,5,7,11],B,p.youngsModulus*p.Iz/L**3);block(K,[2,4,8,10],B.map((r,i)=>r.map((v,j)=>v*sign[i]*sign[j])),p.youngsModulus*p.Iy/L**3);
    block(mass,[0,6],[[2,1],[1,2]],p.density*p.area*L/6);block(mass,[3,9],[[2,1],[1,2]],p.density*(p.Iy+p.Iz)*L/6);
    const BM=[[156,22*L,54,-13*L],[22*L,4*L*L,13*L,-3*L*L],[54,13*L,156,-22*L],[-13*L,-3*L*L,-22*L,4*L*L]];
    block(mass,[1,5,7,11],BM,p.density*p.area*L/420);block(mass,[2,4,8,10],BM.map((r,i)=>r.map((v,j)=>v*sign[i]*sign[j])),p.density*p.area*L/420);return remember(e,{K,M:mass});
  };
  S.recover=function(ds,e,vector,acceleration){
    const u=S.local(e,vector);if(!u||u.some(v=>v===null))return null;
    const stored=ds.raw.analysis.type==='static'?ds.resultsByElementId.get(e.id):null;
    const op=!stored||acceleration?S.operators(e):null;
    let end=stored?stored.localEndForces.slice():mul(op.K,u),coeff=Array.from({length:6},()=>[0,0,0,0]);
    if(stored){end=stored.localEndForces.slice();if(stored.sectionForceCoefficients)coeff=stored.sectionForceCoefficients.map(r=>r.slice());else coeff[0][0]=-end[0];}
    else {if(e.type===1133)coeff[0][0]=-end[0];else {for(let i=0;i<6;i++)coeff[i][0]=-end[i];coeff[4][1]=-e.length*end[2];coeff[5][1]=e.length*end[1];}}
    const localA=S.local(e,acceleration),inertia=localA&&!localA.some(v=>v===null)?mul(op.M,localA):null;
    return {u,end,coeff,inertia,total:inertia?end.map((v,i)=>v+inertia[i]):null};
  };
  S.polynomials=function(e,u){
    const L=e.length,d=e.type===1133?3:6,linear=(a,b)=>[a,b-a,0,0];
    if(d===3)return {u:[0,1,2].map(j=>linear(u[j],u[j+3])),r:[[0,0,0,0],[0,0,0,0],[0,0,0,0]]};
    const hermite=(a,b,c,z)=>[a,L*b,-3*a-2*L*b+3*c-L*z,2*a+L*b-2*c+L*z];
    const y=hermite(u[1],u[5],u[7],u[11]),z=hermite(u[2],-u[4],u[8],-u[10]);
    return {u:[linear(u[0],u[6]),y,z],r:[linear(u[3],u[9]),[ -z[1]/L,-2*z[2]/L,-3*z[3]/L,0],[y[1]/L,2*y[2]/L,3*y[3]/L,0]]};
  };
  S.value=(c,t)=>((c[3]*t+c[2])*t+c[1])*t+c[0];
  S.extrema=function(c){const out=[0,1],a=3*c[3],b=2*c[2],z=c[1],scale=Math.max(Math.abs(a),Math.abs(b),Math.abs(z));if(!scale)return out;
    if(Math.abs(a)<1e-14*scale){if(b){const t=-z/b;if(t>0&&t<1)out.push(t);}}else {const discriminant=(b/scale)**2-4*(a/scale)*(z/scale);if(discriminant>=0){const q=-.5*(b/scale+(b>=0?1:-1)*Math.sqrt(discriminant));for(const t of [q/(a/scale),q?z/scale/q:0])if(t>0&&t<1)out.push(t);}}return out.sort((a,b)=>a-b);
  };
  S.sample=function(e,u,t){const p=S.polynomials(e,u),local=p.u.map(c=>S.value(c,t));return {displacement:S.global(e.localBasis,local),rotation:p.r.map(c=>S.value(c,t))};};
  S.coefficient=function(e,result,quantity){const names=['N','Vy','Vz','Mx','My','Mz'];return quantity==='stress'?result.coeff[0].map(v=>v/e.properties.area):result.coeff[names.indexOf(quantity)];};
  S.bounds=function(e,u){const p=S.polynomials(e,u),bound=rows=>Math.hypot(...rows.map(c=>Math.max(...S.extrema(c).map(t=>Math.abs(S.value(c,t))))));return {translation:bound(p.u),rotation:bound(p.r)};};
  S.spectrum=function(ds,row){const a=ds.raw.analysis;if(!a.spectrumFrequencyHz)return null;const i=ds.physicalIndex.get(row);if(a.nodalDisplacementAmplitudeSpectrum&&i!==undefined)return {x:a.spectrumFrequencyHz,y:a.nodalDisplacementAmplitudeSpectrum[i]};
    const terms=ds.physicalRows[row];if(terms.length!==1||!a.displacementAmplitudeSpectrum)return null;const index=ds.systemRows.get(terms[0][0]);if(index===undefined)return null;return {x:a.spectrumFrequencyHz,y:a.displacementAmplitudeSpectrum[index].map(v=>v*Math.abs(terms[0][1]))};
  };
  S.full=ds=>{const a=ds.raw.analysis;return a.type!=='transient'||!!a.displacements&&a.globalDOFIds.length===ds.raw.model.numberOfDOFs;};
}(window.MKEFPost));
