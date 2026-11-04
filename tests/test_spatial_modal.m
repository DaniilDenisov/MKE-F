function test_spatial_modal()
p=StructFEProblem(exampleCasePath('CaseSpatialTruss.txt')); m=p.GetAnalysisModel(); m.forceBoundaryConditions=zeros(0,9);
r=solveModal(m); assert(abs(r.eigenvalues-3*2e11/(7850*2^2))/r.eigenvalues<1e-12);
p=StructFEProblem(exampleCasePath('CaseSpatialFrame.txt')); m=p.GetAnalysisModel(); m.forceBoundaryConditions=zeros(0,9);
r=solveModal(m); c=buildConstraintTransform(m);
checkModes(m,r,c);
e=m.elementData; U=[eye(3) zeros(3);zeros(3,6);eye(3) zeros(3);zeros(3,6)];
assert(norm(U(:,1:3).'*e.mass*U(:,1:3)-7850*.01*2*eye(3),inf)<1e-10);
t=m; s=t.fixedBoundaryConditions(1); s.node=2; s.fixThetaX=0; t.fixedBoundaryConditions(2)=s;
rt=solveModal(t); expected=3*8e10*1e-6/(7850*(8e-6+2e-6)*2^2);
assert(numel(rt.eigenvalues)==1 && abs(rt.eigenvalues-expected)/expected<1e-12);
% A spatial model restricted to XY matches the current planar mass/stiffness.
map=reshape(1:6,3,2).'; flat=createStructuralElement(113,[0 0;2 0],[1 2],[.01 2e11 7850 2e-6],map);
ids=[1 2 6 7 8 12];
assert(norm(e.stiffness(ids,ids)-flat.stiffness,inf)==0);
assert(norm(e.mass(ids,ids)-flat.mass,inf)<1e-12);
% Consistent bending convergence against the independent cantilever root.
exact=1.875104068711961^2*sqrt(2e11*2e-6/(7850*.01*2^4)); previous=inf;
for count=[1 2 4 8]
    p=beam(count,false); mm=p.GetAnalysisModel(); rr=solveModal(mm);
    err=abs(rr.angularFrequenciesRadPerSec(1)-exact)/exact;
    assert(err<previous); previous=err;
end
assert(previous<3e-6);
% A simply supported beam retains private rotations and converges to pi^2.
exact=pi^2*sqrt(2e11*2e-6/(7850*.01*2^4)); p=beam(8,true); mm=p.GetAnalysisModel(); rr=solveModal(mm);
assert(abs(rr.angularFrequenciesRadPerSec(1)-exact)/exact<2e-5);
checkModes(mm,rr,buildConstraintTransform(mm));
% Equal bending inertias give repeated eigenvalues; compare their subspace.
isotropic=m; er=createStructuralElement(1436,e.nodeCoordinates,[1 2],[.01 2e11 7850 8e10 2e-6 2e-6 1e-6 0 1 0 0],reshape(1:12,6,2).');
isotropic.elementData=er; isotropic.stiffness=sparse(er.stiffness); isotropic.mass=sparse(er.mass);
rr=solveModal(isotropic); assert(abs(rr.eigenvalues(1)-rr.eigenvalues(2))/rr.eigenvalues(1)<1e-10);
P=rr.modeShapes(:,1:2); Q=P*[cos(.3) -sin(.3);sin(.3) cos(.3)];
assert(norm(orth(P)*orth(P).'-orth(Q)*orth(Q).',inf)<1e-12);
end
function checkModes(m,r,c)
P=r.modeShapes; R=c.T.'*(m.stiffness*P-m.mass*P*diag(r.eigenvalues));
assert(norm(R,inf)/max(1,norm(m.stiffness*P,inf))<1e-9);
gram=P.'*m.mass*P; normalized=diag(1./sqrt(diag(gram)))*gram*diag(1./sqrt(diag(gram)));
assert(norm(normalized-eye(size(gram)),inf)<1e-9);
end
function p=beam(n,pinned)
text=sprintf('analysis\nmodal\nnodes\n%d\n',n+1);
for i=0:n, text=[text sprintf('%g,0,0\n',2*i/n)]; end
text=[text sprintf('elems_1436\n%d\n',n)];
for i=1:n, text=[text sprintf('1436,%d,%d,.01,2e11,7850,8e10,8e-6,2e-6,1e-6,0,1,0,0\n',i,i+1)]; end
text=[text sprintf('bcfix_3d\n%d\n',n+1)];
for i=1:n+1
    fy=i==1 || (pinned && i==n+1); rz=(i==1 && ~pinned);
    text=[text sprintf('%d,1,%d,1,1,1,%d\n',i,fy,rz)];
end
if pinned, text=[text sprintf('releases\n2\n1,1,Mz\n%d,2,Mz\n',n)]; end
file=[tempname() '.txt']; fid=fopen(file,'w'); fprintf(fid,'%s',text); fclose(fid); cleanup=onCleanup(@() delete(file)); p=StructFEProblem(file);
end
