function test_spatial_frame()
p=StructFEProblem(exampleCasePath('CaseSpatialFrame.txt')); m=p.GetAnalysisModel(); r=p.RunStatic();
E=2e11; A=.01; L=2; G=8e10; Iy=8e-6; Iz=2e-6; J=1e-6;
expected=[1000*L/(E*A);2000*L^3/(3*E*Iz);3000*L^3/(3*E*Iy);400*L/(G*J);-3000*L^2/(2*E*Iy);2000*L^2/(2*E*Iz)];
assert(norm(r.displacements(7:12)-expected,inf)<1e-12);
assert(norm(r.equilibriumResidual,inf)<1e-8);
e=m.elementData; map=reshape(1:12,6,2).';
% Rotate the whole model (geometry, orientation and load) in 3D.
a=.61; b=.47; Q=[cos(a) -sin(a) 0;sin(a) cos(a) 0;0 0 1]*[cos(b) 0 sin(b);0 1 0;-sin(b) 0 cos(b)];
v=Q*[0;1;0]; er=createStructuralElement(1436,e.nodeCoordinates*Q.',[1 2],[e.properties v.' 0],map);
H=blkdiag(Q,Q,Q,Q); mr=m; mr.nodeCoordinates=m.nodeCoordinates*Q.'; mr.elementData=er;
mr.stiffness=H*m.stiffness*H.'; mr.forceBoundaryConditions(1,3:8)=(blkdiag(Q,Q)*m.forceBoundaryConditions(1,3:8).').';
mr.rigidMotionMap=H*m.rigidMotionMap*blkdiag(Q.',Q.');
rr=solveStatic(mr);
assert(norm(rr.displacements-H*r.displacements,inf)<1e-12);
assert(norm(rr.elementResults.localEndForces-r.elementResults.localEndForces,inf)<1e-8);
assert(norm(er.stiffness-H*e.stiffness*H.',inf)/norm(e.stiffness,inf)<1e-12);
% Six rigid motions, including the lever arm of a rigid rotation.
U=zeros(12,6); U(1:3,1:3)=eye(3); U(7:9,1:3)=eye(3);
for j=1:3
    w=zeros(3,1); w(j)=1; U(4:6,j+3)=w; U(10:12,j+3)=w;
    U(7:9,j+3)=cross(w,[L;0;0]);
end
assert(norm(e.stiffness*U,inf)/norm(e.stiffness,inf)<1e-12);
% 90 degree roll exchanges the two global bending compliances.
er=createStructuralElement(1436,e.nodeCoordinates,[1 2],[e.properties 0 1 0 90],map);
ur=er.stiffness(7:12,7:12)\[0;1;0;0;0;0];
assert(abs(ur(2)-L^3/(3*E*Iy))<1e-15);
failed=false; try, createStructuralElement(1436,e.nodeCoordinates,[1 2],[e.properties 1 0 0 0],map); catch, failed=true; end
assert(failed);
end
