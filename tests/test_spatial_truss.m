function test_spatial_truss()
p=StructFEProblem(exampleCasePath('CaseSpatialTruss.txt')); m=p.GetAnalysisModel(); r=p.RunStatic();
assert(m.dimension==3 && m.dofPerNode==3);
assert(abs(r.displacements(6)-1e-6)<1e-15);
assert(abs(r.reactions(3)+1000)<1e-9);
assert(norm(r.equilibriumResidual,inf)<1e-9);
assert(abs(r.elementResults.axialForce-1000)<1e-9);
map=reshape(1:6,3,2).';
e=createStructuralElement(1133,[1 2 3;3 5 9],[1 2],[.01 2e11 7850],map);
x=e.localBasis(1,:).'; u=[zeros(3,1);x*1e-3];
assert(abs(u.'*e.stiffness*u-.01*2e11/e.length*1e-6)<1e-8);
assert(norm(e.stiffness-e.stiffness.',inf)<1e-6);
assert(norm(e.stiffness*[eye(3);eye(3)],inf)<1e-6);
end
