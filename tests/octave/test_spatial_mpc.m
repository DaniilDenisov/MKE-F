function test_spatial_mpc()
p=StructFEProblem(exampleCasePath('CaseSpatialFrame.txt')); m=p.GetAnalysisModel();
m.forceBoundaryConditions=[10 2 0 1000 0 0 0 0 0];
m.multiPointConstraints=struct('depNode',2,'depDOF',2,'rhs',0,'masters',[2 3 2]);
r=solveStatic(m); ky=3*2e11*2e-6/2^3; kz=3*2e11*8e-6/2^3;
q=2000/(4*ky+kz);
assert(abs(r.displacements(9)-q)<1e-12 && abs(r.displacements(8)-2*q)<1e-12);
assert(norm(r.reactions-r.supportReactions-r.mpcForces,inf)<1e-9);
assert(abs(r.mpcForces(9)+2*r.mpcForces(8))<1e-9);
% Chained translations at a spatial node; third component is uz, not thetaZ.
m.multiPointConstraints(2)=struct('depNode',2,'depDOF',3,'rhs',0,'masters',[2 1 3]);
r=solveStatic(m); assert(abs(r.displacements(8)-6*r.displacements(7))<1e-12);
bad=m; bad.multiPointConstraints(2).masters=[2 2 1]; mustFail(@() solveStatic(bad));
bad=m; bad.multiPointConstraints(2).depDOF=2; mustFail(@() solveStatic(bad));
bad=m; bad.multiPointConstraints(2).depNode=1; mustFail(@() solveStatic(bad));
bad=m; bad.multiPointConstraints(2).masters=[2 7 1]; mustFail(@() solveStatic(bad));
p=StructFEProblem(exampleCasePath('CaseSpatialTruss.txt')); t=p.GetAnalysisModel();
t.multiPointConstraints=struct('depNode',2,'depDOF',3,'rhs',0,'masters',[1 3 1]);
r=solveStatic(t); assert(norm(r.displacements)==0);
bad=t; bad.multiPointConstraints.depDOF=4; mustFail(@() solveStatic(bad));
text=fileread(exampleCasePath('CaseSpatialFrame.txt'));
text=[text sprintf('\nmpc\n1\n2,2,0,1,2,3,2\n')];
file=[tempname() '.txt']; fid=fopen(file,'w'); fprintf(fid,'%s',text); fclose(fid); cleanup=onCleanup(@() delete(file));
p=StructFEProblem(file); r=p.RunStatic(); assert(abs(r.displacements(8)-2*r.displacements(9))<1e-12);
end
function mustFail(f)
failed=false; try, f(); catch err, failed=strcmp(err.identifier,'MKEF:InvalidMPC'); end
assert(failed,'Expected InvalidMPC diagnostic.');
end
