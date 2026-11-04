function test_spatial_loads()
p=StructFEProblem(exampleCasePath('CaseSpatialFrame.txt')); m=p.GetAnalysisModel();
m.forceBoundaryConditions=zeros(0,9); L=2; E=2e11; Iy=8e-6; Iz=2e-6;
m.elementLoads=struct('type',20,'elementId',1,'coordinateSystem',1,'ends',[0 0;1000 1000;2000 2000]);
r=solveStatic(m);
assert(abs(r.displacements(8)-1000*L^4/(8*E*Iz))<1e-12);
assert(abs(r.displacements(9)-2000*L^4/(8*E*Iy))<1e-12);
assert(norm(r.reactions(1:6)-[0;-2000;-4000;0;4000;-2000],inf)<1e-8);
assert(norm(r.elementResults.sectionForceCoefficients*ones(4,1),inf)<1e-8);
m.elementLoads.type=21; m.elementLoads.ends(:,1)=0;
r=solveStatic(m);
assert(abs(r.displacements(8)-11*1000*L^4/(120*E*Iz))<1e-12);
assert(norm(r.equilibriumResidual,inf)<1e-8);
% A non-axis-aligned, rolled member: global and local load inputs agree.
map=reshape(1:12,6,2).';
e=createStructuralElement(1436,[1 2 3;3 5 9],[1 2],[m.elementData.properties 0 1 0 31],map);
m.elementData=e; m.nodeCoordinates=e.nodeCoordinates; m.stiffness=sparse(e.stiffness);
[~,~,~,m.nodalKinematics,m.rigidMotionMap]=buildSpatialDOFs(2,e,struct([]),e.nodeCoordinates);
r=solveStatic(m); m.elementLoads.coordinateSystem=2; m.elementLoads.ends=e.localBasis.'*m.elementLoads.ends;
s=solveStatic(m); assert(norm(r.displacements-s.displacements,inf)<1e-12);
assert(norm(r.equilibriumResidual,inf)<1e-7);
% The TXT reader and direct model produce the same consistent loads.
text=fileread(exampleCasePath('CaseSpatialFrame.txt'));
text=[text sprintf('\neload_linear_3d\n1\n21,1,1,0,0,0,0,1000,2000\n')];
file=[tempname() '.txt']; fid=fopen(file,'w'); fprintf(fid,'%s',text); fclose(fid); cleanup=onCleanup(@() delete(file));
p=StructFEProblem(file); r=p.RunStatic(); assert(norm(r.equilibriumResidual,inf)<1e-8);
end
