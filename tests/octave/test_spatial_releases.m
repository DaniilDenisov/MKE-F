function test_spatial_releases()
% Fixed/pinned in each bending plane. Private rotations remain unknowns.
base=fileread(exampleCasePath('CaseSpatialFrame.txt'));
base=regexprep(base,'bcforce_stat_3d[\s\S]*$','');
for axis={'My','Mz'}
    name=axis{1}; j=2; inertia=2e-6; if strcmp(name,'My'), j=3; inertia=8e-6; end
    load=zeros(1,3); load(j)=1000;
    text=[base sprintf('\nbcfix_3d\n1\n2,0,%d,%d,0,0,0\nreleases\n1\n1,2,%s\neload_uniform_3d\n1\n20,1,1,%g,%g,%g\n',j==2,j==3,name,load)];
    p=fromText(text); m=p.GetAnalysisModel(); r=p.RunStatic();
    moment=6; if j==3, moment=5; end
    assert(abs(r.elementResults.localEndForces(moment+6))<1e-8);
    assert(abs(abs(r.elementResults.localEndForces(moment))-1000*2^2/8)<1e-8);
    u=r.elementResults.localDisplacements;
    assert(abs(abs(u(moment+6))-1000*2^3/(48*2e11*inertia))<1e-12);
    assert(any(strcmp({m.dofRegistry.kind},'elementEnd')));
    assert(norm(m.stiffness*m.rigidMotionMap,inf)/norm(m.stiffness,inf)<1e-12);
    assert(norm(r.equilibriumResidual,inf)<1e-8);
end
% Torsion release at a restrained node: a tip torque has no load path.
p=fromText([base sprintf('\nreleases\n1\n1,1,Mx\nbcforce_stat_3d\n1\n10,2,0,0,0,1,0,0\n')]);
mustFail(@() p.RunStatic(),'MKEF:SingularStiffness');
% Moment on a completely absent joint direction is rejected, not discarded.
mustFail(@() fromText([base sprintf('\nreleases\n1\n1,2,Mz\nbcforce_stat_3d\n1\n10,2,0,0,0,0,0,1\n')]),'MKEF:InactiveRotationLoad');
mustFail(@() fromText([base sprintf('\nreleases\n1\n1,2,Mz\nbcforce_stat_3d\n1\n10,2,0,0,0,0,0,1e-14\n')]),'MKEF:InactiveRotationLoad');
mustFail(@() fromText([base sprintf('\nreleases\n1\n1,2,Mz\nmpc\n1\n2,6,0,1,2,1,1\n')]),'MKEF:InvalidMPC');
% One transferred oblique direction: SVD rank is scalar even for a 1x3 A.
tilted=strrep(base,sprintf('2,0,0'),sprintf('2,3,4'));
tilted=[tilted sprintf('\nreleases\n2\n1,2,My\n1,2,Mz\n')];
p=fromText(tilted); one=p.GetAnalysisModel(); G=full(one.nodalKinematics(10:12,:)); x=[2;3;4]/sqrt(29);
assert(norm(G*G.'-x*x.',inf)<1e-12);
mustFail(@() fromText([tilted sprintf('\nmpc\n1\n2,5,0,1,2,4,1.5\n')]),'MKEF:InvalidMPC');
% Two non-parallel members sharing an all-released rotation plane. The active
% node basis must span transmitted local directions, not global mask indices.
text=sprintf(['analysis\nstatic\nnodes\n3\n0,0,0\n2,3,4\n5,1,6\nelems_1436\n2\n' ...
 '1436,1,2,.01,2e11,7850,8e10,8e-6,2e-6,1e-6,0,1,0,27\n' ...
 '1436,3,2,.01,2e11,7850,8e10,8e-6,2e-6,1e-6,0,1,0,13\n' ...
 'bcfix_3d\n2\n1,1,1,1,1,1,1\n3,1,1,1,1,1,1\n' ...
 'releases\n4\n1,2,My\n1,2,Mz\n2,2,My\n2,2,Mz\n' ...
 'bcforce_stat_3d\n1\n10,2,100,200,300,0,0,0\n']);
p=fromText(text); m=p.GetAnalysisModel(); r=p.RunStatic();
assert(rank(full(m.nodalKinematics(10:12,:)),1e-10)==2);
assert(norm(r.equilibriumResidual,inf)<1e-7);
for e=r.elementResults(:).', assert(norm(e.localEndForces([11 12]),inf)<1e-7); end
assert(norm(m.stiffness*m.rigidMotionMap,inf)/norm(m.stiffness,inf)<1e-12);
if exist('output','dir')~=7, mkdir('output'); end
exportPostprocessorData(m,r,fullfile('output','CaseSpatialObliqueRelease.json'));
% Rotate the complete released model including reference vectors and forces.
a=.47; Q=[cos(a) 0 sin(a);0 1 0;-sin(a) 0 cos(a)];
elements=m.elementData; map=reshape(1:18,6,3).'; coords=m.nodeCoordinates*Q.';
for i=1:2
    e=elements(i); props=[e.properties (Q*e.referenceY.').' e.sectionAngle];
    elements(i)=createStructuralElement(1436,coords(e.nodeNumbers,:),e.nodeNumbers,props,map);
end
[mr.dofMap,mr.dofRegistry,mr.elementData,mr.nodalKinematics,mr.rigidMotionMap]=buildSpatialDOFs(3,elements,m.releases,coords);
mr.numberOfDOFs=numel(mr.dofRegistry); mr.numberOfNodes=3; mr.dimension=3; mr.dofPerNode=6; mr.dofLabels=m.dofLabels;
mr.fixedBoundaryConditions=m.fixedBoundaryConditions; mr.multiPointConstraints=m.multiPointConstraints;
mr.elementLoads=m.elementLoads; mr.nodeCoordinates=coords; mr.forceBoundaryConditions=m.forceBoundaryConditions;
mr.forceBoundaryConditions(1,3:5)=(Q*[100;200;300]).';
[mr.stiffness,mr.mass]=assembleGlobalMatrices(mr.elementData,mr.numberOfDOFs);
s=solveStatic(mr);
for i=1:2, assert(norm(s.elementResults(i).localEndForces-r.elementResults(i).localEndForces,inf)<1e-6); end
end
function p=fromText(text)
file=[tempname() '.txt']; fid=fopen(file,'w'); fprintf(fid,'%s',text); fclose(fid); cleanup=onCleanup(@() delete(file)); p=StructFEProblem(file);
end
function mustFail(f,id)
failed=false; try, f(); catch err, failed=strcmp(err.identifier,id); if ~failed, rethrow(err); end; end
assert(failed,['Expected ' id]);
end
