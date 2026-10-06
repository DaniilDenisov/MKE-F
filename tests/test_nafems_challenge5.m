function test_nafems_challenge5()
%TEST_NAFEMS_CHALLENGE5 Saved inputs, independent analytics and quick profile.
root=fileparts(fileparts(mfilename('fullpath')));
source=fullfile(root,'examples','cases','nafems-challenge-5');
destination=tempname(); mkdir(destination);
cleanup=onCleanup(@() rmdir(destination,'s'));
manifest=generate_nafems_challenge5_cases(destination);
assert(numel(manifest)==20);
assert(strcmp(fileread(fullfile(source,'manifest.csv')),fileread(fullfile(destination,'manifest.csv'))));
for row=manifest
    original=fileread(fullfile(source,row.file));
    assert(~any(original==char(13)) && original(end)==char(10));
    assert(strcmp(original,fileread(fullfile(destination,row.file))));
    % Read every supplied mesh, including the conditional 400 mesh.
    mesh=FEMesh(fullfile(source,row.file));
    assert(mesh.numberOfNodes>0);
    if ~strcmp(row.profile,'quick'), continue; end
    p=StructFEProblem(fullfile(source,row.file),struct('verbose',false,'plotting',false));
    model=p.GetAnalysisModel(); result=p.RunSelected();
    nafems_challenge5_check(model,result,row.kind,row.n);
    if strcmp(row.kind,'frame')
        [points,values]=nafems_challenge5_shapes(model,result.modeShapes(:,1));
        assert(size(points,1)==402 && size(values,1)==804);
        assert(norm(values([1 2 403 404]),inf)<1e-12);
        assert(norm(values(401:402)-values(803:804),inf)<1e-12);
    end
end
% A unit rigid rotation must be reconstructed by cubic interpolation exactly.
p=StructFEProblem(fullfile(source,'frame_n002.txt'),struct('verbose',false,'plotting',false));
m=p.GetAnalysisModel(); rotation=zeros(m.numberOfDOFs,1);
rotation(m.dofMap(:,1))=-m.nodeCoordinates(:,2);
rotation(m.dofMap(:,2))=m.nodeCoordinates(:,1);
for d=1:numel(m.dofRegistry)
    if strcmp(m.dofRegistry(d).component,'thetaZ'), rotation(d)=1; end
end
[points,values]=nafems_challenge5_shapes(m,rotation);
expected=reshape([-points(:,2) points(:,1)].',[],1);
assert(norm(values-expected,inf)<1e-12);
end
