function test_spatial_export()
if exist('output','dir')~=7, mkdir('output'); end
names={'CaseSpatialTruss','CaseSpatialFrame','CaseSpatialTrussModal','CaseSpatialFrameModal','CaseSpatialTrussTransient','CaseSpatialFrameTransient'};
for i=1:numel(names)
    name=names{i}; p=StructFEProblem(exampleCasePath([name '.txt'])); m=p.GetAnalysisModel(); r=p.RunSelected();
    data=exportPostprocessorData(m,r,fullfile('output',[name '.json']));
    assert(data.version==7 && data.model.dimension==3);
    assert(~isfield(data.model,'stiffness') && ~isfield(data.model.elements,'localMass'));
    q=jsondecode(fileread(fullfile('output',[name '.json'])));
    entries=q.model.nodeKinematics.entries;
    G=sparse(entries(:,1),entries(:,2),entries(:,3),q.model.nodeKinematics.rows,q.model.nodeKinematics.columns);
    if isfield(m,'nodalKinematics'), assert(norm(G-m.nodalKinematics,inf)<1e-12); else, assert(isequal(G,speye(m.numberOfDOFs))); end
    if strcmp(r.analysisType,'static')
        for j=1:numel(q.model.elements)
            e=q.model.elements(j); local=e.localFromSystem*r.displacements(e.globalDOFs);
            assert(norm(local-r.elementResults(j).localDisplacements,inf)<1e-12);
        end
    end
    if strcmp(r.analysisType,'transient')
        assert(data.analysis.timeStep==r.timeStep);
        [~,physicalSpectrum]=computeSingleSidedSpectrum(full(G*r.displacements),r.timeStep);
        assert(norm(data.analysis.nodalDisplacementAmplitudeSpectrum-physicalSpectrum,inf)<1e-12);
        options=struct('selectedGlobalDOFs',m.numberOfDOFs,'timeStride',3,'transientFields',{{'displacements','mpcMultipliers'}});
        partial=exportPostprocessorData(m,r,fullfile('output',[name 'Partial.json']),options);
        assert(size(partial.analysis.displacements,1)==1 && partial.analysis.time(end)==r.time(end));
        assert(~isfield(partial.analysis,'velocities'));
    end
end
p=StructFEProblem(exampleCasePath('CaseSpatialFrame.txt')); m=p.GetAnalysisModel(); r=p.RunStatic();
bad=r; bad.displacements(1)=NaN; mustFail(@() createPostprocessorData(m,bad));
bad=m; bad.nodeCoordinates(1,1)=inf; mustFail(@() createPostprocessorData(bad,r));
bad=m; bad.elementData.transformation(1,1)=2; mustFail(@() createPostprocessorData(bad,r));
bad=r; bad.elementResults.localEndForces=1; mustFail(@() createPostprocessorData(m,bad));
% Combined projected rotations, MPC, linear member loads and nodal loads.
text=fileread(exampleCasePath('CaseSpatialFrame.txt'));
text=[text sprintf('\nreleases\n1\n1,2,Mz\nmpc\n1\n2,2,0,1,2,3,2\neload_linear_3d\n1\n21,1,1,0,0,0,0,1000,2000\n')];
file=[tempname() '.txt']; fid=fopen(file,'w'); fprintf(fid,'%s',text); fclose(fid); cleanup=onCleanup(@() delete(file));
p=StructFEProblem(file); m=p.GetAnalysisModel(); r=p.RunStatic();
data=exportPostprocessorData(m,r,fullfile('output','CaseSpatialCombined.json'));
assert(numel(data.model.mpcs)==1 && numel(data.model.releases)==1 && numel(data.model.elementLoads)==1);
assert(norm(r.equilibriumResidual,inf)<1e-8);
end
function mustFail(f)
failed=false; try, f(); catch, failed=true; end
assert(failed);
end
