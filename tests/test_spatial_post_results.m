function test_spatial_post_results()
% Browser parity fixtures: inclined frame, rolled section, private rotation,
% projected physical rotations and a translational MPC. No global K/M in JSON.
text=sprintf(['analysis\nmodal\nnodes\n2\n0,0,0\n2,1,3\n' ...
 'elems_1436\n1\n1436,1,2,.01,2e11,7850,8e10,8e-6,2e-6,1e-6,0,1,0,23\n' ...
 'bcfix_3d\n1\n1,1,1,1,1,1,1\nreleases\n1\n1,2,Mz\nmpc\n1\n2,2,0,1,2,3,2\n']);
file=[tempname() '.txt'];fid=fopen(file,'w');fprintf(fid,'%s',text);fclose(fid);cleanup=onCleanup(@() delete(file));
p=StructFEProblem(file);m=p.GetAnalysisModel();modal=solveModal(m);
exportPostprocessorData(m,modal,fullfile('output','CaseSpatialTiltedModal.json'));
u0=modal.modeShapes*(1:size(modal.modeShapes,2)).';
r=solveTransient(m,struct('timeStep',.0002,'duration',.03,'initialDisplacement',u0));
full=exportPostprocessorData(m,r,fullfile('output','CaseSpatialTiltedTransient.json'),struct('timeStride',3));
[~,spectrum]=computeSingleSidedSpectrum(m.nodalKinematics*r.displacements,r.timeStep);
assert(norm(full.analysis.nodalDisplacementAmplitudeSpectrum-spectrum,inf)<1e-12);
mixed=find(sum(spones(m.nodalKinematics),2)>1);
assert(~isempty(mixed));
wrong=abs(m.nodalKinematics)*r.displacementAmplitudeSpectrum;
assert(norm(wrong(mixed,:)-spectrum(mixed,:),inf)>1e-6*norm(spectrum(mixed,:),inf));
private=find(strcmp({m.dofRegistry.kind},'elementEnd'));
partial=exportPostprocessorData(m,r,fullfile('output','CaseSpatialTiltedPartial.json'),struct('selectedGlobalDOFs',private,'timeStride',3));
assert(isempty(partial.analysis.physicalDOFIds));
e=m.elementData;B=e.transformation;ids=e.dofs;indices=[1 2 4 7 10];
reference=struct('K',e.localStiffness,'M',e.localMass,'frames',indices-1, ...
 'physical',m.nodalKinematics*r.displacements(:,1:3:end), ...
 'u',B*full.analysis.displacements(ids,indices), ...
 'elastic',e.localStiffness*B*full.analysis.displacements(ids,indices), ...
 'inertia',e.localMass*B*full.analysis.accelerations(ids,indices));
assert(norm(reference.elastic(12,:)+reference.inertia(12,:),inf)<1e-8*norm(reference.elastic,inf));
assert(norm(reference.elastic(12,:),inf)>1e-6);
fid=fopen(fullfile('output','spatial-post-parity.json'),'w');fprintf(fid,'%s',jsonencode(reference));fclose(fid);
end
