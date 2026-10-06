function test_releases()
% Independent closed-form, kinematic and dynamic release acceptance checks.
if exist('output','dir') ~= 7, mkdir('output'); end
p = StructFEProblem(exampleCasePath('CaseReleaseStatic.txt'));
m = p.GetAnalysisModel(); r = p.RunStatic();
assert(isequal(m.dofMap,[1 2 0;3 4 0]));
assert(isequal(m.elementData.dofs,[1;2;5;3;4;6]));
assert(m.numberOfDOFs == 6);
assert(norm(r.elementResults.localEndForces([3 6]),inf)<1e-8);
assert(norm(r.reactions([2 4])-[2000;2000],inf)<1e-8);
assert(norm(r.displacements([5 6])-[-1;1]*1000*4^3/(24*2e11*8e-6),inf)<1e-12);
assert(norm(r.equilibriumResidual,inf)<1e-8);
assert(issparse(m.stiffness) && issparse(m.mass));
assert(norm(m.stiffness-m.stiffness.',inf)<1e-8 && norm(m.mass-m.mass.',inf)<1e-12);
modal=solveModal(withoutElementLoads(m));
expected=sort([120;2520]*2e11*8e-6/(7850*.01*4^4));
assert(norm(modal.eigenvalues(1:2)-expected,inf)<1e-7);
assert(all(modal.eigenvalues>0));
phi=modal.modeShapes(:,1); omega=modal.angularFrequenciesRadPerSec(1);
t=solveTransient(withoutElementLoads(m),struct('timeStep',1e-4,'duration',.01,'initialDisplacement',phi));
energy=sum(t.velocities.*(m.mass*t.velocities))/2+sum(t.displacements.*(m.stiffness*t.displacements))/2;
assert(max(abs(energy-energy(1)))/energy(1)<1e-9);
assert(norm(t.displacements(:,end)-phi*cos(omega*t.time(end)),inf)<1e-4*norm(phi,inf));
assert(norm(t.reactions([5 6],:),inf)/norm(m.stiffness*t.displacements,inf)<1e-9);

% One released end: fixed-pinned beam under uniform load.
text=fileread(exampleCasePath('CaseReleaseStatic.txt'));
one=strrep(text,sprintf('2\n1,1,Mz\n1,2,Mz'),sprintf('1\n1,2,Mz'));
one=strrep(one,'4,1,0,0,0','1,1,0,0,0');
q=fromText(one); a=q.RunStatic();
assert(abs(a.elementResults.localEndForces(3)-1000*4^2/8)<1e-7);
assert(abs(a.elementResults.localEndForces(6))<1e-8);
assert(abs(a.reactions(2)-2500)<1e-8);

% Rotate a fixed/fixed translation model: local end actions are invariant.
horizontal=strrep(text,'6,2,0,0,0','4,2,0,0,0');
q=fromText(horizontal); hr=q.RunStatic();
q=fromText(strrep(horizontal,'4,0,0', '2.4,3.2,0')); ir=q.RunStatic();
assert(norm(hr.elementResults.localDisplacements-ir.elementResults.localDisplacements,inf)<1e-12);
assert(norm(hr.elementResults.localEndForces-ir.elementResults.localEndForces,inf)<1e-8);
assert(norm(ir.equilibriumResidual,inf)<1e-8);

% All rotational restraints are ignored without restraining private ends.
q=fromText(strrep(text,'4,1,0,0,0','1,1,0,0,0')); a=q.RunStatic();
assert(numel(q.mesh.warnings)==1 && norm(a.displacements-r.displacements)<1e-12);
exportPostprocessorData(q.GetAnalysisModel(),a,fullfile('output','CaseReleaseWarnings.json'));

% Equivalent disconnected ends tied by translational MPCs.
joint=sprintf(['analysis\nmodal\nnodes\n3\n0,0,0\n2,0,0\n4,0,0\nelems_113\n2\n' ...
 '113,1,2,.01,2e11,7850,8e-6\n113,2,3,.01,2e11,7850,8e-6\n' ...
 'releases\n2\n1,2,Mz\n2,1,Mz\nbcfix\n2\n1,1,0,0,0\n1,3,0,0,0\n']);
p=fromText(joint); jm=p.GetAnalysisModel(); jr=p.RunModal();
reference=strrep(joint,sprintf('nodes\n3\n0,0,0\n2,0,0\n4,0,0'),sprintf('nodes\n4\n0,0,0\n2,0,0\n4,0,0\n2,0,0'));
reference=strrep(reference,'113,2,3,','113,4,3,');
reference=strrep(reference,sprintf('releases\n2\n1,2,Mz\n2,1,Mz'),sprintf('mpc\n2\n4,1,0,1,2,1,1\n4,2,0,1,2,2,1'));
q=fromText(reference); qr=q.RunModal();
assert(norm(jr.eigenvalues-qr.eigenvalues,inf)/norm(qr.eigenvalues,inf)<1e-10);
jm.forceBoundaryConditions=[10 2 0 -1000 0 0]; a=solveStatic(jm);
assert(a.elementResults(1).globalDOFs(4)==a.elementResults(2).globalDOFs(1));
assert(abs(a.elementResults(1).localDisplacements(6)-a.elementResults(2).localDisplacements(3))>1e-6);
assert(abs(a.elementResults(1).localEndForces(6))<1e-8 && abs(a.elementResults(2).localEndForces(3))<1e-8);

% A mixed rigid/released joint keeps its shared rotation plus one private DOF.
mixed=strrep(joint,sprintf('releases\n2\n1,2,Mz\n2,1,Mz'),sprintf('releases\n1\n1,2,Mz'));
p=fromText(mixed); mm=p.GetAnalysisModel();
assert(mm.dofMap(2,3)>0 && mm.numberOfDOFs==10);
assert(mm.elementData(1).dofs(6)~=mm.dofMap(2,3));
assert(mm.elementData(2).dofs(3)==mm.dofMap(2,3));

% Repeatable sections and the empty-release compatibility path.
p=fromText(strrep(text,sprintf('releases\n2\n1,1,Mz\n1,2,Mz'),sprintf('releases\n1\n1,1,Mz\nreleases\n1\n1,2,Mz')));
assert(isequal(p.mesh.iMnod,m.dofMap));
legacy=fileread(exampleCasePath('CaseBeam.txt'));
p=fromText(legacy); q=fromText([legacy sprintf('\nreleases\n0\n')]);
assert(isequal(p.K,q.K) && isequal(p.M,q.M) && isequal(p.mesh.iMnod,q.mesh.iMnod));

% Combined v4 + MPC + element loads, including one mode and one selected DOF.
combined=[text sprintf('\nmpc\n1\n2,1,0,1,1,1,1\n')];
p=fromText(combined); cm=p.GetAnalysisModel(); cr=p.RunStatic();
exportPostprocessorData(cm,cr,fullfile('output','CaseReleaseMPCStatic.json'));
cm=withoutElementLoads(cm); mr=solveModal(cm);
assert(numel(mr.eigenvalues)==2);
exportPostprocessorData(cm,mr,fullfile('output','CaseReleaseMPCModal.json'));
tr=solveTransient(cm,struct('timeStep',1e-4,'duration',.001,'initialDisplacement',mr.modeShapes(:,1)));
exportPostprocessorData(cm,tr,fullfile('output','CaseReleaseMPCTransient.json'));
exportPostprocessorData(cm,tr,fullfile('output','CaseReleaseMPCPartial.json'),struct('selectedGlobalDOFs',cm.numberOfDOFs,'timeStride',3));
single=mr; single.eigenvalues=mr.eigenvalues(1); single.frequenciesHz=mr.frequenciesHz(1);
single.angularFrequenciesRadPerSec=mr.angularFrequenciesRadPerSec(1); single.modeShapes=mr.modeShapes(:,1);
single.supportReactions=mr.supportReactions(:,1); single.mpcForces=mr.mpcForces(:,1); single.mpcMultipliers=mr.mpcMultipliers(:,1);
exportPostprocessorData(cm,single,fullfile('output','CaseReleaseMPCSingle.json'));

bad={strrep(text,'1,2,Mz','1,1,Mz'),strrep(text,'1,2,Mz','99,2,Mz'), ...
 strrep(text,'1,2,Mz','1,3,Mz'),strrep(text,'1,2,Mz','1,2,V'), ...
 [text sprintf('\nbcforce_stat\n1\n10,1,0,0,1\n')], ...
 [text sprintf('\nmpc\n1\n1,3,0,1,2,1,1\n')]};
for i=1:numel(bad), mustFail(@() fromText(bad{i})); end
mustFail(@() fromText(strrep(fileread(exampleCasePath('CaseReleaseTransient.txt')),'0.02,2,2','0.02,1,3')));
mustFail(@() fromText([fileread(exampleCasePath('CaseMPCTruss.txt')) sprintf('\nreleases\n1\n1,1,Mz\n')]));
mustFail(@() fromText([text sprintf('\nmpc\n1\n2,1,0,1,1,3,1\n')]));
mustFail(@() fromText([text sprintf('\nreleases\n1\n1,NaN,Mz\n')]));
mechanism=strrep(text,sprintf('bcfix\n2\n4,1,0,0,0\n6,2,0,0,0'),sprintf('bcfix\n1\n4,1,0,0,0'));
q=fromText(mechanism); mustFail(@() q.RunStatic());

% v4 roundtrip, all analyses and sparse selection including internal rotations.
for name={'CaseReleaseStatic','CaseReleaseModal','CaseReleaseTransient'}
    p=StructFEProblem(exampleCasePath([name{1} '.txt'])); result=p.RunSelected(); model=p.GetAnalysisModel();
    data=exportPostprocessorData(model,result,fullfile('output',[name{1} '.json']));
    assert(data.version==4);
    decoded=jsondecode(fileread(fullfile('output',[name{1} '.json'])));
    assert(numel(decoded.model.dofRegistry)==model.numberOfDOFs);
    if strcmp(result.analysisType,'transient')
        exportPostprocessorData(model,result,fullfile('output','CaseReleasePartial.json'), ...
            struct('selectedGlobalDOFs',model.numberOfDOFs,'timeStride',3));
    end
end
end

function m=withoutElementLoads(m)
m.elementLoads=struct([]);
end
function p=fromText(text)
file=[tempname() '.txt']; fid=fopen(file,'w'); fprintf(fid,'%s',text); fclose(fid);
cleanup=onCleanup(@() delete(file)); p=StructFEProblem(file);
end
function mustFail(action)
failed=false; try, action(); catch, failed=true; end
assert(failed,'Expected release validation failure.');
end
