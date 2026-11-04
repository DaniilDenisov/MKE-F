function test_planar_recovery()
% Fresh planar v8 cases and references from solver-owned operators.
directory=fullfile('output','planar-recovery');if ~exist(directory,'dir'),mkdir(directory);end
names={'CasePreprocessorStatic','CaseUniformFrame','CaseLinearFrame','CaseLinearMPC', ...
    'CaseLinearRelease','CaseMPCTruss','CaseMPCFrame','CaseSectionBeam','CaseReleaseStatic', ...
    'CasePreprocessorModal','CasePreprocessorTransient','CaseMPCModal','CaseMPCTransient', ...
    'CaseReleaseModal','CaseReleaseTransient'};
manifest={};
for i=1:numel(names)
    p=StructFEProblem(exampleCasePath([names{i} '.txt']));m=p.GetAnalysisModel();r=p.RunSelected();
    manifest{end+1}=writeCase(directory,names{i},m,r,struct());
end
% Inclined released frame with MPC and a global linear load.
text=fileread(exampleCasePath('CaseReleaseStatic.txt'));
text=strrep(text,sprintf('4,0,0'),sprintf('3,2,0'));
text=regexprep(text,'eload_uniform[\s\S]*','');
text=[text sprintf('\nmpc\n1\n2,1,0,1,1,1,1\n')];
file=fullfile(directory,'inclined.txt');fid=fopen(file,'w');fprintf(fid,'%s',text);fclose(fid);
p=StructFEProblem(file);m=p.GetAnalysisModel();modal=solveModal(m);
manifest{end+1}=writeCase(directory,'inclined-modal',m,modal,struct());
r=solveTransient(m,struct('timeStep',.0001,'duration',.002,'initialDisplacement',modal.modeShapes(:,1)*.001));
manifest{end+1}=writeCase(directory,'inclined-transient',m,r,struct('timeStride',3));
manifest{end+1}=writeCase(directory,'inclined-partial',m,r,struct('selectedGlobalDOFs',m.numberOfDOFs,'timeStride',3,'transientFields',{{'displacements'}}));
manifest{end+1}=writeCase(directory,'inclined-acceleration',m,r,struct('transientFields',{{'accelerations'}}));
single=modal;
for f={'eigenvalues','frequenciesHz','angularFrequenciesRadPerSec'},single.(f{1})=single.(f{1})(1);end
for f={'modeShapes','supportReactions','mpcForces','mpcMultipliers'},single.(f{1})=single.(f{1})(:,1);end
manifest{end+1}=writeCase(directory,'inclined-single',m,single,struct());
text=[text sprintf('eload_linear\n1\n21,1,2,200,-1000,-100,-2000\n')];
fid=fopen(file,'w');fprintf(fid,'%s',text);fclose(fid);
p=StructFEProblem(file);manifest{end+1}=writeCase(directory,'inclined-static',p.GetAnalysisModel(),p.RunStatic(),struct());
% Truss dynamics use 4x4 K/M and carry no fictitious shear/bending forces.
p=StructFEProblem(exampleCasePath('CaseMPCTruss.txt'));m=p.GetAnalysisModel();modal=solveModal(m);
manifest{end+1}=writeCase(directory,'truss-modal',m,modal,struct());
r=solveTransient(m,struct('timeStep',.0001,'duration',.002,'initialDisplacement',modal.modeShapes(:,1)*.001));
manifest{end+1}=writeCase(directory,'truss-transient',m,r,struct());
bad=m;bad.elementData(1).localMass(1,1)=NaN;
mustFail(@()createPostprocessorData(bad,modal,struct('formatVersion',8)));
bad=m;bad.elementData(1).localStiffness(1,2)=1e30;
mustFail(@()createPostprocessorData(bad,modal,struct('formatVersion',8)));
% Independent closed form: pinned beam, L=4, q=-1000 => midpoint M=2000.
p=StructFEProblem(exampleCasePath('CaseReleaseStatic.txt'));r=p.RunStatic();c=r.elementResults.sectionForceCoefficients;
assert(abs(polyval(fliplr(c(3,:)),.5)-2000)<1e-8);
assert(abs(polyval(fliplr(c(2,:)),.5))<1e-8);
fid=fopen(fullfile(directory,'manifest.json'),'w');fprintf(fid,'%s',jsonencode(manifest));fclose(fid);
fprintf('  Planar v8: %d cases, solver K/M, coefficients, releases/MPC and partial data PASS\n',numel(manifest));
end

function entry=writeCase(directory,name,m,r,options)
options.formatVersion=8;reference={};
for mode={'inline','table','auto'}
    options.operatorStorage=mode{1};
    data=exportPostprocessorData(m,r,fullfile(directory,[name '-' mode{1} '.json']),options);
    assert(data.version==8 && data.model.dimension==2);
    if strcmp(mode{1},'inline'),original=data;end
    assert(isequal(data.analysis,original.analysis));
    for i=1:numel(m.elementData)
        e=m.elementData(i);out=data.model.elements(i);
        assert(isequal(out.globalDOFs,e.dofs(:).') && isequal(out.localFromSystem,e.transformation));
        assert(norm(e.mass-e.transformation.'*e.localMass*e.transformation,inf)<1e-12*norm(e.mass,inf));
        if strcmp(r.analysisType,'static')
            assert(~isfield(out,'recoveryOperators') && ~isfield(out,'recoveryOperatorId'));
            assert(isequal(data.analysis.elementResults(i).sectionForceCoefficients,r.elementResults(i).sectionForceCoefficients));
        else
            if data.model.recovery.version==2,op=data.model.recoveryOperators(out.recoveryOperatorId);else,op=out.recoveryOperators;end
            assert(isequal(op.K,e.localStiffness) && isequal(op.M,e.localMass));
        end
    end
end
for i=1:numel(m.elementData)
    e=m.elementData(i);
    if strcmp(r.analysisType,'static')
        ref=struct('end',r.elementResults(i).localEndForces(:).', ...
            'coeff',{rows(r.elementResults(i).sectionForceCoefficients)});
    else
        if strcmp(r.analysisType,'modal')
            u=r.modeShapes;acc=-u.*repmat(r.eigenvalues(:).',size(u,1),1);
        else
            indices=1:data.analysis.sampling.timeStride:numel(r.time);if indices(end)~=numel(r.time),indices(end+1)=numel(r.time);end
            u=r.displacements(:,indices);acc=r.accelerations(:,indices);
        end
        T=e.transformation;
        ref=struct('local',{rows(T*u(e.dofs,:))},'elastic',{rows(e.localStiffness*T*u(e.dofs,:))},'inertia',{rows(e.localMass*T*acc(e.dofs,:))});
    end
    reference{end+1}=ref;
end
fid=fopen(fullfile(directory,[name '-reference.json']),'w');fprintf(fid,'%s',encodeResultJson(reference));fclose(fid);
entry=struct('name',name,'analysis',r.analysisType,'elements',numel(m.elementData));
end
function value=rows(A)
value=cell(size(A,1),1);for i=1:size(A,1),value{i}=num2cell(A(i,:));end
end
function mustFail(f)
try,f();catch err,assert(strncmp(err.identifier,'MKEF:Postprocessor',18));return;end
error('Expected invalid planar export to fail.');
end
