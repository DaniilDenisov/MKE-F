function test_recovery_export()
% Export v8 results and independent solver recovery references.
directory=fullfile('output','recovery'); if ~exist(directory,'dir'),mkdir(directory);end
names={'CaseSpatialTruss','CaseSpatialFrame','CaseSpatialTrussModal', ...
    'CaseSpatialFrameModal','CaseSpatialTrussTransient','CaseSpatialFrameTransient'};
manifest={};
for i=1:numel(names)
    p=StructFEProblem(exampleCasePath([names{i} '.txt']));
    m=p.GetAnalysisModel();r=p.RunSelected();
    manifest{end+1}=writePair(directory,names{i},m,r,struct());
end
% Unreleased frame baselines can be repeated for reader/storage measurements
% without incorrectly sharing a private released-end coordinate between members.
for kind={'modal','transient'}
    source=['CaseSpatialFrame' upper(kind{1}(1)) kind{1}(2:end)];
    text=regexprep(fileread(exampleCasePath([source '.txt'])),'releases\s+1\s+1,2,Mz\s*','');
    file=fullfile(directory,['benchmark-' kind{1} '.txt']);
    fid=fopen(file,'w');fprintf(fid,'%s',text);fclose(fid);
    p=StructFEProblem(file);m=p.GetAnalysisModel();r=p.RunSelected();
    manifest{end+1}=writePair(directory,['benchmark-' kind{1}],m,r,struct());
end
% Inclined, rolled catalog frame, released end and a projected MPC.
text=sprintf(['analysis\nmodal\nnodes\n2\n0,0,0\n2,1,3\n' ...
 'elems_1436\n1\n1436,1,2,.01,2e11,7850,8e10,8e-6,2e-6,1e-6,0,1,0,23\n' ...
 'bcfix_3d\n1\n1,1,1,1,1,1,1\nreleases\n1\n1,2,Mz\n' ...
 'mpc\n1\n2,2,0,1,2,3,2\nsections\n1\n' ...
 '{"id":1,"name":"Rectangle","type":"rectangle","dimensions":{"b":0.1,"h":0.2}}\n' ...
 'element_sections\n1\n1,1,23\n']);
file=fullfile(directory,'released.txt');fid=fopen(file,'w');fprintf(fid,'%s',text);fclose(fid);
p=StructFEProblem(file);m=p.GetAnalysisModel();modal=solveModal(m);
manifest{end+1}=writePair(directory,'released-modal',m,modal,struct());
r=solveTransient(m,struct('timeStep',.0002,'duration',.01, ...
    'initialDisplacement',modal.modeShapes*(1:size(modal.modeShapes,2)).'));
manifest{end+1}=writePair(directory,'released-transient',m,r,struct('timeStride',3));
manifest{end+1}=writePair(directory,'released-partial',m,r, ...
    struct('selectedGlobalDOFs',m.numberOfDOFs,'timeStride',3,'transientFields',{{'displacements'}}));
manifest{end+1}=writePair(directory,'released-no-displacements',m,r, ...
    struct('transientFields',{{'accelerations'}}));
% Reject invalid operators/options before writing a result file.
bad=m;bad.elementData(1).localMass(1,1)=NaN;
mustFail(@() createPostprocessorData(bad,modal,struct('formatVersion',8)),'MKEF:PostprocessorInvalidResult');
bad=m;bad.elementData(1).localStiffness(1,2)=1e9;
mustFail(@() createPostprocessorData(bad,modal,struct('formatVersion',8)),'MKEF:PostprocessorInvalidModel');
for field={'localStiffness','localMass'}
    bad=m;bad.elementData(1).(field{1})=bad.elementData(1).(field{1})*1e-40;
    small=createPostprocessorData(bad,modal,struct('formatVersion',8));
    key='K';if strcmp(field{1},'localMass'),key='M';end
    assert(isequal(small.model.elements(1).recoveryOperators.(key),bad.elementData(1).(field{1})));
end
for value={9,'8',true,[8 8]}
    mustFail(@() createPostprocessorData(m,modal,struct('formatVersion',value{1})),'MKEF:PostprocessorInvalidOptions');
end
% The mixed 2D/3D browser integration must not depend on a stale manual export.
planar=StructFEProblem(exampleCasePath('CasePreprocessorStatic.txt'));
exportPostprocessorData(planar.GetAnalysisModel(),planar.RunSelected(),fullfile('output','static-frame.json'));
fid=fopen(fullfile(directory,'manifest.json'),'w');fprintf(fid,'%s',jsonencode(manifest));fclose(fid);
fprintf('  Recovery export: %d fresh v8 cases, exact solver K/M, partial data and invalid inputs PASS\n',numel(manifest));
end

function entry=writePair(directory,name,m,r,options)
options.formatVersion=8;
options.operatorStorage='inline';
data=exportPostprocessorData(m,r,fullfile(directory,[name '-v8.json']),options);
tableOptions=options;tableOptions.operatorStorage='table';
exportPostprocessorData(m,r,fullfile(directory,[name '-table.json']),tableOptions);
assert(data.version==8);
assert(~isfield(data.model,'stiffness') && ~isfield(data.model,'mass'));
reference=cell(1,numel(m.elementData));
for i=1:numel(m.elementData)
    e=m.elementData(i);reference{i}=struct('K',e.localStiffness,'M',e.localMass);
    if strcmp(r.analysisType,'static')
        assert(~isfield(data.model.elements(i),'recoveryOperators'));
        reference{i}.elastic=r.elementResults(i).localEndForces;
    else
        op=data.model.elements(i).recoveryOperators;
        assert(isequal(op.K,e.localStiffness) && isequal(op.M,e.localMass));
        % Compare the unmodified solver state; the exporter preserves every double.
        if strcmp(r.analysisType,'modal')
            u=r.modeShapes;a=-u.*repmat(r.eigenvalues(:).',size(u,1),1);
        else
            indices=1:data.analysis.sampling.timeStride:numel(r.time);
            if indices(end)~=numel(r.time),indices(end+1)=numel(r.time);end
            u=r.displacements(:,indices);
            a=r.accelerations(:,indices);
        end
        transformation=e.transformation;
        local=transformation*u(e.dofs,:);
        reference{i}.elastic=jsonMatrix(e.localStiffness*local);
        reference{i}.inertia=jsonMatrix(e.localMass*(transformation*a(e.dofs,:)));
        reference{i}.local=jsonMatrix(local);
    end
end
fid=fopen(fullfile(directory,[name '-reference.json']),'w');fprintf(fid,'%s',encodeResultJson(reference));fclose(fid);
entry=struct('name',name,'analysis',r.analysisType,'elements',numel(m.elementData), ...
    'operatorBytes',dir(fullfile(directory,[name '-v8.json'])).bytes);
end

function mustFail(operation,identifier)
try,operation();catch err,assert(strcmp(err.identifier,identifier),err.message);return;end
error('Expected failure: %s',identifier);
end

function rows=jsonMatrix(values)
rows=cell(size(values,1),1);
for i=1:size(values,1),rows{i}=num2cell(values(i,:));end
end
