function test_spatial_sections()
if ~exist('output','dir'), mkdir('output'); end
types={'rectangle','circle','square_tube','rectangular_tube','circular_tube','channel','tee','ibeam'};
rows={};
for k=1:numel(types)
    s=struct('id',k,'name',types{k},'type',types{k},'dimensions',struct('b',.1,'h',.2,'d',.1,'t',.005,'tw',.008,'tf',.012));
    p=spatialSectionProperties(s); bigger=s; fields=fieldnames(s.dimensions);
    for j=1:numel(fields), bigger.dimensions.(fields{j})=3*s.dimensions.(fields{j}); end
    q=spatialSectionProperties(bigger); assert(abs(q.torsionConstant/p.torsionConstant-81)<1e-9);
    manual=s; manual.torsion=struct('mode','manual','value',.123); q=spatialSectionProperties(manual); assert(q.torsionConstant==.123);
    coeff=[10 20 -30 0;zeros(3,4);7 -5 3 9;-3 4 10 -7];
    stress=recoverSpatialSectionStress(p,coeff);
    rows{end+1}=struct('section',s,'properties',p,'coefficients',coeff,'stress',stress);
end
% Rounded sections exercise both runtimes, interior extrema and the actual solve/export path.
for kind={'square_tube','rectangular_tube'}
    s=struct('id',1,'name',['Rounded ' kind{1}],'type',kind{1},'dimensions',struct('b',.1,'t',.005,'inner_radius',.005));
    if strcmp(kind{1},'rectangular_tube'),s.dimensions.h=.2;end
    p=spatialSectionProperties(s);
    coeff=[10 20 -30 0;zeros(3,4);7 -35 30 9;-3 40 -100 70];
    rows{end+1}=struct('section',s,'properties',p,'coefficients',coeff,'stress',recoverSpatialSectionStress(p,coeff));
    sharp=s;sharp.dimensions.inner_radius=0;sp=spatialSectionProperties(sharp);
    assert(p.area<sp.area && p.Iy<sp.Iy && p.Iz<sp.Iz && p.torsionConstant~=sp.torsionConstant);
    for radius=[-.001 .045 NaN Inf]
        bad=s;bad.dimensions.inner_radius=radius;
        try,sectionProperties(bad);error('Expected invalid radius');catch err,assert(strcmp(err.identifier,'MKEF:InvalidSection'));end
    end
    for dim=[2 3]
        if dim==2
            base=fileread(exampleCasePath('CaseSectionBeam.txt'));
            text=regexprep(base,'\{"id":1[^\r\n]+',jsonencode(s));
        else
            base=fileread(exampleCasePath('CaseSpatialFrame.txt'));
            text=[base sprintf('\nsections\n1\n%s\nelement_sections\n1\n1,1,0\n',jsonencode(s))];
        end
        name=sprintf('Rounded%dD_%s',dim,kind{1});file=fullfile('output',[name '.txt']);
        fid=fopen(file,'w');fprintf(fid,'%s',text);fclose(fid);
        problem=StructFEProblem(file);m=problem.GetAnalysisModel();r=problem.RunSelected();
        data=exportPostprocessorData(m,r,fullfile('output',[name '.json']));
        assert(data.model.sections(1).dimensions.inner_radius==.005);
        assert(abs(data.model.elements(1).properties.area/p.area-1)<1e-12);
        if dim==3
            assert(abs(data.model.elements(1).properties.torsionConstant/p.torsionConstant-1)<1e-12);
        else
            assert(abs(data.model.elements(1).properties.momentOfInertia/p.Iz-1)<1e-12);
            assert(abs(r.elementResults.sectionStress.maximum.value/(500*p.etaMax/p.Iz)-1)<1e-10);
        end
    end
end
s=struct('type','rectangle','dimensions',struct('b',1,'h',1)); p=spatialSectionProperties(s);
assert(abs(p.torsionConstant-0.140577014955)<2e-12);
if ~exist('output','dir'), mkdir('output'); end
fid=fopen('output/spatial-section-parity.json','w'); fprintf(fid,'%s',jsonencode(rows)); fclose(fid);
fprintf('Spatial section properties PASS\n');
for name={'CaseSpatialFrame','CaseSpatialFrameModal','CaseSpatialFrameTransient','CaseSpatialTruss'}
    k=1; if strcmp(name{1},'CaseSpatialTruss'), k=6; end
    s=rows{k}.section; s.torsion=struct('mode','manual','value',1e-6);
    text=[fileread(exampleCasePath([name{1} '.txt'])) sprintf('\nsections\n1\n%s\nelement_sections\n1\n1,%d,-37.5\n',jsonencode(s),s.id)];
    file=fullfile('output',[name{1} 'Section.txt']);fid=fopen(file,'w');fprintf(fid,'%s',text);fclose(fid);
    problem=StructFEProblem(file); m=problem.GetAnalysisModel();r=problem.RunSelected();
    data=exportPostprocessorData(m,r,fullfile('output',[name{1} 'Section.json']));
    assert(data.model.elements(1).sectionAngle==-37.5);
    assert(data.model.elements(1).sectionId==s.id);
    if strcmp(r.analysisType,'static'),assert(isstruct(data.analysis.elementResults(1).sectionStress));end
end
fprintf('Spatial section TXT / solver / v7 PASS\n');
% Legacy open-profile models still build and export, but every solver rejects them.
base=fileread(exampleCasePath('CaseSpatialFrame.txt'));
for k=1:numel(types)
    s=rows{k}.section;
    text=[base sprintf('\nsections\n1\n%s\nelement_sections\n1\n1,%d,0\n',jsonencode(s),s.id)];
    file=fullfile('output',['LegacyFrame-' s.type '.txt']);fid=fopen(file,'w');fprintf(fid,'%s',text);fclose(fid);
    problem=StructFEProblem(file);m=problem.GetAnalysisModel();
    if k<=5
        r=solveStatic(m); assert(norm(r.equilibriumResidual,inf)<1e-7);
    else
        mustReject(@() solveStatic(m)); mustReject(@() solveModal(m));
        mustReject(@() solveTransient(m,struct('timeStep',.001,'duration',.002)));
        mustReject(@() problem.RunSelected());
        % Reconstruct a pre-policy result with the same numerical matrices.
        numerical=m; numerical.elementData(1).sectionId=0;
        r=solveStatic(numerical);r.elementResults=recoverElementResults(m,r.displacements);
        exportPostprocessorData(m,r,fullfile('output',['CaseSpatialLegacy-' s.type '.json']));
        % Unassigned old catalog entries and manual numerical properties remain valid.
        validateSpatialFrameSections(numerical);
    end
end
fprintf('Spatial frame section policy and legacy exports PASS\n');
% Released inclined member: compare the actual elastic forces, not condensed ends.
text=fileread(fullfile('output','CaseSpatialFrameModalSection.txt'));
text=strrep(text,sprintf('2,0,0'),sprintf('2,1,3'));
text=[text sprintf('\nmpc\n1\n2,2,0,1,2,3,2\n')];
file=fullfile('output','spatial-section-released.txt');fid=fopen(file,'w');fprintf(fid,'%s',text);fclose(fid);
problem=StructFEProblem(file);m=problem.GetAnalysisModel();modal=solveModal(m);
exportPostprocessorData(m,modal,'output/CaseSpatialSectionReleasedModal.json');
transient=solveTransient(m,struct('timeStep',.0002,'duration',.01,'initialDisplacement',modal.modeShapes*(1:size(modal.modeShapes,2)).'));
exportPostprocessorData(m,transient,'output/CaseSpatialSectionReleasedTransient.json',struct('timeStride',3));
exportPostprocessorData(m,transient,'output/CaseSpatialSectionReleasedPartial.json',struct('selectedGlobalDOFs',m.numberOfDOFs,'timeStride',3));
reference={};
for frame=[1 4 10]
    er=recoverElementResults(m,transient.displacements(:,frame));
    reference{end+1}=struct('frame',(frame-1)/3,'coefficients',er.sectionForceCoefficients,'stress',er.sectionStress);
end
fid=fopen('output/spatial-section-dynamic-parity.json','w');fprintf(fid,'%s',jsonencode(reference));fclose(fid);
fprintf('Released spatial section dynamic fixtures PASS\n');
end
function mustReject(operation)
try, operation(); catch err
    assert(strcmp(err.identifier,'MKEF:UnsupportedFrameSection'));
    assert(~isempty(strfind(err.message,'members 1 (')));
    assert(~isempty(strfind(err.message,'shear-center')));return;
end
error('Expected MKEF:UnsupportedFrameSection');
end
