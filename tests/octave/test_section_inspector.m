function test_section_inspector()
% Fresh static cases for the paired section maps (including zero and compression).
if ~exist('output','dir'),mkdir('output');end
circle=struct('id',1,'name','Circle 100','type','circle','dimensions',struct('d',.1));
tube=struct('id',1,'name','Square tube 100 x 5','type','square_tube','dimensions',struct('b',.1,'t',.005));
frame=fileread(exampleCasePath('CaseSpatialFrame.txt'));
truss=fileread(exampleCasePath('CaseSpatialTruss.txt'));
cases={
    'SectionInspectorTorsion',strrep(frame,'10,2,1000,2000,3000,400,0,0','10,2,0,0,0,400,0,0'),circle;
    'SectionInspectorZero',strrep(frame,'10,2,1000,2000,3000,400,0,0','10,2,0,0,0,0,0,0'),circle;
    'SectionInspectorCompression',strrep(truss,'10,2,0,0,1000,0,0,0','10,2,0,0,-1000,0,0,0'),tube
};
for i=1:size(cases,1)
    source=[cases{i,2} sprintf('\nsections\n1\n%s\nelement_sections\n1\n1,1,0\n',jsonencode(cases{i,3}))];
    file=fullfile('output',[cases{i,1} '.txt']);fid=fopen(file,'w');fprintf(fid,'%s',source);fclose(fid);
    problem=StructFEProblem(file);model=problem.GetAnalysisModel();result=problem.RunSelected();
    data=exportPostprocessorData(model,result,fullfile('output',[cases{i,1} '.json']));
    if i<=2,assert(abs(data.analysis.elementResults(1).sectionStress.maximum.value)<1e-8);end
    if i==3,assert(data.analysis.elementResults(1).sectionStress.maximum.value<0);end
end
fprintf('Section inspector zero, circular torsion and compression fixtures PASS\n');
end
