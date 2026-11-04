function test_result_numbers()
% Node verifies IEEE bits independently: Octave jsondecode can itself lose an ULP.
directory=fullfile('output','result-numbers');if ~exist(directory,'dir'),mkdir(directory);end
x=[0 -0.0 eps(0) -eps(0) realmin realmin/2 1e-300 1e-40 1e-20 eps/2 eps 1-eps/2 1 1+eps pi realmax];
payload=struct('values',{num2cell(x)},'singleton',{{{1e-40}}},'empty',{{}}, ...
    'text',sprintf('Quote " slash \\ newline\nUnicode: тест'),'flags',[true false]);
write(directory,'compact.json',encodeResultJson(payload));
write(directory,'pretty.json',encodeResultJson(payload,true));
reference=struct('values',{cellstr(num2hex(x))});
quiet=struct('verbose',false,'plotting',false);
problem=StructFEProblem(exampleCasePath('CaseSpatialFrame.txt'),quiet);m=problem.GetAnalysisModel();r=solveModal(m);
m.elementData(1).localStiffness=m.elementData(1).localStiffness*1e-40;
m.elementData(1).localMass=m.elementData(1).localMass*1e-40;
r.modeShapes=r.modeShapes*1e-40;
data=exportPostprocessorData(m,r,fullfile(directory,'tiny-modal.json'),struct('operatorStorage','inline'));
reference.K=hexRows(m.elementData(1).localStiffness);reference.M=hexRows(m.elementData(1).localMass);reference.modeShapes=hexRows(r.modeShapes);
assert(isequal(data.analysis.modeShapes,r.modeShapes));
problem=StructFEProblem(exampleCasePath('CasePreprocessorTransient.txt'),quiet);m=problem.GetAnalysisModel();r=problem.RunSelected();
for name={'displacements','velocities','accelerations','loadHistory','reactions'}
    r.(name{1})=r.(name{1})*1e-40;reference.(name{1})=hexRows(r.(name{1}));
end
exportPostprocessorData(m,r,fullfile(directory,'tiny-transient.json'));
write(directory,'reference.json',jsonencode(reference));
% Malformed numbers must never silently become JSON null.
for bad={NaN,Inf,1i,{Inf}}
    try,encodeResultJson(bad{1});error('Expected failure');catch err,assert(strcmp(err.identifier,'MKEF:PostprocessorEncodingFailed'));end
end
fprintf('  Result numbers: finite doubles, subnormals, matrices and histories exported for bit checks PASS\n');
end
function rows=hexRows(a)
rows=cell(size(a,1),1);for i=1:size(a,1),rows{i}=cellstr(num2hex(a(i,:)));end
end
function write(directory,name,text)
fid=fopen(fullfile(directory,name),'wb');fwrite(fid,unicode2native(text,'UTF-8'),'uint8');fclose(fid);
end
