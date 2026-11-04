function test_result_resources()
% Actual solved 200-member planar beam, not duplicated result records.
directory=fullfile('output','result-resources');if ~exist(directory,'dir'),mkdir(directory);end
n=200;nodes=cell(n+1,1);elements=cell(n,1);assignments=cell(n,1);
for i=1:n+1,nodes{i}=sprintf('%g,0,0',(i-1)/4);end
for i=1:n,elements{i}=sprintf('113,%d,%d,.02,2e11,7800,.00006666666666666667',i,i+1);assignments{i}=sprintf('%d,1,0',i);end
section=struct('id',1,'name','Rectangle','type','rectangle','dimensions',struct('b',.1,'h',.2));
source=sprintf('analysis\nstatic\nnodes\n%d\n%s\nelems_113\n%d\n%s\nbcfix\n1\n1,1,0,0,0\nbcforce_stat\n1\n10,%d,0,-1000,0\nsections\n1\n%s\nelement_sections\n%d\n%s\n', ...
    n+1,strjoin(nodes,'\n'),n,strjoin(elements,'\n'),n+1,jsonencode(section),n,strjoin(assignments,'\n'));
filename=fullfile(directory,'beam.txt');fid=fopen(filename,'w');fprintf(fid,'%s',source);fclose(fid);
p=StructFEProblem(filename,struct('verbose',false,'plotting',false));m=p.GetAnalysisModel();
static=solveStatic(m);modal=solveModal(m);transient=solveTransient(m,struct('timeStep',.001,'duration',.02,'initialDisplacement',modal.modeShapes(:,1)*.001));
rows={};
for result={static,modal,transient}
    r=result{1};name=['beam-2d-' r.analysisType];started=tic();
    exportPostprocessorData(m,r,fullfile(directory,[name '.json']));
    rows{end+1}=struct('name',name,'elements',n,'exportSeconds',toc(started),'bytes',dir(fullfile(directory,[name '.json'])).bytes);
end
p=StructFEProblem(fullfile('output','recovery-pool','grid.txt'),struct('verbose',false,'plotting',false));m=p.GetAnalysisModel();
m.forceBoundaryConditions=[10 m.numberOfNodes 1000 -3000 -1000 0 0 0 0];r=solveStatic(m);started=tic();
name='grid-3d-static';exportPostprocessorData(m,r,fullfile(directory,[name '.json']));
rows{end+1}=struct('name',name,'elements',numel(m.elementData),'exportSeconds',toc(started),'bytes',dir(fullfile(directory,[name '.json'])).bytes);
fid=fopen(fullfile(directory,'exports.json'),'w');fprintf(fid,'%s',encodeResultJson(rows));fclose(fid);
fprintf('  Physical 200-member beam: static/modal/transient resource fixtures PASS\n');
end
