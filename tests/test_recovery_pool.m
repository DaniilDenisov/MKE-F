function test_recovery_pool()
% Exact pooling, real solved grids and clearly labelled storage benchmarks.
directory=fullfile('output','recovery-pool');if ~exist(directory,'dir'),mkdir(directory);end
a=struct('K',diag([1 2]),'M',diag([3 4]));b=a;b.K(1,1)=1+eps(1);
c=a;c.M(1,1)=3+eps(3);
sample=struct('recovery',struct('version',1),'elements', ...
    struct('recoveryOperators',{a,b,a,c,b}));
pooled=poolRecoveryOperators(sample,'table');
assert(isequal([pooled.elements.recoveryOperatorId],[1 2 1 3 2]));
assert(numel(pooled.recoveryOperators)==3);
assert(isequal(pooled.recoveryOperators(2).K,b.K));
assert(isequal(pooled.recoveryOperators(3).M,c.M));
single=sample;single.elements=single.elements(1);
assert(isequal(poolRecoveryOperators(single,'auto'),single));
assert(isequal(poolRecoveryOperators(sample,'inline'),sample));
rare=single;rare.elements=repmat(single.elements,1,20);
for i=1:20,rare.elements(i).recoveryOperators.K(1,1)=i;end
rare.elements(20)=rare.elements(19);
assert(isequal(poolRecoveryOperators(rare,'auto'),rare),'One small repeat cannot pay for all the references');

% Real model: 3x3 bays, two storeys, 48 nodes and 80 frame elements.
filename=fullfile(directory,'grid.txt');fid=fopen(filename,'w');fprintf(fid,'%s',gridText());fclose(fid);
p=StructFEProblem(filename);grid=p.GetAnalysisModel();modal=solveModal(grid);
transient=solveTransient(grid,struct('timeStep',.001,'duration',.02, ...
    'initialDisplacement',modal.modeShapes(:,1)*.001));
report={};
report{end+1}=measure(directory,'solved-grid-modal',grid,modal,false);
report{end+1}=measure(directory,'solved-grid-transient',grid,transient,false);
% Same operators, different orientation/kinematics. The actual force recovery
% for every member is also compared in JavaScript for both storage layouts.
assert(report{1}.uniqueOperators==2);

p=StructFEProblem(exampleCasePath('CaseSpatialFrame.txt'));base=p.GetAnalysisModel();r=solveModal(base);
report{end+1}=measure(directory,'single-modal',base,r,false);
for pattern={'repeated','mixed','unique'}
    m=base;m.elementData=repmat(base.elementData,1,2000);
    for i=1:numel(m.elementData)
        factor=1;
        if strcmp(pattern{1},'mixed'),factor=1+mod(i-1,8)/10;end
        if strcmp(pattern{1},'unique'),factor=1+i/10000;end
        m.elementData(i).properties(2)=base.elementData.properties(2)*factor;
        m.elementData(i).localStiffness=base.elementData.localStiffness*factor;
    end
    % Histories belong to the base one-element solve. Repeated topology is a
    % storage/reader benchmark, NOT a solved 2000-member physical structure.
    report{end+1}=measure(directory,['storage-' pattern{1}],m,r,true);
end
% Static storage is unchanged; invalid option combinations are rejected.
r=solveStatic(base);
inline=createPostprocessorData(base,r,struct('formatVersion',8,'operatorStorage','inline'));
table=createPostprocessorData(base,r,struct('formatVersion',8,'operatorStorage','table'));
assert(isequal(inline.model,table.model));
try,createPostprocessorData(base,r,struct('operatorStorage','bad'));error('Test:ExpectedFailure','Invalid storage accepted');
catch err,assert(strcmp(err.identifier,'MKEF:PostprocessorInvalidOptions'));end
fid=fopen(fullfile(directory,'measurements.json'),'w');fprintf(fid,'%s',jsonencode(report));fclose(fid);
fprintf('  Operator pooling: exact K/M pairs, adjacent doubles, actual grids and six storage measurements PASS\n');
end

function row=measure(directory,name,m,r,synthetic)
row=struct('name',name,'elements',numel(m.elementData),'synthetic',synthetic,'variants',struct());
for mode={'inline','table','auto'}
    options=struct('formatVersion',8,'operatorStorage',mode{1},'title',name);
    file=fullfile(directory,[name '-' mode{1} '.json']);started=tic();
    data=exportPostprocessorData(m,r,file,options);elapsed=toc(started);
    bytes=dir(file).bytes;
    row.variants.(mode{1})=struct('bytes',bytes,'exportSeconds',elapsed,'recoveryVersion',data.model.recovery.version);
    if strcmp(mode{1},'inline'),original=data;end
    if strcmp(mode{1},'table')
        row.uniqueOperators=numel(data.model.recoveryOperators);
        for i=1:numel(m.elementData)
            op=data.model.recoveryOperators(data.model.elements(i).recoveryOperatorId);
            assert(isequal(op.K,m.elementData(i).localStiffness));
            assert(isequal(op.M,m.elementData(i).localMass));
        end
    end
    assert(isequal(data.analysis,original.analysis));
end
assert(row.variants.auto.bytes<=row.variants.inline.bytes);
expected=1;if row.variants.table.bytes<row.variants.inline.bytes,expected=2;end
assert(row.variants.auto.recoveryVersion==expected);
if expected==2,assert(row.variants.auto.bytes==row.variants.table.bytes);end
fprintf('  %s: %d elements, %d pairs; %d -> %d bytes (%.1f%% saved)\n', ...
    name,row.elements,row.uniqueOperators,row.variants.inline.bytes,row.variants.auto.bytes, ...
    100*(1-row.variants.auto.bytes/row.variants.inline.bytes));
end

function text=gridText()
n=3;id=@(x,y,z) 1+x+(n+1)*y+(n+1)^2*z;
nodes={};elements={};supports={};
for z=0:2
    for y=0:n
        for x=0:n
            node=id(x,y,z);nodes{end+1}=sprintf('%g,%g,%g',2*x,2*y,3*z);
            if z==0,supports{end+1}=sprintf('%d,1,1,1,1,1,1',node);continue;end
            if z>0,elements{end+1}=member(id(x,y,z-1),node,[0 1 0],0);end
            if x<n,elements{end+1}=member(node,id(x+1,y,z),[0 0 1],15);end
            if y<n,elements{end+1}=member(node,id(x,y+1,z),[0 0 1],30);end
        end
    end
end
text=sprintf('analysis\nmodal\nnodes\n%d\n%s\nelems_1436\n%d\n%s\nbcfix_3d\n%d\n%s\n', ...
    numel(nodes),strjoin(nodes,'\n'),numel(elements),strjoin(elements,'\n'),numel(supports),strjoin(supports,'\n'));
text=[text sprintf('releases\n1\n2,2,Mz\nmpc\n1\n18,1,0,1,19,1,0.5\n')];
end
function text=member(a,b,reference,angle)
text=sprintf('1436,%d,%d,.01,2e11,7850,8e10,8e-6,2e-6,1e-6,%g,%g,%g,%g',a,b,reference,angle);
end
