function test_mpc()
%TEST_MPC Independent KKT/nullspace verification, histories and input errors.
p = StructFEProblem(exampleCasePath('CaseMPCFrame.txt'));
m = p.GetAnalysisModel(); c = buildConstraintTransform(m);
assert(issparse(c.T) && issparse(c.C));
assert(norm(c.C*c.T,inf) < 1e-14);
assert(isequal(c.dependentDOFs,8));
s = solveStatic(m);
S = speye(m.numberOfDOFs); A = [S(c.fixedDOFs,:); c.C];
augmented = [m.stiffness -A.'; A sparse(size(A,1),size(A,1))];
reference = augmented \ [s.loadVector; zeros(size(A,1),1)];
assert(norm(reference(1:m.numberOfDOFs)-s.displacements,inf) < 1e-10);
assert(norm(reference(m.numberOfDOFs+numel(c.fixedDOFs)+1:end)-s.mpcMultipliers,inf) < 1e-7);
assert(norm(s.reactions-s.supportReactions-s.mpcForces,inf) < 1e-7);
assert(abs(s.mpcMultipliers) > 1);
assert(abs(s.supportReactions(2)-s.reactions(2)) > 1);
assert(norm(c.T.'*s.reactions,inf) < 1e-7);
assert(norm(s.equilibriumResidual,inf) < 1e-7);
assert(norm(c.T.'*m.stiffness*c.T-(c.T.'*m.stiffness*c.T).',inf) < 1e-8);
modal = solveModal(m);
Z = null(full(A));
eigenvalues = sort(eig(Z.'*m.stiffness*Z,Z.'*m.mass*Z));
assert(norm((modal.eigenvalues-eigenvalues)./eigenvalues,inf) < 1e-9);
assert(norm(c.C*modal.modeShapes,inf) < 1e-12);
assert(size(modal.modeShapes,2) == numel(c.independentDOFs));
assert(norm(c.T.'*m.mass*c.T-(c.T.'*m.mass*c.T).',inf) < 1e-12);
dynamic = m; dynamic.forceBoundaryConditions(:,1) = 13;
t = solveTransient(dynamic,struct('timeStep',1e-4,'duration',.001));
for name = {'displacements','velocities','accelerations'}
    assert(norm(c.C*t.(name{1}),inf) < 1e-8);
end
assert(norm(t.reactions-t.supportReactions-t.mpcForces,inf) < 1e-6);
bad = zeros(m.numberOfDOFs,1); bad(8) = 1;
expect('MKEF:InvalidInitialConditions',@() solveTransient(m,struct('timeStep',1e-4,'duration',.001,'initialDisplacement',bad)));
expect('MKEF:InvalidInitialConditions',@() solveTransient(m,struct('timeStep',1e-4,'duration',.001,'initialVelocity',bad)));
chain = m;
chain.multiPointConstraints(2) = equation(2,2,[1 2 .25;3 1 .1]);
chain.multiPointConstraints(3) = equation(3,1,[2 1 .8]);
cc = buildConstraintTransform(chain); result = solveStatic(chain);
assert(norm(cc.C*cc.T,inf) < 1e-14);
assert(norm(result.reactions-result.supportReactions-result.mpcForces,inf) < 1e-7);
S = speye(chain.numberOfDOFs); AA = [S(cc.fixedDOFs,:); cc.C];
reference = [chain.stiffness -AA.';AA sparse(size(AA,1),size(AA,1))] \ [result.loadVector;zeros(size(AA,1),1)];
assert(norm(reference(chain.numberOfDOFs+numel(cc.fixedDOFs)+1:end)-result.mpcMultipliers,inf)<1e-7);
% A fully eliminated model is valid in statics and has no dynamic modes.
locked=m; locked.fixedBoundaryConditions=[1 1 0 0 0;1 2 0 0 0;3 3 0 0 0];
lockedResult=solveStatic(locked); assert(all(lockedResult.displacements==0));
expect('MKEF:NoFreeDOFs',@() solveModal(locked));
% Input validation and unsupported topology.
badModel = m; badModel.multiPointConstraints(2) = m.multiPointConstraints(1); expect('MKEF:InvalidMPC',@() buildConstraintTransform(badModel));
badModel = m; badModel.multiPointConstraints(2) = equation(2,2,[3 2 1]); expect('MKEF:InvalidMPC',@() buildConstraintTransform(badModel));
badModel = m; badModel.fixedBoundaryConditions = [1 1 0 0 0;6 3 0 0 0]; expect('MKEF:InvalidMPC',@() buildConstraintTransform(badModel));
badRecords = {equation(4,1,[1 1 1]),equation(3,4,[1 1 1]),equation(3,2,[3 2 1]),equation(3,2,[1 2 0]),equation(3,2,[1 2 1;1 2 2]),equation(3,2,[1 2 Inf])};
for i=1:numel(badRecords)
    badModel=m; badModel.multiPointConstraints=badRecords{i}; expect('MKEF:InvalidMPC',@() buildConstraintTransform(badModel));
end
badModel=m; badModel.multiPointConstraints.rhs=1; expect('MKEF:InvalidMPC',@() buildConstraintTransform(badModel));
% Export shape, selected DOFs and time decimation.
data = createPostprocessorData(m,s); assert(data.version == 3 && numel(data.model.mpcs) == 1);
data = createPostprocessorData(dynamic,t,struct('timeStride',3,'selectedGlobalDOFs',[2 8]));
assert(size(data.analysis.mpcForces,1) == 2);
assert(size(data.analysis.mpcMultipliers,1) == 1);
assert(size(data.analysis.mpcMultipliers,2) == numel(data.analysis.time));
filename = [tempname() '.json']; cleanup = onCleanup(@() delete(filename));
exportPostprocessorData(m,s,filename);
text = fileread(filename); assert(~isempty(strfind(text,'"mpcMultipliers":[')));
clear cleanup;
testRefinedTruss();
testParserFailures();
end

function testRefinedTruss()
filename = [tempname() '.txt']; cleanup = onCleanup(@() delete(filename));
previous = Inf;
for count = [1 2 4 8 16]
    fid = fopen(filename,'w');
    fprintf(fid,'analysis\nmodal\nnodes\n%d\n',count+1);
    for i=0:count, fprintf(fid,'%g,0,0\n',i/count); end
    fprintf(fid,'elems_112\n%d\n',count);
    for i=1:count, fprintf(fid,'112,%d,%d,1,100,1\n',i,i+1); end
    fprintf(fid,'bcfix\n2\n1,1,0,0,0\n2,%d,0,0,0\nmpc\n%d\n',count+1,count-1);
    for i=2:count
        xi=(i-1)/count;
        fprintf(fid,'%d,2,0,2,1,2,%.17g,%d,2,%.17g\n',i,1-xi,count+1,xi);
    end
    fclose(fid);
    p=StructFEProblem(filename); r=p.RunModal();
    assert(numel(r.frequenciesHz)==count && all(r.frequenciesHz>0));
    assert(r.frequenciesHz(1)<previous); previous=r.frequenciesHz(1);
end
% Fixed/free axial rod: f1 = sqrt(E/rho)/(4 L) = 2.5.
assert(abs(previous/2.5-1)<.001);
clear cleanup;
end

function testParserFailures()
base=fileread(exampleCasePath('CaseMPCFrame.txt'));
filename=[tempname() '.txt']; cleanup=onCleanup(@() delete(filename));
records={'3,2,0,2,1,2,0.5','3,2,0,1,1,2,NaN','3,2,0,0'};
for i=1:numel(records)
    fid=fopen(filename,'w'); fprintf(fid,'%s',strrep(base,'3,2,0,2,1,2,0.5,2,2,0.5',records{i})); fclose(fid);
    expect('MKEF:MalformedInput',@() FEMesh(filename));
end
% Repeated MPC sections append in source order.
fid=fopen(filename,'w'); fprintf(fid,'%s\nmpc\n1\n3,1,0,1,2,1,1\n',base); fclose(fid);
mesh=FEMesh(filename); assert(numel(mesh.multiPointConstraints)==2);
assert(mesh.multiPointConstraints(2).sourceLine>mesh.multiPointConstraints(1).sourceLine);
clear cleanup;
end

function m = equation(node,dof,masters)
m = struct('depNode',node,'depDOF',dof,'rhs',0,'masters',masters,'sourceLine',0);
end
function expect(id,operation)
try, operation(); catch e, assert(strcmp(e.identifier,id),e.message); return; end
error('Test:MissingError','Expected %s.',id);
end
