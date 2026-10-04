function test_element_loads()
%TEST_ELEMENT_LOADS Independent beam solutions and load-work invariants.
p = StructFEProblem(exampleCasePath('CaseUniformFrame.txt'));
m = p.GetAnalysisModel();
r = solveStatic(m);
L = 2; E = 2e11; I = 1e-5; q = -1000;
near(r.reactions(1:3), [0; -q*L; -q*L^2/2]);
near(r.displacements(5:6), [q*L^4/(8*E*I); q*L^3/(6*E*I)]);
near(r.elementResults.localEndForces(4:6), zeros(3, 1));
[~, recovered] = StressCalc(m, r.displacements);
near(recovered.localEndForces, r.elementResults.localEndForces);
near(r.equilibriumResidual, zeros(3, 1));

% Simple supports: constrain translations, leave both rotations free.
s = m; s.fixedBoundaryConditions = [4 1 0 0 0; 4 2 0 0 0];
rs = solveStatic(s);
near(rs.reactions([2 5]), [-q*L/2; -q*L/2]);
near(rs.elementResults.localEndForces([3 6]), [0; 0]);
b = m; b.fixedBoundaryConditions = [1 1 0 0 0; 1 2 0 0 0];
rb = solveStatic(b);
near(rb.elementResults.localEndForces([3 6]), [-q*L^2/12; q*L^2/12]);

% Pure end moment and superposition, including multiple element records.
z = m; z.elementLoads = m.elementLoads([]);
z.forceBoundaryConditions = [10 2 0 0 120 0];
rz = solveStatic(z); near(rz.reactions(1:3), [0; 0; -120]);
s = m; s.forceBoundaryConditions = [10 2 7 19 120 0];
s.elementLoads(2) = s.elementLoads(1);
rs = solveStatic(s);
z.forceBoundaryConditions = s.forceBoundaryConditions;
rz = solveStatic(z);
near(rs.displacements, 2*r.displacements + rz.displacements);
s.elementLoads(2).qy = -s.elementLoads(1).qy;
near(buildStaticLoad(s), buildNodalLoadVector(s));

% Axial intensity gives the resultant and the mean axial force.
a = m; a.elementLoads.qx = 300; a.elementLoads.qy = 0;
ra = solveStatic(a); near(ra.reactions(1), -300*L);
near(ra.elementResults.axialForce, 300*L/2);

% Rotated/translated geometry, local/global equivalent records.
coords = [3 4; 4.2 5.6];
e = createStructuralElement(113, coords, [1 2], [0.01 E 7800 I], [1 2 3; 4 5 6]);
t = m; t.elementData = e; t.nodeCoordinates = coords;
t.stiffness = sparse(e.stiffness); t.mass = sparse(e.mass);
t.elementLoads.coordinateSystem = 1; t.elementLoads.qx = 70; t.elementLoads.qy = -90;
fLocal = buildStaticLoad(t);
globalQ = e.transformation(1:2, 1:2).' * [70; -90];
t.elementLoads.coordinateSystem = 2; t.elementLoads.qx = globalQ(1); t.elementLoads.qy = globalQ(2);
near(buildStaticLoad(t), fLocal);
[vectors, ~, resultant] = getElementLoadData(t);
near(calculateStaticEquilibrium(t, fLocal, zeros(6, 1)), resultant);
% Independent Simpson integration of load work against cubic Hermite shape.
u = [0.1; 0.2; -0.3; 0.4; -0.1; 0.6];
xis = [0 .5 1]; weights = [1 4 1]/6; work = 0;
for j = 1:3
    x = xis(j);
    axial = (1-x)*u(1) + x*u(4);
    transverse = (1-3*x^2+2*x^3)*u(2) + L*(x-2*x^2+x^3)*u(3) + ...
        (3*x^2-2*x^3)*u(5) + L*(-x^2+x^3)*u(6);
    work = work + L*weights(j)*(70*axial - 90*transverse);
end
near(vectors{1}.'*u, work);

mustFail(@() solveModal(m), 'MKEF:AnalysisLoadMismatch');
mustFail(@() solveTransient(m, struct()), 'MKEF:AnalysisLoadMismatch');
mustFail(@() buildTransientLoad(m, 0.1, 2), 'MKEF:AnalysisLoadMismatch');
bad = m; bad.elementLoads.coordinateSystem = 3;
mustFail(@() buildStaticLoad(bad), 'MKEF:InvalidElementLoad');
bad = m; bad.elementLoads.qy = 0;
mustFail(@() buildStaticLoad(bad), 'MKEF:InvalidElementLoad');
bad = m; bad.elementLoads.elementId = 3;
mustFail(@() buildStaticLoad(bad), 'MKEF:InvalidElementLoad');

data = createPostprocessorData(m, r);
assert(data.version == 2 && numel(data.model.elementLoads) == 1);
near(data.analysis.elementResults.equivalentLocalLoadVector(:), [0; -1000; -1000/3; 0; -1000; 1000/3]);
assert(createPostprocessorData(z, rz).version == 1);
filename = [tempname() '.json']; cleanup = onCleanup(@() delete(filename));
exportPostprocessorData(m, r, filename);
json = fileread(filename);
assert(~isempty(strfind(json, '"elementLoads":[')));
assert(~isempty(strfind(json, '"nodalLoads":[]')));
clear cleanup;

% Both parser and direct-kernel validation cover invalid input.
source = fileread(exampleCasePath('CaseUniformFrame.txt'));
record = '20,1,2,0,-1000';
invalid = {'20,1,2,,5', '20,1,2,0,NaN', '20,1,2,0,Inf', ...
    '20,1,2,0,0', '20,0,2,0,1', '20,1.5,2,0,1', '20,1,3,0,1', ...
    '21,1,2,0,1', '20,1,2,1'};
for i = 1:numel(invalid)
    parseText(strrep(source, record, invalid{i}), false);
end
parseText(strrep(source, sprintf('analysis\nstatic'), sprintf('analysis\nmodal')), false);
parseText(strrep(source, sprintf('analysis\nstatic'), sprintf('analysis\ntransient,0.01,0.1,2,1')), false);
truss = strrep(source, 'elems_113', 'elems_112');
truss = strrep(truss, '113,1,2,0.01,200000000000,7800,0.00001', '112,1,2,0.01,200000000000,7800');
parseText(truss, false);
repeat = parseText([source sprintf('\neload_uniform\n1\n20,1,1,2,0\n')], true);
assert(numel(repeat.elementLoads) == 2);
legacy = parseText(strrep(source, sprintf('analysis\nstatic\n'), ''), true);
legacyModel = createAnalysisModel(m.stiffness, m.mass, legacy);
near(solveStatic(legacyModel).displacements, r.displacements);
mustFail(@() solveModal(legacyModel), 'MKEF:AnalysisLoadMismatch');
end

function near(actual, expected)
assert(norm(actual(:)-expected(:), inf) <= 1e-8*max(1, norm(expected(:), inf)));
end
function mustFail(operation, identifier)
try
    operation();
catch e
    assert(strcmp(e.identifier, identifier), e.message); return;
end
error('Expected error %s.', identifier);
end

function mesh = parseText(text, valid)
filename = [tempname() '.txt'];
fid = fopen(filename, 'wt'); fprintf(fid, '%s', text); fclose(fid);
cleanup = onCleanup(@() delete(filename));
mesh = [];
try
    mesh = FEMesh(filename);
catch e
    if valid, rethrow(e); end
    assert(~isempty(strfind(e.message, filename)), 'Parser error must identify source file.');
    return;
end
assert(valid, 'Invalid element-load input was accepted.');
end
