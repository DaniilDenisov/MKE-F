function run_verification_tests()
%RUN_VERIFICATION_TESTS Self-contained numerical checks for GNU Octave.

testsDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(testsDir);
previousDir = pwd;
cleanup = onCleanup(@() cd(previousDir));
cd(rootDir);
addpath(rootDir);
addpath(testsDir);

requireOctave();
fprintf('Running self-contained verification tests...\n');

runNamedTest('assembled matrix invariants', @testMatrixInvariants);
runNamedTest('single axial truss', @testSingleAxialTruss);
runNamedTest('cantilever beam stiffness', @testCantileverBeam);
runNamedTest('transient load histories', @testTransientLoadHistories);

fprintf('All verification tests passed.\n');
clear cleanup;
end

function testMatrixInvariants()
options = quietOptions();
caseFiles = {
    'ANSYSBeamStatic01.txt'
    'Case1ElementBeam.txt'
    'CaseATransSite.txt'
    'CaseBeam.txt'
    'CaseBeamDyn.txt'
    'CaseBeamFreq.txt'
    'CaseFig11.7p363 MarioPaz.txt'
};

for i = 1:numel(caseFiles)
    problem = StructFEProblem(caseFiles{i}, options);
    expectedDOFs = problem.mesh.numberOfNodes * problem.mesh.dofPerNode;
    assert(isequal(size(problem.K), [expectedDOFs expectedDOFs]));
    assert(isequal(size(problem.M), [expectedDOFs expectedDOFs]));
    assert(all(isfinite(problem.K(:))));
    assert(all(isfinite(problem.M(:))));
    assertRelativeSmall(problem.K - problem.K.', problem.K, 1e-12, ...
        ['Stiffness matrix is not symmetric: ' caseFiles{i}]);
end

% Beam mass symmetry is a known defect scheduled for Commit 3.
trussProblem = StructFEProblem('CaseATransSite.txt', options);
assertRelativeSmall(trussProblem.M - trussProblem.M.', trussProblem.M, ...
    1e-12, 'Truss mass matrix is not symmetric.');
end

function testSingleAxialTruss()
options = quietOptions();
caseFile = fullfile('tests', 'fixtures', 'CaseSingleTruss.txt');
problem = StructFEProblem(caseFile, options);

elementLength = 2.0;
area = 0.01;
youngsModulus = 2e11;
density = 7850;
force = 1000;

loadVector = zeros(4, 1);
loadVector(3) = force;
displacement = zeros(4, 1);
displacement(3) = problem.K(3, 3) \ loadVector(3);

expectedDisplacement = force * elementLength / (area * youngsModulus);
assertClose(displacement(3), expectedDisplacement, 1e-12, 1e-15, ...
    'Axial displacement does not match FL/(EA).');

reaction = problem.K * displacement - loadVector;
assertClose(reaction(1), -force, 1e-12, 1e-9, ...
    'Axial support reaction is incorrect.');
assertClose(reaction(1) + force, 0, 0, 1e-9, ...
    'Axial forces are not in equilibrium.');

xDOFs = [1 3];
assembledMassX = sum(sum(problem.M(xDOFs, xDOFs)));
expectedMass = density * area * elementLength;
assertClose(assembledMassX, expectedMass, 1e-12, 1e-12, ...
    'Truss translational mass is incorrect.');
end

function testCantileverBeam()
options = quietOptions();
problem = StructFEProblem('Case1ElementBeam.txt', options);

elementLength = 0.5;
youngsModulus = 2e11;
momentOfInertia = 0.33e-8;
force = 100;

loadVector = zeros(6, 1);
loadVector(5) = force;
freeDOFs = 4:6;
displacement = zeros(6, 1);
displacement(freeDOFs) = ...
    problem.K(freeDOFs, freeDOFs) \ loadVector(freeDOFs);

expectedDeflection = force * elementLength^3 / ...
    (3 * youngsModulus * momentOfInertia);
expectedRotation = force * elementLength^2 / ...
    (2 * youngsModulus * momentOfInertia);
assertClose(displacement(5), expectedDeflection, 1e-11, 1e-14, ...
    'Cantilever tip deflection does not match PL^3/(3EI).');
assertClose(displacement(6), expectedRotation, 1e-11, 1e-14, ...
    'Cantilever tip rotation does not match PL^2/(2EI).');

reaction = problem.K * displacement - loadVector;
assertClose(reaction(2), -force, 1e-11, 1e-9, ...
    'Cantilever shear reaction is incorrect.');
assertClose(reaction(3), -force * elementLength, 1e-11, 1e-9, ...
    'Cantilever moment reaction is incorrect.');
end

function testTransientLoadHistories()
options = quietOptions();
timeStep = 1e-4;
stepCount = 4;

pulseProblem = StructFEProblem('CaseBeamDyn.txt', options);
pulseProblem.ts = timeStep;
pulseProblem.tsNum = stepCount;
pulseProblem.F = zeros(size(pulseProblem.K, 1), stepCount);
pulseProblem.ApplyForceBC();
pulseExpected = zeros(size(pulseProblem.F));
pulseDOF = pulseProblem.mesh.iMnod(3, 2);
pulseExpected(pulseDOF, 1) = -1000;
assertClose(pulseProblem.F, pulseExpected, 0, 1e-12, ...
    'Type 10 transient load is not a one-step rectangular pulse.');

harmonicProblem = StructFEProblem('CaseBeamFreq.txt', options);
harmonicProblem.ts = timeStep;
harmonicProblem.tsNum = stepCount;
harmonicProblem.F = zeros(size(harmonicProblem.K, 1), stepCount);
harmonicProblem.ApplyForceBC();
harmonicExpected = zeros(size(harmonicProblem.F));
harmonicDOF = harmonicProblem.mesh.iMnod(3, 2);
time = (1:stepCount) * timeStep;
harmonicExpected(harmonicDOF, :) = -1000 * sin(2*pi*135*time);
assertClose(harmonicProblem.F, harmonicExpected, 1e-12, 1e-12, ...
    'Type 11 harmonic load history is incorrect.');
end

function options = quietOptions()
options = struct('verbose', false, 'plotting', false);
end

function requireOctave()
if ~exist('OCTAVE_VERSION', 'builtin')
    error('MKEF:OctaveRequired', ...
        'GNU Octave is the only supported runtime for this test suite.');
end
end

function runNamedTest(name, testFunction)
fprintf('  %-36s', [name ':']);
testFunction();
fprintf(' PASS\n');
end

function assertRelativeSmall(value, reference, relativeTolerance, message)
valueNorm = norm(value(:), inf);
referenceNorm = max(norm(reference(:), inf), 1);
if valueNorm > relativeTolerance * referenceNorm
    error('MKEF:VerificationFailed', '%s Error=%g, scale=%g.', ...
        message, valueNorm, referenceNorm);
end
end

function assertClose(actual, expected, relativeTolerance, absoluteTolerance, message)
difference = max(abs(actual(:) - expected(:)));
scale = max(abs(expected(:)));
if isempty(difference)
    difference = 0;
end
if isempty(scale)
    scale = 0;
end
limit = absoluteTolerance + relativeTolerance * scale;
if difference > limit
    error('MKEF:VerificationFailed', '%s Error=%g, tolerance=%g.', ...
        message, difference, limit);
end
end
