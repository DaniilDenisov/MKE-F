function run_verification_tests()
%RUN_VERIFICATION_TESTS Самодостаточные численные проверки для GNU Octave.
% Используемые формулы и инварианты не требуют MATLAB, ANSYS или доступа
% к учебникам, упомянутым в примерах.

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
runNamedTest('element matrices', @test_element_matrices);
runNamedTest('single axial truss', @testSingleAxialTruss);
runNamedTest('cantilever beam stiffness', @testCantileverBeam);
runNamedTest('transient load histories', @testTransientLoadHistories);
runNamedTest('functional analysis core', @testFunctionalAnalysisCore);

fprintf('All verification tests passed.\n');
clear cleanup;
end

function testMatrixInvariants()
options = quietOptions();

% Загружаем все старые примеры, чтобы проверить чтение входных данных и глобальную
% сборку матриц рамных и ферменных моделей. Эти тесты не воспроизводят результаты примеров.
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

% Ошибка согласованной матрицы масс балки, исправленная в коммите 496569, проявлялась
% на глобальном уровне как нарушение симметрии. Поэтому явно проверяем каждую M.
for i = 1:numel(caseFiles)
    problem = StructFEProblem(caseFiles{i}, options);
    assertRelativeSmall(problem.M - problem.M.', problem.M, 1e-12, ...
        ['Mass matrix is not symmetric: ' caseFiles{i}]);
end
end

function testSingleAxialTruss()
options = quietOptions();
caseFile = fullfile('tests', 'fixtures', 'CaseSingleTruss.txt');
problem = StructFEProblem(caseFile, options);

% Для горизонтального стержня из одного элемента известно точное перемещение
% u = FL/(EA). Также проверяем равновесие реакций и собранную матрицу масс.
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

% Теория Эйлера-Бернулли даёт точные прогиб и угол поворота конца этой консоли
% из одного элемента под действием поперечной узловой силы.
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

% Нагрузка типа 10 в динамическом расчёте намеренно действует как прямоугольный
% импульс длительностью в один шаг, а не как постоянная ступенчатая нагрузка.
pulseProblem = StructFEProblem('CaseBeamDyn.txt', options);
pulseLoads = buildTransientLoad(pulseProblem.GetAnalysisModel(), ...
    timeStep, stepCount);
pulseExpected = zeros(size(pulseLoads));
pulseDOF = pulseProblem.mesh.iMnod(3, 2);
pulseExpected(pulseDOF, 1) = -1000;
assertClose(pulseLoads, pulseExpected, 0, 1e-12, ...
    'Type 10 transient load is not a one-step rectangular pulse.');

% Для типа 11 в столбцах времени должны вычисляться значения F0*sin(2*pi*f*t).
harmonicProblem = StructFEProblem('CaseBeamFreq.txt', options);
harmonicLoads = buildTransientLoad(harmonicProblem.GetAnalysisModel(), ...
    timeStep, stepCount);
harmonicExpected = zeros(size(harmonicLoads));
harmonicDOF = harmonicProblem.mesh.iMnod(3, 2);
time = (1:stepCount) * timeStep;
harmonicExpected(harmonicDOF, :) = -1000 * sin(2*pi*135*time);
assertClose(harmonicLoads, harmonicExpected, 1e-12, 1e-12, ...
    'Type 11 harmonic load history is incorrect.');
end

function testFunctionalAnalysisCore()
options = quietOptions();
problem = StructFEProblem('CaseBeamDyn.txt', options);
originalK = problem.K;
originalM = problem.M;
sentinelF = (1:size(problem.K, 1)).';
problem.F = sentinelF;

% Численное ядро принимает только структуру данных модели и возвращает
% именованные результаты, не изменяя фасад или переданную структуру.
model = problem.GetAnalysisModel();
modelBefore = model;
directStatic = solveStatic(model);
directModal = solveModal(model);
directTransient = solveTransient(model, ...
    struct('timeStep', 1e-4, 'duration', 4e-4));
assert(isequal(model, modelBefore));
assert(isfield(directStatic, 'displacements'));
assert(isfield(directStatic, 'reactions'));
assert(isfield(directStatic, 'loadVector'));
assert(isfield(directModal, 'frequenciesHz'));
assert(isfield(directModal, 'modeShapes'));
assert(isfield(directTransient, 'time'));
assert(isfield(directTransient, 'displacements'));
assert(isfield(directTransient, 'velocities'));
assert(isfield(directTransient, 'accelerations'));
assert(isfield(directTransient, 'loadHistory'));

% Повторные вызовы и другой порядок видов анализа должны давать те же
% результаты. Поля K, M и старое совместимое поле F остаются нетронутыми.
staticFirst = problem.RunStatic();
modalFirst = problem.RunModal();
transientFirst = problem.RunTransient(1e-4, 4e-4, 3, 2);

transientSecond = problem.RunTransient(1e-4, 4e-4, 3, 2);
staticSecond = problem.RunStatic();
modalSecond = problem.RunModal();

assertClose(staticSecond.displacements, staticFirst.displacements, ...
    1e-12, 1e-14, 'Static result depends on call order.');
assertClose(staticSecond.reactions, staticFirst.reactions, ...
    1e-12, 1e-10, 'Static reactions depend on call order.');
assertClose(modalSecond.frequenciesHz, modalFirst.frequenciesHz, ...
    1e-12, 1e-12, 'Modal result depends on call order.');
assertClose(transientSecond.displacements, transientFirst.displacements, ...
    1e-12, 1e-14, 'Transient result depends on call order.');
assertClose(transientSecond.loadHistory, transientFirst.loadHistory, ...
    0, 1e-12, 'Transient loads depend on call order.');

assertClose(staticFirst.displacements, directStatic.displacements, ...
    1e-12, 1e-14, 'Static facade and core results differ.');
assertClose(modalFirst.frequenciesHz, directModal.frequenciesHz, ...
    1e-12, 1e-12, 'Modal facade and core results differ.');
assertClose(transientFirst.displacements, directTransient.displacements, ...
    1e-12, 1e-14, 'Transient facade and core results differ.');

assert(isequal(problem.K, originalK));
assert(isequal(problem.M, originalM));
assert(isequal(problem.F, sentinelF));
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
% Масштабируем ошибку инварианта по величине проверяемой матрицы.
valueNorm = norm(value(:), inf);
referenceNorm = max(norm(reference(:), inf), 1);
if valueNorm > relativeTolerance * referenceNorm
    error('MKEF:VerificationFailed', '%s Error=%g, scale=%g.', ...
        message, valueNorm, referenceNorm);
end
end

function assertClose(actual, expected, relativeTolerance, absoluteTolerance, message)
% Около нуля используем абсолютный допуск, для больших значений — относительный.
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
