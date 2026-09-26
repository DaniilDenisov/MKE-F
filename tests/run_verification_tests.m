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
runNamedTest('nodal load semantics', @testNodalLoadSemantics);
runNamedTest('functional analysis core', @testFunctionalAnalysisCore);
runNamedTest('free-DOF reduction', @testFreeDOFReduction);
runNamedTest('constraint validation', @testConstraintValidation);
runNamedTest('Newmark transient analysis', @test_newmark_transient);

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

result = problem.RunStatic();
loadVector = result.loadVector;
displacement = result.displacements;

expectedDisplacement = force * elementLength / (area * youngsModulus);
assertClose(displacement(3), expectedDisplacement, 1e-12, 1e-15, ...
    'Axial displacement does not match FL/(EA).');

reaction = result.reactions;
assertClose(reaction(1), -force, 1e-12, 1e-9, ...
    'Axial support reaction is incorrect.');
assertClose(reaction(1) + force, 0, 0, 1e-9, ...
    'Axial forces are not in equilibrium.');
assertClose(result.equilibriumResidual, zeros(3, 1), 0, 1e-9, ...
    'The axial bar is not in global equilibrium.');

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

result = problem.RunStatic();
loadVector = result.loadVector;
displacement = result.displacements;

expectedDeflection = force * elementLength^3 / ...
    (3 * youngsModulus * momentOfInertia);
expectedRotation = force * elementLength^2 / ...
    (2 * youngsModulus * momentOfInertia);
assertClose(displacement(5), expectedDeflection, 1e-11, 1e-14, ...
    'Cantilever tip deflection does not match PL^3/(3EI).');
assertClose(displacement(6), expectedRotation, 1e-11, 1e-14, ...
    'Cantilever tip rotation does not match PL^2/(2EI).');

reaction = result.reactions;
assertClose(reaction(2), -force, 1e-11, 1e-9, ...
    'Cantilever shear reaction is incorrect.');
assertClose(reaction(3), -force * elementLength, 1e-11, 1e-9, ...
    'Cantilever moment reaction is incorrect.');
assertClose(result.equilibriumResidual, zeros(3, 1), 0, 1e-9, ...
    'The cantilever is not in global equilibrium.');
end

function testNodalLoadSemantics()
options = quietOptions();
timeStep = 1e-4;
stepCount = 4;

% Рамный узел принимает [Fx, Fy, Mz]. Нагрузки в одном узле и в разных
% узлах должны суммироваться через общую карту степеней свободы.
staticProblem = StructFEProblem('Case1ElementBeam.txt', options);
staticModel = staticProblem.GetAnalysisModel();
staticModel.forceBoundaryConditions = ...
    [10, 1, -4,  6,  8; ...
     10, 2, 10, 20, 30; ...
     10, 2,  5, -2,  7];
staticLoads = buildStaticLoad(staticModel);
assertClose(staticLoads, [-4; 6; 8; 15; 18; 37], 0, 1e-12, ...
    'Static Fx, Fy, and Mz loads were not assembled by nodal DOF.');
staticResult = solveStatic(staticModel);
assertClose(staticResult.loadVector, staticLoads, 0, 1e-12, ...
    'Static analysis did not preserve the assembled nodal loads.');
assertClose(staticResult.equilibriumResidual, zeros(3, 1), 0, 1e-8, ...
    'Static nodal loads and reactions are not in equilibrium.');

% Нагрузка типа 10 в динамическом расчёте намеренно действует как прямоугольный
% импульс длительностью в один шаг для обратной совместимости.
pulseProblem = StructFEProblem('CaseBeamDyn.txt', options);
pulseLoads = buildTransientLoad(pulseProblem.GetAnalysisModel(), ...
    timeStep, stepCount);
pulseExpected = zeros(size(pulseLoads));
pulseDOF = pulseProblem.mesh.iMnod(3, 2);
pulseExpected(pulseDOF, 2) = -1000;
assertClose(pulseLoads, pulseExpected, 0, 1e-12, ...
    'Legacy type 10 transient load is not a one-step rectangular pulse.');

% Тип 12 и маркер bcforce_pulse задают тот же импульс явно, включая Mz.
explicitPulseProblem = StructFEProblem(fullfile('tests', 'fixtures', ...
    'CaseExplicitPulseBeam.txt'), options);
explicitPulseLoads = buildTransientLoad(...
    explicitPulseProblem.GetAnalysisModel(), timeStep, stepCount);
explicitPulseExpected = zeros(size(explicitPulseLoads));
explicitPulseDOFs = explicitPulseProblem.mesh.iMnod(2, :);
explicitPulseExpected(explicitPulseDOFs, 2) = [0; -1000; 25];
assertClose(explicitPulseLoads, explicitPulseExpected, 0, 1e-12, ...
    'Explicit type 12 pulse history is incorrect.');
explicitPulseProblem.ts = timeStep;
explicitPulseProblem.tsNum = stepCount;
explicitPulseProblem.ApplyForceBC();
assertClose(explicitPulseProblem.F, explicitPulseExpected, 0, 1e-12, ...
    'ApplyForceBC does not delegate to the complete nodal load builder.');

% Тип 13 и маркер bcforce_step сохраняют нагрузку на всех шагах.
stepProblem = StructFEProblem(fullfile('tests', 'fixtures', ...
    'CaseStepBeam.txt'), options);
stepModel = stepProblem.GetAnalysisModel();
stepLoads = buildTransientLoad(stepModel, timeStep, stepCount);
stepExpected = zeros(size(stepLoads));
stepDOFs = stepProblem.mesh.iMnod(2, :);
stepExpected(stepDOFs, :) = ...
    repmat([0; -1000; 25], 1, stepCount + 1);
assertClose(stepLoads, stepExpected, 0, 1e-12, ...
    'Type 13 persistent step history is incorrect.');
stepResult = stepProblem.RunTransient(timeStep, ...
    stepCount * timeStep, 2, 2);
assertClose(stepResult.loadHistory, stepExpected, 0, 1e-12, ...
    'RunTransient did not use the persistent step history.');
assertThrows('MKEF:TimeDependentLoadInStaticAnalysis', ...
    @() buildStaticLoad(stepModel));

% Для типа 11 все компоненты, включая Mz, равны F0*sin(2*pi*f*t).
harmonicProblem = StructFEProblem('CaseBeamFreq.txt', options);
harmonicModel = harmonicProblem.GetAnalysisModel();
harmonicModel.forceBoundaryConditions(1, 5) = 25;
harmonicLoads = buildTransientLoad(harmonicModel, ...
    timeStep, stepCount);
harmonicExpected = zeros(size(harmonicLoads));
harmonicDOFs = harmonicProblem.mesh.iMnod(3, :);
time = (0:stepCount) * timeStep;
harmonicExpected(harmonicDOFs, :) = ...
    [0; -1000; 25] * sin(2*pi*135*time);
assertClose(harmonicLoads, harmonicExpected, 1e-12, 1e-12, ...
    'Type 11 harmonic load history is incorrect.');

% У ферменного узла нет вращательной СС, поэтому ненулевой Mz отклоняется,
% а не игнорируется как несуществующая сила Fz.
trussProblem = StructFEProblem(fullfile('tests', 'fixtures', ...
    'CaseSingleTruss.txt'), options);
trussModel = trussProblem.GetAnalysisModel();
trussModel.forceBoundaryConditions(1, 5) = 1;
assertThrows('MKEF:UnsupportedLoadComponent', ...
    @() buildStaticLoad(trussModel));
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

function testFreeDOFReduction()
options = quietOptions();
problem = StructFEProblem('CaseBeam.txt', options);
originalK = problem.K;
originalM = problem.M;
model = problem.GetAnalysisModel();

[fixedDOFs, freeDOFs] = partitionDOFs(model);
assert(isequal(fixedDOFs, 1:5));
assert(isequal(freeDOFs, 6:9));

% Совместимый метод теперь только возвращает разбиение и не редактирует матрицы.
[facadeFixedDOFs, facadeFreeDOFs] = problem.ApplyFixBC();
assert(isequal(facadeFixedDOFs, fixedDOFs));
assert(isequal(facadeFreeDOFs, freeDOFs));
assert(isequal(problem.K, originalK));
assert(isequal(problem.M, originalM));

modalResult = problem.RunModal();
assert(numel(modalResult.frequenciesHz) == numel(freeDOFs));
assert(all(modalResult.frequenciesHz > 0));
assert(all(all(modalResult.modeShapes(fixedDOFs, :) == 0)));

transientResult = problem.RunTransient(1e-4, 4e-4, 3, 2);
assert(all(all(transientResult.displacements(fixedDOFs, :) == 0)));
assert(all(all(transientResult.velocities(fixedDOFs, :) == 0)));
assert(all(all(transientResult.accelerations(fixedDOFs, :) == 0)));

% Нагрузка на заделанную СС не входит в редуцированную систему, но должна
% сохраниться в полном векторе нагрузки и учитываться в реакции.
cantilever = StructFEProblem('Case1ElementBeam.txt', options);
loadedSupportModel = cantilever.GetAnalysisModel();
loadedSupportModel.forceBoundaryConditions(end + 1, :) = ...
    [10, 1, 25, -30, 0];
staticResult = solveStatic(loadedSupportModel);
assertClose(staticResult.loadVector(1:2), [25; -30], 0, 1e-12, ...
    'Loads on restrained DOFs were not preserved.');
assert(all(staticResult.displacements(staticResult.fixedDOFs) == 0));
assertClose(staticResult.reactions(staticResult.freeDOFs), ...
    zeros(numel(staticResult.freeDOFs), 1), 0, 1e-9, ...
    'Free DOFs contain non-zero reactions.');
assertClose(staticResult.equilibriumResidual, zeros(3, 1), ...
    0, 1e-9, 'Static loads and reactions are not in equilibrium.');
end

function testConstraintValidation()
options = quietOptions();
problem = StructFEProblem(fullfile('tests', 'fixtures', ...
    'CaseSingleTruss.txt'), options);
model = problem.GetAnalysisModel();

duplicateModel = model;
duplicateModel.fixedBoundaryConditions(end + 1, :) = [2, 1, 0, 0, 0];
assertThrows('MKEF:DuplicateConstraint', ...
    @() partitionDOFs(duplicateModel));

invalidTypeModel = model;
invalidTypeModel.fixedBoundaryConditions(1, 1) = 99;
assertThrows('MKEF:InvalidConstraint', ...
    @() partitionDOFs(invalidTypeModel));

invalidNodeModel = model;
invalidNodeModel.fixedBoundaryConditions(1, 2) = 3;
assertThrows('MKEF:InvalidConstraint', ...
    @() partitionDOFs(invalidNodeModel));

mechanismModel = model;
mechanismModel.fixedBoundaryConditions = [1, 1, 0, 0, 0];
assertThrows('MKEF:SingularStiffness', ...
    @() solveStatic(mechanismModel));

fullyFixedModel = model;
fullyFixedModel.fixedBoundaryConditions = ...
    [1, 1, 0, 0, 0; 1, 2, 0, 0, 0];
assertThrows('MKEF:NoFreeDOFs', @() solveModal(fullyFixedModel));
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

function assertThrows(expectedIdentifier, operation)
try
    operation();
catch exception
    if strcmp(exception.identifier, expectedIdentifier)
        return;
    end
    error('MKEF:VerificationFailed', ...
        'Expected error %s, received %s.', ...
        expectedIdentifier, exception.identifier);
end
error('MKEF:VerificationFailed', 'Expected error %s.', expectedIdentifier);
end
