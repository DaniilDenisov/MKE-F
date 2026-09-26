function result = example_04_load_histories(verbose)
%EXAMPLE_04_LOAD_HISTORIES Build static, pulse, step, and harmonic loads.
%   RESULT = EXAMPLE_04_LOAD_HISTORIES() prints exact nodal load histories
%   on a short time grid and checks their semantics. Pass false for quiet use.

if nargin < 1
    verbose = true;
end
validateVerbose(verbose);
addRepositoryRoot();

model = createOneElementFrameModel();
loadedNode = 2;
loadedDOFs = model.dofMap(loadedNode, :);
components = [10; -20; 3];
timeStep = 0.25;
stepCount = 4;
time = (0:stepCount) * timeStep;

staticModel = model;
staticModel.forceBoundaryConditions = ...
    [10 loadedNode 10 -20  3 0; ...
     10 loadedNode  5   2 -1 0];
staticLoad = buildStaticLoad(staticModel);
expectedStaticLoad = zeros(model.numberOfDOFs, 1);
expectedStaticLoad(loadedDOFs) = [15; -18; 2];
assertClose(staticLoad, expectedStaticLoad, 0, 1e-12, ...
    'Multiple static nodal loads were not summed by component.');

pulseModel = model;
pulseModel.forceBoundaryConditions = ...
    [12 loadedNode components.' 0];
pulseHistory = buildTransientLoad(pulseModel, timeStep, stepCount);
expectedPulseHistory = zeros(model.numberOfDOFs, stepCount + 1);
expectedPulseHistory(loadedDOFs, 2) = components;
assertClose(pulseHistory, expectedPulseHistory, 0, 1e-12, ...
    'The explicit one-step pulse is incorrect.');

legacyPulseModel = model;
legacyPulseModel.forceBoundaryConditions = ...
    [10 loadedNode components.' 0];
legacyPulseHistory = buildTransientLoad( ...
    legacyPulseModel, timeStep, stepCount);
assertClose(legacyPulseHistory, expectedPulseHistory, 0, 1e-12, ...
    'The backward-compatible type 10 pulse is incorrect.');

stepModel = model;
stepModel.forceBoundaryConditions = ...
    [13 loadedNode components.' 0];
stepHistory = buildTransientLoad(stepModel, timeStep, stepCount);
expectedStepHistory = zeros(model.numberOfDOFs, stepCount + 1);
expectedStepHistory(loadedDOFs, :) = repmat(components, 1, stepCount + 1);
assertClose(stepHistory, expectedStepHistory, 0, 1e-12, ...
    'The persistent step history is incorrect.');

harmonicModel = model;
harmonicFrequencyHz = 1;
harmonicModel.forceBoundaryConditions = ...
    [11 loadedNode components.' harmonicFrequencyHz];
harmonicHistory = buildTransientLoad( ...
    harmonicModel, timeStep, stepCount);
expectedHarmonicHistory = zeros(model.numberOfDOFs, stepCount + 1);
expectedHarmonicHistory(loadedDOFs, :) = ...
    components * [0 1 0 -1 0];
assertClose(harmonicHistory, expectedHarmonicHistory, ...
    0, 1e-12, 'The harmonic load samples are incorrect.');

trussModel = createOneElementTrussModel();
trussModel.forceBoundaryConditions = [10 2 0 -20 3 0];
assertThrows('MKEF:UnsupportedLoadComponent', ...
    @() buildStaticLoad(trussModel));

if verbose
    fprintf('\nExample 04: nodal load histories\n');
    fprintf('Time grid [s]:\n');
    disp(time);
    fprintf('Loaded frame DOFs [Fx Fy Mz] -> [ux uy thetaZ]:\n');
    disp(loadedDOFs);
    fprintf('Summed static load vector:\n');
    disp(staticLoad);
    fprintf('Explicit pulse at loaded DOFs (nonzero only at t=dt):\n');
    disp(pulseHistory(loadedDOFs, :));
    fprintf('Persistent step at loaded DOFs (active from t=0):\n');
    disp(stepHistory(loadedDOFs, :));
    fprintf('One-hertz harmonic load at loaded DOFs:\n');
    disp(harmonicHistory(loadedDOFs, :));
    fprintf('Discrete pulse impulse components F0*dt:\n');
    disp(components * timeStep);
    fprintf('Checks: summation, pulse, step, harmonic, and truss Mz rejection PASS.\n');
end

result = struct();
result.time = time;
result.loadedDOFs = loadedDOFs;
result.staticLoad = staticLoad;
result.pulseHistory = pulseHistory;
result.legacyPulseHistory = legacyPulseHistory;
result.stepHistory = stepHistory;
result.harmonicHistory = harmonicHistory;
result.pulseImpulse = components * timeStep;
end

function model = createOneElementFrameModel()
nodeCoordinates = [0 0; 1 0];
dofMap = [1 2 3; 4 5 6];
element = createStructuralElement(113, nodeCoordinates, [1 2], ...
    [0.0004 200e9 7800 5.33e-8], dofMap);
[globalK, globalM] = assembleGlobalMatrices(element, 6);
model = commonModel(globalK, globalM, dofMap, nodeCoordinates, element, 3);
end

function model = createOneElementTrussModel()
nodeCoordinates = [0 0; 1 0];
dofMap = [1 2; 3 4];
element = createStructuralElement(112, nodeCoordinates, [1 2], ...
    [0.01 200e9 7850], dofMap);
[globalK, globalM] = assembleGlobalMatrices(element, 4);
model = commonModel(globalK, globalM, dofMap, nodeCoordinates, element, 2);
end

function model = commonModel(globalK, globalM, dofMap, ...
        nodeCoordinates, element, dofPerNode)
model = struct();
model.stiffness = globalK;
model.mass = globalM;
model.dofMap = dofMap;
model.fixedBoundaryConditions = zeros(0, 5);
model.forceBoundaryConditions = zeros(0, 6);
model.nodeCoordinates = nodeCoordinates;
model.elementData = element;
model.numberOfNodes = 2;
model.dofPerNode = dofPerNode;
model.numberOfDOFs = 2 * dofPerNode;
end

function addRepositoryRoot()
exampleDirectory = fileparts(mfilename('fullpath'));
referenceDirectory = fileparts(exampleDirectory);
repositoryRoot = fileparts(referenceDirectory);
addpath(repositoryRoot);
end

function validateVerbose(verbose)
if ~isscalar(verbose) || ...
        ~(islogical(verbose) || (isnumeric(verbose) && isfinite(verbose))) || ...
        ~ismember(double(verbose), [0 1])
    error('MKEF:InvalidExampleOption', ...
        'verbose must be a scalar logical value.');
end
end

function assertThrows(expectedIdentifier, operation)
try
    operation();
catch exception
    if strcmp(exception.identifier, expectedIdentifier)
        return;
    end
    error('MKEF:ReferenceExampleFailed', ...
        'Expected error %s, received %s.', ...
        expectedIdentifier, exception.identifier);
end
error('MKEF:ReferenceExampleFailed', ...
    'Expected error %s.', expectedIdentifier);
end

function assertClose(actual, expected, relativeTolerance, ...
        absoluteTolerance, message)
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
    error('MKEF:ReferenceExampleFailed', ...
        '%s Error=%g, tolerance=%g.', message, difference, limit);
end
end
