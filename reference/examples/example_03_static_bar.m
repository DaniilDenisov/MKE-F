function result = example_03_static_bar(verbose)
%EXAMPLE_03_STATIC_BAR Solve one axial bar and check the analytical result.
%   RESULT = EXAMPLE_03_STATIC_BAR() prints the reduced system, expands the
%   solution to all DOFs, and checks displacement, reaction, and equilibrium.
%   Pass false to run the example without diagnostic output.

if nargin < 1
    verbose = true;
end
validateVerbose(verbose);
addRepositoryRoot();

elementLength = 2.0;
area = 0.01;
youngsModulus = 200e9;
density = 7850;
appliedForce = 1000;

nodeCoordinates = [0 0; elementLength 0];
dofMap = [1 2; 3 4];
element = createStructuralElement(112, nodeCoordinates, [1 2], ...
    [area youngsModulus density], dofMap);
[globalK, globalM] = assembleGlobalMatrices(element, 4);

model = struct();
model.stiffness = globalK;
model.mass = globalM;
model.dofMap = dofMap;
model.fixedBoundaryConditions = [1 1 0 0 0; 2 2 0 0 0];
model.forceBoundaryConditions = [10 2 appliedForce 0 0 0];
model.nodeCoordinates = nodeCoordinates;
model.elementData = element;
model.numberOfNodes = 2;
model.dofPerNode = 2;
model.numberOfDOFs = 4;

analysis = solveStatic(model);
expectedDisplacement = appliedForce * elementLength / ...
    (area * youngsModulus);
expectedDisplacements = [0; 0; expectedDisplacement; 0];
expectedReactions = [-appliedForce; 0; 0; 0];

assert(isequal(analysis.fixedDOFs, [1 2 4]), ...
    'MKEF:ReferenceExampleFailed: unexpected fixed DOFs.');
assert(isequal(analysis.freeDOFs, 3), ...
    'MKEF:ReferenceExampleFailed: unexpected free DOFs.');
assertClose(analysis.displacements, expectedDisplacements, ...
    1e-12, 1e-15, 'The complete displacement vector is incorrect.');
assertClose(analysis.reactions, expectedReactions, ...
    1e-12, 1e-9, 'The support reaction is incorrect.');
assertClose(analysis.equilibriumResidual, zeros(3, 1), ...
    0, 1e-9, 'Global force or moment equilibrium failed.');
assertClose(globalK * analysis.displacements - ...
    analysis.loadVector - analysis.reactions, zeros(4, 1), ...
    0, 1e-9, 'The complete static system is not in equilibrium.');
assertClose(analysis.reactions(analysis.freeDOFs), 0, ...
    0, 1e-9, 'A free DOF contains a support reaction.');

freeDOFs = analysis.freeDOFs;
reducedK = globalK(freeDOFs, freeDOFs);
reducedLoad = analysis.loadVector(freeDOFs);
assertClose(reducedK * analysis.displacements(freeDOFs), reducedLoad, ...
    1e-12, 1e-9, 'The reduced static equation is not satisfied.');

if verbose
    fprintf('\nExample 03: static axial bar\n');
    fprintf('Global DOF map [ux uy]:\n');
    disp(dofMap);
    fprintf('Fixed DOFs: ');
    disp(analysis.fixedDOFs);
    fprintf('Free DOFs: ');
    disp(analysis.freeDOFs);
    fprintf('Reduced stiffness Kff [N/m]:\n');
    disp(full(reducedK));
    fprintf('Reduced load Ff [N]:\n');
    disp(reducedLoad);
    fprintf('Complete displacement vector [m]:\n');
    disp(analysis.displacements);
    fprintf('Reaction vector [N]:\n');
    disp(analysis.reactions);
    fprintf('Analytical displacement FL/(EA) = %.12g m.\n', ...
        expectedDisplacement);
    fprintf('Global equilibrium residual [Fx; Fy; Mz]:\n');
    disp(analysis.equilibriumResidual);
    fprintf('Checks: reduction, expansion, FL/(EA), reaction, and equilibrium PASS.\n');
end

result = analysis;
result.model = model;
result.reducedStiffness = reducedK;
result.reducedLoad = reducedLoad;
result.expectedDisplacement = expectedDisplacement;
end

function addRepositoryRoot()
exampleDirectory = fileparts(mfilename('fullpath'));
referenceDirectory = fileparts(exampleDirectory);
repositoryRoot = fileparts(referenceDirectory);
addpath(repositoryRoot);
setup();
end

function validateVerbose(verbose)
if ~isscalar(verbose) || ...
        ~(islogical(verbose) || (isnumeric(verbose) && isfinite(verbose))) || ...
        ~ismember(double(verbose), [0 1])
    error('MKEF:InvalidExampleOption', ...
        'verbose must be a scalar logical value.');
end
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
