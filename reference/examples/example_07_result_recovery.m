function result = example_07_result_recovery(verbose)
%EXAMPLE_07_RESULT_RECOVERY Recover truss and frame element results.
%   RESULT = EXAMPLE_07_RESULT_RECOVERY() checks truss strain, stress, and
%   axial force, then checks frame local end forces against analytical
%   values, support reactions, and the applied nodal load. Pass false to
%   run the example without diagnostic output.

if nargin < 1
    verbose = true;
end
validateVerbose(verbose);
repositoryRoot = addRepositoryRoot();

options = struct('verbose', false, 'plotting', false);

trussProblem = StructFEProblem(fullfile(repositoryRoot, 'tests', ...
    'fixtures', 'CaseSingleTruss.txt'), options);
trussModel = trussProblem.GetAnalysisModel();
trussAnalysis = solveStatic(trussModel);
truss = trussAnalysis.elementResults(1);

trussForce = 1000;
trussArea = 0.01;
trussYoungsModulus = 2e11;
expectedTrussStress = trussForce / trussArea;
expectedTrussStrain = expectedTrussStress / trussYoungsModulus;
expectedTrussEndForces = [-trussForce; 0; trussForce; 0];

assert(truss.type == 112, ...
    'MKEF:ReferenceExampleFailed: the truss has an unexpected type.');
assertClose(truss.localDisplacements, ...
    trussModel.elementData(1).transformation * ...
    trussAnalysis.displacements(truss.globalDOFs), ...
    1e-12, 1e-15, 'The truss displacement transformation is incorrect.');
assertClose(truss.localEndForces, expectedTrussEndForces, ...
    1e-12, 1e-9, 'The truss local end forces are incorrect.');
assertClose(truss.axialStrain, expectedTrussStrain, ...
    1e-12, 1e-15, 'The truss axial strain is incorrect.');
assertClose(truss.axialStress, expectedTrussStress, ...
    1e-12, 1e-9, 'The truss axial stress is incorrect.');
assertClose(truss.axialForce, trussForce, ...
    1e-12, 1e-9, 'The truss axial force is incorrect.');
assertClose(truss.axialForce, trussArea * truss.axialStress, ...
    1e-12, 1e-9, 'The truss force and stress are inconsistent.');
assertClose(truss.localEndForces(1), trussAnalysis.reactions(1), ...
    1e-12, 1e-9, 'The truss support-end force disagrees with the reaction.');
assertClose(truss.localEndForces(3), trussAnalysis.loadVector(3), ...
    1e-12, 1e-9, 'The truss free-end force disagrees with the nodal load.');

frameProblem = StructFEProblem(fullfile(repositoryRoot, ...
    'examples', 'cases', 'Case1ElementBeam.txt'), options);
frameModel = frameProblem.GetAnalysisModel();
frameAnalysis = solveStatic(frameModel);
frame = frameAnalysis.elementResults(1);

frameForce = 100;
frameLength = 0.5;
expectedFrameEndForces = ...
    [0; -frameForce; -frameForce * frameLength; 0; frameForce; 0];

assert(frame.type == 113, ...
    'MKEF:ReferenceExampleFailed: the frame has an unexpected type.');
assertClose(frame.localDisplacements, ...
    frameModel.elementData(1).transformation * ...
    frameAnalysis.displacements(frame.globalDOFs), ...
    1e-12, 1e-15, 'The frame displacement transformation is incorrect.');
assertClose(frame.localEndForces, expectedFrameEndForces, ...
    1e-11, 1e-8, 'The frame local end forces are incorrect.');
assertClose(frame.localEndForces(1:3), frameAnalysis.reactions(1:3), ...
    1e-11, 1e-8, 'The frame support-end forces disagree with reactions.');
assertClose(frame.localEndForces(4:6), frameAnalysis.loadVector(4:6), ...
    1e-11, 1e-8, 'The frame free-end forces disagree with the nodal load.');
assertClose(frame.axialForce, 0, 0, 1e-9, ...
    'The transversely loaded frame has a spurious axial force.');

if verbose
    fprintf('\nExample 07: element-result recovery\n');
    fprintf('Truss local displacements [u1; v1; u2; v2]:\n');
    disp(truss.localDisplacements);
    fprintf('Truss local end forces [N1; V1; N2; V2] [N]:\n');
    disp(truss.localEndForces);
    fprintf('Truss strain = %.12g, stress = %.12g Pa, axial force = %.12g N.\n', ...
        truss.axialStrain, truss.axialStress, truss.axialForce);
    fprintf('Frame local end forces [N1; V1; M1; N2; V2; M2]:\n');
    disp(frame.localEndForces);
    fprintf('Frame support reactions [Fx; Fy; Mz]:\n');
    disp(frameAnalysis.reactions(1:3));
    fprintf(['Checks: DOF gathering, local transformation, truss stress, ' ...
        'axial force, frame end forces, reactions, and analytical values PASS.\n']);
end

result = struct();
result.truss = trussAnalysis;
result.trussModel = trussModel;
result.frame = frameAnalysis;
result.frameModel = frameModel;
result.expectedTrussStrain = expectedTrussStrain;
result.expectedTrussStress = expectedTrussStress;
result.expectedTrussAxialForce = trussForce;
result.expectedFrameEndForces = expectedFrameEndForces;
end

function repositoryRoot = addRepositoryRoot()
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
