function result = example_05_modal_frame(verbose)
%EXAMPLE_05_MODAL_FRAME Solve a frame eigenproblem and verify its modes.
%   RESULT = EXAMPLE_05_MODAL_FRAME() prints the reduced modal system,
%   frequencies, normalized modes, orthogonality, and residuals. Pass false
%   to run the example without diagnostic output.

if nargin < 1
    verbose = true;
end
validateVerbose(verbose);
repositoryRoot = addRepositoryRoot();

options = struct('verbose', false, 'plotting', false);
problem = StructFEProblem(fullfile(repositoryRoot, 'examples', 'cases', ...
    'CaseBeam.txt'), options);
model = problem.GetAnalysisModel();
analysis = solveModal(model);

freeDOFs = analysis.freeDOFs;
fixedDOFs = analysis.fixedDOFs;
reducedK = model.stiffness(freeDOFs, freeDOFs);
reducedM = model.mass(freeDOFs, freeDOFs);
reducedModes = analysis.modeShapes(freeDOFs, :);
eigenvalues = analysis.eigenvalues(:);
modeCount = numel(eigenvalues);

massNormalizedModes = zeros(size(reducedModes));
displayModes = zeros(size(analysis.modeShapes));
relativeResiduals = zeros(modeCount, 1);

for modeNumber = 1:modeCount
    mode = reducedModes(:, modeNumber);
    modalMass = mode.' * reducedM * mode;
    assert(modalMass > 0, ...
        'MKEF:ReferenceExampleFailed: modal mass must be positive.');
    massNormalizedModes(:, modeNumber) = mode / sqrt(modalMass);

    fullMode = analysis.modeShapes(:, modeNumber);
    [maximumComponent, maximumIndex] = max(abs(fullMode));
    assert(maximumComponent > 0, ...
        'MKEF:ReferenceExampleFailed: a mode shape is identically zero.');
    signMultiplier = sign(fullMode(maximumIndex));
    displayModes(:, modeNumber) = ...
        signMultiplier * fullMode / maximumComponent;

    residual = reducedK * mode - eigenvalues(modeNumber) * reducedM * mode;
    residualScale = max([norm(reducedK * mode), ...
        abs(eigenvalues(modeNumber)) * norm(reducedM * mode), 1]);
    relativeResiduals(modeNumber) = norm(residual) / residualScale;
end

modalMassMatrix = massNormalizedModes.' * reducedM * massNormalizedModes;
modalStiffnessMatrix = ...
    massNormalizedModes.' * reducedK * massNormalizedModes;

assert(isequal(fixedDOFs, 1:5), ...
    'MKEF:ReferenceExampleFailed: unexpected fixed DOFs.');
assert(isequal(freeDOFs, 6:9), ...
    'MKEF:ReferenceExampleFailed: unexpected free DOFs.');
assert(all(all(analysis.modeShapes(fixedDOFs, :) == 0)), ...
    'MKEF:ReferenceExampleFailed: restrained mode components are not zero.');
assert(all(eigenvalues > 0) && all(diff(eigenvalues) >= 0), ...
    'MKEF:ReferenceExampleFailed: eigenvalues are not positive and sorted.');
assert(max(relativeResiduals) < 1e-10, ...
    'MKEF:ReferenceExampleFailed: a modal residual is too large.');
assertClose(modalMassMatrix, eye(modeCount), ...
    1e-10, 1e-12, 'Mass-normalized modes are not M-orthonormal.');
assertClose(modalStiffnessMatrix, diag(eigenvalues), ...
    1e-10, 1e-7, 'Modal stiffness is not diagonal.');
assert(all(all(displayModes(fixedDOFs, :) == 0)), ...
    'MKEF:ReferenceExampleFailed: display scaling changed restrained DOFs.');
assertClose(max(abs(displayModes), [], 1), ones(1, modeCount), ...
    0, 1e-12, 'Display modes are not scaled to unit maximum.');

if verbose
    fprintf('\nExample 05: modal frame\n');
    fprintf('Fixed DOFs: ');
    disp(fixedDOFs);
    fprintf('Free DOFs: ');
    disp(freeDOFs);
    fprintf('Natural frequencies [Hz]:\n');
    disp(analysis.frequenciesHz);
    fprintf('Maximum relative residual by mode: %.6g\n', ...
        max(relativeResiduals));
    fprintf('Phi'' * Mff * Phi after mass normalization:\n');
    disp(modalMassMatrix);
    fprintf('Phi'' * Kff * Phi after mass normalization:\n');
    disp(modalStiffnessMatrix);
    fprintf('Full modes scaled to max(abs(component)) = 1:\n');
    disp(displayModes);
    fprintf(['Checks: reduction, expansion, residuals, sorting, ' ...
        'mass orthogonality, and display scaling PASS.\n']);
end

result = analysis;
result.model = model;
result.reducedStiffness = reducedK;
result.reducedMass = reducedM;
result.massNormalizedReducedModes = massNormalizedModes;
result.modalMassMatrix = modalMassMatrix;
result.modalStiffnessMatrix = modalStiffnessMatrix;
result.displayModeShapes = displayModes;
result.relativeResiduals = relativeResiduals;
end

function repositoryRoot = addRepositoryRoot()
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
