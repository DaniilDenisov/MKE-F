function result = example_06_newmark_sdof(verbose)
%EXAMPLE_06_NEWMARK_SDOF Compare Newmark response with free vibration.
%   RESULT = EXAMPLE_06_NEWMARK_SDOF() integrates an undamped oscillator
%   on three time grids, checks second-order convergence, energy, initial
%   acceleration, and its spectrum. Pass false for quiet use.

if nargin < 1
    verbose = true;
end
validateVerbose(verbose);
addRepositoryRoot();

mass = 2;
stiffness = 18;
initialDisplacement = 0.02;
initialVelocity = -0.01;
omega = sqrt(stiffness / mass);
naturalFrequencyHz = omega / (2 * pi);
period = 2 * pi / omega;
periodCount = 3;
stepsPerPeriod = [50 100 200];
analyses = cell(size(stepsPerPeriod));
maximumDisplacementErrors = zeros(size(stepsPerPeriod));
maximumVelocityErrors = zeros(size(stepsPerPeriod));

model = makeSDOFModel(mass, stiffness);

for gridNumber = 1:numel(stepsPerPeriod)
    timeStep = period / stepsPerPeriod(gridNumber);
    duration = periodCount * period;
    options = struct(...
        'timeStep', timeStep, ...
        'duration', duration, ...
        'initialDisplacement', [initialDisplacement; 0], ...
        'initialVelocity', [initialVelocity; 0]);
    analysis = solveTransient(model, options);
    exactDisplacement = initialDisplacement * cos(omega * analysis.time) + ...
        (initialVelocity / omega) * sin(omega * analysis.time);
    exactVelocity = -initialDisplacement * omega * sin(omega * analysis.time) + ...
        initialVelocity * cos(omega * analysis.time);

    maximumDisplacementErrors(gridNumber) = ...
        max(abs(analysis.displacements(1, :) - exactDisplacement));
    maximumVelocityErrors(gridNumber) = ...
        max(abs(analysis.velocities(1, :) - exactVelocity));
    analyses{gridNumber} = analysis;
end

analysis = analyses{end};
timeStep = analysis.timeStep;
exactDisplacement = initialDisplacement * cos(omega * analysis.time) + ...
    (initialVelocity / omega) * sin(omega * analysis.time);
exactVelocity = -initialDisplacement * omega * sin(omega * analysis.time) + ...
    initialVelocity * cos(omega * analysis.time);
exactAcceleration = -omega^2 * exactDisplacement;
energy = 0.5 * mass * analysis.velocities(1, :).^2 + ...
    0.5 * stiffness * analysis.displacements(1, :).^2;
relativeEnergyDrift = max(abs(energy - energy(1))) / energy(1);
effectiveStiffness = stiffness + 4 * mass / timeStep^2;

[~, peakIndex] = max(analysis.displacementAmplitudeSpectrum(1, 2:end));
peakIndex = peakIndex + 1;
spectrumPeakHz = analysis.spectrumFrequencyHz(peakIndex);
frequencyResolutionHz = 1 / (numel(analysis.time) * timeStep);

assert(analysis.newmarkBeta == 0.25 && analysis.newmarkGamma == 0.5, ...
    'MKEF:ReferenceExampleFailed: unexpected Newmark parameters.');
assert(analysis.time(1) == 0 && ...
    abs(analysis.time(end) - periodCount * period) < 1e-12, ...
    'MKEF:ReferenceExampleFailed: the time grid endpoints are incorrect.');
assert(numel(analysis.time) == analysis.stepCount + 1, ...
    'MKEF:ReferenceExampleFailed: the time grid length is incorrect.');
assertClose(analysis.accelerations(1, 1), ...
    -stiffness * initialDisplacement / mass, ...
    1e-12, 1e-12, 'Initial acceleration does not satisfy equilibrium.');
assert(all(all(analysis.displacements(analysis.fixedDOFs, :) == 0)), ...
    'MKEF:ReferenceExampleFailed: restrained displacements are not zero.');
assert(all(all(analysis.velocities(analysis.fixedDOFs, :) == 0)), ...
    'MKEF:ReferenceExampleFailed: restrained velocities are not zero.');
assert(all(all(analysis.accelerations(analysis.fixedDOFs, :) == 0)), ...
    'MKEF:ReferenceExampleFailed: restrained accelerations are not zero.');
assert(maximumDisplacementErrors(2) < 0.35 * maximumDisplacementErrors(1) && ...
    maximumDisplacementErrors(3) < 0.35 * maximumDisplacementErrors(2), ...
    'MKEF:ReferenceExampleFailed: displacement does not converge quadratically.');
assert(maximumVelocityErrors(2) < 0.35 * maximumVelocityErrors(1) && ...
    maximumVelocityErrors(3) < 0.35 * maximumVelocityErrors(2), ...
    'MKEF:ReferenceExampleFailed: velocity does not converge quadratically.');
assertClose(analysis.displacements(1, :), exactDisplacement, ...
    2e-3, 1e-8, 'The refined displacement history is inaccurate.');
assertClose(analysis.velocities(1, :), exactVelocity, ...
    2e-3, 1e-8, 'The refined velocity history is inaccurate.');
assertClose(analysis.accelerations(1, :), exactAcceleration, ...
    2e-3, 1e-8, 'The refined acceleration history is inaccurate.');
assert(relativeEnergyDrift < 1e-10, ...
    'MKEF:ReferenceExampleFailed: undamped energy drift is too large.');
assert(abs(spectrumPeakHz - naturalFrequencyHz) <= frequencyResolutionHz, ...
    'MKEF:ReferenceExampleFailed: the spectrum peak misses the natural frequency.');

if verbose
    fprintf('\nExample 06: Newmark undamped SDOF\n');
    fprintf('Natural frequency: %.9g Hz; period: %.9g s.\n', ...
        naturalFrequencyHz, period);
    fprintf('Steps per period and maximum displacement error:\n');
    disp([stepsPerPeriod(:), maximumDisplacementErrors(:)]);
    fprintf('Steps per period and maximum velocity error:\n');
    disp([stepsPerPeriod(:), maximumVelocityErrors(:)]);
    fprintf('Initial acceleration: %.9g m/s^2.\n', ...
        analysis.accelerations(1, 1));
    fprintf('Effective stiffness K + 4M/dt^2: %.9g.\n', ...
        effectiveStiffness);
    fprintf('Relative energy drift: %.6g.\n', relativeEnergyDrift);
    fprintf('Spectrum peak: %.9g Hz; bin width: %.9g Hz.\n', ...
        spectrumPeakHz, frequencyResolutionHz);
    fprintf(['Checks: initial state, time grid, analytical history, ' ...
        'quadratic convergence, energy, and spectrum PASS.\n']);
end

result = analysis;
result.model = model;
result.mass = mass;
result.stiffness = stiffness;
result.omega = omega;
result.naturalFrequencyHz = naturalFrequencyHz;
result.period = period;
result.stepsPerPeriod = stepsPerPeriod;
result.maximumDisplacementErrors = maximumDisplacementErrors;
result.maximumVelocityErrors = maximumVelocityErrors;
result.exactDisplacement = exactDisplacement;
result.exactVelocity = exactVelocity;
result.exactAcceleration = exactAcceleration;
result.energy = energy;
result.relativeEnergyDrift = relativeEnergyDrift;
result.effectiveStiffness = effectiveStiffness;
result.spectrumPeakHz = spectrumPeakHz;
result.frequencyResolutionHz = frequencyResolutionHz;
end

function model = makeSDOFModel(mass, stiffness)
model = struct();
model.stiffness = diag([stiffness, 1]);
model.mass = diag([mass, 1]);
model.dofMap = [1 2];
model.fixedBoundaryConditions = [2 1 0 0 0];
model.forceBoundaryConditions = zeros(0, 6);
model.nodeCoordinates = [0 0];
model.numberOfNodes = 1;
model.dofPerNode = 2;
model.numberOfDOFs = 2;
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
