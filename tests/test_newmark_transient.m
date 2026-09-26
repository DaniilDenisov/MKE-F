function test_newmark_transient()
%TEST_NEWMARK_TRANSIENT Independent checks for average-acceleration Newmark.

testFreeVibrationAndInitialState();
testImpulseResponseAndConvergence();
testBeamImpulseRefinement();
testSingleSidedSpectrum();
testInitialConditionValidation();
end

function testFreeVibrationAndInitialState()
mass = 2;
stiffness = 18;
omega = sqrt(stiffness / mass);
period = 2 * pi / omega;
timeStep = period / 200;
stepCount = 1000;
initialDisplacement = 0.02;
initialVelocity = -0.01;
model = makeSDOFModel(mass, stiffness, zeros(0, 6));
options = struct(...
    'timeStep', timeStep, ...
    'duration', stepCount * timeStep, ...
    'initialDisplacement', [initialDisplacement; 0], ...
    'initialVelocity', [initialVelocity; 0]);

result = solveTransient(model, options);
expectedDisplacement = initialDisplacement * cos(omega * result.time) + ...
    (initialVelocity / omega) * sin(omega * result.time);
expectedVelocity = -initialDisplacement * omega * sin(omega * result.time) + ...
    initialVelocity * cos(omega * result.time);

assert(numel(result.time) == stepCount + 1);
assert(size(result.loadHistory, 2) == stepCount + 1);
assertClose(result.time(1), 0, 0, 0, ...
    'The first transient column is not t=0.');
assertClose(result.time(end), options.duration, 1e-12, 1e-12, ...
    'The final transient time is incorrect.');
assertClose(result.accelerations(1, 1), ...
    -stiffness * initialDisplacement / mass, 1e-12, 1e-12, ...
    'Initial acceleration was not calculated from equilibrium.');
assertClose(result.displacements(1, :), expectedDisplacement, ...
    5e-3, 1e-7, 'Free-vibration displacement is inaccurate.');
assertClose(result.velocities(1, :), expectedVelocity, ...
    5e-3, 1e-7, 'Free-vibration velocity is inaccurate.');
assert(all(result.displacements(2, :) == 0));
assert(all(result.velocities(2, :) == 0));
assert(all(result.accelerations(2, :) == 0));
assert(norm(result.equilibriumResidual, inf) < 1e-10);

energy = 0.5 * mass * result.velocities(1, :).^2 + ...
    0.5 * stiffness * result.displacements(1, :).^2;
relativeEnergyDrift = max(abs(energy - energy(1))) / energy(1);
if relativeEnergyDrift > 1e-10
    error('MKEF:VerificationFailed', ...
        'Undamped Newmark energy drift is %g.', relativeEnergyDrift);
end
end

function testImpulseResponseAndConvergence()
% With F0=I/dt, the one-sample pulse has fixed discrete impulse I. As dt
% decreases its response converges to the analytical ideal-impulse response.
mass = 2;
stiffness = 18;
omega = sqrt(stiffness / mass);
period = 2 * pi / omega;
impulse = 4;
duration = 2 * period;

[coarseError, responseScale] = impulseError(...
    mass, stiffness, impulse, period / 40, duration);
[fineError, ~] = impulseError(...
    mass, stiffness, impulse, period / 80, duration);

if fineError >= 0.4 * coarseError
    error('MKEF:VerificationFailed', ...
        'SDOF impulse response does not converge under time-step refinement.');
end
if fineError > 0.01 * responseScale
    error('MKEF:VerificationFailed', ...
        'Refined SDOF impulse response is not close to the analytical solution.');
end
end

function [errorValue, responseScale] = impulseError(...
    mass, stiffness, impulse, timeStep, duration)
omega = sqrt(stiffness / mass);
forceAmplitude = impulse / timeStep;
model = makeSDOFModel(mass, stiffness, ...
    [12, 1, forceAmplitude, 0, 0, 0]);
result = solveTransient(model, ...
    struct('timeStep', timeStep, 'duration', duration));

assertClose(sum(result.loadHistory(1, :)) * timeStep, impulse, ...
    1e-12, 1e-12, 'The discrete pulse impulse is not F0*dt.');
assert(result.loadHistory(1, 1) == 0);
assert(result.loadHistory(1, 2) == forceAmplitude);

% The single discrete force sample is centered at t=dt. Ignore its two-step
% numerical representation when comparing the subsequent free response.
comparison = result.time >= 2 * timeStep;
analytical = impulse / (mass * omega) * ...
    sin(omega * (result.time(comparison) - timeStep));
numerical = result.displacements(1, comparison);
errorValue = max(abs(numerical - analytical));
responseScale = impulse / (mass * omega);
end

function testBeamImpulseRefinement()
options = struct('verbose', false, 'plotting', false);
problem = StructFEProblem('CaseBeamDyn.txt', options);
baseModel = problem.GetAnalysisModel();
impulse = -0.1;
timeSteps = [2e-4, 1e-4, 5e-5];
duration = 0.01;
results = cell(1, numel(timeSteps));

for i = 1:numel(timeSteps)
    model = baseModel;
    model.forceBoundaryConditions = ...
        [12, 3, 0, impulse / timeSteps(i), 0];
    results{i} = solveTransient(model, struct(...
        'timeStep', timeSteps(i), ...
        'duration', duration + timeSteps(i)));
end

selectedDOF = baseModel.dofMap(3, 2);
coarseIndices = 4:40;
coarseResponse = results{1}.displacements(...
    selectedDOF, coarseIndices + 2);
mediumResponse = results{2}.displacements(...
    selectedDOF, 2 * coarseIndices + 2);
fineResponse = results{3}.displacements(...
    selectedDOF, 4 * coarseIndices + 2);
coarseDifference = max(abs(coarseResponse - mediumResponse));
fineDifference = max(abs(mediumResponse - fineResponse));

if fineDifference >= 0.5 * coarseDifference
    error('MKEF:VerificationFailed', ...
        'Beam impulse response does not converge under time-step refinement.');
end
end

function testSingleSidedSpectrum()
timeStep = 0.01;
sampleCount = 200;
frequency = 5;
time = (0:(sampleCount - 1)) * timeStep;
signal = 2 + 3 * sin(2 * pi * frequency * time);
[frequencyHz, amplitude] = ...
    computeSingleSidedSpectrum(signal, timeStep);

[~, frequencyIndex] = min(abs(frequencyHz - frequency));
assertClose(frequencyHz(frequencyIndex), frequency, 0, 1e-12, ...
    'The FFT frequency bins are incorrect.');
assertClose(amplitude(1), 2, 1e-12, 1e-12, ...
    'The FFT DC amplitude is incorrectly normalized.');
assertClose(amplitude(frequencyIndex), 3, 1e-12, 1e-12, ...
    'The single-sided FFT amplitude is incorrectly normalized.');
end

function testInitialConditionValidation()
model = makeSDOFModel(2, 18, zeros(0, 6));
baseOptions = struct('timeStep', 0.01, 'duration', 0.1);

invalidSize = baseOptions;
invalidSize.initialDisplacement = 1;
assertThrows('MKEF:InvalidInitialConditions', ...
    @() solveTransient(model, invalidSize));

restrainedMotion = baseOptions;
restrainedMotion.initialVelocity = [0; 1];
assertThrows('MKEF:InvalidInitialConditions', ...
    @() solveTransient(model, restrainedMotion));

facadeOptions = struct('verbose', false, 'plotting', false);
problem = StructFEProblem('Case1ElementBeam.txt', facadeOptions);
facadeInitialDisplacement = zeros(6, 1);
facadeInitialDisplacement(5) = 1e-5;
facadeResult = problem.RunTransient(1e-4, 2e-4, 2, 2, ...
    struct('initialDisplacement', facadeInitialDisplacement));
assertClose(facadeResult.displacements(:, 1), ...
    facadeInitialDisplacement, 0, 0, ...
    'RunTransient did not pass initial conditions to the numerical core.');
end

function model = makeSDOFModel(mass, stiffness, forceBoundaryConditions)
% DOF 1 is the oscillator; DOF 2 is restrained by support type 2.
model = struct();
model.stiffness = diag([stiffness, 1]);
model.mass = diag([mass, 1]);
model.dofMap = [1, 2];
model.fixedBoundaryConditions = [2, 1, 0, 0, 0];
model.forceBoundaryConditions = forceBoundaryConditions;
model.nodeCoordinates = [0, 0];
model.numberOfNodes = 1;
model.dofPerNode = 2;
model.numberOfDOFs = 2;
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
