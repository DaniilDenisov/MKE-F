function result = solveTransient(model, options)
%SOLVETRANSIENT Integrate an undamped model with Newmark average acceleration.
% All state and load histories are allocated for this call and returned in a
% result struct. The model is never modified.

if ~isstruct(options) || ~isfield(options, 'timeStep') || ...
        ~isfield(options, 'duration')
    error('MKEF:InvalidTransientOptions', ...
        'Transient options require timeStep and duration fields.');
end

timeStep = options.timeStep;
duration = options.duration;
if ~isscalar(timeStep) || ~isfinite(timeStep) || timeStep <= 0 || ...
        ~isscalar(duration) || ~isfinite(duration) || duration < timeStep
    error('MKEF:InvalidTransientOptions', ...
        'timeStep must be positive and duration must cover at least one step.');
end

stepCount = fix(duration / timeStep);
loads = buildTransientLoad(model, timeStep, stepCount);
[fixedDOFs, freeDOFs] = partitionDOFs(model);
reducedK = model.stiffness(freeDOFs, freeDOFs);
reducedM = model.mass(freeDOFs, freeDOFs);
validateReducedSystem(reducedK, reducedM, 'transient');

reducedAccelerations = zeros(numel(freeDOFs), stepCount + 1);
reducedVelocities = zeros(numel(freeDOFs), stepCount + 1);
reducedDisplacements = zeros(numel(freeDOFs), stepCount + 1);

gamma = 0.5;
beta = 0.25;
a0 = 1 / (beta * timeStep^2);
a2 = 1 / (beta * timeStep);
a3 = 1 / (2 * beta) - 1;
a6 = timeStep * (1 - gamma);
a7 = gamma * timeStep;
effectiveK = reducedK + a0 * reducedM;

for step = 1:stepCount
    effectiveLoad = loads(freeDOFs, step) + reducedM * ...
        (a0 * reducedDisplacements(:, step) + ...
         a2 * reducedVelocities(:, step) + ...
         a3 * reducedAccelerations(:, step));
    reducedDisplacements(:, step + 1) = effectiveK \ effectiveLoad;
    reducedAccelerations(:, step + 1) = ...
        a0 * (reducedDisplacements(:, step + 1) - ...
        reducedDisplacements(:, step)) - ...
        a2 * reducedVelocities(:, step) - ...
        a3 * reducedAccelerations(:, step);
    reducedVelocities(:, step + 1) = reducedVelocities(:, step) + ...
        a6 * reducedAccelerations(:, step) + ...
        a7 * reducedAccelerations(:, step + 1);
end

accelerations = zeros(model.numberOfDOFs, stepCount + 1);
velocities = zeros(model.numberOfDOFs, stepCount + 1);
displacements = zeros(model.numberOfDOFs, stepCount + 1);
accelerations(freeDOFs, :) = reducedAccelerations;
velocities(freeDOFs, :) = reducedVelocities;
displacements(freeDOFs, :) = reducedDisplacements;

result = struct();
result.analysisType = 'transient';
result.time = (0:stepCount) * timeStep;
result.displacements = displacements;
result.velocities = velocities;
result.accelerations = accelerations;
result.loadHistory = loads;
result.timeStep = timeStep;
result.stepCount = stepCount;
result.duration = stepCount * timeStep;
result.requestedDuration = duration;
result.fixedDOFs = fixedDOFs;
result.freeDOFs = freeDOFs;
end
