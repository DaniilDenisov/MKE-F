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
[constrainedK, constrainedM, fixedDOFs] = ...
    applyFixedBoundaryConditions(model, 'transient');

accelerations = zeros(model.numberOfDOFs, stepCount + 1);
velocities = zeros(model.numberOfDOFs, stepCount + 1);
displacements = zeros(model.numberOfDOFs, stepCount + 1);

gamma = 0.5;
beta = 0.25;
a0 = 1 / (beta * timeStep^2);
a2 = 1 / (beta * timeStep);
a3 = 1 / (2 * beta) - 1;
a6 = timeStep * (1 - gamma);
a7 = gamma * timeStep;
effectiveK = constrainedK + a0 * constrainedM;

for step = 1:stepCount
    effectiveLoad = loads(:, step) + constrainedM * ...
        (a0 * displacements(:, step) + ...
         a2 * velocities(:, step) + a3 * accelerations(:, step));
    displacements(:, step + 1) = effectiveK \ effectiveLoad;
    accelerations(:, step + 1) = ...
        a0 * (displacements(:, step + 1) - displacements(:, step)) - ...
        a2 * velocities(:, step) - a3 * accelerations(:, step);
    velocities(:, step + 1) = velocities(:, step) + ...
        a6 * accelerations(:, step) + a7 * accelerations(:, step + 1);
end

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
end
