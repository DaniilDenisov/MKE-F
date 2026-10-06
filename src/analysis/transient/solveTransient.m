function result = solveTransient(model, options)
%SOLVETRANSIENT Integrate an undamped model with Newmark average acceleration.
% All state and load histories are allocated for this call and returned in a
% result struct. The model is never modified.

rejectElementLoads(model);
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
constraints = buildConstraintTransform(model);
fixedDOFs = constraints.fixedDOFs;
freeDOFs = constraints.freeDOFs;
independent = constraints.independentDOFs;
T = constraints.T;
reducedK = model.stiffness(freeDOFs, freeDOFs);
reducedM = model.mass(freeDOFs, freeDOFs);
reducedLoads = loads(freeDOFs,:);
if constraints.hasMPC
    reducedK = T.' * model.stiffness * T;
    reducedM = T.' * model.mass * T;
    reducedLoads = T.' * loads;
end
validateReducedSystem(reducedK, reducedM, 'transient');

reducedAccelerations = zeros(numel(independent), stepCount + 1);
reducedVelocities = zeros(numel(independent), stepCount + 1);
reducedDisplacements = zeros(numel(independent), stepCount + 1);

initialDisplacement = getInitialCondition(options, ...
    'initialDisplacement', model.numberOfDOFs);
initialVelocity = getInitialCondition(options, ...
    'initialVelocity', model.numberOfDOFs);
if any(initialDisplacement(fixedDOFs) ~= 0) || ...
        any(initialVelocity(fixedDOFs) ~= 0)
    error('MKEF:InvalidInitialConditions', ...
        'Initial displacement and velocity must be zero at restrained DOFs.');
end
if constraints.hasMPC
    if norm(initialDisplacement-T*initialDisplacement(independent),inf) > 1e-10*max(1,norm(initialDisplacement,inf)) || ...
            norm(initialVelocity-T*initialVelocity(independent),inf) > 1e-10*max(1,norm(initialVelocity,inf))
        error('MKEF:InvalidInitialConditions','Initial displacement and velocity must satisfy all MPC equations.');
    end
end
reducedDisplacements(:, 1) = initialDisplacement(independent);
reducedVelocities(:, 1) = initialVelocity(independent);
reducedAccelerations(:, 1) = reducedM \ ...
    (reducedLoads(:, 1) - reducedK * reducedDisplacements(:, 1));

% Newmark-beta average-acceleration method, equations (22)-(23) in
% H.P. Gavin, Numerical Integration in Structural Dynamics, Duke University:
% https://people.duke.edu/~hpgavin/StructuralDynamics/NumericalIntegration.pdf
gamma = 0.5;
beta = 0.25;
newmarkA0 = 1 / (beta * timeStep^2);
newmarkA2 = 1 / (beta * timeStep);
newmarkA3 = 1 / (2 * beta) - 1;
newmarkA6 = timeStep * (1 - gamma);
newmarkA7 = gamma * timeStep;
effectiveK = reducedK + newmarkA0 * reducedM;
effectiveFactor = chol((effectiveK + effectiveK.') / 2);

for step = 1:stepCount
    effectiveLoad = reducedLoads(:, step + 1) + reducedM * ...
        (newmarkA0 * reducedDisplacements(:, step) + ...
         newmarkA2 * reducedVelocities(:, step) + ...
         newmarkA3 * reducedAccelerations(:, step));
    reducedDisplacements(:, step + 1) = effectiveFactor \ ...
        (effectiveFactor.' \ effectiveLoad);
    reducedAccelerations(:, step + 1) = ...
        newmarkA0 * (reducedDisplacements(:, step + 1) - ...
        reducedDisplacements(:, step)) - ...
        newmarkA2 * reducedVelocities(:, step) - ...
        newmarkA3 * reducedAccelerations(:, step);
    reducedVelocities(:, step + 1) = reducedVelocities(:, step) + ...
        newmarkA6 * reducedAccelerations(:, step) + ...
        newmarkA7 * reducedAccelerations(:, step + 1);
end

accelerations = full(T * reducedAccelerations);
velocities = full(T * reducedVelocities);
displacements = full(T * reducedDisplacements);

dynamicResidual = model.mass * accelerations + ...
    model.stiffness * displacements - loads;
[spectrumFrequencyHz, displacementAmplitudeSpectrum] = ...
    computeSingleSidedSpectrum(displacements, timeStep);

result = struct();
result.analysisType = 'transient';
result.time = (0:stepCount) * timeStep;
result.displacements = displacements;
result.velocities = velocities;
result.accelerations = accelerations;
result.loadHistory = loads;
result.reactions = dynamicResidual;
result.equilibriumResidual = T.' * dynamicResidual;
result.spectrumFrequencyHz = spectrumFrequencyHz;
result.displacementAmplitudeSpectrum = displacementAmplitudeSpectrum;
result.timeStep = timeStep;
result.stepCount = stepCount;
result.duration = stepCount * timeStep;
result.requestedDuration = duration;
result.fixedDOFs = fixedDOFs;
result.freeDOFs = freeDOFs;
result.newmarkBeta = beta;
result.newmarkGamma = gamma;
result = recoverConstraintForces(result,constraints,dynamicResidual, ...
    max([norm(model.mass*accelerations,inf),norm(model.stiffness*displacements,inf),norm(loads,inf)]));
end

function values = getInitialCondition(options, fieldName, numberOfDOFs)
if ~isfield(options, fieldName)
    values = zeros(numberOfDOFs, 1);
    return;
end

values = options.(fieldName);
if ~isnumeric(values) || ~isvector(values) || ...
        numel(values) ~= numberOfDOFs || any(~isfinite(values(:)))
    error('MKEF:InvalidInitialConditions', ...
        '%s must be a finite vector with one value per model DOF.', ...
        fieldName);
end
values = values(:);
end
