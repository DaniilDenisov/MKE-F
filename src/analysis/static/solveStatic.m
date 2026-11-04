function result = solveStatic(model)
%SOLVESTATIC Solve a static model without modifying the supplied model.

validateSpatialFrameSections(model);
loads = buildStaticLoad(model);
constraints = buildConstraintTransform(model);
fixedDOFs = constraints.fixedDOFs;
freeDOFs = constraints.freeDOFs;
T = constraints.T;
reducedK = model.stiffness(freeDOFs, freeDOFs);
if constraints.hasTransform, reducedK = T.' * model.stiffness * T; end
validateReducedSystem(reducedK, [], 'static');

displacements = zeros(model.numberOfDOFs, 1);
if constraints.hasTransform
    displacements = full(T * (reducedK \ (T.' * loads)));
elseif ~isempty(freeDOFs)
    displacements(freeDOFs) = reducedK \ loads(freeDOFs);
end
reactions = model.stiffness * displacements - loads;
elementResults = recoverElementResults(model, displacements);
equilibriumResidual = calculateStaticEquilibrium(model, loads, reactions);
if isfield(model, 'elementLoads') && ~isempty(model.elementLoads)
    [~, ~, sourceResultant] = getElementLoadData(model);
    equilibriumResidual = calculateStaticEquilibrium(model, ...
        buildNodalLoadVector(model), reactions) + sourceResultant;
end
equilibriumScale = max([norm(loads, 1), norm(reactions, 1), 1]);
lengthScale = max([abs(model.nodeCoordinates(:)); 1]);
equilibriumTolerance = 1e-9 * ...
    [equilibriumScale; equilibriumScale; equilibriumScale * lengthScale];
if isfield(model,'dimension') && model.dimension==3
    equilibriumTolerance=1e-9*equilibriumScale*[ones(3,1);lengthScale*ones(3,1)];
end
if any(abs(equilibriumResidual) > equilibriumTolerance)
    error('MKEF:EquilibriumFailure', ...
        'Static force or moment equilibrium failed; residual norm is %g.', ...
        norm(equilibriumResidual, inf));
end

result = struct();
result.analysisType = 'static';
result.displacements = displacements;
if isfield(model,'nodalKinematics'), result.nodalDisplacements=full(model.nodalKinematics*displacements); end
result.reactions = reactions;
result.loadVector = loads;
result.elementResults = elementResults;
result.fixedDOFs = fixedDOFs;
result.freeDOFs = freeDOFs;
result.equilibriumResidual = equilibriumResidual;
result = recoverConstraintForces(result, constraints, reactions, max(norm(model.stiffness*displacements,inf),norm(loads,inf)));
end
