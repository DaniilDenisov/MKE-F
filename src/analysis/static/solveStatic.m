function result = solveStatic(model)
%SOLVESTATIC Solve a static model without modifying the supplied model.

loads = buildStaticLoad(model);
constraints = buildConstraintTransform(model);
fixedDOFs = constraints.fixedDOFs;
freeDOFs = constraints.freeDOFs;
T = constraints.T;
reducedK = model.stiffness(freeDOFs, freeDOFs);
if constraints.hasMPC, reducedK = T.' * model.stiffness * T; end
validateReducedSystem(reducedK, [], 'static');

displacements = zeros(model.numberOfDOFs, 1);
if constraints.hasMPC
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
if any(abs(equilibriumResidual) > equilibriumTolerance)
    error('MKEF:EquilibriumFailure', ...
        'Static force or moment equilibrium failed; residual norm is %g.', ...
        norm(equilibriumResidual, inf));
end

result = struct();
result.analysisType = 'static';
result.displacements = displacements;
result.reactions = reactions;
result.loadVector = loads;
result.elementResults = elementResults;
result.fixedDOFs = fixedDOFs;
result.freeDOFs = freeDOFs;
result.equilibriumResidual = equilibriumResidual;
result = recoverConstraintForces(result, constraints, reactions, max(norm(model.stiffness*displacements,inf),norm(loads,inf)));
end
