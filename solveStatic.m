function result = solveStatic(model)
%SOLVESTATIC Solve a static model without modifying the supplied model.

loads = buildStaticLoad(model);
[constrainedK, ~, fixedDOFs] = ...
    applyFixedBoundaryConditions(model, 'static');
displacements = constrainedK \ loads;

result = struct();
result.analysisType = 'static';
result.displacements = displacements;
result.reactions = model.stiffness * displacements - loads;
result.loadVector = loads;
result.fixedDOFs = fixedDOFs;
end
