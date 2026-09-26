function loadVector = buildStaticLoad(model)
%BUILDSTATICLOAD Build a new static nodal-load vector from model data.

loadVector = zeros(model.numberOfDOFs, 1);
for i = 1:size(model.forceBoundaryConditions, 1)
    boundaryCondition = model.forceBoundaryConditions(i, :);
    boundaryType = boundaryCondition(1);
    if boundaryType == 11
        error('MKEF:HarmonicLoadRequiresTimeStep', ...
            'A harmonic load cannot be used in a static analysis.');
    end
    if boundaryType ~= 10
        continue;
    end

    nodeNumber = boundaryCondition(2);
    nodeDOFs = model.dofMap(nodeNumber, :);
    componentCount = min(2, model.dofPerNode);
    for component = 1:componentCount
        loadVector(nodeDOFs(component)) = ...
            loadVector(nodeDOFs(component)) + boundaryCondition(component + 2);
    end
end
end
