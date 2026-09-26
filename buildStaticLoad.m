function loadVector = buildStaticLoad(model)
%BUILDSTATICLOAD Build a new static nodal-load vector from model data.

loadVector = zeros(model.numberOfDOFs, 1);
for i = 1:size(model.forceBoundaryConditions, 1)
    boundaryCondition = model.forceBoundaryConditions(i, :);
    boundaryType = boundaryCondition(1);
    if boundaryType ~= 10
        if ismember(boundaryType, [11, 12, 13])
            error('MKEF:TimeDependentLoadInStaticAnalysis', ...
                'Load type %d is time-dependent and cannot be used in a static analysis.', ...
                boundaryType);
        end
        error('MKEF:UnsupportedLoadType', ...
            'Unsupported nodal load type %g.', boundaryType);
    end

    nodeNumber = boundaryCondition(2);
    nodeDOFs = model.dofMap(nodeNumber, :);
    components = getNodalLoadComponents(model, boundaryCondition);
    loadVector(nodeDOFs) = loadVector(nodeDOFs) + components;
end
end
