function [fixedDOFs, freeDOFs] = partitionDOFs(model)
%PARTITIONDOFS Convert support definitions to validated DOF index vectors.

fixedDOFs = [];
for i = 1:size(model.fixedBoundaryConditions, 1)
    boundaryCondition = model.fixedBoundaryConditions(i, :);
    boundaryType = boundaryCondition(1);
    nodeNumber = boundaryCondition(2);

    if ~isfinite(boundaryType) || boundaryType ~= fix(boundaryType) || ...
            ~ismember(boundaryType, 1:4)
        error('MKEF:InvalidConstraint', ...
            'Constraint %d has unsupported type %g.', i, boundaryType);
    end
    if ~isfinite(nodeNumber) || nodeNumber ~= fix(nodeNumber) || ...
            nodeNumber < 1 || nodeNumber > model.numberOfNodes
        error('MKEF:InvalidConstraint', ...
            'Constraint %d refers to invalid node %g.', i, nodeNumber);
    end

    switch boundaryType
        case 1
            localFixedDOFs = 1:model.dofPerNode;
        case 2
            localFixedDOFs = 2:model.dofPerNode;
        case 3
            localFixedDOFs = [1, 3:model.dofPerNode];
        case 4
            localFixedDOFs = [1:2, 4:model.dofPerNode];
    end

    newFixedDOFs = model.dofMap(nodeNumber, localFixedDOFs);
    duplicateDOFs = intersect(fixedDOFs, newFixedDOFs);
    if ~isempty(duplicateDOFs)
        error('MKEF:DuplicateConstraint', ...
            'Constraint %d repeats restrained global DOF %d.', ...
            i, duplicateDOFs(1));
    end
    fixedDOFs = [fixedDOFs, newFixedDOFs]; %#ok<AGROW>
end

fixedDOFs = sort(fixedDOFs);
freeDOFs = setdiff(1:model.numberOfDOFs, fixedDOFs);
end
