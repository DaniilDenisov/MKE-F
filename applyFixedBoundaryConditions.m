function [constrainedK, constrainedM, fixedDOFs] = ...
    applyFixedBoundaryConditions(model, analysisType)
%APPLYFIXEDBOUNDARYCONDITIONS Apply the legacy constraint transformation.
% This helper edits local matrix copies only. Commit 5 will replace this
% transformation with an explicit free-DOF reduction.

constrainedK = model.stiffness;
constrainedM = model.mass;
fixedDOFs = [];
isStatic = strcmp(analysisType, 'static');

for i = 1:size(model.fixedBoundaryConditions, 1)
    boundaryCondition = model.fixedBoundaryConditions(i, :);
    boundaryType = boundaryCondition(1);
    nodeNumber = boundaryCondition(2);
    nodeDOFs = model.dofMap(nodeNumber, :);

    switch boundaryType
        case 1
            localFixedDOFs = 1:model.dofPerNode;
        case 2
            localFixedDOFs = 2:model.dofPerNode;
        case 3
            localFixedDOFs = [1, 3:model.dofPerNode];
        case 4
            localFixedDOFs = [1:2, 4:model.dofPerNode];
        otherwise
            localFixedDOFs = [];
    end

    for localDOF = localFixedDOFs
        globalDOF = nodeDOFs(localDOF);
        fixedDOFs(end + 1) = globalDOF; %#ok<AGROW>
        constrainedM(globalDOF, :) = 0;
        constrainedM(:, globalDOF) = 0;
        constrainedM(globalDOF, globalDOF) = 1;
        constrainedK(globalDOF, :) = 0;
        constrainedK(:, globalDOF) = 0;
        if isStatic
            constrainedK(globalDOF, globalDOF) = 1;
        end
    end
end

fixedDOFs = unique(fixedDOFs);
end
