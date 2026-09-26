function loadHistory = buildTransientLoad(model, timeStep, stepCount)
%BUILDTRANSIENTLOAD Build a fresh load history for one transient run.
% Type 10 intentionally acts in the first integration step only. Type 11
% follows F0*sin(2*pi*f*t), sampled at t = dt, 2*dt, ..., stepCount*dt.

loadHistory = zeros(model.numberOfDOFs, stepCount);
for i = 1:size(model.forceBoundaryConditions, 1)
    boundaryCondition = model.forceBoundaryConditions(i, :);
    boundaryType = boundaryCondition(1);
    nodeNumber = boundaryCondition(2);
    nodeDOFs = model.dofMap(nodeNumber, :);
    componentCount = min(2, model.dofPerNode);

    if boundaryType == 10
        for component = 1:componentCount
            loadHistory(nodeDOFs(component), 1) = ...
                loadHistory(nodeDOFs(component), 1) + ...
                boundaryCondition(component + 2);
        end
    elseif boundaryType == 11
        frequency = boundaryCondition(6);
        time = (1:stepCount) * timeStep;
        for component = 1:componentCount
            loadHistory(nodeDOFs(component), :) = ...
                loadHistory(nodeDOFs(component), :) + ...
                boundaryCondition(component + 2) * ...
                sin(2 * pi * frequency * time);
        end
    end
end
end
