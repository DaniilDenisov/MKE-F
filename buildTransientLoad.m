function loadHistory = buildTransientLoad(model, timeStep, stepCount)
%BUILDTRANSIENTLOAD Build a fresh load history for one transient run.
% Type 10 is the backward-compatible one-step pulse. Type 11 follows
% F0*sin(2*pi*f*t). Type 12 is an explicit one-step pulse, and type 13 is a
% persistent step. Columns correspond to t = 0, dt, ..., stepCount*dt.

loadHistory = zeros(model.numberOfDOFs, stepCount + 1);
time = (0:stepCount) * timeStep;
for i = 1:size(model.forceBoundaryConditions, 1)
    boundaryCondition = model.forceBoundaryConditions(i, :);
    boundaryType = boundaryCondition(1);
    nodeNumber = boundaryCondition(2);
    nodeDOFs = model.dofMap(nodeNumber, :);
    components = getNodalLoadComponents(model, boundaryCondition);

    if boundaryType == 10 || boundaryType == 12
        loadHistory(nodeDOFs, 2) = ...
            loadHistory(nodeDOFs, 2) + components;
    elseif boundaryType == 11
        frequency = boundaryCondition(6);
        loadHistory(nodeDOFs, :) = loadHistory(nodeDOFs, :) + ...
            components * sin(2 * pi * frequency * time);
    elseif boundaryType == 13
        loadHistory(nodeDOFs, :) = loadHistory(nodeDOFs, :) + ...
            repmat(components, 1, stepCount + 1);
    else
        error('MKEF:UnsupportedLoadType', ...
            'Unsupported nodal load type %g.', boundaryType);
    end
end
end
