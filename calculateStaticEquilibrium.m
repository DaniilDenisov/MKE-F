function residual = calculateStaticEquilibrium(model, loads, reactions)
%CALCULATESTATICEQUILIBRIUM Return global Fx, Fy, and Mz imbalance.

externalForces = loads + reactions;
xDOFs = model.dofMap(:, 1);
yDOFs = model.dofMap(:, 2);
forceX = externalForces(xDOFs);
forceY = externalForces(yDOFs);

momentZ = sum(model.nodeCoordinates(:, 1) .* forceY - ...
    model.nodeCoordinates(:, 2) .* forceX);
if model.dofPerNode >= 3
    momentDOFs = model.dofMap(:, 3);
    momentZ = momentZ + sum(externalForces(momentDOFs));
end

residual = [sum(forceX); sum(forceY); momentZ];
end
