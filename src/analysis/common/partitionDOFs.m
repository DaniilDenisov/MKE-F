function [fixedDOFs, freeDOFs] = partitionDOFs(model)
%PARTITIONDOFS Partition masks; legacy numeric records are adapted at entry.
supports = normalizeSupports(model.fixedBoundaryConditions, model.dofPerNode);
restrained = false(1, model.numberOfDOFs);
labels = {'ux', 'uy', 'thetaZ'};
for i = 1:numel(supports)
    s = supports(i);
    if s.node ~= fix(s.node) || s.node < 1 || s.node > model.numberOfNodes
        error('MKEF:InvalidConstraint', 'Constraint %d refers to invalid node %g.', i, s.node);
    end
    local = find([s.fixUx s.fixUy s.fixThetaZ]);
    for dof = local
        id = model.dofMap(s.node, dof);
        if restrained(id)
            error('MKEF:DuplicateConstraint', 'Constraint %d repeats node %d DOF %s.', i, s.node, labels{dof});
        end
        restrained(id) = true;
    end
end
fixedDOFs = find(restrained);
freeDOFs = find(~restrained);
end
