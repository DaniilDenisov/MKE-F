function [fixedDOFs, freeDOFs] = partitionDOFs(model)
%PARTITIONDOFS Partition masks; legacy numeric records are adapted at entry.
if isfield(model,'nodalKinematics')
    c=buildSpatialConstraints(model); fixedDOFs=c.fixedDOFs; freeDOFs=c.freeDOFs; return;
end
dimension=2; if isfield(model,'dimension'), dimension=model.dimension; end
supports = normalizeSupports(model.fixedBoundaryConditions, model.dofPerNode,dimension);
restrained = false(1, model.numberOfDOFs);
labels = {'ux', 'uy', 'thetaZ'};
if dimension==3, labels=model.dofLabels; end
for i = 1:numel(supports)
    s = supports(i);
    if s.node ~= fix(s.node) || s.node < 1 || s.node > model.numberOfNodes
        error('MKEF:InvalidConstraint', 'Constraint %d refers to invalid node %g.', i, s.node);
    end
    if dimension==3
        local=find([s.fixUx s.fixUy s.fixUz s.fixThetaX s.fixThetaY s.fixThetaZ]);
    else
        local = find([s.fixUx s.fixUy s.fixThetaZ]);
    end
    for dof = local
        id = model.dofMap(s.node, dof);
        if id == 0, continue; end
        if restrained(id)
            error('MKEF:DuplicateConstraint', 'Constraint %d repeats node %d DOF %s.', i, s.node, labels{dof});
        end
        restrained(id) = true;
    end
end
fixedDOFs = find(restrained);
freeDOFs = find(~restrained);
end
