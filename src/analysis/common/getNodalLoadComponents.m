function components = getNodalLoadComponents(model, boundaryCondition)
%GETNODALLOADCOMPONENTS Return loads conjugate to the model nodal DOFs.
% A 2D truss node uses [Fx; Fy]. A 2D frame node uses [Fx; Fy; Mz].

requiredColumns = 2 + model.dofPerNode;
if numel(boundaryCondition) < requiredColumns
    error('MKEF:InvalidNodalLoad', ...
        'The nodal load does not provide all %d model components.', ...
        model.dofPerNode);
end

if model.dofPerNode == 2 && numel(boundaryCondition) >= 5 && ...
        boundaryCondition(5) ~= 0
    error('MKEF:UnsupportedLoadComponent', ...
        ['A 2D truss node accepts only Fx and Fy. It has no rotational ' ...
         'degree of freedom for a nodal moment Mz.']);
end

components = boundaryCondition(3:requiredColumns).';
end
