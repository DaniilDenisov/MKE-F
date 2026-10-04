function [loadVector, localVectors] = buildElementLoadVector(model)
%BUILDELEMENTLOADVECTOR Assemble consistent element loads by the DOF map.
loadVector = zeros(model.numberOfDOFs, 1);
if ~isfield(model, 'elementLoads') || isempty(model.elementLoads)
    localVectors = {};
    return;
end
[localVectors, ~, ~] = getElementLoadData(model);
for i = 1:numel(model.elementData)
    element = model.elementData(i);
    loadVector(element.dofs) = loadVector(element.dofs) + ...
        element.transformation.' * localVectors{i};
end
end
