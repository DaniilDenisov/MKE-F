function [globalK, globalM] = assembleGlobalMatrices(elements, numberOfDOFs)
%ASSEMBLEGLOBALMATRICES Assemble sparse global matrices from triplet arrays.

if ~isstruct(elements) || ~isscalar(numberOfDOFs) || ...
        numberOfDOFs < 1 || numberOfDOFs ~= fix(numberOfDOFs)
    error('MKEF:InvalidAssemblyInput', ...
        'Assembly requires an element struct array and a positive DOF count.');
end

entryCount = 0;
for elementNumber = 1:numel(elements)
    entryCount = entryCount + numel(elements(elementNumber).dofs)^2;
end

rows = zeros(entryCount, 1);
columns = zeros(entryCount, 1);
stiffnessValues = zeros(entryCount, 1);
massValues = zeros(entryCount, 1);
cursor = 1;

for elementNumber = 1:numel(elements)
    element = elements(elementNumber);
    validateElementForAssembly(element, numberOfDOFs, elementNumber);
    dofs = element.dofs(:);
    localEntryCount = numel(dofs)^2;
    indices = cursor:(cursor + localEntryCount - 1);
    [elementRows, elementColumns] = ndgrid(dofs, dofs);
    rows(indices) = elementRows(:);
    columns(indices) = elementColumns(:);
    stiffnessValues(indices) = element.stiffness(:);
    massValues(indices) = element.mass(:);
    cursor = cursor + localEntryCount;
end

globalK = sparse(rows, columns, stiffnessValues, ...
    numberOfDOFs, numberOfDOFs);
globalM = sparse(rows, columns, massValues, ...
    numberOfDOFs, numberOfDOFs);
end

function validateElementForAssembly(element, numberOfDOFs, elementNumber)
requiredFields = {'dofs', 'stiffness', 'mass'};
for i = 1:numel(requiredFields)
    if ~isfield(element, requiredFields{i})
        error('MKEF:InvalidAssemblyInput', ...
            'Element %d is missing field %s.', ...
            elementNumber, requiredFields{i});
    end
end
dofs = element.dofs(:);
if any(dofs < 1) || any(dofs > numberOfDOFs) || ...
        any(dofs ~= fix(dofs)) || numel(unique(dofs)) ~= numel(dofs) || ...
        ~isequal(size(element.stiffness), [numel(dofs), numel(dofs)]) || ...
        ~isequal(size(element.mass), [numel(dofs), numel(dofs)]) || ...
        any(~isfinite(element.stiffness(:))) || any(~isfinite(element.mass(:)))
    error('MKEF:InvalidAssemblyInput', ...
        'Element %d has invalid DOFs or matrix dimensions.', elementNumber);
end
end
