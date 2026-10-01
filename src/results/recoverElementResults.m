function elementResults = recoverElementResults(model, displacements)
%RECOVERELEMENTRESULTS Recover local element deformation and force results.
% Positive axial strain, stress, and force denote tension. localEndForces
% use the local element DOF order and are forces exerted on the element.

if ~isfield(model, 'elementData') || ~isstruct(model.elementData)
    error('MKEF:InvalidRecoveryModel', ...
        'The analysis model does not contain element recovery data.');
end
if ~isnumeric(displacements) || ~isreal(displacements) || ...
        ~isvector(displacements) || numel(displacements) ~= model.numberOfDOFs || ...
        any(~isfinite(displacements(:)))
    error('MKEF:InvalidRecoveryDisplacements', ...
        'Displacements must be a finite vector with one value per model DOF.');
end
displacements = displacements(:);

emptyResult = struct( ...
    'elementNumber', [], ...
    'type', [], ...
    'nodeNumbers', [], ...
    'globalDOFs', [], ...
    'length', [], ...
    'localDisplacements', [], ...
    'localEndForces', [], ...
    'axialStrain', [], ...
    'axialStress', [], ...
    'axialForce', []);
elementResults = repmat(emptyResult, numel(model.elementData), 1);

for elementNumber = 1:numel(model.elementData)
    data = model.elementData(elementNumber);
    validateElementData(data, model.numberOfDOFs, elementNumber);
    globalDisplacements = displacements(data.dofs);
    localDisplacements = data.transformation * globalDisplacements;
    localEndForces = data.localStiffness * localDisplacements;

    switch data.type
        case 112
            endAxialDOF = 3;
        case 113
            endAxialDOF = 4;
        otherwise
            error('MKEF:UnsupportedElementRecovery', ...
                'Element %d has unsupported type %g.', elementNumber, data.type);
    end

    axialStrain = ...
        (localDisplacements(endAxialDOF) - localDisplacements(1)) / data.length;
    axialStress = data.properties(2) * axialStrain;

    elementResults(elementNumber).elementNumber = elementNumber;
    elementResults(elementNumber).type = data.type;
    elementResults(elementNumber).nodeNumbers = data.nodeNumbers;
    elementResults(elementNumber).globalDOFs = data.dofs;
    elementResults(elementNumber).length = data.length;
    elementResults(elementNumber).localDisplacements = localDisplacements;
    elementResults(elementNumber).localEndForces = localEndForces;
    elementResults(elementNumber).axialStrain = axialStrain;
    elementResults(elementNumber).axialStress = axialStress;
    elementResults(elementNumber).axialForce = data.properties(1) * axialStress;
end
end

function validateElementData(data, numberOfDOFs, elementNumber)
requiredFields = {'type', 'nodeNumbers', 'dofs', 'length', ...
    'transformation', 'localStiffness', 'properties'};
for i = 1:numel(requiredFields)
    if ~isfield(data, requiredFields{i})
        error('MKEF:InvalidRecoveryModel', ...
            'Element %d is missing recovery field %s.', ...
            elementNumber, requiredFields{i});
    end
end
if data.length <= 0 || ~isfinite(data.length) || ...
        any(data.dofs < 1) || any(data.dofs > numberOfDOFs) || ...
        any(data.dofs ~= fix(data.dofs)) || ...
        size(data.transformation, 2) ~= numel(data.dofs) || ...
        size(data.localStiffness, 1) ~= size(data.transformation, 1) || ...
        size(data.localStiffness, 2) ~= size(data.transformation, 1)
    error('MKEF:InvalidRecoveryModel', ...
        'Element %d contains inconsistent recovery data.', elementNumber);
end
end
