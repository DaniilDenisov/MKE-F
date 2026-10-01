function element = createStructuralElement(type, nodeCoordinates, ...
        nodeNumbers, properties, dofMap)
%CREATESTRUCTURALELEMENT Build the common immutable element-data interface.
% Fields used by assembly and result recovery are:
%   type, nodeNumbers, nodeCoordinates, properties, dofs, length,
%   transformation, localStiffness, stiffness, and mass.

validateGeometry(nodeCoordinates, type);
validateTopology(nodeNumbers, dofMap, type);
validateProperties(properties, type);

nodeCoordinates = nodeCoordinates(:, 1:2);
nodeNumbers = nodeNumbers(:).';
properties = properties(:).';
delta = nodeCoordinates(2, :) - nodeCoordinates(1, :);
length = norm(delta);
c = delta(1) / length;
s = delta(2) / length;

switch type
    case 112
        transformation = [c s 0 0; -s c 0 0; ...
                          0 0 c s; 0 0 -s c];
        [localStiffness, localMass] = ...
            trussMatrices(properties, length);
    case 113
        transformation = [c s 0 0 0 0; ...
                         -s c 0 0 0 0; ...
                          0 0 1 0 0 0; ...
                          0 0 0 c s 0; ...
                          0 0 0 -s c 0; ...
                          0 0 0 0 0 1];
        [localStiffness, localMass] = ...
            frameMatrices(properties, length);
    otherwise
        error('MKEF:UnsupportedElementType', ...
            'Unsupported structural element type %g.', type);
end

element = struct();
element.type = type;
element.nodeNumbers = nodeNumbers;
element.nodeCoordinates = nodeCoordinates;
element.properties = properties;
element.dofs = reshape(dofMap(nodeNumbers, :).', [], 1);
element.length = length;
element.transformation = transformation;
element.localStiffness = localStiffness;
element.stiffness = transformation.' * localStiffness * transformation;
element.mass = transformation.' * localMass * transformation;
end

function validateGeometry(nodeCoordinates, type)
if ~isnumeric(nodeCoordinates) || ~isreal(nodeCoordinates) || ...
        size(nodeCoordinates, 1) ~= 2 || size(nodeCoordinates, 2) < 2 || ...
        any(~isfinite(nodeCoordinates(:, 1:2))(:))
    error('MKEF:InvalidElementGeometry', ...
        'Element type %g requires two finite XY node coordinates.', type);
end
if norm(nodeCoordinates(2, 1:2) - nodeCoordinates(1, 1:2)) == 0
    error('MKEF:InvalidElementGeometry', ...
        'Element length must be greater than zero.');
end
end

function validateTopology(nodeNumbers, dofMap, type)
expectedDOFs = 2 + (type == 113);
if ~ismember(type, [112 113])
    error('MKEF:UnsupportedElementType', ...
        'Unsupported structural element type %g.', type);
end
if ~isnumeric(nodeNumbers) || numel(nodeNumbers) ~= 2 || ...
        any(nodeNumbers ~= fix(nodeNumbers)) || any(nodeNumbers < 1) || ...
        any(nodeNumbers > size(dofMap, 1)) || size(dofMap, 2) ~= expectedDOFs
    error('MKEF:InvalidElementTopology', ...
        'Element type %g has invalid nodes or DOF mapping.', type);
end
end

function validateProperties(properties, type)
propertyCount = 3 + (type == 113);
if ~isnumeric(properties) || ~isreal(properties) || ...
        numel(properties) < propertyCount || ...
        any(~isfinite(properties(1:propertyCount))) || ...
        any(properties(1:propertyCount) <= 0)
    if type == 112
        description = 'Truss properties A, E, and rho';
    else
        description = 'Frame properties A, E, rho, and I';
    end
    error('MKEF:InvalidElementProperties', ...
        '%s must be positive finite values.', description);
end
end

function [K, M] = trussMatrices(properties, length)
area = properties(1);
youngsModulus = properties(2);
density = properties(3);
axial = area * youngsModulus / length;
K = [ axial 0 -axial 0; ...
          0 0      0 0; ...
     -axial 0  axial 0; ...
          0 0      0 0];
M = density * area * length / 6 * ...
    [2 0 1 0; 0 2 0 1; 1 0 2 0; 0 1 0 2];
end

function [K, M] = frameMatrices(properties, length)
area = properties(1);
youngsModulus = properties(2);
density = properties(3);
momentOfInertia = properties(4);
axial = youngsModulus * area / length;
bending12 = 12 * youngsModulus * momentOfInertia / length^3;
bending6 = 6 * youngsModulus * momentOfInertia / length^2;
bending4 = 4 * youngsModulus * momentOfInertia / length;
bending2 = 2 * youngsModulus * momentOfInertia / length;
K = [ axial       0          0 -axial        0          0; ...
           0 bending12   bending6      0 -bending12   bending6; ...
           0  bending6   bending4      0  -bending6   bending2; ...
      -axial       0          0  axial        0          0; ...
           0 -bending12 -bending6      0  bending12  -bending6; ...
           0  bending6   bending2      0  -bending6   bending4];
M = density * area * length / 420 * ...
    [140          0               0  70          0               0; ...
       0        156       22*length   0         54      -13*length; ...
       0  22*length     4*length^2   0  13*length     -3*length^2; ...
      70          0               0 140          0               0; ...
       0         54       13*length   0        156      -22*length; ...
       0 -13*length    -3*length^2   0 -22*length      4*length^2];
end
