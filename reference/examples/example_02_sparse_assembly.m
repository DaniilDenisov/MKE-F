function result = example_02_sparse_assembly(verbose)
%EXAMPLE_02_SPARSE_ASSEMBLY Assemble two truss elements into sparse matrices.
%   RESULT = EXAMPLE_02_SPARSE_ASSEMBLY() prints the DOF maps, element
%   contributions, and assembled matrices. Pass false for a quiet run.

if nargin < 1
    verbose = true;
end
validateVerbose(verbose);
addRepositoryRoot();

nodeCoordinates = [0 0; 1 0; 2 0];
dofMap = reshape(1:6, 2, 3).';
properties = [0.01 200e9 7850];
connectivity = [1 2; 2 3];

elementCells = cell(2, 1);
for elementNumber = 1:2
    nodes = connectivity(elementNumber, :);
    elementCells{elementNumber} = createStructuralElement(112, ...
        nodeCoordinates(nodes, :), nodes, properties, dofMap);
end
elements = vertcat(elementCells{:});
[globalK, globalM] = assembleGlobalMatrices(elements, 6);

% This dense scatter sum is an independent, transparent check of the
% production sparse triplet assembler. Element matrices still come from
% createStructuralElement; no educational solver is reimplemented here.
expectedK = zeros(6, 6);
expectedM = zeros(6, 6);
for elementNumber = 1:numel(elements)
    dofs = elements(elementNumber).dofs;
    expectedK(dofs, dofs) = expectedK(dofs, dofs) + ...
        elements(elementNumber).stiffness;
    expectedM(dofs, dofs) = expectedM(dofs, dofs) + ...
        elements(elementNumber).mass;
end

assert(issparse(globalK), ...
    'MKEF:ReferenceExampleFailed: stiffness matrix is not sparse.');
assert(issparse(globalM), ...
    'MKEF:ReferenceExampleFailed: mass matrix is not sparse.');
assertClose(full(globalK), expectedK, 1e-12, 1e-8, ...
    'Sparse stiffness assembly differs from the scatter sum.');
assertClose(full(globalM), expectedM, 1e-12, 1e-12, ...
    'Sparse mass assembly differs from the scatter sum.');
assertSymmetric(globalK, 'The assembled stiffness matrix is not symmetric.');
assertSymmetric(globalM, 'The assembled mass matrix is not symmetric.');
assertNumericalRank(globalK, 2, ...
    'The two-bar chain must have two independent axial deformation modes.');

sharedXDOF = dofMap(2, 1);
sharedYDOF = dofMap(2, 2);
expectedSharedK = elements(1).stiffness(3, 3) + ...
    elements(2).stiffness(1, 1);
expectedSharedM = elements(1).mass(3, 3) + ...
    elements(2).mass(1, 1);
assertClose(globalK(sharedXDOF, sharedXDOF), expectedSharedK, ...
    1e-12, 1e-8, 'Shared-node stiffness contributions were not added.');
assertClose(globalM(sharedXDOF, sharedXDOF), expectedSharedM, ...
    1e-12, 1e-12, 'Shared-node mass contributions were not added.');
assert(globalK(sharedYDOF, sharedYDOF) == 0, ...
    'MKEF:ReferenceExampleFailed: a horizontal truss has transverse stiffness.');

rigidX = [1; 0; 1; 0; 1; 0];
rigidY = [0; 1; 0; 1; 0; 1];
totalPhysicalMass = properties(3) * properties(1) * 2;
assertClose(rigidX.' * globalM * rigidX, totalPhysicalMass, ...
    1e-12, 1e-12, 'The chain mass is incorrect in X translation.');
assertClose(rigidY.' * globalM * rigidY, totalPhysicalMass, ...
    1e-12, 1e-12, 'The chain mass is incorrect in Y translation.');

if verbose
    fprintf('\nExample 02: sparse global assembly\n');
    fprintf('Node coordinates [m]:\n');
    disp(nodeCoordinates);
    fprintf('Global DOF map [ux uy]:\n');
    disp(dofMap);
    for elementNumber = 1:numel(elements)
        fprintf('Element %d, nodes [%d %d], global DOFs:\n', ...
            elementNumber, connectivity(elementNumber, 1), ...
            connectivity(elementNumber, 2));
        disp(elements(elementNumber).dofs.');
        fprintf('Element %d global stiffness contribution:\n', elementNumber);
        disp(elements(elementNumber).stiffness);
    end
    fprintf('Assembled sparse K (shown as a full matrix) [N/m]:\n');
    disp(full(globalK));
    fprintf('Assembled sparse M (shown as a full matrix) [kg]:\n');
    disp(full(globalM));
    fprintf('Stored nonzeros: nnz(K) = %d, nnz(M) = %d.\n', ...
        nnz(globalK), nnz(globalM));
    fprintf('At shared DOF %d: K(%d,%d) = %.6g, M(%d,%d) = %.6g.\n', ...
        sharedXDOF, sharedXDOF, sharedXDOF, ...
        globalK(sharedXDOF, sharedXDOF), sharedXDOF, sharedXDOF, ...
        globalM(sharedXDOF, sharedXDOF));
    fprintf('Checks: sparse storage, scatter equality, symmetry, rank, and mass PASS.\n');
end

result = struct();
result.nodeCoordinates = nodeCoordinates;
result.connectivity = connectivity;
result.dofMap = dofMap;
result.elements = elements;
result.stiffness = globalK;
result.mass = globalM;
result.sharedDOFs = dofMap(2, :);
result.totalPhysicalMass = totalPhysicalMass;
end

function addRepositoryRoot()
exampleDirectory = fileparts(mfilename('fullpath'));
referenceDirectory = fileparts(exampleDirectory);
repositoryRoot = fileparts(referenceDirectory);
addpath(repositoryRoot);
end

function validateVerbose(verbose)
if ~isscalar(verbose) || ...
        ~(islogical(verbose) || (isnumeric(verbose) && isfinite(verbose))) || ...
        ~ismember(double(verbose), [0 1])
    error('MKEF:InvalidExampleOption', ...
        'verbose must be a scalar logical value.');
end
end

function assertSymmetric(matrix, message)
errorNorm = norm(matrix - matrix.', inf);
scale = max(norm(matrix, inf), 1);
if errorNorm > 1e-12 * scale
    error('MKEF:ReferenceExampleFailed', ...
        '%s Error=%g, scale=%g.', message, errorNorm, scale);
end
end

function assertNumericalRank(matrix, expectedRank, message)
singularValues = svd(full(matrix));
tolerance = max(singularValues) * 1e-10;
actualRank = sum(singularValues > tolerance);
if actualRank ~= expectedRank
    error('MKEF:ReferenceExampleFailed', ...
        '%s Rank=%d, expected=%d.', message, actualRank, expectedRank);
end
end

function assertClose(actual, expected, relativeTolerance, ...
        absoluteTolerance, message)
difference = max(abs(actual(:) - expected(:)));
scale = max(abs(expected(:)));
if isempty(difference)
    difference = 0;
end
if isempty(scale)
    scale = 0;
end
limit = absoluteTolerance + relativeTolerance * scale;
if difference > limit
    error('MKEF:ReferenceExampleFailed', ...
        '%s Error=%g, tolerance=%g.', message, difference, limit);
end
end
