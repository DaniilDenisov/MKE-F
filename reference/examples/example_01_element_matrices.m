function result = example_01_element_matrices(verbose)
%EXAMPLE_01_ELEMENT_MATRICES Inspect real truss and frame element matrices.
%   RESULT = EXAMPLE_01_ELEMENT_MATRICES() prints the important matrices,
%   checks their invariants, and returns the element data for experiments.
%   Pass false to run the example without diagnostic output.

if nargin < 1
    verbose = true;
end
validateVerbose(verbose);
addRepositoryRoot();

coordinates = [0 0; 3 4];
nodeNumbers = [1 2];
trussDOFMap = [1 2; 3 4];
frameDOFMap = [1 2 3; 4 5 6];

trussProperties = [0.01 200e9 7850];
frameProperties = [0.02 210e9 7850 8e-6];

truss = createStructuralElement(112, coordinates, nodeNumbers, ...
    trussProperties, trussDOFMap);
frame = createStructuralElement(113, coordinates, nodeNumbers, ...
    frameProperties, frameDOFMap);

assertClose(truss.length, 5, 0, 1e-14, ...
    'The truss length must equal five metres.');
assertClose(frame.length, 5, 0, 1e-14, ...
    'The frame length must equal five metres.');
assertSymmetric(truss.localStiffness, ...
    'The local truss stiffness matrix is not symmetric.');
assertSymmetric(truss.stiffness, ...
    'The global truss stiffness matrix is not symmetric.');
assertSymmetric(truss.mass, ...
    'The global truss mass matrix is not symmetric.');
assertSymmetric(frame.localStiffness, ...
    'The local frame stiffness matrix is not symmetric.');
assertSymmetric(frame.stiffness, ...
    'The global frame stiffness matrix is not symmetric.');
assertSymmetric(frame.mass, ...
    'The global frame mass matrix is not symmetric.');
assertClose(truss.transformation * truss.transformation.', eye(4), ...
    1e-12, 1e-12, 'The truss transformation must be orthogonal.');
assertClose(frame.transformation * frame.transformation.', eye(6), ...
    1e-12, 1e-12, 'The frame transformation must be orthogonal.');
assertNumericalRank(truss.stiffness, 1, ...
    'A free planar truss element must have one deformational mode.');
assertNumericalRank(frame.stiffness, 3, ...
    'A free planar frame element must have three deformational modes.');

trussRigidX = [1; 0; 1; 0];
trussRigidY = [0; 1; 0; 1];
frameRigidX = [1; 0; 0; 1; 0; 0];
frameRigidY = [0; 1; 0; 0; 1; 0];
trussMass = trussProperties(3) * trussProperties(1) * truss.length;
frameMass = frameProperties(3) * frameProperties(1) * frame.length;
assertClose(trussRigidX.' * truss.mass * trussRigidX, trussMass, ...
    1e-12, 1e-12, 'The truss translational mass is incorrect in X.');
assertClose(trussRigidY.' * truss.mass * trussRigidY, trussMass, ...
    1e-12, 1e-12, 'The truss translational mass is incorrect in Y.');
assertClose(frameRigidX.' * frame.mass * frameRigidX, frameMass, ...
    1e-12, 1e-12, 'The frame translational mass is incorrect in X.');
assertClose(frameRigidY.' * frame.mass * frameRigidY, frameMass, ...
    1e-12, 1e-12, 'The frame translational mass is incorrect in Y.');

if verbose
    fprintf('\nExample 01: element matrices\n');
    fprintf('Coordinates [m]:\n');
    disp(coordinates);
    fprintf('Length = %.6g m, direction cosines c = %.6g, s = %.6g\n', ...
        truss.length, truss.transformation(1, 1), ...
        truss.transformation(1, 2));
    fprintf('\nTruss DOFs [u1x u1y u2x u2y]:\n');
    disp(truss.dofs.');
    fprintf('Truss transformation T:\n');
    disp(truss.transformation);
    fprintf('Truss local stiffness K_local [N/m]:\n');
    disp(truss.localStiffness);
    fprintf('Truss global stiffness T''*K_local*T [N/m]:\n');
    disp(truss.stiffness);
    fprintf('Truss consistent mass [kg]:\n');
    disp(truss.mass);
    fprintf('\nFrame DOFs [u1x u1y theta1 u2x u2y theta2]:\n');
    disp(frame.dofs.');
    fprintf('Frame transformation T:\n');
    disp(frame.transformation);
    fprintf('Frame local stiffness K_local:\n');
    disp(frame.localStiffness);
    fprintf('Frame global stiffness T''*K_local*T:\n');
    disp(frame.stiffness);
    fprintf('Frame consistent mass:\n');
    disp(frame.mass);
    fprintf('Checks: symmetry, orthogonality, rank, and total mass PASS.\n');
end

result = struct();
result.coordinates = coordinates;
result.truss = truss;
result.frame = frame;
result.expectedRanks = struct('truss', 1, 'frame', 3);
result.physicalMasses = struct('truss', trussMass, 'frame', frameMass);
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
