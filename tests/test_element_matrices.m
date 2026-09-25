function test_element_matrices()
%TEST_ELEMENT_MATRICES Verify 2D truss and frame element kernels.

testTrussOrientations();
testBeamOrientations();
testInvalidElementData();
end

function testTrussOrientations()
area = 0.02;
youngsModulus = 210e9;
density = 7850;
properties = [area youngsModulus density];
coordinates = {
    [0 0; 5 0]
    [0 0; 0 5]
    [1 2; 4 6]
};
referenceStiffnessSpectrum = [];
referenceMassSpectrum = [];

for i = 1:numel(coordinates)
    coords = coordinates{i};
    [K, M] = assembleTruss(coords, properties);
    [expectedK, expectedM, length] = expectedTrussMatrices(coords, properties);

    assertClose(K, expectedK, 1e-12, 1e-8, ...
        'Truss stiffness matrix has incorrect coefficients.');
    assertClose(M, expectedM, 1e-12, 1e-12, ...
        'Truss mass matrix has incorrect coefficients.');
    assertSymmetric(K, 'Truss stiffness matrix is not symmetric.');
    assertSymmetric(M, 'Truss mass matrix is not symmetric.');
    assertNumericalRank(K, 1, ...
        'A free truss element must have three zero-energy modes.');

    rigidX = [1; 0; 1; 0];
    rigidY = [0; 1; 0; 1];
    expectedMass = density * area * length;
    assertClose(rigidX.' * M * rigidX, expectedMass, 1e-12, 1e-12, ...
        'Truss total mass in global X is incorrect.');
    assertClose(rigidY.' * M * rigidY, expectedMass, 1e-12, 1e-12, ...
        'Truss total mass in global Y is incorrect.');

    stiffnessSpectrum = sort(eig((K + K.') / 2));
    massSpectrum = sort(eig((M + M.') / 2));
    if isempty(referenceStiffnessSpectrum)
        referenceStiffnessSpectrum = stiffnessSpectrum;
        referenceMassSpectrum = massSpectrum;
    else
        assertClose(stiffnessSpectrum, referenceStiffnessSpectrum, ...
            1e-12, 1e-6, 'Truss stiffness is not rotation invariant.');
        assertClose(massSpectrum, referenceMassSpectrum, ...
            1e-12, 1e-12, 'Truss mass is not rotation invariant.');
    end
end
end

function testBeamOrientations()
area = 0.02;
youngsModulus = 210e9;
density = 7850;
momentOfInertia = 8e-6;
properties = [area youngsModulus density momentOfInertia];
coordinates = {
    [0 0; 5 0]
    [0 0; 0 5]
    [1 2; 4 6]
};
referenceStiffnessSpectrum = [];
referenceMassSpectrum = [];

for i = 1:numel(coordinates)
    coords = coordinates{i};
    [K, M] = assembleBeam(coords, properties);
    [expectedK, expectedM, length] = expectedBeamMatrices(coords, properties);

    assertClose(K, expectedK, 1e-12, 1e-8, ...
        'Beam stiffness matrix has incorrect coefficients.');
    assertClose(M, expectedM, 1e-12, 1e-12, ...
        'Beam mass matrix has incorrect coefficients.');
    assertSymmetric(K, 'Beam stiffness matrix is not symmetric.');
    assertSymmetric(M, 'Beam mass matrix is not symmetric.');
    assertNumericalRank(K, 3, ...
        'A free beam element must have three rigid-body modes.');

    rigidX = [1; 0; 0; 1; 0; 0];
    rigidY = [0; 1; 0; 0; 1; 0];
    expectedMass = density * area * length;
    assertClose(rigidX.' * M * rigidX, expectedMass, 1e-12, 1e-12, ...
        'Beam total mass in global X is incorrect.');
    assertClose(rigidY.' * M * rigidY, expectedMass, 1e-12, 1e-12, ...
        'Beam total mass in global Y is incorrect.');

    stiffnessSpectrum = sort(eig((K + K.') / 2));
    massSpectrum = sort(eig((M + M.') / 2));
    if isempty(referenceStiffnessSpectrum)
        referenceStiffnessSpectrum = stiffnessSpectrum;
        referenceMassSpectrum = massSpectrum;
    else
        assertClose(stiffnessSpectrum, referenceStiffnessSpectrum, ...
            1e-12, 1e-5, 'Beam stiffness is not rotation invariant.');
        assertClose(massSpectrum, referenceMassSpectrum, ...
            1e-12, 1e-10, 'Beam mass is not rotation invariant.');
    end
end
end

function testInvalidElementData()
validCoords = [0 0; 1 0];

assertThrows('MKEF:InvalidElementGeometry', ...
    @() setupTruss([0 0; 0 0], [1 2e11 7850]));
assertThrows('MKEF:InvalidElementGeometry', ...
    @() setupBeam([0 0; 0 0], [1 2e11 7850 1]));
assertThrows('MKEF:InvalidElementGeometry', ...
    @() setupTruss([0 0; Inf 0], [1 2e11 7850]));
assertThrows('MKEF:InvalidElementGeometry', ...
    @() setupBeam([0 0; NaN 1], [1 2e11 7850 1]));

assertThrows('MKEF:InvalidElementProperties', ...
    @() setupTruss(validCoords, [0 2e11 7850]));
assertThrows('MKEF:InvalidElementProperties', ...
    @() setupTruss(validCoords, [1 -2e11 7850]));
assertThrows('MKEF:InvalidElementProperties', ...
    @() setupTruss(validCoords, [1 2e11 0]));
assertThrows('MKEF:InvalidElementProperties', ...
    @() setupBeam(validCoords, [1 2e11 7850 0]));
assertThrows('MKEF:InvalidElementProperties', ...
    @() setupBeam(validCoords, [1 2e11 Inf 1]));
end

function [K, M] = assembleTruss(coords, properties)
element = Truss2DElement();
element.SetupElement(coords, [1 2], properties);
[K, M] = element.Assembler(zeros(4), zeros(4), [1 2; 3 4]);
end

function [K, M] = assembleBeam(coords, properties)
element = Beam2DElement();
element.SetupElement(coords, [1 2], properties);
[K, M] = element.Assembler(zeros(6), zeros(6), [1 2 3; 4 5 6]);
end

function setupTruss(coords, properties)
element = Truss2DElement();
element.SetupElement(coords, [1 2], properties);
end

function setupBeam(coords, properties)
element = Beam2DElement();
element.SetupElement(coords, [1 2], properties);
end

function [K, M, length] = expectedTrussMatrices(coords, properties)
delta = coords(2,:) - coords(1,:);
length = norm(delta);
c = delta(1) / length;
s = delta(2) / length;
area = properties(1);
youngsModulus = properties(2);
density = properties(3);

K = area * youngsModulus / length * ...
    [ c*c  c*s -c*c -c*s; ...
      c*s  s*s -c*s -s*s; ...
     -c*c -c*s  c*c  c*s; ...
     -c*s -s*s  c*s  s*s];
M = density * area * length / 6 * ...
    [2 0 1 0; 0 2 0 1; 1 0 2 0; 0 1 0 2];
end

function [K, M, length] = expectedBeamMatrices(coords, properties)
delta = coords(2,:) - coords(1,:);
length = norm(delta);
c = delta(1) / length;
s = delta(2) / length;
area = properties(1);
youngsModulus = properties(2);
density = properties(3);
momentOfInertia = properties(4);

axial = youngsModulus * area / length;
bending12 = 12 * youngsModulus * momentOfInertia / length^3;
bending6 = 6 * youngsModulus * momentOfInertia / length^2;
bending4 = 4 * youngsModulus * momentOfInertia / length;
bending2 = 2 * youngsModulus * momentOfInertia / length;
localK = [ axial       0          0 -axial        0          0; ...
               0 bending12   bending6      0 -bending12   bending6; ...
               0  bending6   bending4      0  -bending6   bending2; ...
          -axial       0          0  axial        0          0; ...
               0 -bending12 -bending6      0  bending12  -bending6; ...
               0  bending6   bending2      0  -bending6   bending4];
localM = density * area * length / 420 * ...
    [140       0            0  70       0            0; ...
       0     156    22*length   0      54   -13*length; ...
       0 22*length 4*length^2   0 13*length -3*length^2; ...
      70       0            0 140       0            0; ...
       0      54    13*length   0     156   -22*length; ...
       0 -13*length -3*length^2 0 -22*length 4*length^2];
transform = [ c s 0 0 0 0; -s c 0 0 0 0; 0 0 1 0 0 0; ...
              0 0 0 c s 0; 0 0 0 -s c 0; 0 0 0 0 0 1];
K = transform.' * localK * transform;
M = transform.' * localM * transform;
end

function assertSymmetric(matrix, message)
errorNorm = norm(matrix - matrix.', inf);
scale = max(norm(matrix, inf), 1);
if errorNorm > 1e-12 * scale
    error('MKEF:VerificationFailed', '%s Error=%g, scale=%g.', ...
        message, errorNorm, scale);
end
end

function assertNumericalRank(matrix, expectedRank, message)
singularValues = svd(matrix);
tolerance = max(singularValues) * 1e-10;
actualRank = sum(singularValues > tolerance);
if actualRank ~= expectedRank
    error('MKEF:VerificationFailed', '%s Rank=%d, expected=%d.', ...
        message, actualRank, expectedRank);
end
end

function assertThrows(expectedIdentifier, operation)
try
    operation();
catch exception
    if strcmp(exception.identifier, expectedIdentifier)
        return;
    end
    error('MKEF:VerificationFailed', ...
        'Expected error %s, received %s.', ...
        expectedIdentifier, exception.identifier);
end
error('MKEF:VerificationFailed', 'Expected error %s.', expectedIdentifier);
end

function assertClose(actual, expected, relativeTolerance, absoluteTolerance, message)
difference = max(abs(actual(:) - expected(:)));
scale = max(abs(expected(:)));
limit = absoluteTolerance + relativeTolerance * scale;
if difference > limit
    error('MKEF:VerificationFailed', '%s Error=%g, tolerance=%g.', ...
        message, difference, limit);
end
end
