function test_result_recovery()
%TEST_RESULT_RECOVERY Verify member results against elementary solutions.

testSingleTrussRecovery();
testRotatedTrussRecovery();
testCantileverEndForces();
end

function testSingleTrussRecovery()
options = struct('verbose', false, 'plotting', false);
problem = StructFEProblem(fullfile('tests', 'fixtures', ...
    'CaseSingleTruss.txt'), options);
model = problem.GetAnalysisModel();
result = solveStatic(model);
element = result.elementResults(1);

force = 1000;
area = 0.01;
youngsModulus = 2e11;
expectedStress = force / area;
expectedStrain = expectedStress / youngsModulus;

assert(element.type == 112);
assert(isequal(element.nodeNumbers, [1 2]));
assert(isequal(element.globalDOFs, [1; 2; 3; 4]));
assertClose(element.axialStrain, expectedStrain, 1e-12, 1e-15, ...
    'Recovered truss strain is incorrect.');
assertClose(element.axialStress, expectedStress, 1e-12, 1e-9, ...
    'Recovered truss stress is incorrect.');
assertClose(element.axialForce, force, 1e-12, 1e-9, ...
    'Recovered truss axial force is incorrect.');
assertClose(element.localEndForces, [-force; 0; force; 0], ...
    1e-12, 1e-9, 'Recovered truss local end forces are incorrect.');

[stressFromCompatibilityFunction, compatibilityResults] = ...
    StressCalc(model, result.displacements);
assertClose(stressFromCompatibilityFunction, expectedStress, ...
    1e-12, 1e-9, 'StressCalc does not use the repaired recovery path.');
assertClose(compatibilityResults(1).axialForce, force, ...
    1e-12, 1e-9, 'StressCalc returned inconsistent element results.');
end

function testRotatedTrussRecovery()
coordinates = [1 2; 4 6];
area = 0.02;
youngsModulus = 210e9;
element = createStructuralElement(112, coordinates, [1 2], ...
    [area youngsModulus 7850], [1 2; 3 4]);

length = 5;
extension = 1e-3;
direction = (coordinates(2,:) - coordinates(1,:)) / length;
globalDisplacements = [0; 0; direction(1)*extension; direction(2)*extension];
model = struct();
model.numberOfDOFs = 4;
model.elementData = element;

result = recoverElementResults(model, globalDisplacements);
expectedStrain = extension / length;
assertClose(result.localDisplacements, [0; 0; extension; 0], ...
    1e-12, 1e-15, 'Rotated truss local displacements are incorrect.');
assertClose(result.axialStrain, expectedStrain, 1e-12, 1e-15, ...
    'Rotated truss strain does not use the element transformation.');
assertClose(result.axialForce, area * youngsModulus * expectedStrain, ...
    1e-12, 1e-8, 'Rotated truss axial force is incorrect.');
end

function testCantileverEndForces()
options = struct('verbose', false, 'plotting', false);
problem = StructFEProblem('Case1ElementBeam.txt', options);
result = problem.RunStatic();
element = result.elementResults(1);

force = 100;
length = 0.5;
expectedEndForces = [0; -force; -force*length; 0; force; 0];
assert(element.type == 113);
assert(isequal(element.globalDOFs, (1:6).'));
assertClose(element.localEndForces, expectedEndForces, ...
    1e-11, 1e-8, 'Recovered frame local end forces are incorrect.');
assertClose(element.localEndForces(1:3), result.reactions(1:3), ...
    1e-11, 1e-8, 'Frame fixed-end forces disagree with support reactions.');
assertClose(element.localEndForces(4:6), result.loadVector(4:6), ...
    1e-11, 1e-8, 'Frame free-end forces disagree with the applied nodal load.');
assertClose(element.axialForce, 0, 0, 1e-9, ...
    'Transversely loaded cantilever has a spurious axial force.');
end

function assertClose(actual, expected, relativeTolerance, absoluteTolerance, message)
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
    error('MKEF:VerificationFailed', '%s Error=%g, tolerance=%g.', ...
        message, difference, limit);
end
end
