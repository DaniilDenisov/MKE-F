function test_input_validation()
%TEST_INPUT_VALIDATION Verify actionable parsing and presentation checks.

fixtures = fullfile('tests', 'fixtures');
testRepeatedSections(fixtures);
testMalformedRecords(fixtures);
testMixedMeshRejection(fixtures);
testUnknownMarker(fixtures);
testFileOpenFailure(fixtures);
testPlotSelectionAndDefaults();
end

function testRepeatedSections(fixtures)
filename = fullfile(fixtures, 'CaseRepeatedSections.txt');
mesh = FEMesh(filename);
assert(mesh.numberOfFixBCs == 2);
assert(mesh.numberOfForceBCs == 7);
assert(size(mesh.allFixBCs, 1) == 2);
assert(isequal(size(mesh.allForceBCs), [7 6]));
assert(sum(mesh.allForceBCs(:, 3)) == 70);

problem = StructFEProblem(filename);
result = problem.RunStatic();
assert(abs(result.loadVector(4) - 70) < 1e-12);
assert(abs(result.reactions(1) + 70) < 1e-9);
end

function testMalformedRecords(fixtures)
assertThrows('MKEF:MalformedInput', ...
    @() FEMesh(fullfile(fixtures, 'CaseMalformedProperties.txt')), ':7:');
assertThrows('MKEF:MalformedInput', ...
    @() FEMesh(fullfile(fixtures, 'CaseInvalidNode.txt')), ...
    'invalid node ID 3');
end

function testMixedMeshRejection(fixtures)
assertThrows('MKEF:UnsupportedMixedMesh', ...
    @() FEMesh(fullfile(fixtures, 'CaseMixedElements.txt')), ':9:');
end

function testUnknownMarker(fixtures)
assertThrows('MKEF:UnknownInputMarker', ...
    @() FEMesh(fullfile(fixtures, 'CaseUnknownMarker.txt')), ':8:');
end

function testFileOpenFailure(fixtures)
assertThrows('MKEF:InputFileOpenFailed', ...
    @() FEMesh(fullfile(fixtures, 'file-that-does-not-exist.txt')), ...
    'Cannot open input file');
end

function testPlotSelectionAndDefaults()
problem = StructFEProblem('CaseBeamDyn.txt');
assert(~problem.verbose);
assert(~problem.plotting);
assertThrows('MKEF:InvalidPlotSelection', ...
    @() problem.RunTransient(1e-4, 2e-4, 0, 1), ...
    'integer indices within the model');
assertThrows('MKEF:InvalidPlotSelection', ...
    @() problem.RunTransient(1e-4, 2e-4, 1, 4), ...
    'integer indices within the model');
end

function assertThrows(expectedIdentifier, operation, messageFragment)
try
    operation();
catch exception
    if ~strcmp(exception.identifier, expectedIdentifier)
        error('MKEF:VerificationFailed', ...
            'Expected error %s, received %s.', ...
            expectedIdentifier, exception.identifier);
    end
    if isempty(strfind(exception.message, messageFragment)) %#ok<STREMP>
        error('MKEF:VerificationFailed', ...
            'Error message "%s" does not contain "%s".', ...
            exception.message, messageFragment);
    end
    return;
end
error('MKEF:VerificationFailed', 'Expected error %s.', expectedIdentifier);
end
