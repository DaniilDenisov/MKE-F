function test_postprocessor_export()
%TEST_POSTPROCESSOR_EXPORT Verify schema-v1 conversion and safe file output.

options = struct('verbose', false, 'plotting', false);
problem = StructFEProblem(exampleCasePath('Case1ElementBeam.txt'), options);
model = problem.GetAnalysisModel();
modelBefore = model;

staticResult = problem.RunStatic();
staticBefore = staticResult;
data = createPostprocessorData(model, staticResult, ...
    struct('title', 'Cantilever', 'lengthUnit', 'm', ...
           'forceUnit', 'N', 'momentUnit', 'N*m', 'stressUnit', 'Pa'));
assert(strcmp(data.format, 'mkef-postprocessor'));
assert(data.version == 1);
assert(strcmp(data.analysis.type, 'static'));
assert(strcmp(data.metadata.title, 'Cantilever'));
assert(strcmp(data.metadata.units.length, 'm'));
assert(isequal(data.model.dofMap, model.dofMap));
assert(isequal([data.model.nodes.id], 1:model.numberOfNodes));
assert(data.model.elements(1).type == 113);
assert(~isfield(data.model, 'stiffness'));
assert(~isfield(data.model.elements(1), 'localStiffness'));
assertClose(data.analysis.displacements, staticResult.displacements);
assertClose(data.analysis.elementResults(1).localEndForces, ...
    staticResult.elementResults(1).localEndForces);
assert(isequal(model, modelBefore));
assert(isequal(staticResult, staticBefore));

trussProblem = StructFEProblem(fullfile('tests', 'fixtures', ...
    'CaseSingleTruss.txt'), options);
trussData = createPostprocessorData(trussProblem.GetAnalysisModel(), ...
    trussProblem.RunStatic());
assert(trussData.model.elements(1).type == 112);
assert(isfield(trussData.model.elements(1).properties, 'area'));
assert(~isfield(trussData.model.elements(1).properties, 'momentOfInertia'));
assertClose(trussData.analysis.elementResults(1).axialForce, 1000);

frameProblem = StructFEProblem( ...
    exampleCasePath('ANSYSBeamStatic01.txt'), options);
frameData = createPostprocessorData(frameProblem.GetAnalysisModel(), ...
    frameProblem.RunStatic());
assert(numel(frameData.model.nodes) == frameProblem.mesh.numberOfNodes);
assert(numel(frameData.model.elements) == frameProblem.mesh.numberOfElems);

modalResult = problem.RunModal();
modalData = createPostprocessorData(model, modalResult);
assert(strcmp(modalData.analysis.type, 'modal'));
assert(isequal(size(modalData.analysis.modeShapes), ...
    size(modalResult.modeShapes)));
assertClose(modalData.analysis.modeShapes, modalResult.modeShapes);

transientResult = problem.RunTransient(1e-4, 5e-4, 2, 2);
transientOptions = struct('transientFields', ...
    {{'displacements', 'reactions', 'spectrum'}}, ...
    'selectedGlobalDOFs', [6 5], 'timeStride', 2);
transientData = createPostprocessorData(model, transientResult, transientOptions);
expectedIndices = [1 3 5 6];
assert(isequal(transientData.analysis.globalDOFIds, [6 5]));
assert(isequal(transientData.analysis.time, ...
    transientResult.time(expectedIndices)));
assertClose(transientData.analysis.displacements, ...
    transientResult.displacements([6 5], expectedIndices));
assertClose(transientData.analysis.reactions, ...
    transientResult.reactions([6 5], expectedIndices));
assertClose(transientData.analysis.displacementAmplitudeSpectrum, ...
    transientResult.displacementAmplitudeSpectrum([6 5], :));
assert(~isfield(transientData.analysis, 'velocities'));
assert(transientData.analysis.sampling.originalSampleCount == 6);
assert(transientData.analysis.sampling.exportedSampleCount == 4);
assert(transientData.analysis.sampling.timeStride == 2);

temporaryDirectory = tempname();
mkdir(temporaryDirectory);
cleanup = onCleanup(@() rmdir(temporaryDirectory, 's'));
filename = fullfile(temporaryDirectory, 'result.json');
writtenData = exportPostprocessorData(model, staticResult, filename, ...
    struct('prettyPrint', true));
decoded = jsondecode(fileread(filename));
assert(strcmp(decoded.format, 'mkef-postprocessor'));
assert(decoded.version == 1);
assertClose(decoded.analysis.displacements, staticResult.displacements);
assert(isequal(writtenData.analysis.displacements, data.analysis.displacements));
jsonText = fileread(filename);
assertJsonArray(jsonText, 'nodes');
assertJsonArray(jsonText, 'elements');
assertJsonArray(jsonText, 'nodeIds');
assertJsonArray(jsonText, 'dofMap');
assertJsonArray(jsonText, 'supports');
assertJsonArray(jsonText, 'displacements');
assertJsonArray(jsonText, 'elementResults');
assertJsonArray(jsonText, 'localEndForces');

frameFilename = fullfile(temporaryDirectory, 'static-frame.json');
exportPostprocessorData(frameProblem.GetAnalysisModel(), ...
    frameProblem.RunStatic(), frameFilename);
frameJson = fileread(frameFilename);
assertJsonArray(frameJson, 'supports');

singleModeResult = modalResult;
singleModeResult.frequenciesHz = modalResult.frequenciesHz(1);
singleModeResult.angularFrequenciesRadPerSec = ...
    modalResult.angularFrequenciesRadPerSec(1);
singleModeResult.modeShapes = modalResult.modeShapes(:, 1);
modalFilename = fullfile(temporaryDirectory, 'single-mode.json');
exportPostprocessorData(model, singleModeResult, modalFilename);
modalJson = fileread(modalFilename);
assertJsonArray(modalJson, 'frequenciesHz');
assertJsonArray(modalJson, 'angularFrequenciesRadPerSec');
assertJsonMatrix(modalJson, 'modeShapes');

singleDofOptions = transientOptions;
singleDofOptions.selectedGlobalDOFs = 6;
transientFilename = fullfile(temporaryDirectory, 'single-dof.json');
exportPostprocessorData(model, transientResult, transientFilename, ...
    singleDofOptions);
transientJson = fileread(transientFilename);
assertJsonArray(transientJson, 'time');
assertJsonArray(transientJson, 'globalDOFIds');
assertJsonMatrix(transientJson, 'displacements');
assertJsonMatrix(transientJson, 'reactions');
assertJsonMatrix(transientJson, 'displacementAmplitudeSpectrum');

fid = fopen(filename, 'wb');
fwrite(fid, uint8('preserve-me'), 'uint8');
fclose(fid);
invalidResult = staticResult;
invalidResult.displacements(1) = Inf;
assertThrows('MKEF:PostprocessorInvalidResult', ...
    @() exportPostprocessorData(model, invalidResult, filename));
assert(strcmp(fileread(filename), 'preserve-me'));

invalidModel = model;
invalidModel.nodeCoordinates(1, 1) = NaN;
assertThrows('MKEF:PostprocessorInvalidModel', ...
    @() createPostprocessorData(invalidModel, staticResult));
assertThrows('MKEF:PostprocessorInvalidOptions', ...
    @() createPostprocessorData(model, staticResult, struct('timeStride', 2)));
assertThrows('MKEF:PostprocessorInvalidOptions', ...
    @() createPostprocessorData(model, transientResult, ...
        struct('selectedGlobalDOFs', model.numberOfDOFs + 1)));
assertThrows('MKEF:PostprocessorInvalidOptions', ...
    @() createPostprocessorData(model, staticResult, struct('typo', true)));

clear cleanup;
end

function assertClose(actual, expected)
difference = max(abs(actual(:) - expected(:)));
if isempty(difference)
    difference = 0;
end
assert(difference <= 1e-12 * max([max(abs(expected(:))), 1]));
end

function assertThrows(identifier, operation)
try
    operation();
catch exception
    assert(strcmp(exception.identifier, identifier), ...
        ['Expected ' identifier ', received ' exception.identifier '.']);
    return;
end
error('MKEF:VerificationFailed', 'Expected error %s.', identifier);
end

function assertJsonArray(json, fieldName)
pattern = ['"' fieldName '"\s*:\s*\['];
assert(~isempty(regexp(json, pattern, 'once')), ...
    ['Expected JSON array for field ' fieldName '.']);
end

function assertJsonMatrix(json, fieldName)
pattern = ['"' fieldName '"\s*:\s*\[\s*\['];
assert(~isempty(regexp(json, pattern, 'once')), ...
    ['Expected nested JSON arrays for field ' fieldName '.']);
end
