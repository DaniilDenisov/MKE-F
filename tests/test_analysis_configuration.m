function test_analysis_configuration()
%TEST_ANALYSIS_CONFIGURATION Verify configured and legacy case execution.

staticProblem = StructFEProblem(exampleCasePath( ...
    'CasePreprocessorStatic.txt'));
assert(strcmp(staticProblem.mesh.analysisConfiguration.type, 'static'));
selectedStatic = staticProblem.RunSelected();
explicitStatic = staticProblem.RunStatic();
assertClose(selectedStatic.displacements, explicitStatic.displacements, ...
    'Configured static result differs from RunStatic.');

modalProblem = StructFEProblem(exampleCasePath( ...
    'CasePreprocessorModal.txt'));
assert(strcmp(modalProblem.mesh.analysisConfiguration.type, 'modal'));
selectedModal = modalProblem.RunSelected();
explicitModal = modalProblem.RunModal();
assertClose(selectedModal.frequenciesHz, explicitModal.frequenciesHz, ...
    'Configured modal result differs from RunModal.');

transientProblem = StructFEProblem(exampleCasePath( ...
    'CasePreprocessorTransient.txt'));
configuration = transientProblem.mesh.analysisConfiguration;
assert(strcmp(configuration.type, 'transient'));
assert(configuration.timeStep == 1e-4);
assert(configuration.duration == 1e-3);
assert(configuration.monitorNode == 3);
assert(configuration.monitorDOF == 2);
selectedTransient = transientProblem.RunSelected();
explicitTransient = transientProblem.RunTransient(1e-4, 1e-3, 3, 2);
assertClose(selectedTransient.displacements, ...
    explicitTransient.displacements, ...
    'Configured transient result differs from RunTransient.');

legacyProblem = StructFEProblem(exampleCasePath('Case1ElementBeam.txt'));
assert(isempty(fieldnames(legacyProblem.mesh.analysisConfiguration)));
assertThrows('MKEF:MissingAnalysisConfiguration', ...
    @() legacyProblem.RunSelected());

fixtures = fullfile('tests', 'fixtures');
assertThrows('MKEF:DuplicateInputSection', ...
    @() FEMesh(fullfile(fixtures, 'CaseDuplicateAnalysis.txt')));
assertThrows('MKEF:MalformedInput', ...
    @() FEMesh(fullfile(fixtures, 'CaseInvalidAnalysis.txt')));
assertThrows('MKEF:MalformedInput', ...
    @() FEMesh(fullfile(fixtures, 'CaseInvalidMonitor.txt')));
assertThrows('MKEF:MalformedInput', ...
    @() FEMesh(fullfile(fixtures, 'CaseInvalidMonitorDOF.txt')));
assertThrows('MKEF:AnalysisLoadMismatch', ...
    @() FEMesh(fullfile(fixtures, 'CaseModalWithLoad.txt')));
assertThrows('MKEF:AnalysisLoadMismatch', ...
    @() FEMesh(fullfile(fixtures, 'CaseStaticWithHarmonicLoad.txt')));
end

function assertClose(actual, expected, message)
difference = max(abs(actual(:) - expected(:)));
scale = max([abs(expected(:)); 1]);
if difference > 1e-11 * scale
    error('MKEF:VerificationFailed', '%s', message);
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
