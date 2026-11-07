function prepare_release_fixtures()
%PREPARE_RELEASE_FIXTURES Fresh inputs for browser suites outside the core tests.
setup();
addpath(fullfile(pwd, 'tests'));
test_section_inspector();
export_mpc_examples();
problem = StructFEProblem(fullfile('examples', 'cases', 'CaseUniformFrame.txt'), ...
    struct('verbose', false, 'plotting', false));
exportPostprocessorData(problem.GetAnalysisModel(), problem.RunStatic(), ...
    fullfile('output', 'uniform-result.json'));
end
