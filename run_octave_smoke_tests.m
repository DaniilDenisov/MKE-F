function run_octave_smoke_tests()
%RUN_OCTAVE_SMOKE_TESTS Portable headless smoke tests for MATLAB and Octave.

rootDir = fileparts(mfilename('fullpath'));
previousDir = pwd;
cleanup = onCleanup(@() cd(previousDir));
cd(rootDir);
addpath(rootDir);

options = struct('verbose', false, 'plotting', false);
caseFiles = {
    'ANSYSBeamStatic01.txt'
    'Case1ElementBeam.txt'
    'CaseATransSite.txt'
    'CaseBeam.txt'
    'CaseBeamDyn.txt'
    'CaseBeamFreq.txt'
    'CaseFig11.7p363 MarioPaz.txt'
};

fprintf('Runtime: %s\n', runtimeName());
fprintf('Loading %d input cases...\n', numel(caseFiles));

for i = 1:numel(caseFiles)
    problem = StructFEProblem(caseFiles{i}, options);
    expectedDOFs = problem.mesh.numberOfNodes * problem.mesh.dofPerNode;
    assert(isequal(size(problem.K), [expectedDOFs expectedDOFs]));
    assert(isequal(size(problem.M), [expectedDOFs expectedDOFs]));
    assert(all(isfinite(problem.K(:))));
    assert(all(isfinite(problem.M(:))));
    fprintf('  OK: %s\n', caseFiles{i});
end

% Exercise every analysis path with deliberately short transient histories.
problem = StructFEProblem('Case1ElementBeam.txt', options);
problem.RunStatic();

problem = StructFEProblem('CaseBeam.txt', options);
problem.RunModal();

problem = StructFEProblem('CaseBeamDyn.txt', options);
problem.RunTransient(1e-4, 3e-4, 3, 2);

problem = StructFEProblem('CaseBeamFreq.txt', options);
problem.RunTransient(1e-4, 3e-4, 3, 2);

problem = StructFEProblem('CaseATransSite.txt', options);
problem.RunTransient(1e-4, 3e-4, 5, 2);

fprintf('All smoke tests passed.\n');

% Keep the cleanup object alive until the function exits.
clear cleanup;
end

function name = runtimeName()
if exist('OCTAVE_VERSION', 'builtin')
    name = ['GNU Octave ' OCTAVE_VERSION];
else
    name = ['MATLAB ' version];
end
end
