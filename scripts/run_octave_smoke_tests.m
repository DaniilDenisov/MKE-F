function run_octave_smoke_tests()
%RUN_OCTAVE_SMOKE_TESTS Дымовые тесты GNU Octave без графического интерфейса.

scriptsDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptsDir);
previousDir = pwd;
cleanup = onCleanup(@() cd(previousDir));
cd(rootDir);
addpath(rootDir);
setup();
caseDirectory = fullfile(rootDir, 'examples', 'cases');

options = struct('verbose', false, 'plotting', false);

% Примеры из репозитория - интеграционные тесты, а не решения. Охватывают прямые и наклонные балки,
% непрямолинейную раму, ферму, статические, импульсные и гармонические нагрузки.
% В smoke-наборе примеры ANSYS и Mario Paz проверяют чтение и сборку. Матрицы,
% частоты и формы Mario Paz отдельно сверяются с книгой в проверочном наборе.
caseFiles = {
    fullfile(caseDirectory, 'ANSYSBeamStatic01.txt')
    fullfile(caseDirectory, 'Case1ElementBeam.txt')
    fullfile(caseDirectory, 'CaseATransSite.txt')
    fullfile(caseDirectory, 'CaseBeam.txt')
    fullfile(caseDirectory, 'CaseBeamDyn.txt')
    fullfile(caseDirectory, 'CaseBeamFreq.txt')
    fullfile(caseDirectory, 'CaseFig11.7p363 MarioPaz.txt')
    fullfile(caseDirectory, 'CasePreprocessorStatic.txt')
    fullfile(caseDirectory, 'CasePreprocessorModal.txt')
    fullfile(caseDirectory, 'CasePreprocessorTransient.txt')
};

fprintf('Runtime: %s\n', runtimeName());
fprintf('Loading %d input cases...\n', numel(caseFiles));

for i = 1:numel(caseFiles)
    % Конструктор полностью читает входной файл, создаёт сетку и структуры
    % элементов, нумерует степени свободы и собирает разреженные матрицы K и M.
    problem = StructFEProblem(caseFiles{i}, options);
    expectedDOFs = problem.mesh.numberOfNodes * problem.mesh.dofPerNode;
    assert(isequal(size(problem.K), [expectedDOFs expectedDOFs]));
    assert(isequal(size(problem.M), [expectedDOFs expectedDOFs]));
    assert(issparse(problem.K));
    assert(issparse(problem.M));
    assert(all(isfinite(problem.K(:))));
    assert(all(isfinite(problem.M(:))));
    fprintf('  OK: %s\n', caseFiles{i});
end

% Проверка всех видов анализа с умышленно короткими временными историями.
% Это интеграционные проверки отсутствия сбоев; аналитические значения
% проверяются отдельно в наборе проверочных тестов.
problem = StructFEProblem(fullfile(caseDirectory, ...
    'Case1ElementBeam.txt'), options);
problem.RunStatic();

problem = StructFEProblem(fullfile(caseDirectory, 'CaseBeam.txt'), options);
problem.RunModal();

problem = StructFEProblem(fullfile(caseDirectory, 'CaseBeamDyn.txt'), options);
problem.RunTransient(1e-4, 3e-4, 3, 2);

problem = StructFEProblem(fullfile(caseDirectory, 'CaseBeamFreq.txt'), options);
problem.RunTransient(1e-4, 3e-4, 3, 2);

problem = StructFEProblem(fullfile(caseDirectory, 'CaseATransSite.txt'), options);
problem.RunTransient(1e-4, 3e-4, 5, 2);

problem = StructFEProblem(fullfile(caseDirectory, ...
    'CasePreprocessorStatic.txt'), options);
problem.RunSelected();

problem = StructFEProblem(fullfile(caseDirectory, ...
    'CasePreprocessorModal.txt'), options);
problem.RunSelected();

problem = StructFEProblem(fullfile(caseDirectory, ...
    'CasePreprocessorTransient.txt'), options);
problem.RunSelected();

fprintf('All smoke tests passed.\n');

% Объект очистки должен существовать до выхода из функции.
clear cleanup;
end

function name = runtimeName()
if ~exist('OCTAVE_VERSION', 'builtin')
    error('MKEF:OctaveRequired', ...
        'GNU Octave is the only supported runtime for this project.');
end
name = ['GNU Octave ' OCTAVE_VERSION];
end
