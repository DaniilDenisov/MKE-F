function run_octave_smoke_tests()
%RUN_OCTAVE_SMOKE_TESTS Дымовые тесты GNU Octave без графического интерфейса.

rootDir = fileparts(mfilename('fullpath'));
previousDir = pwd;
cleanup = onCleanup(@() cd(previousDir));
cd(rootDir);
addpath(rootDir);

options = struct('verbose', false, 'plotting', false);

% Примеры из репозитория - интеграционные тесты, а не решения. Охватывают прямые и наклонные балки,
% непрямолинейную раму, ферму, статические, импульсные и гармонические нагрузки.
% Примеры ANSYS и Mario Paz пока проверяют только успешное чтение и сборку;
% сравнения с опубликованными результатами не выполняются.
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
    % Конструктор полностью читает входной файл, создаёт сетку и объекты
    % элементов, нумерует степени свободы и собирает глобальные матрицы K и M.
    problem = StructFEProblem(caseFiles{i}, options);
    expectedDOFs = problem.mesh.numberOfNodes * problem.mesh.dofPerNode;
    assert(isequal(size(problem.K), [expectedDOFs expectedDOFs]));
    assert(isequal(size(problem.M), [expectedDOFs expectedDOFs]));
    assert(all(isfinite(problem.K(:))));
    assert(all(isfinite(problem.M(:))));
    fprintf('  OK: %s\n', caseFiles{i});
end

% Проверка всех видов анализа с умышленно короткими временными историями.
% Это интеграционные проверки отсутствия сбоев; аналитические значения
% проверяются отдельно в наборе проверочных тестов.
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
