function run_octave_tests()
%RUN_OCTAVE_TESTS Запуск полного набора тестов GNU Octave.
% Это одна точка входа для локального запуска и GitHub Actions.

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
addpath(fullfile(rootDir, 'tests'));

if ~exist('OCTAVE_VERSION', 'builtin')
    error('MKEF:OctaveRequired', ...
        'GNU Octave is the only supported runtime for this project.');
end

% Смоук-тесты проверяют чтение данных, сборку и все виды анализа.
% Затем проверочные тесты сравнивают выбранные результаты с независимыми формулами.
run_octave_smoke_tests();
run_verification_tests();
end
