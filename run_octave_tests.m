function run_octave_tests()
%RUN_OCTAVE_TESTS Run the complete GNU Octave test suite.

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
addpath(fullfile(rootDir, 'tests'));

if ~exist('OCTAVE_VERSION', 'builtin')
    error('MKEF:OctaveRequired', ...
        'GNU Octave is the only supported runtime for this project.');
end

run_octave_smoke_tests();
run_verification_tests();
end
