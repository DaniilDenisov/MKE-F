function filename = exampleCasePath(name)
%EXAMPLECASEPATH Return the absolute path of a supplied example case.

testsDirectory = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(testsDirectory);
filename = fullfile(repositoryRoot, 'examples', 'cases', name);
end
