function setup()
%SETUP Add the MKE-F source and launcher directories to the Octave path.

repositoryRoot = fileparts(mfilename('fullpath'));
sourceDirectories = {
    fullfile(repositoryRoot, 'scripts')
    fullfile(repositoryRoot, 'src', 'model')
    fullfile(repositoryRoot, 'src', 'analysis', 'common')
    fullfile(repositoryRoot, 'src', 'analysis', 'static')
    fullfile(repositoryRoot, 'src', 'analysis', 'modal')
    fullfile(repositoryRoot, 'src', 'analysis', 'transient')
    fullfile(repositoryRoot, 'src', 'results')
};

for i = 1:numel(sourceDirectories)
    addpath(sourceDirectories{i});
end
end
