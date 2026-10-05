function results = run_reference_examples(verbose)
%RUN_REFERENCE_EXAMPLES Run every executable reference example in order.
%   RESULTS = RUN_REFERENCE_EXAMPLES() runs all examples with diagnostic
%   output and returns each result under a field matching its function name.
%   Pass false for a quiet assertion-only run. The runner locates the
%   repository from its own file, so it is independent of the current
%   working directory once reference/examples is on the Octave path.

if nargin < 1
    verbose = true;
end
validateVerbose(verbose);

exampleDirectory = fileparts(mfilename('fullpath'));
referenceDirectory = fileparts(exampleDirectory);
repositoryRoot = fileparts(referenceDirectory);
addpath(exampleDirectory);
addpath(repositoryRoot);
setup();

exampleFunctions = { ...
    @example_01_element_matrices, ...
    @example_02_sparse_assembly, ...
    @example_03_static_bar, ...
    @example_04_load_histories, ...
    @example_05_modal_frame, ...
    @example_06_newmark_sdof, ...
    @example_07_result_recovery, ...
    @example_08_pinned_portal_uniform};

results = struct();
if verbose
    fprintf('\nMKE-F reference examples\n');
end

for exampleNumber = 1:numel(exampleFunctions)
    exampleFunction = exampleFunctions{exampleNumber};
    exampleName = func2str(exampleFunction);
    if verbose
        fprintf('\n[%d/%d] %s\n', ...
            exampleNumber, numel(exampleFunctions), exampleName);
    end
    results.(exampleName) = exampleFunction(verbose);
end

if verbose
    fprintf('\nAll %d reference examples PASS.\n', numel(exampleFunctions));
end
end

function validateVerbose(verbose)
if ~isscalar(verbose) || ...
        ~(islogical(verbose) || (isnumeric(verbose) && isfinite(verbose))) || ...
        ~ismember(double(verbose), [0 1])
    error('MKEF:InvalidExampleOption', ...
        'verbose must be a scalar logical value.');
end
end
