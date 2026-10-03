%RUN_SOLVER_JOB Fixed non-interactive entry point for the container runner.

arguments = argv();
if numel(arguments) ~= 3
    fprintf(2, ['MKEF_ERROR\tMKEF:InvalidRunnerArguments\t' ...
        'Expected input path, output path, and result title.\n']);
    exit(2);
end

inputPath = arguments{1};
outputPath = arguments{2};
resultTitle = arguments{3};
repositoryRoot = fileparts(fileparts(mfilename('fullpath')));

try
    addpath(repositoryRoot);
    setup();

    problemOptions = struct('verbose', false, 'plotting', false);
    problem = StructFEProblem(inputPath, problemOptions);
    model = problem.GetAnalysisModel();
    result = problem.RunSelected();

    exportOptions = struct('title', resultTitle, ...
        'lengthUnit', '', 'forceUnit', '', 'momentUnit', '', ...
        'stressUnit', '', 'timeUnit', '', 'prettyPrint', false);
    if strcmp(result.analysisType, 'transient')
        exportOptions.transientFields = {'displacements', 'velocities', ...
            'accelerations', 'loadHistory', 'reactions', 'spectrum'};
        exportOptions.timeStride = 1;
    end
    exportPostprocessorData(model, result, outputPath, exportOptions);
catch exception
    identifier = exception.identifier;
    if isempty(identifier)
        identifier = 'MKEF:SolverFailure';
    end
    message = strrep(exception.message, sprintf('\n'), ' ');
    message = strrep(message, sprintf('\r'), ' ');
    fprintf(2, 'MKEF_ERROR\t%s\t%s\n', identifier, message);
    for frameIndex = 1:numel(exception.stack)
        frame = exception.stack(frameIndex);
        fprintf(2, '  at %s (%s:%d)\n', ...
            frame.name, frame.file, frame.line);
    end
    exit(1);
end

exit(0);
