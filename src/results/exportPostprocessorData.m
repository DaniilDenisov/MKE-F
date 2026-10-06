function data = exportPostprocessorData(model, result, filename, options)
%EXPORTPOSTPROCESSORDATA Validate and atomically write postprocessor JSON.

if nargin < 4
    options = struct();
end
if isstring(filename) && isscalar(filename)
    filename = char(filename);
end
if ~ischar(filename) || isempty(filename) || size(filename, 1) ~= 1
    error('MKEF:PostprocessorInvalidFilename', ...
        'filename must be a nonempty scalar string.');
end

data = createPostprocessorData(model, result, options);
prettyPrint = false;
if isfield(options, 'prettyPrint')
    prettyPrint = logical(options.prettyPrint);
end
try
    % A scalar numeric value and a one-element vector have the same shape in
    % MATLAB/Octave, as do a scalar struct and a one-element struct array.
    % jsonencode therefore emits JSON scalars/objects for schema fields that
    % must remain arrays.  Convert only the JSON-facing copy to cells so the
    % array and matrix dimensions are explicit even when a dimension is one.
    jsonData = prepareJsonData(data);
    json = jsonencode(jsonData, 'PrettyPrint', prettyPrint);
catch exception
    error('MKEF:PostprocessorEncodingFailed', ...
        'Could not encode postprocessor JSON: %s', exception.message);
end

targetDirectory = fileparts(filename);
if isempty(targetDirectory)
    targetDirectory = pwd;
end
if exist(targetDirectory, 'dir') ~= 7
    error('MKEF:PostprocessorWriteFailed', ...
        'Target directory does not exist: %s', targetDirectory);
end

temporaryPath = tempname(targetDirectory);
temporaryExists = false;
cleanup = onCleanup(@removeTemporary);
fid = fopen(temporaryPath, 'wb');
if fid < 0
    error('MKEF:PostprocessorWriteFailed', ...
        'Could not create temporary output in %s.', targetDirectory);
end
temporaryExists = true;

try
    bytes = unicode2native(json, 'UTF-8');
    written = fwrite(fid, bytes, 'uint8');
    closeStatus = fclose(fid);
    fid = -1;
    if written ~= numel(bytes) || closeStatus ~= 0
        error('MKEF:PostprocessorWriteFailed', ...
            'Could not completely write %s.', filename);
    end
catch exception
    if fid >= 0
        fclose(fid);
        fid = -1;
    end
    if strcmp(exception.identifier, 'MKEF:PostprocessorWriteFailed')
        rethrow(exception);
    end
    error('MKEF:PostprocessorWriteFailed', ...
        'Could not write %s: %s', filename, exception.message);
end

backupPath = '';
originalMoved = false;
if exist(filename, 'file') == 2
    backupPath = tempname(targetDirectory);
    [status, message] = rename(filename, backupPath);
    if status ~= 0
        error('MKEF:PostprocessorWriteFailed', ...
            'Could not prepare %s for replacement: %s', filename, message);
    end
    originalMoved = true;
end

[status, message] = rename(temporaryPath, filename);
if status ~= 0
    if originalMoved
        rename(backupPath, filename);
    end
    error('MKEF:PostprocessorWriteFailed', ...
        'Could not replace %s: %s', filename, message);
end
temporaryExists = false;
if originalMoved && exist(backupPath, 'file') == 2
    delete(backupPath);
end
clear cleanup;

    function removeTemporary()
        if temporaryExists && exist(temporaryPath, 'file') == 2
            delete(temporaryPath);
        end
    end
end

function output = prepareJsonData(data)
output = data;

output.model.nodes = recordArray(data.model.nodes);
output.model.elements = recordArray(data.model.elements);
for i = 1:numel(output.model.elements)
    output.model.elements{i}.nodeIds = ...
        numericArray(output.model.elements{i}.nodeIds);
end
output.model.dofMap = numericMatrix(data.model.dofMap);
output.model.supports = recordArray(data.model.supports);
if isfield(data.model, 'elementLoads')
    output.model.nodalLoads = recordArray(data.model.nodalLoads);
    output.model.elementLoads = recordArray(data.model.elementLoads);
    for i = 1:numel(output.model.elementLoads)
        load = output.model.elementLoads{i};
        if load.type == 20, fields = {'qx1','qy1','qx2','qy2'};
        else, fields = {'qx','qy'}; end
        fields = fields(isfield(load, fields));
        output.model.elementLoads{i} = rmfield(load, fields);
    end
end

switch data.analysis.type
    case 'static'
        vectorFields = {'displacements', 'loadVector', 'reactions', ...
            'equilibriumResidual'};
        for i = 1:numel(vectorFields)
            name = vectorFields{i};
            output.analysis.(name) = numericArray(data.analysis.(name));
        end
        output.analysis.elementResults = ...
            recordArray(data.analysis.elementResults);
        for i = 1:numel(output.analysis.elementResults)
            output.analysis.elementResults{i}.localEndForces = numericArray( ...
                output.analysis.elementResults{i}.localEndForces);
            if isfield(data.model, 'elementLoads')
                output.analysis.elementResults{i}.equivalentLocalLoadVector = ...
                    numericArray(output.analysis.elementResults{i}.equivalentLocalLoadVector);
            end
        end

    case 'modal'
        output.analysis.frequenciesHz = ...
            numericArray(data.analysis.frequenciesHz);
        output.analysis.angularFrequenciesRadPerSec = ...
            numericArray(data.analysis.angularFrequenciesRadPerSec);
        output.analysis.modeShapes = numericMatrix(data.analysis.modeShapes);

    case 'transient'
        output.analysis.time = numericArray(data.analysis.time);
        output.analysis.globalDOFIds = ...
            numericArray(data.analysis.globalDOFIds);
        historyFields = {'displacements', 'velocities', 'accelerations', ...
            'loadHistory', 'reactions'};
        for i = 1:numel(historyFields)
            name = historyFields{i};
            if isfield(data.analysis, name)
                output.analysis.(name) = ...
                    numericMatrix(data.analysis.(name));
            end
        end
        if isfield(data.analysis, 'spectrumFrequencyHz')
            output.analysis.spectrumFrequencyHz = ...
                numericArray(data.analysis.spectrumFrequencyHz);
            output.analysis.displacementAmplitudeSpectrum = numericMatrix( ...
                data.analysis.displacementAmplitudeSpectrum);
        end
end
if isfield(data.model,'mpcs')
    output.model.mpcs = recordArray(data.model.mpcs);
    for i = 1:numel(output.model.mpcs)
        output.model.mpcs{i}.masters = recordArray(data.model.mpcs(i).masters);
    end
    output.model.dependentDOFs = numericArray(data.model.dependentDOFs);
    output.model.independentDOFs = numericArray(data.model.independentDOFs);
    names = {'supportReactions','mpcForces','mpcMultipliers'};
    for i = 1:numel(names)
        name = names{i};
        if isfield(data.analysis,name)
            if strcmp(data.analysis.type,'static')
                output.analysis.(name) = numericArray(data.analysis.(name));
            else
                output.analysis.(name) = numericMatrix(data.analysis.(name));
            end
        end
    end
end
if isfield(data.model, 'dofRegistry')
    output.model.dofRegistry = recordArray(data.model.dofRegistry);
    output.model.releases = recordArray(data.model.releases);
    output.model.warnings = data.model.warnings;
    for i = 1:numel(output.model.elements)
        output.model.elements{i}.globalDOFs = numericArray(data.model.elements(i).globalDOFs);
    end
end

end

function output = recordArray(value)
output = cell(1, numel(value));
for i = 1:numel(value)
    output{i} = value(i);
end
end

function output = numericArray(value)
output = num2cell(double(value(:).'));
end

function output = numericMatrix(value)
output = cell(size(value, 1), 1);
for row = 1:size(value, 1)
    output{row} = numericArray(value(row, :));
end
end
