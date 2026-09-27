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
    json = jsonencode(data, 'PrettyPrint', prettyPrint);
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
