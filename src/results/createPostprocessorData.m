function data = createPostprocessorData(model, result, options)
%CREATEPOSTPROCESSORDATA Convert an analysis model/result to schema version 1 or 2.
% The returned struct contains only postprocessing data; assembled and
% element matrices are deliberately excluded.

if nargin < 3
    options = struct();
end

analysisType = validateAnalysisType(result);
[settings, ~] = validateOptions(options, analysisType);
modelInfo = validateModel(model);

data = struct();
data.format = 'mkef-postprocessor';
data.version = 1;
data.metadata = createMetadata(settings, analysisType);
data.model = createModelData(model, modelInfo);
if isfield(model, 'elementLoads') && ~isempty(model.elementLoads)
    if ~strcmp(analysisType, 'static')
        rejectElementLoads(model);
    end
    getElementLoadData(model);
    data.version = 2;
    data.model.elementLoads = model.elementLoads;
    records = model.forceBoundaryConditions;
    validateNumeric(records, 'model.forceBoundaryConditions');
    if size(records, 2) < 5 || any(records(:, 1) ~= 10) || ...
            any(records(:, 2) ~= fix(records(:, 2))) || ...
            any(records(:, 2) < 1 | records(:, 2) > model.numberOfNodes)
        invalidModel('Version 2 requires finite static nodal loads on existing nodes.');
    end
    nodalLoads = repmat(struct('type', 10, 'nodeId', 0, 'fx', 0, 'fy', 0, 'mz', 0), size(records, 1), 1);
    for i = 1:size(records, 1)
        nodalLoads(i) = struct('type', records(i, 1), 'nodeId', records(i, 2), ...
            'fx', records(i, 3), 'fy', records(i, 4), 'mz', records(i, 5));
    end
    data.model.nodalLoads = nodalLoads;
end

switch analysisType
    case 'static'
        data.analysis = createStaticData(result, modelInfo);
    case 'modal'
        data.analysis = createModalData(result, modelInfo);
    case 'transient'
        data.analysis = createTransientData(result, modelInfo, settings);
    otherwise
        error('MKEF:PostprocessorUnsupportedAnalysis', ...
            'Unsupported analysis type %s.', analysisType);
end

if data.version == 2
    for i = 1:numel(data.analysis.elementResults)
        source = result.elementResults(i);
        requireFields(source, {'equivalentLocalLoadVector'}, 'element result');
        validateVector(source.equivalentLocalLoadVector, 6, 'equivalentLocalLoadVector');
        data.analysis.elementResults(i).equivalentLocalLoadVector = ...
            double(source.equivalentLocalLoadVector(:).');
    end
end
end

function analysisType = validateAnalysisType(result)
if ~isstruct(result) || ~isscalar(result) || ...
        ~isfield(result, 'analysisType') || ...
        ~(ischar(result.analysisType) || ...
          (isstring(result.analysisType) && isscalar(result.analysisType)))
    error('MKEF:PostprocessorInvalidResult', ...
        'The result must contain a scalar analysisType string.');
end
analysisType = char(result.analysisType);
if ~ismember(analysisType, {'static', 'modal', 'transient'})
    error('MKEF:PostprocessorUnsupportedAnalysis', ...
        'Unsupported analysis type %s.', analysisType);
end
end

function [settings, supplied] = validateOptions(options, analysisType)
if ~isstruct(options) || ~isscalar(options)
    error('MKEF:PostprocessorInvalidOptions', ...
        'Exporter options must be a scalar struct.');
end

allowed = {'title', 'lengthUnit', 'forceUnit', 'momentUnit', ...
    'stressUnit', 'timeUnit', 'transientFields', 'timeStride', ...
    'selectedGlobalDOFs', 'prettyPrint'};
supplied = fieldnames(options).';
unknown = setdiff(supplied, allowed);
if ~isempty(unknown)
    error('MKEF:PostprocessorInvalidOptions', ...
        'Unknown exporter option %s.', unknown{1});
end

transientOnly = {'transientFields', 'timeStride', 'selectedGlobalDOFs'};
if ~strcmp(analysisType, 'transient')
    invalid = intersect(supplied, transientOnly);
    if ~isempty(invalid)
        error('MKEF:PostprocessorInvalidOptions', ...
            'Option %s is valid only for transient results.', invalid{1});
    end
end

settings = struct();
settings.title = [upper(analysisType(1)) analysisType(2:end) ' analysis'];
settings.lengthUnit = '';
settings.forceUnit = '';
settings.momentUnit = '';
settings.stressUnit = '';
settings.timeUnit = '';
settings.transientFields = {'displacements', 'velocities', ...
    'accelerations', 'loadHistory', 'reactions', 'spectrum'};
settings.timeStride = 1;
settings.selectedGlobalDOFs = [];
settings.prettyPrint = false;

textFields = {'title', 'lengthUnit', 'forceUnit', 'momentUnit', ...
    'stressUnit', 'timeUnit'};
for i = 1:numel(textFields)
    name = textFields{i};
    if isfield(options, name)
        settings.(name) = scalarText(options.(name), name);
    end
end

if isfield(options, 'transientFields')
    fields = options.transientFields;
    if ischar(fields) || (isstring(fields) && isscalar(fields))
        fields = {char(fields)};
    elseif isstring(fields)
        fields = cellstr(fields(:));
    end
    if ~iscell(fields)
        error('MKEF:PostprocessorInvalidOptions', ...
            'transientFields must be a cell array of field names.');
    end
    fields = fields(:).';
    validFields = {'displacements', 'velocities', 'accelerations', ...
        'loadHistory', 'reactions', 'spectrum'};
    for i = 1:numel(fields)
        fields{i} = scalarText(fields{i}, 'transientFields');
        if ~ismember(fields{i}, validFields)
            error('MKEF:PostprocessorInvalidOptions', ...
                'Unsupported transient field %s.', fields{i});
        end
    end
    if numel(unique(fields)) ~= numel(fields)
        error('MKEF:PostprocessorInvalidOptions', ...
            'transientFields must not contain duplicates.');
    end
    settings.transientFields = fields;
end

if isfield(options, 'timeStride')
    value = options.timeStride;
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
            ~isfinite(value) || value < 1 || value ~= fix(value)
        error('MKEF:PostprocessorInvalidOptions', ...
            'timeStride must be a positive integer.');
    end
    settings.timeStride = double(value);
end

if isfield(options, 'selectedGlobalDOFs')
    value = options.selectedGlobalDOFs;
    if ~isnumeric(value) || ~isvector(value) || ~isreal(value) || ...
            isempty(value) || any(~isfinite(value(:))) || ...
            any(value(:) < 1) || any(value(:) ~= fix(value(:))) || ...
            numel(unique(value(:))) ~= numel(value)
        error('MKEF:PostprocessorInvalidOptions', ...
            'selectedGlobalDOFs must contain unique positive integers.');
    end
    settings.selectedGlobalDOFs = double(value(:).');
end

if isfield(options, 'prettyPrint')
    value = options.prettyPrint;
    if ~(islogical(value) && isscalar(value)) && ...
            ~(isnumeric(value) && isscalar(value) && ismember(value, [0 1]))
        error('MKEF:PostprocessorInvalidOptions', ...
            'prettyPrint must be a scalar logical value.');
    end
    settings.prettyPrint = logical(value);
end
end

function value = scalarText(value, fieldName)
if isstring(value) && isscalar(value)
    value = char(value);
end
if ~ischar(value) || size(value, 1) > 1
    error('MKEF:PostprocessorInvalidOptions', ...
        '%s must be a scalar string.', fieldName);
end
end

function info = validateModel(model)
required = {'numberOfNodes', 'dofPerNode', 'numberOfDOFs', 'nodeCoordinates', ...
    'dofMap', 'elementData', 'fixedBoundaryConditions'};
requireFields(model, required, 'model');

validatePositiveInteger(model.numberOfNodes, 'model.numberOfNodes');
validatePositiveInteger(model.numberOfDOFs, 'model.numberOfDOFs');
if ~isnumeric(model.dofPerNode) || ~isscalar(model.dofPerNode) || ...
        ~ismember(model.dofPerNode, [2 3])
    invalidModel('model.dofPerNode must be 2 or 3.');
end
if model.numberOfDOFs ~= model.numberOfNodes * model.dofPerNode
    invalidModel('The model node and DOF counts are inconsistent.');
end

validateNumeric(model.nodeCoordinates, 'model.nodeCoordinates');
if size(model.nodeCoordinates, 1) ~= model.numberOfNodes || ...
        size(model.nodeCoordinates, 2) < 2
    invalidModel('model.nodeCoordinates has inconsistent dimensions.');
end

validateNumeric(model.dofMap, 'model.dofMap');
if ~isequal(size(model.dofMap), ...
        [model.numberOfNodes model.dofPerNode]) || ...
        any(model.dofMap(:) ~= fix(model.dofMap(:))) || ...
        ~isequal(sort(model.dofMap(:)).', 1:model.numberOfDOFs)
    invalidModel('model.dofMap must contain every global DOF exactly once.');
end

if ~isstruct(model.elementData) || isempty(model.elementData)
    invalidModel('model.elementData must be a nonempty struct array.');
end
for i = 1:numel(model.elementData)
    element = model.elementData(i);
    requireFields(element, {'type', 'nodeNumbers', 'properties'}, ...
        sprintf('model.elementData(%d)', i));
    if ~isnumeric(element.type) || ~isscalar(element.type) || ...
            ~ismember(element.type, [112 113])
        invalidModel(sprintf('Element %d has an unsupported type.', i));
    end
    expectedDOFs = 2 + (element.type == 113);
    if expectedDOFs ~= model.dofPerNode
        invalidModel(sprintf('Element %d does not match dofPerNode.', i));
    end
    validateNumeric(element.nodeNumbers, ...
        sprintf('model.elementData(%d).nodeNumbers', i));
    if numel(element.nodeNumbers) ~= 2 || ...
            any(element.nodeNumbers(:) ~= fix(element.nodeNumbers(:))) || ...
            any(element.nodeNumbers(:) < 1) || ...
            any(element.nodeNumbers(:) > model.numberOfNodes) || ...
            element.nodeNumbers(1) == element.nodeNumbers(2)
        invalidModel(sprintf('Element %d has invalid node references.', i));
    end
    validateNumeric(element.properties, ...
        sprintf('model.elementData(%d).properties', i));
    propertyCount = 3 + (element.type == 113);
    if numel(element.properties) < propertyCount || ...
            any(element.properties(1:propertyCount) <= 0)
        invalidModel(sprintf('Element %d has invalid properties.', i));
    end
end

validateNumeric(model.fixedBoundaryConditions, ...
    'model.fixedBoundaryConditions');
if size(model.fixedBoundaryConditions, 2) ~= 5
    invalidModel('model.fixedBoundaryConditions must have five columns.');
end
for i = 1:size(model.fixedBoundaryConditions, 1)
    support = model.fixedBoundaryConditions(i, :);
    if support(1) ~= fix(support(1)) || ~ismember(support(1), 1:4) || ...
            support(2) ~= fix(support(2)) || support(2) < 1 || ...
            support(2) > model.numberOfNodes || any(support(3:5) ~= 0)
        invalidModel(sprintf('Support %d is invalid.', i));
    end
end

% Reuse the production constraint rules, including duplicate-DOF rejection.
try
    partitionDOFs(model);
catch exception
    error('MKEF:PostprocessorInvalidModel', ...
        'Invalid supports: %s', exception.message);
end

info = struct('numberOfNodes', double(model.numberOfNodes), ...
    'dofPerNode', double(model.dofPerNode), ...
    'numberOfDOFs', double(model.numberOfDOFs), ...
    'numberOfElements', double(numel(model.elementData)));
end

function metadata = createMetadata(settings, analysisType)
metadata = struct();
metadata.title = settings.title;
metadata.createdUtc = strftime('%Y-%m-%dT%H:%M:%SZ', gmtime(time()));
metadata.generator = 'MKE-F';
metadata.octaveVersion = OCTAVE_VERSION;
metadata.analysisType = analysisType;
metadata.units = struct('length', settings.lengthUnit, ...
    'force', settings.forceUnit, 'moment', settings.momentUnit, ...
    'stress', settings.stressUnit, 'time', settings.timeUnit);
end

function output = createModelData(model, info)
output = struct();
output.dimension = 2;
output.dofPerNode = info.dofPerNode;
if info.dofPerNode == 2
    output.dofLabels = {'ux', 'uy'};
else
    output.dofLabels = {'ux', 'uy', 'rz'};
end

nodes = repmat(struct('id', 0, 'x', 0, 'y', 0), info.numberOfNodes, 1);
for i = 1:info.numberOfNodes
    nodes(i).id = i;
    nodes(i).x = model.nodeCoordinates(i, 1);
    nodes(i).y = model.nodeCoordinates(i, 2);
end
output.nodes = nodes;

emptyElement = struct('id', 0, 'type', 0, 'nodeIds', [], ...
    'properties', struct());
elements = repmat(emptyElement, info.numberOfElements, 1);
for i = 1:info.numberOfElements
    source = model.elementData(i);
    properties = struct('area', source.properties(1), ...
        'youngsModulus', source.properties(2), ...
        'density', source.properties(3));
    if source.type == 113
        properties.momentOfInertia = source.properties(4);
    end
    elements(i).id = i;
    elements(i).type = source.type;
    elements(i).nodeIds = double(source.nodeNumbers(:).');
    elements(i).properties = properties;
end
output.elements = elements;
output.dofMap = double(model.dofMap);

supports = repmat(struct('type', 0, 'nodeId', 0), ...
    size(model.fixedBoundaryConditions, 1), 1);
for i = 1:numel(supports)
    supports(i).type = model.fixedBoundaryConditions(i, 1);
    supports(i).nodeId = model.fixedBoundaryConditions(i, 2);
end
output.supports = supports;
end

function output = createStaticData(result, info)
required = {'displacements', 'loadVector', 'reactions', ...
    'equilibriumResidual', 'elementResults'};
requireFields(result, required, 'result');
validateVector(result.displacements, info.numberOfDOFs, ...
    'result.displacements');
validateVector(result.loadVector, info.numberOfDOFs, 'result.loadVector');
validateVector(result.reactions, info.numberOfDOFs, 'result.reactions');
validateVector(result.equilibriumResidual, 3, 'result.equilibriumResidual');
if ~isstruct(result.elementResults) || ...
        numel(result.elementResults) ~= info.numberOfElements
    invalidResult('result.elementResults has inconsistent dimensions.');
end

items = repmat(struct('elementId', 0, 'type', 0, ...
    'localEndForces', [], 'axialStrain', 0, 'axialStress', 0, ...
    'axialForce', 0), info.numberOfElements, 1);
for i = 1:info.numberOfElements
    source = result.elementResults(i);
    requireFields(source, {'type', 'localEndForces', 'axialStrain', ...
        'axialStress', 'axialForce'}, sprintf('result.elementResults(%d)', i));
    expectedSize = 4 + 2 * (source.type == 113);
    if ~ismember(source.type, [112 113])
        invalidResult(sprintf('Element result %d has an unsupported type.', i));
    end
    validateVector(source.localEndForces, expectedSize, ...
        sprintf('result.elementResults(%d).localEndForces', i));
    validateScalar(source.axialStrain, sprintf('elementResults(%d).axialStrain', i));
    validateScalar(source.axialStress, sprintf('elementResults(%d).axialStress', i));
    validateScalar(source.axialForce, sprintf('elementResults(%d).axialForce', i));
    items(i).elementId = i;
    items(i).type = source.type;
    items(i).localEndForces = double(source.localEndForces(:).');
    items(i).axialStrain = source.axialStrain;
    items(i).axialStress = source.axialStress;
    items(i).axialForce = source.axialForce;
end

output = struct('type', 'static', ...
    'displacements', double(result.displacements(:).'), ...
    'loadVector', double(result.loadVector(:).'), ...
    'reactions', double(result.reactions(:).'), ...
    'equilibriumResidual', double(result.equilibriumResidual(:).'), ...
    'elementResults', items);
end

function output = createModalData(result, info)
required = {'frequenciesHz', 'angularFrequenciesRadPerSec', 'modeShapes'};
requireFields(result, required, 'result');
validateNumeric(result.frequenciesHz, 'result.frequenciesHz');
validateNumeric(result.angularFrequenciesRadPerSec, ...
    'result.angularFrequenciesRadPerSec');
validateNumeric(result.modeShapes, 'result.modeShapes');
modeCount = numel(result.frequenciesHz);
if modeCount < 1 || numel(result.angularFrequenciesRadPerSec) ~= modeCount || ...
        ~isequal(size(result.modeShapes), [info.numberOfDOFs modeCount])
    invalidResult('Modal result dimensions are inconsistent.');
end
output = struct('type', 'modal', ...
    'frequenciesHz', double(result.frequenciesHz(:).'), ...
    'angularFrequenciesRadPerSec', ...
        double(result.angularFrequenciesRadPerSec(:).'), ...
    'modeShapes', double(result.modeShapes));
end

function output = createTransientData(result, info, settings)
requireFields(result, {'time'}, 'result');
validateNumeric(result.time, 'result.time');
if ~isvector(result.time) || isempty(result.time)
    invalidResult('result.time must be a nonempty vector.');
end
time = double(result.time(:).');
if any(diff(time) <= 0)
    invalidResult('result.time must be strictly increasing.');
end

if isempty(settings.selectedGlobalDOFs)
    selectedDOFs = 1:info.numberOfDOFs;
else
    selectedDOFs = settings.selectedGlobalDOFs;
    if any(selectedDOFs > info.numberOfDOFs)
        error('MKEF:PostprocessorInvalidOptions', ...
            'selectedGlobalDOFs exceeds the model DOF count.');
    end
end

indices = 1:settings.timeStride:numel(time);
if indices(end) ~= numel(time)
    indices(end + 1) = numel(time);
end
output = struct();
output.type = 'transient';
output.time = time(indices);
output.globalDOFIds = double(selectedDOFs);

historyFields = {'displacements', 'velocities', 'accelerations', ...
    'loadHistory', 'reactions'};
for i = 1:numel(historyFields)
    name = historyFields{i};
    if ismember(name, settings.transientFields)
        requireFields(result, {name}, 'result');
        validateNumeric(result.(name), ['result.' name]);
        if ~isequal(size(result.(name)), [info.numberOfDOFs numel(time)])
            invalidResult(sprintf('result.%s has inconsistent dimensions.', name));
        end
        output.(name) = double(result.(name)(selectedDOFs, indices));
    end
end

if ismember('spectrum', settings.transientFields)
    requireFields(result, {'spectrumFrequencyHz', ...
        'displacementAmplitudeSpectrum'}, 'result');
    validateNumeric(result.spectrumFrequencyHz, 'result.spectrumFrequencyHz');
    validateNumeric(result.displacementAmplitudeSpectrum, ...
        'result.displacementAmplitudeSpectrum');
    frequency = result.spectrumFrequencyHz(:).';
    if ~isequal(size(result.displacementAmplitudeSpectrum), ...
            [info.numberOfDOFs numel(frequency)])
        invalidResult('The displacement spectrum dimensions are inconsistent.');
    end
    output.spectrumFrequencyHz = double(frequency);
    output.displacementAmplitudeSpectrum = ...
        double(result.displacementAmplitudeSpectrum(selectedDOFs, :));
end

output.sampling = struct('originalSampleCount', numel(time), ...
    'exportedSampleCount', numel(indices), ...
    'timeStride', settings.timeStride);
end

function requireFields(value, fields, context)
if ~isstruct(value) || ~isscalar(value)
    if strcmp(context, 'model')
        invalidModel('The model must be a scalar struct.');
    else
        invalidResult([context ' must be a scalar struct.']);
    end
end
for i = 1:numel(fields)
    if ~isfield(value, fields{i})
        if strcmp(context, 'model')
            invalidModel(['Missing field model.' fields{i} '.']);
        else
            invalidResult(['Missing field ' context '.' fields{i} '.']);
        end
    end
end
end

function validateNumeric(value, path)
if ~isnumeric(value) || ~isreal(value) || issparse(value) || ...
        any(~isfinite(value(:)))
    if strncmp(path, 'model', 5)
        invalidModel([path ' must contain finite, real, dense numeric values.']);
    else
        invalidResult([path ' must contain finite, real, dense numeric values.']);
    end
end
end

function validateVector(value, count, path)
validateNumeric(value, path);
if ~isvector(value) || numel(value) ~= count
    invalidResult([path ' has inconsistent dimensions.']);
end
end

function validateScalar(value, path)
validateNumeric(value, ['result.' path]);
if ~isscalar(value)
    invalidResult(['result.' path ' must be scalar.']);
end
end

function validatePositiveInteger(value, path)
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value < 1 || value ~= fix(value)
    invalidModel([path ' must be a positive integer.']);
end
end

function invalidModel(message)
error('MKEF:PostprocessorInvalidModel', '%s', message);
end

function invalidResult(message)
error('MKEF:PostprocessorInvalidResult', '%s', message);
end
