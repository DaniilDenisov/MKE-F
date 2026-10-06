function [localVectors, intensities, resultant] = getElementLoadData(model)
%GETELEMENTLOADDATA Integrate full-span uniform/linear frame loads exactly.
% intensities contains summed local [qx1 qy1 qx2 qy2] for each element.
count = numel(model.elementData);
localVectors = cell(count, 1);
intensities = zeros(count, 4);
resultant = zeros(3, 1);
for i = 1:count
    localVectors{i} = zeros(numel(model.elementData(i).dofs), 1);
end
if ~isfield(model, 'elementLoads') || isempty(model.elementLoads), return; end
loads = model.elementLoads;
baseFields = {'type', 'elementId', 'coordinateSystem'};
if ~isstruct(loads) || ~all(isfield(loads, baseFields))
    error('MKEF:InvalidElementLoad', 'Element loads require type, elementId and coordinateSystem.');
end
for i = 1:numel(loads)
    load = loads(i);
    validateScalars(load, baseFields);
    if load.type == 20
        fields = {'qx', 'qy'};
    elseif load.type == 21
        fields = {'qx1', 'qy1', 'qx2', 'qy2'};
    else
        error('MKEF:InvalidElementLoad', 'Unsupported element load type.');
    end
    validateScalars(load, fields);
    if load.type == 20
        ends = [load.qx load.qx; load.qy load.qy];
    else
        ends = [load.qx1 load.qx2; load.qy1 load.qy2];
    end
    id = load.elementId;
    if id ~= fix(id) || id < 1 || id > count || ...
            ~ismember(load.coordinateSystem, [1 2]) || all(ends(:) == 0)
        error('MKEF:InvalidElementLoad', 'Invalid element load %d.', i);
    end
    element = model.elementData(id);
    if element.type ~= 113
        error('MKEF:UnsupportedElementLoad', 'Distributed loads require frame element 113.');
    end
    rotation = element.transformation(1:2, 1:2);
    if load.coordinateSystem == 2
        globalEnds = ends; ends = rotation * ends;
    else
        globalEnds = rotation.' * ends;
    end
    L = element.length;
    a = ends(:, 1); b = ends(:, 2);
    localVectors{id} = localVectors{id} + ...
        [L*(2*a(1)+b(1))/6; L*(7*a(2)+3*b(2))/20; L^2*(3*a(2)+2*b(2))/60; ...
         L*(a(1)+2*b(1))/6; L*(3*a(2)+7*b(2))/20; -L^2*(2*a(2)+3*b(2))/60];
    intensities(id, :) = intensities(id, :) + ends(:).';
    % Integrate r(x) cross q(x) independently of the FE load vector.
    origin = element.nodeCoordinates(1, 1:2).';
    delta = element.nodeCoordinates(2, 1:2).' - origin;
    force = L * (globalEnds(:, 1) + globalEnds(:, 2))/2;
    firstMoment = L * (globalEnds(:, 1) + 2*globalEnds(:, 2))/6;
    moment = origin(1)*force(2)-origin(2)*force(1) + ...
        delta(1)*firstMoment(2)-delta(2)*firstMoment(1);
    resultant = resultant + [force; moment];
end
end

function validateScalars(load, fields)
for j = 1:numel(fields)
    if ~isfield(load, fields{j})
        error('MKEF:InvalidElementLoad', 'Missing component %s.', fields{j});
    end
    value = load.(fields{j});
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ~isfinite(value)
        error('MKEF:InvalidElementLoad', 'Element load components and identifiers must be finite scalars.');
    end
end
end
