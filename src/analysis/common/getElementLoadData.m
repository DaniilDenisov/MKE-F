function [localVectors, intensities, resultant] = getElementLoadData(model)
%GETELEMENTLOADDATA Validate uniform loads and integrate in the original basis.
count = numel(model.elementData);
localVectors = cell(count, 1);
intensities = zeros(count, 2);
resultant = zeros(3, 1);
for i = 1:count
    localVectors{i} = zeros(numel(model.elementData(i).dofs), 1);
end
if ~isfield(model, 'elementLoads') || isempty(model.elementLoads)
    return;
end
loads = model.elementLoads;
fields = {'type', 'elementId', 'coordinateSystem', 'qx', 'qy'};
if ~isstruct(loads) || ~all(isfield(loads, fields))
    error('MKEF:InvalidElementLoad', 'Element loads require type, elementId, coordinateSystem, qx and qy.');
end
for i = 1:numel(loads)
    load = loads(i);
    for j = 1:numel(fields)
        value = load.(fields{j});
        if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ~isfinite(value)
            error('MKEF:InvalidElementLoad', 'Element load components and identifiers must be finite scalars.');
        end
    end
    id = load.elementId;
    if load.type ~= 20 || id ~= fix(id) || id < 1 || id > count || ...
            ~ismember(load.coordinateSystem, [1 2]) || (load.qx == 0 && load.qy == 0)
        error('MKEF:InvalidElementLoad', 'Invalid uniform element load %d.', i);
    end
    element = model.elementData(id);
    if element.type ~= 113
        error('MKEF:UnsupportedElementLoad', 'Uniform element loads require frame element 113.');
    end
    rotation = element.transformation(1:2, 1:2);
    q = [load.qx; load.qy];
    if load.coordinateSystem == 2
        globalQ = q;
        q = rotation * q;
    else
        globalQ = rotation.' * q;
    end
    L = element.length;
    localVectors{id} = localVectors{id} + ...
        [q(1)*L/2; q(2)*L/2; q(2)*L^2/12; ...
         q(1)*L/2; q(2)*L/2; -q(2)*L^2/12];
    intensities(id, :) = intensities(id, :) + q.';
    % Independent integration of the source distribution's force and moment.
    center = mean(element.nodeCoordinates(:, 1:2), 1);
    force = globalQ * L;
    resultant = resultant + [force; center(1)*force(2)-center(2)*force(1)];
end
end
