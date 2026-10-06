function supports = normalizeSupports(records, dofPerNode)
%NORMALIZESUPPORTS Canonical masks with a legacy numeric input adapter.
persistent catalog;
if isempty(catalog)
    root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    catalog = jsondecode(fileread(fullfile(root, 'shared', 'support-definitions.json')));
end
if dofPerNode == 2
    table = catalog.truss;
else
    table = catalog.frame;
end
empty = struct('node', {}, 'fixUx', {}, 'fixUy', {}, 'fixThetaZ', {});
if isnumeric(records)
    if size(records, 2) ~= 5 || any(~isfinite(records(:)))
        error('MKEF:InvalidConstraint', 'Constraints require finite [type,node,0,0,0] records.');
    end
    supports = empty;
    for i = 1:size(records, 1)
        type = records(i, 1);
        if type ~= fix(type) || type < 1 || type > size(table, 1)
            error('MKEF:InvalidConstraint', 'Constraint %d has unsupported type %g.', i, type);
        end
        if any(records(i, 3:5) ~= 0)
            error('MKEF:UnsupportedConstraintValue', 'Only homogeneous zero constraints are supported.');
        end
        supports(i, 1) = struct('node', records(i, 2), 'fixUx', logical(table(type, 1)), ...
            'fixUy', logical(table(type, 2)), 'fixThetaZ', logical(table(type, 3)));
    end
elseif isstruct(records) && all(isfield(records, fieldnames(empty)))
    supports = records(:);
else
    error('MKEF:InvalidConstraint', 'Expected support masks or legacy numeric constraints.');
end
for i = 1:numel(supports)
    s = supports(i);
    values = {s.node, s.fixUx, s.fixUy, s.fixThetaZ};
    if ~all(cellfun(@(x) (isnumeric(x) || islogical(x)) && isscalar(x) && isfinite(x), values))
        error('MKEF:InvalidConstraint', 'Constraint %d has invalid scalar fields.', i);
    end
    mask = [s.fixUx s.fixUy s.fixThetaZ];
    if any(~ismember(mask, [0 1])) || ~any(mask) || (dofPerNode == 2 && mask(3))
        error('MKEF:InvalidConstraint', 'Invalid fixed DOFs at node %g.', s.node);
    end
end
end
