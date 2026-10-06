function [dofMap, registry, elements] = buildDOFRegistry(nodeCount, dofPerNode, elements, releases)
%BUILDDOFREGISTRY Assign shared nodal and private released-end coordinates.
released = false(numel(elements), 2);
for i = 1:numel(releases)
    r = releases(i);
    line = 0;
    if isfield(r, 'sourceLine'), line = r.sourceLine; end
    if ~isscalar(r.elementId) || ~isfinite(r.elementId) || r.elementId ~= fix(r.elementId) || ...
            r.elementId < 1 || r.elementId > numel(elements) || ...
            ~isscalar(r.end) || ~ismember(r.end, [1 2]) || ~strcmp(r.component, 'Mz')
        error('MKEF:InvalidRelease', 'Release at source line %d: expected existing element,end (1 or 2),Mz.', line);
    end
    if elements(r.elementId).type ~= 113 || released(r.elementId,r.end)
        error('MKEF:InvalidRelease', 'Release Mz at source line %d: element %d end %d is not a frame or is duplicated.', line,r.elementId,r.end);
    end
    released(r.elementId,r.end) = true;
end
active = true(nodeCount,dofPerNode);
if any(released(:))
    attached = false(nodeCount,1); rigid = attached;
    for i = 1:numel(elements)
        for j = 1:2
            node = elements(i).nodeNumbers(j);
            attached(node) = true;
            if ~released(i,j), rigid(node) = true; end
        end
    end
    active(attached & ~rigid,3) = false;
end
dofMap = zeros(nodeCount,dofPerNode);
registry = struct('id',{},'kind',{},'nodeId',{},'component',{},'elementId',{},'end',{});
labels = {'ux','uy','thetaZ'};
for node = 1:nodeCount
    for d = 1:dofPerNode
        if active(node,d)
            id = numel(registry)+1;
            dofMap(node,d) = id;
            registry(id) = struct('id',id,'kind','node','nodeId',node, ...
                'component',labels{d},'elementId',0,'end',0);
        end
    end
end
for i = 1:numel(elements)
    ids = reshape(dofMap(elements(i).nodeNumbers,:).',[],1);
    for j = 1:2
        if released(i,j)
            id = numel(registry)+1;
            ids(3*j) = id;
            registry(id) = struct('id',id,'kind','elementEnd', ...
                'nodeId',elements(i).nodeNumbers(j),'component','thetaZ','elementId',i,'end',j);
        end
    end
    elements(i).dofs = ids;
end
end
