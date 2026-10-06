function constraints = buildConstraintTransform(model)
%BUILDCONSTRAINTTRANSFORM Validate SPC/MPC and build sparse u = T*q.
[fixed, free] = partitionDOFs(model);
n = model.numberOfDOFs;
mpcs = struct([]);
if isfield(model, 'multiPointConstraints'), mpcs = model.multiPointConstraints; end
if ~isstruct(mpcs), error('MKEF:InvalidMPC', 'MPC records must be structures.'); end
count = numel(mpcs);
dependent = zeros(1, count);
masterDOFs = cell(1, count);
for i = 1:count
    m = mpcs(i);
    if ~all(isfield(m, {'depNode','depDOF','rhs','masters'}))
        fail(m, i, 'Missing MPC fields.');
    end
    dependent(i) = dofID(model, m.depNode, m.depDOF, m, i);
    if ~isnumeric(m.rhs) || ~isscalar(m.rhs) || ~isfinite(m.rhs) || m.rhs ~= 0
        fail(m, i, 'Only rhs = 0 is supported.');
    end
    a = m.masters;
    if ~isnumeric(a) || ~isreal(a) || isempty(a) || size(a,2) ~= 3 || any(~isfinite(a(:))) || ~any(a(:,3))
        fail(m, i, 'Expected nonzero masters [node,dof,coefficient].');
    end
    ids = zeros(1,size(a,1));
    for j = 1:size(a,1), ids(j) = dofID(model,a(j,1),a(j,2),m,i); end
    if numel(unique(ids)) ~= numel(ids), fail(m,i,'Repeated master DOF.'); end
    if any(ids == dependent(i)), fail(m,i,'Dependent DOF cannot be its own master.'); end
    if any(fixed == dependent(i)), fail(m,i,'Dependent DOF is also restrained by an SPC.'); end
    masterDOFs{i} = ids;
end
if numel(unique(dependent)) ~= count
    [~, first] = unique(dependent, 'stable');
    duplicate = setdiff(1:count, first);
    fail(mpcs(duplicate(1)),duplicate(1),'Dependent DOF occurs in more than one MPC.');
end
owner = zeros(1,n); owner(dependent) = 1:count;
indegree = zeros(1,count); children = cell(1,count);
for i = 1:count
    parents = owner(masterDOFs{i}(mpcs(i).masters(:,3).' ~= 0));
    parents = parents(parents > 0);
    indegree(i) = numel(parents);
    for parent = parents, children{parent}(end+1) = i; end
end
order = zeros(1,count); queue = find(indegree == 0); head = 1; done = 0;
while head <= numel(queue)
    i = queue(head); head = head + 1; done = done + 1; order(done) = i;
    for child = children{i}
        indegree(child) = indegree(child)-1;
        if indegree(child) == 0, queue(end+1) = child; end
    end
end
if done ~= count
    i = find(indegree > 0,1); fail(mpcs(i),i,'Cycle in MPC dependency graph.');
end
independent = setdiff(free, dependent);
T = sparse(independent,1:numel(independent),1,n,numel(independent));
entryCount = count + sum(cellfun(@numel,masterDOFs));
row = zeros(1,entryCount); col = row; val = row; cursor = 0;
for i = 1:count
    ids = masterDOFs{i};
    indices = cursor + (1:(1+numel(ids))); cursor = indices(end);
    row(indices) = i;
    col(indices) = [dependent(i) ids];
    val(indices) = [1 -mpcs(i).masters(:,3).'];
end
C = sparse(row,col,val,count,n);
for i = order
    T(dependent(i),:) = mpcs(i).masters(:,3).' * T(masterDOFs{i},:);
end
if any(~isfinite(nonzeros(T))), error('MKEF:InvalidMPC','MPC chain overflowed during elimination.'); end
constraints = struct('T',T,'C',C,'fixedDOFs',fixed,'freeDOFs',free, ...
    'dependentDOFs',dependent,'independentDOFs',independent,'order',order,'hasMPC',count > 0);
end

function id = dofID(model,node,dof,m,i)
if ~isnumeric(node) || ~isscalar(node) || ~isfinite(node) || node ~= fix(node) || node < 1 || node > model.numberOfNodes
    fail(m,i,'Invalid node reference.');
end
if ~isnumeric(dof) || ~isscalar(dof) || ~isfinite(dof) || dof ~= fix(dof) || dof < 1 || dof > model.dofPerNode
    fail(m,i,sprintf('Invalid DOF at node %d (thetaZ is unavailable for trusses).',node));
end
id = model.dofMap(node,dof);
if id == 0, fail(m,i,sprintf('Node %d: thetaZ is absent because all connected ends release Mz.',node)); end
end

function fail(m,i,message)
line = '';
if isfield(m,'sourceLine'), line = sprintf(' (source line %d)',m.sourceLine); end
error('MKEF:InvalidMPC','MPC %d%s: %s',i,line,message);
end
