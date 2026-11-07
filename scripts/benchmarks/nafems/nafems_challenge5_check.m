function metrics = nafems_challenge5_check(model, result, kind, n)
%NAFEMS_CHALLENGE5_CHECK Physical and algebraic checks, independent of generation.
K = model.stiffness; M = model.mass;
c = buildConstraintTransform(model); T = c.T;
Kr = T.'*K*T; Mr = T.'*M*T;
metrics = struct('dofs', size(T,2), 'massKg', 0, 'residual', 0, 'error1', NaN, 'error2', NaN);
check(norm(K-K.',inf)/norm(K,inf) <= 1e-10, 'Stiffness symmetry');
check(norm(M-M.',inf)/norm(M,inf) <= 1e-10, 'Mass symmetry');
E = 210e9; rho = 7800; A = pi*.01^2; I = pi*.01^4/4;
isFrame = strncmp(kind,'frame',5);
lengthSum = 1; if isFrame, lengthSum = 1+sqrt(2); end
for d = 1:2
    translation = zeros(model.numberOfDOFs,1);
    translation(model.dofMap(:,d)) = 1;
    mass = translation.'*M*translation;
    check(abs(mass/(rho*A*lengthSum)-1) < 1e-12, 'Total translational mass');
    metrics.massKg = mass;
end
properties = vertcat(model.elementData.properties);
check(max(max(abs(properties./repmat([A E rho I],size(properties,1),1)-1))) < 1e-12, 'Section/material');
if isFrame
    check(model.numberOfNodes == 2*n+1 && numel(model.elementData) == 2*n, 'Frame mesh size');
    check(norm(model.nodeCoordinates(1:3,:)-[0 0;0 1;1 1],inf) < 1e-14, 'Physical joints');
    check(metrics.dofs == 6*n && model.numberOfDOFs == 6*n+4, 'Frame DOF count');
    check(isequal(c.fixedDOFs(:).', [1 2 3 4]), 'Frame supports');
    check(numel(model.releases) == 4 && all(model.dofMap(1:3,3) == 0), 'Physical hinges');
    actual = [[model.releases.elementId].' [model.releases.end].'];
    check(isequal(sortrows(actual),sortrows([1 1;n 2;n+1 1;2*n 2])), 'Release placement');
    check(sum(strcmp({model.dofRegistry.kind},'elementEnd')) == 4, 'Private rotations');
    for member = 1:2
        start = [0 1]; L = 1;
        if member == 2, start = [0 0]; L = sqrt(2); end
        for j = 1:n
            e = model.elementData((member-1)*n+j);
            expected = [start+(j-1)/n*([1 1]-start); start+j/n*([1 1]-start)];
            check(norm(e.nodeCoordinates-expected,inf) < 1e-13 && abs(e.length-L/n)<1e-13, 'Member subdivision');
        end
    end
end
if strcmp(result.analysisType,'static')
    for release = model.releases
        moment = result.elementResults(release.elementId).localEndForces(3*release.end);
        check(abs(moment) <= 1e-8, 'Static released-end moment');
    end
    check(norm(result.equilibriumResidual,inf) <= 1e-8, 'Static equilibrium');
    return;
end
lambda = result.eigenvalues;
check(all(isfinite(lambda)) && all(lambda>0), 'Positive spectrum');
check(numel(lambda) == size(T,2), 'Complete modal spectrum');
% All used constraints are SPC: independent coordinates are rows of phi.
q = result.modeShapes(c.independentDOFs,:);
residual = Kr*q-Mr*q*diag(lambda);
scale = (norm(Kr,inf)+abs(lambda.')*norm(Mr,inf)).*max(abs(q),[],1);
metrics.residual = max(max(abs(residual),[],1)./scale);
check(metrics.residual <= 1e-8, 'Eigenpair residual');
if strcmp(kind,'axial') || strcmp(kind,'bending')
    k = (1:min(2,numel(lambda))).';
    if strcmp(kind,'axial')
        expected = (2*k-1)/4*sqrt(E/rho);
        limits = [.0005;.004];
    else
        expected = k.^2*pi/2*sqrt(E*I/(rho*A));
        limits = [.001;.001];
    end
    errors = abs(result.frequenciesHz(k)./expected-1);
    metrics.error1 = errors(1);
    if numel(k)>1, metrics.error2 = errors(2); end
    if n==16, check(all(errors < limits(k)), 'Analytical frequency convergence'); end
end
end

function check(condition, description)
if ~condition, error('NAFEMS5:VerificationFailed','%s failed.',description); end
end
