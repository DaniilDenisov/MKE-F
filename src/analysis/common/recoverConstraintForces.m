function result = recoverConstraintForces(result, constraints, residual, scale)
%RECOVERCONSTRAINTFORCES r = supportReactions + C'*lambda.
if isfield(constraints,'spatial')
    multipliers=zeros(size(constraints.A,1),size(residual,2));
    if ~isempty(constraints.pivots)
        multipliers=constraints.A(:,constraints.pivots).'\residual(constraints.pivots,:);
    end
    ns=size(constraints.S,1); sp=multipliers(1:ns,:); mp=multipliers(ns+1:end,:);
    result.supportReactions=full(constraints.S.'*sp);
    result.mpcForces=full(constraints.C.'*mp); result.mpcMultipliers=full(mp);
    result.nodalSupportReactions=full(constraints.physicalS.'*sp);
    result.nodalMPCForces=full(constraints.physicalC.'*mp);
    result.dependentDOFs=constraints.dependentDOFs; result.independentDOFs=constraints.independentDOFs;
    if norm(residual-result.supportReactions-result.mpcForces,inf)>1e-8*max(1,scale) || ...
            norm(constraints.T.'*residual,inf)>1e-8*max(1,scale)*max(1,norm(constraints.T,1))
        error('MKEF:EquilibriumFailure','Spatial constraint reaction recovery failed.');
    end
    return;
end
order = constraints.order;
dep = constraints.dependentDOFs(order);
C = constraints.C(order,:);
multipliers = zeros(size(C,1),size(residual,2));
if ~isempty(dep)
    % Acyclic ordering makes this system triangular, including chained MPCs.
    multipliers(order,:) = C(:,dep).' \ residual(dep,:);
end
mpcForces = full(constraints.C.' * multipliers);
supportReactions = zeros(size(residual));
fixed = constraints.fixedDOFs;
supportReactions(fixed,:) = residual(fixed,:) - mpcForces(fixed,:);
remainder = residual - supportReactions - mpcForces;
tolerance = 1e-8 * max(1,scale);
if norm(remainder,inf) > tolerance || norm(constraints.T.'*residual,inf) > tolerance*max(1,norm(constraints.T,1))
    error('MKEF:EquilibriumFailure','Constraint force recovery failed.');
end
result.supportReactions = supportReactions;
result.mpcForces = mpcForces;
result.mpcMultipliers = multipliers;
result.dependentDOFs = constraints.dependentDOFs;
result.independentDOFs = constraints.independentDOFs;
end
