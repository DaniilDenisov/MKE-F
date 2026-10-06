function result = recoverConstraintForces(result, constraints, residual, scale)
%RECOVERCONSTRAINTFORCES r = supportReactions + C'*lambda.
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
