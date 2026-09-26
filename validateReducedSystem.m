function validateReducedSystem(stiffness, mass, analysisType)
%VALIDATEREDUCEDSYSTEM Reject mechanisms and invalid reduced matrices.

if isempty(stiffness)
    if strcmp(analysisType, 'static')
        return;
    end
    error('MKEF:NoFreeDOFs', ...
        'The model has no free degrees of freedom for %s analysis.', ...
        analysisType);
end

if any(~isfinite(stiffness(:)))
    error('MKEF:InvalidReducedSystem', ...
        'The reduced stiffness matrix contains non-finite entries.');
end

symStiffness = (stiffness + stiffness.') / 2;
[~, stiffnessFailure] = chol(symStiffness);
if stiffnessFailure ~= 0
    error('MKEF:SingularStiffness', ...
        ['The reduced stiffness matrix is singular or not positive definite. ' ...
         'Check for insufficient supports or a structural mechanism.']);
end

if strcmp(analysisType, 'static')
    return;
end
if isempty(mass) || any(~isfinite(mass(:)))
    error('MKEF:InvalidReducedSystem', ...
        'The reduced mass matrix is empty or contains non-finite entries.');
end

symMass = (mass + mass.') / 2;
[~, massFailure] = chol(symMass);
if massFailure ~= 0
    error('MKEF:SingularMass', ...
        ['The reduced mass matrix is singular or not positive definite. ' ...
         'Check element connectivity and material density.']);
end
end
