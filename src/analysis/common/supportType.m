function type = supportType(support, dofPerNode)
%SUPPORTTYPE Canonical file code for a validated mask.
mask = [support.fixUx support.fixUy support.fixThetaZ];
for candidate = 1:(6 + (dofPerNode == 3))
    s = normalizeSupports([candidate support.node 0 0 0], dofPerNode);
    if isequal(logical(mask), [s.fixUx s.fixUy s.fixThetaZ])
        type = candidate;
        return;
    end
end
error('MKEF:InvalidConstraint', 'No support type for this mask.');
end
