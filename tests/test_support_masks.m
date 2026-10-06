function test_support_masks()
%TEST_SUPPORT_MASKS Complete catalog, canonicalization and diagnostics.
expected = {[1 2 3], [2 3], [1 3], [1 2], 1, 2, 3};
model = struct('dofPerNode',3,'numberOfNodes',1,'numberOfDOFs',3,'dofMap',[1 2 3]);
for type = 1:7
    model.fixedBoundaryConditions = normalizeSupports([type 1 0 0 0], 3);
    [fixed, free] = partitionDOFs(model);
    assert(isequal(fixed, expected{type}));
    assert(isequal(free, setdiff(1:3, expected{type})));
    assert(supportType(model.fixedBoundaryConditions,3) == type);
end
for type = 1:6
    s = normalizeSupports([type 1 0 0 0],2);
    canonical = [1 2 3 1 3 2];
    assert(supportType(s,2) == canonical(type));
end
try
    normalizeSupports([7 1 0 0 0],2);
    error('Test:MissingError','Expected invalid constraint.');
catch e
    assert(strcmp(e.identifier,'MKEF:InvalidConstraint'));
end
model.fixedBoundaryConditions = [4 1 0 0 0; 5 1 0 0 0];
try
    partitionDOFs(model);
    error('Test:MissingError','Expected duplicate constraint.');
catch e
    assert(strcmp(e.identifier,'MKEF:DuplicateConstraint'));
    assert(~isempty(strfind(e.message,'node 1 DOF ux')));
end
end
