function test_element_registry()
r = elementTypeInfo();
assert(isequal([r.code],[112 113 1133 1436]));
assert(isequal({r.implementation},{'truss2d','frame2d','truss3d','frame3d'}));
assert(isequal([r.dimension],[2 2 3 3]));
for code = [112 113 1133 1436]
    info=elementTypeInfo(code);
    p=ones(1,numel(info.properties));
    coords=[zeros(1,info.dimension);1 zeros(1,info.dimension-1)];
    if code==1436, p=[p 0 1 0 0]; end
    map=reshape(1:2*info.dofPerNode,info.dofPerNode,2).';
    e=createStructuralElement(code,coords,[1 2],p,map);
    assert(isequal(size(e.stiffness),[numel(map) numel(map)]));
end
assert(isequal(elementTypeInfo(1133).dofLabels,{'ux';'uy';'uz'}));
assert(isequal(elementTypeInfo(1436).dofLabels,{'ux';'uy';'uz';'thetaX';'thetaY';'thetaZ'}));
failed=false; try, elementTypeInfo(999); catch err, failed=strcmp(err.identifier,'MKEF:UnsupportedElementType'); end
assert(failed);
end
