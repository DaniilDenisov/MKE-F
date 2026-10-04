function loadVector = buildStaticLoad(model)
%BUILDSTATICLOAD Combine fresh nodal and consistent element load vectors.
loadVector = buildNodalLoadVector(model) + buildElementLoadVector(model);
end
