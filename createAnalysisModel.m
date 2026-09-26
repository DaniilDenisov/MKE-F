function model = createAnalysisModel(K, M, mesh)
%CREATEANALYSISMODEL Copy assembled model data for side-effect-free solvers.
% The returned value is a plain struct. Numerical solvers receive this value
% instead of the mutable StructFEProblem facade or its FEMesh handle.

if ~ismatrix(K) || size(K, 1) ~= size(K, 2)
    error('MKEF:InvalidModel', 'The stiffness matrix must be square.');
end
if ~isequal(size(M), size(K))
    error('MKEF:InvalidModel', ...
        'The mass and stiffness matrices must have the same dimensions.');
end

numberOfDOFs = mesh.numberOfNodes * mesh.dofPerNode;
if size(K, 1) ~= numberOfDOFs
    error('MKEF:InvalidModel', ...
        'Matrix dimensions do not match the mesh degrees of freedom.');
end

model = struct();
model.stiffness = K;
model.mass = M;
model.dofMap = mesh.iMnod;
model.fixedBoundaryConditions = mesh.allFixBCs;
model.forceBoundaryConditions = mesh.allForceBCs;
model.numberOfNodes = mesh.numberOfNodes;
model.dofPerNode = mesh.dofPerNode;
model.numberOfDOFs = numberOfDOFs;
end
