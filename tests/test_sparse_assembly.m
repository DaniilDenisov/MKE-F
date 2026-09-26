function test_sparse_assembly()
%TEST_SPARSE_ASSEMBLY Verify indexed triplet assembly and linear storage.

testSuppliedCasesAreSparse();
testLargeTrussChain();
end

function testSuppliedCasesAreSparse()
options = struct('verbose', false, 'plotting', false);
caseFiles = {
    'ANSYSBeamStatic01.txt'
    'Case1ElementBeam.txt'
    'CaseATransSite.txt'
    'CaseBeam.txt'
    'CaseBeamDyn.txt'
    'CaseBeamFreq.txt'
    'CaseFig11.7p363 MarioPaz.txt'
};

for i = 1:numel(caseFiles)
    problem = StructFEProblem(caseFiles{i}, options);
    assert(issparse(problem.K));
    assert(issparse(problem.M));
    elements = problem.mesh.allMeshElems;
    for elementNumber = 1:numel(elements)
        expectedDOFs = reshape(problem.mesh.iMnod( ...
            elements(elementNumber).nodeNumbers, :).', [], 1);
        assert(isequal(elements(elementNumber).dofs, expectedDOFs));
    end
end
end

function testLargeTrussChain()
elementCount = 250;
nodeCount = elementCount + 1;
dofMap = reshape(1:(2*nodeCount), 2, nodeCount).';
elements = cell(elementCount, 1);
properties = [0.01 2e11 7850];

for elementNumber = 1:elementCount
    nodes = [elementNumber elementNumber + 1];
    coordinates = [elementNumber - 1, 0; elementNumber, 0];
    elements{elementNumber} = createStructuralElement(112, ...
        coordinates, nodes, properties, dofMap);
end
elements = vertcat(elements{:});
[K, M] = assembleGlobalMatrices(elements, 2*nodeCount);

assert(issparse(K));
assert(issparse(M));
assert(isequal(size(K), [2*nodeCount 2*nodeCount]));
assert(isequal(size(M), size(K)));
assert(nnz(K) <= 4*nodeCount);
assert(nnz(M) <= 6*nodeCount);
assert(norm(K - K.', inf) == 0);
assert(norm(M - M.', inf) == 0);
end
