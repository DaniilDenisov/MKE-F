function result = solveModal(model)
%SOLVEMODAL Solve the undamped generalized eigenproblem without side effects.

[constrainedK, constrainedM, fixedDOFs] = ...
    applyFixedBoundaryConditions(model, 'modal');
[modeShapes, eigenvalueMatrix] = eig(constrainedK, constrainedM);
eigenvalues = diag(eigenvalueMatrix);
[eigenvalues, order] = sort(eigenvalues);
modeShapes = modeShapes(:, order);
angularFrequencies = sqrt(eigenvalues);

result = struct();
result.analysisType = 'modal';
result.frequenciesHz = angularFrequencies / (2 * pi);
result.angularFrequenciesRadPerSec = angularFrequencies;
result.eigenvalues = eigenvalues;
result.modeShapes = modeShapes;
result.fixedDOFs = fixedDOFs;
end
