function result = solveModal(model)
%SOLVEMODAL Solve the undamped generalized eigenproblem without side effects.

[fixedDOFs, freeDOFs] = partitionDOFs(model);
reducedK = model.stiffness(freeDOFs, freeDOFs);
reducedM = model.mass(freeDOFs, freeDOFs);
validateReducedSystem(reducedK, reducedM, 'modal');

[reducedModeShapes, eigenvalueMatrix] = eig(reducedK, reducedM);
eigenvalues = diag(eigenvalueMatrix);
[eigenvalues, order] = sort(eigenvalues);
reducedModeShapes = reducedModeShapes(:, order);
modeShapes = zeros(model.numberOfDOFs, numel(freeDOFs));
modeShapes(freeDOFs, :) = reducedModeShapes;
angularFrequencies = sqrt(eigenvalues);

result = struct();
result.analysisType = 'modal';
result.frequenciesHz = angularFrequencies / (2 * pi);
result.angularFrequenciesRadPerSec = angularFrequencies;
result.eigenvalues = eigenvalues;
result.modeShapes = modeShapes;
result.fixedDOFs = fixedDOFs;
result.freeDOFs = freeDOFs;
end
