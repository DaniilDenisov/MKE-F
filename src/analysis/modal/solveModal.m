function result = solveModal(model)
%SOLVEMODAL Solve the undamped generalized eigenproblem without side effects.

rejectElementLoads(model);
constraints = buildConstraintTransform(model);
fixedDOFs = constraints.fixedDOFs;
freeDOFs = constraints.freeDOFs;
T = constraints.T;
reducedK = model.stiffness(freeDOFs, freeDOFs);
reducedM = model.mass(freeDOFs, freeDOFs);
if constraints.hasMPC
    reducedK = T.' * model.stiffness * T;
    reducedM = T.' * model.mass * T;
end
validateReducedSystem(reducedK, reducedM, 'modal');

[reducedModeShapes, eigenvalueMatrix] = eig(reducedK, reducedM);
eigenvalues = diag(eigenvalueMatrix);
[eigenvalues, order] = sort(eigenvalues);
reducedModeShapes = reducedModeShapes(:, order);
modeShapes = full(T * reducedModeShapes);
angularFrequencies = sqrt(eigenvalues);

result = struct();
result.analysisType = 'modal';
result.frequenciesHz = angularFrequencies / (2 * pi);
result.angularFrequenciesRadPerSec = angularFrequencies;
result.eigenvalues = eigenvalues;
result.modeShapes = modeShapes;
result.fixedDOFs = fixedDOFs;
result.freeDOFs = freeDOFs;
residual = model.stiffness * modeShapes - model.mass * modeShapes * diag(eigenvalues);
result = recoverConstraintForces(result,constraints,residual,norm(model.stiffness*modeShapes,inf));
end
