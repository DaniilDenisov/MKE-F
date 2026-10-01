function [sigma, elementResults] = StressCalc(model, displacements)
%STRESSCALC Recover axial stresses using the current analysis-model format.
% This legacy-named function replaces the obsolete five-argument version.
% Use result.elementResults from solveStatic for complete element results.

elementResults = recoverElementResults(model, displacements);
if any([elementResults.type] ~= 112)
    error('MKEF:UnsupportedStressRecovery', ...
        'StressCalc supports truss elements only; use elementResults for frames.');
end
sigma = reshape([elementResults.axialStress], [], 1);
end
