function [sigma, elementResults] = StressCalc(model, displacements)
%STRESSCALC Recover axial stresses using the current analysis-model format.
% This legacy-named function replaces the obsolete five-argument version.
% Use result.elementResults from solveStatic for complete element results.

elementResults = recoverElementResults(model, displacements);
% For frames this is the mean axial stress, not the extreme bending stress.
sigma = reshape([elementResults.axialStress], [], 1);
end
