function plotTransientResult(result, globalDOF)
%PLOTTRANSIENTRESULT Plot a returned transient result as a convenience.
% Plotting deliberately lives outside solveTransient so the numerical core
% remains usable in headless processes.

response = result.displacements(globalDOF, :);

subplot(2, 1, 1);
plot(result.time, response);
title('Displacement (Selected DOF)');

subplot(2, 1, 2);
plot(result.spectrumFrequencyHz, ...
    result.displacementAmplitudeSpectrum(globalDOF, :));
title('Single Sided Amplitude Spectrum');
hold off;
end
