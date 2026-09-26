function plotTransientResult(result, globalDOF)
%PLOTTRANSIENTRESULT Plot a returned transient result as a convenience.
% Plotting deliberately lives outside solveTransient so the numerical core
% remains usable in headless processes.

response = result.displacements(globalDOF, :);
sampleCount = numel(response);

subplot(2, 1, 1);
plot(result.time, response);
title('Displacement (Selected DOF)');

transform = fft(response);
samplingFrequency = 1 / result.timeStep;
halfIndex = floor(sampleCount / 2) + 1;
frequency = samplingFrequency * (0:(halfIndex - 1)) / sampleCount;
amplitude = abs(transform / sampleCount);
amplitude = amplitude(1:halfIndex);
if sampleCount > 2
    lastDoubledIndex = halfIndex;
    if mod(sampleCount, 2) == 0
        lastDoubledIndex = halfIndex - 1;
    end
    amplitude(2:lastDoubledIndex) = 2 * amplitude(2:lastDoubledIndex);
end

subplot(2, 1, 2);
plot(frequency, amplitude);
title('Single Sided Amplitude Spectrum');
hold off;
end
