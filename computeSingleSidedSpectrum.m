function [frequencyHz, amplitude] = ...
    computeSingleSidedSpectrum(signals, timeStep)
%COMPUTESINGLESIDEDSPECTRUM Return a correctly normalized one-sided FFT.
% Each row of signals is treated as one real-valued sampled response.

sampleCount = size(signals, 2);
samplingFrequency = 1 / timeStep;
frequencyCount = floor(sampleCount / 2) + 1;
frequencyHz = samplingFrequency * (0:(frequencyCount - 1)) / sampleCount;

twoSidedAmplitude = abs(fft(signals, [], 2) / sampleCount);
amplitude = twoSidedAmplitude(:, 1:frequencyCount);
if sampleCount > 2
    if mod(sampleCount, 2) == 0
        doubledColumns = 2:(frequencyCount - 1);
    else
        doubledColumns = 2:frequencyCount;
    end
    amplitude(:, doubledColumns) = 2 * amplitude(:, doubledColumns);
end
end
