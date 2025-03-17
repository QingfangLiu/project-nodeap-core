
function [avgPowerFiltered, powerRatio] = computePowerMeasures(data, sampleInterval, lowFreq, highFreq)
    % This function calculates alff and falff
    %   alff: average power in the specified frequency range 
    %   falff: the fraction of this power in the frequency range to the power of the whole frequency range
    %
    % Parameters:
    %   data - the time series data vector
    %   sampleInterval - the time interval between samples (in seconds)
    %   lowFreq - the lower bound of the frequency range (in Hz)
    %   highFreq - the upper bound of the frequency range (in Hz)
    %
    % Returns:
    %   avgPowerFiltered - the average power in the specified frequency range
    %                      
    %   powerRatio - the ratio of avgPowerFiltered to the average power without
    %                bandpass filtering

    N = length(data);    % Number of data points
    fs = 1 / sampleInterval;   % Sampling frequency
    [powerSpectrum, freq] = periodogram(data, [], N, fs); % use matlab function to do fft & compute power spectrum

    % Plot the periodogram
%     figure;
%     plot(freq, powerSpectrum);
%     title('Periodogram Power Spectral Density Estimate');
%     xlabel('Frequency (Hz)');
%     ylabel('Power/Frequency (dB/Hz)');

    freqIdx = freq >= lowFreq & freq <= highFreq;  % Find indices of frequencies within the desired range
    allPowerFiltered = sqrt(powerSpectrum(freqIdx));  % sqrt of the power spectrum within the given freq range
    sumPowerFiltered = sum(allPowerFiltered);
    avgPowerFiltered = mean(allPowerFiltered);

    allfreqIdx = freq >= 0 & freq <= 0.25;  % find indices of freq [0,0.25]
    allPower = sqrt(powerSpectrum(allfreqIdx)); % sqrt of the whole power spectrum
    sumPower = sum(allPower);       
    powerRatio = sumPowerFiltered / sumPower; % ratio of the spectrum within the freq range relative to the whole freq spectrum

end


%%
% to test and display this function
% data = dat(:,1);
% TR = 1.5;
% sampleInterval = TR; 
% lowFreq = 0.01; % Lower bound of frequency range in Hz
% highFreq = 0.08; % Upper bound of frequency range in Hz
% t = 1.5:1.5:(N*1.5); % time vector 

