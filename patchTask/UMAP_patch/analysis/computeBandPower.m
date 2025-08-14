function [bandPower, timeVec] = computeBandPower(lfp, fs, band, winSize)
    % lfp: vector
    % fs: sampling rate (Hz)
    % band: [low high] frequency
    % winSize: sliding window size in seconds (e.g., 2)

    windowSamples = round(winSize * fs);
    stepSamples = windowSamples / 2;


    % Bandpass filter
    bpFilt = designfilt('bandpassiir','FilterOrder',4, ...
         'HalfPowerFrequency1',band(1),'HalfPowerFrequency2',band(2), ...
         'SampleRate',fs);
    lfpFilt = filtfilt(bpFilt, lfp);

    % Compute power
    powerSignal = lfpFilt .^ 2;

    % Moving average (sliding window)
    bandPower = movmean(powerSignal, windowSamples);

    % Downsample
    bandPower = bandPower(1:stepSamples:end);
    timeVec = (0:length(bandPower)-1) * stepSamples / fs;
end





