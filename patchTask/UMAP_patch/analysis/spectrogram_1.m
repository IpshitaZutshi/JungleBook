%% Analyze the sleep photometry data for mice with two simultaneous recording sites
%
% USAGE
%   attempts to make a spectrogram with photometry data and ephys data
%
% INPUTS 
%    
%
%    =========================================================================

%{

basepath = pwd;
[~, currentFolderName] = fileparts(basepath);

photometry_file = dir(fullfile(basepath, '*PhotometryBehavHPC.mat'));
load(photometry_file.name)

photometry_file = dir(fullfile(basepath, '*PhotometryBehavStriatum.mat'));
load(photometry_file.name)

signal = photometry_hpc.grabDA_z;      % DA signal vector
signal = photometry_striatum.grabDA_z;      % DA signal vector

Fs = photometry_striatum.sampling_rate;       % sampling frequency in Hz
N = length(signal);                       % number of samples

% Remove mean to avoid DC offset affecting FFT
signal_detrended = signal - mean(signal);

% Compute FFT
Y = fft(signal_detrended);

% Compute two-sided spectrum then single-sided spectrum
P2 = abs(Y / N);          % two-sided spectrum magnitude
P1 = P2(1:floor(N/2)+1); % single-sided spectrum

% Double amplitudes except DC and Nyquist (if even number of points)
P1(2:end-1) = 2*P1(2:end-1);

% Frequency vector for plotting
Fs = double(Fs);
f = Fs*(0:floor(N/2))/N;

% Plot power spectrum
figure;
plot(f, P1);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Single-Sided Amplitude Spectrum of Hippocampus DA Signal');
xlim([0 10]); % Limit x-axis to 10 Hz or adjust depending on expected frequency range
grid on;

t = (0:N-1)/Fs;
plot(t, signal);
xlabel('Time (s)');
ylabel('Signal amplitude');
title('Time domain signal');

%}






% Load LFP from a specific channel
channels = [0:41 43:47 49:55 60 62:69 73:127]; %N17
channels = [70]; % N14
%lfp = bz_GetLFP(channels);
lfp = bz_GetLFP(channels, 'restrict', [0 3600]);

% channels = [5 10 12];  % Example list
% for ch = channels
%     lfp = bz_GetLFP(ch, 'restrict', [0 300]);
%     lfp_signal = double(lfp.data);
%     fs = lfp.samplingRate;
% 
%     figure;
%     spectrogram(lfp_signal, 2*fs, fs, 2^nextpow2(2*fs), fs, 'yaxis');
%     title(['Spectrogram - Channel ' num2str(ch)]);
%     ylim([0 100]);
% end

% Extract signal and sampling rate
lfp_signal = double(lfp.data);  % convert to double for spectrogram
fs = lfp.samplingRate;

% Parameters for spectrogram
window = 2 * fs;        % 2-second window
noverlap = window / 2;  % 50% overlap
nfft = 2^nextpow2(window);  % FFT length

% Design notch filter parameters
f0 = 60;                  % Frequency to remove (Hz)
Q = 35;                   % Quality factor, higher means narrower notch
wo = f0/(fs/2);           % Normalized frequency
bw = wo/Q;

% Design the notch filter
[b, a] = iirnotch(wo, bw);

% Apply the notch filter to your LFP signal
lfp_filtered = filtfilt(b, a, lfp_signal);

% Now generate the spectrogram from the filtered signal
figure;
spectrogram(lfp_filtered, window, noverlap, nfft, fs, 'yaxis');
title('LFP Spectrogram (60 Hz Notch Filter Applied)');
ylim([0 100]);
colormap jet;
colorbar;

% Generate spectrogram
figure;
spectrogram(lfp_signal, window, noverlap, nfft, fs, 'yaxis');
title('LFP Spectrogram');
ylim([0 100]);  % Adjust based on frequency range of interest
colormap jet;
colorbar;



basepath = pwd;
if ~isempty (dir(fullfile(basepath, '*PhotometryBehavHPC.mat'))) 
    hpc_photometry_file = dir(fullfile(basepath, '*PhotometryBehavHPC.mat'));
    load(hpc_photometry_file.name);
    photometry_data.hpc = photometry_hpc;
end

if ~isempty (dir(fullfile(basepath, '*PhotometryBehavStriatum.mat'))) 
    striatum_photometry_file = dir(fullfile(basepath, '*PhotometryBehavStriatum.mat'));
    load(striatum_photometry_file.name);
    photometry_data.striatum = photometry_striatum;
end


% Load photometry signal (z-scored)
hpc_signal = photometry_data.hpc.grabDA_z;
fs_phot = double(photometry_data.hpc.sampling_rate);

% Generate spectrogram
figure;
spectrogram(hpc_signal, 10*fs_phot, 5*fs_phot, [], fs_phot, 'yaxis');  % Use longer windows for slow data
title('Photometry Spectrogram (HPC)');
ylim([0 5]);  % Typically photometry is slow (<1 Hz oscillations)
colormap jet;
colorbar;


striatum_signal = photometry_data.striatum.grabDA_z;
fs_phot = double(photometry_data.striatum.sampling_rate);

% Generate spectrogram
figure;
spectrogram(striatum_signal, 10*fs_phot, 5*fs_phot, [], fs_phot, 'yaxis');  % Use longer windows for slow data
title('Photometry Spectrogram (Striatum)');
ylim([0 5]);  % Typically photometry is slow (<1 Hz oscillations)
colormap jet;
colorbar;













