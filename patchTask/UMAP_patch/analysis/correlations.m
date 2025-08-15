
basepath = pwd;
sampling_rate = 130;

%{

if ~isempty (dir(fullfile(basepath, '*PhotometryBehavHPC.mat'))) 
    hpc_photometry_file = dir(fullfile(basepath, '*PhotometryBehavHPC.mat'));
    load(hpc_photometry_file.name);

    hpc_signal_raw = photometry_hpc.grabDA_z;
    hpc_time = photometry_hpc.timestamps;
end

if ~isempty (dir(fullfile(basepath, '*PhotometryBehavStriatum.mat'))) 
    striatum_photometry_file = dir(fullfile(basepath, '*PhotometryBehavStriatum.mat'));
    load(striatum_photometry_file.name);

    str_signal_raw = photometry_striatum.grabDA_z;
    str_time = photometry_striatum.timestamps;
end

if ~isempty (dir(fullfile(basepath, '*HPCSleepPhotomSynced.mat'))) 
    hpc_photometry_file = dir(fullfile(basepath, '*HPCSleepPhotomSynced.mat'));
    load(hpc_photometry_file.name);

    hpc_signal_raw = hpc_sleep_sync.grabDA_z;
    hpc_time = hpc_sleep_sync.timestamps;
end

if ~isempty (dir(fullfile(basepath, '*StriatumSleepPhotomSynced.mat'))) 
    striatum_photometry_file = dir(fullfile(basepath, '*StriatumSleepPhotomSynced.mat'));
    load(striatum_photometry_file.name);

    str_signal_raw = striatum_sleep_sync.grabDA_z;
    str_time = striatum_sleep_sync.timestamps;
end



% ---- Interpolate both signals onto a common timebase ----
t_start = max(hpc_time(1), str_time(1));
t_end = min(hpc_time(end), str_time(end));
N = min(length(hpc_time), length(str_time));

t_common = linspace(t_start, t_end, N);
hpc_interp = interp1(hpc_time, hpc_signal_raw, t_common, 'linear');
str_interp = interp1(str_time, str_signal_raw, t_common, 'linear');

% ---- Optional: z-score the signals ----
hpc_z = zscore(hpc_interp);
str_z = zscore(str_interp);

% ---- Coherence parameters ----
window = round(sampling_rate * 60);        % 60 sec window
noverlap = round(window * 0.9); % high overlap
nfft = 4 * sampling_rate;

% ---- Compute true spectral coherence ----
[coh_true, freq] = mscohere(hpc_z, str_z, window, noverlap, nfft, sampling_rate);

% ---- Bootstrap 95% CI using shuffled signals ----
n_shuffles = 100;
coh_shuffled = zeros(n_shuffles, length(freq));

for i = 1:n_shuffles
    shift = randi(length(hpc_z));
    str_shuff = circshift(str_z, shift); % circular shift
    coh_shuffled(i,:) = mscohere(hpc_z, str_shuff, window, noverlap, nfft, sampling_rate);
end

coh_ci_lower = prctile(coh_shuffled, 2.5);
coh_ci_upper = prctile(coh_shuffled, 97.5);

% ---- Cross-correlation ----
maxLagSec = 100;
maxLagSamples = maxLagSec * sampling_rate;

[xcorr_vals, lags] = xcorr(hpc_z, str_z, maxLagSamples, 'coeff');
%[xcorr_vals, lags] = xcorr(hpc_signal_raw, str_signal_raw, maxLagSamples);
lag_times = lags / sampling_rate;

% ---- Plotting ----
figure('color','white', 'Position', [100 100 1000 500]);

% -- Coherence plot --
subplot(1,2,1);
plot(freq, coh_true, 'b', 'LineWidth', 2); hold on;
freq=freq';
fill([freq fliplr(freq)], ...
     [coh_ci_upper fliplr(coh_ci_lower)], ...
     [0.5 0.5 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
xlim([0 60]);
xlabel('Frequency (Hz)');
ylabel('Coherence');
title('HPC-Striatum Coherence');

% -- Cross-correlation plot --
subplot(1,2,2);
plot(lag_times, xcorr_vals, 'b', 'LineWidth', 2); hold on;
plot([0 0], ylim, 'r--', 'LineWidth', 1.5);
xlabel('Lags (sec)');
ylabel('Cross-correlation');
title('HPC vs Striatum Cross-Correlation');

sgtitle('Dopamine Signal Coherence and Cross-Correlation');

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% correlations at ripples

if ~isempty (dir(fullfile(basepath, '*ripples.events.mat'))) 
    ripple_file = dir(fullfile(basepath, '*ripples.events.mat'));
    load(ripple_file.name);
end
load(dir(fullfile(basepath,'*.MergePoints.events.mat')));

if pre_post == 1 
    sleep_start = MergePoints.timestamps(1,1); % start of sleep session
    sleep_end = MergePoints.timestamps(1,2); % time stamp of end of sleep session 
elseif pre_post == 2
    sleep_start = MergePoints.timestamps(3,1); 
    sleep_end = MergePoints.timestamps(3,2);  
end

%% Parameters
sampling_rate = 130;
win_sec = 10;                      % ±10 seconds window
win_samples = win_sec * sampling_rate * 2;   % total window in samples
window = round(sampling_rate * 5);           % window for coherence (5s)
noverlap = round(0.8 * window);
nfft = 2^nextpow2(window);
f_range = [0 0.5];                % Plot range for photometry signals

% change to proper start time
adjusted_ts_hpc = hpc_sleep_sync.timestamps + sleep_start;
adjusted_ts_str = striatum_sleep_sync.timestamps + sleep_start;

%% Interpolate to common timebase
N = min(length(adjusted_ts_hpc), length(adjusted_ts_str));
t_common = linspace(max(adjusted_ts_hpc(1), adjusted_ts_str(1)), ...
                    min(adjusted_ts_hpc(end), adjusted_ts_str(end)), N);

hpc_interp = interp1(adjusted_ts_hpc, hpc_sleep_sync.grabDA_z, t_common);
str_interp = interp1(adjusted_ts_str, striatum_sleep_sync.grabDA_z, t_common);

%% Loop over ripples
ripple_peaks = ripples.peaks; % Assumes ripple peaks are in seconds
good_ripples = ripple_peaks(ripple_peaks > t_common(1) + win_sec & ripple_peaks < t_common(end) - win_sec);

n_ripples = length(good_ripples);
all_coh = [];

for i = 1:n_ripples
    ripple_time = good_ripples(i);
   
    % Get indices around ripple
    idx_center = find(abs(t_common - ripple_time) == min(abs(t_common - ripple_time)), 1);
    idx_start = idx_center - win_sec*sampling_rate;
    idx_end   = idx_center + win_sec*sampling_rate - 1;
   
    if idx_start < 1 || idx_end > length(hpc_interp)
        continue
    end
   
    hpc_snip = zscore(hpc_interp(idx_start:idx_end));
    str_snip = zscore(str_interp(idx_start:idx_end));
   
    [cxy,f] = mscohere(hpc_snip, str_snip, window, noverlap, nfft, sampling_rate);
    all_coh = [all_coh; cxy'];
end

%% Average Coherence
mean_coh = mean(all_coh, 1);
sem_coh = std(all_coh, [], 1) / sqrt(size(all_coh, 1));

%% Plot
figure('color','white');
fill([f, fliplr(f)], [mean_coh + sem_coh, fliplr(mean_coh - sem_coh)], ...
    [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5); % shaded area
plot(f, mean_coh, 'b', 'LineWidth', 2); % mean coherence
% xlim(f_range);
xlabel('Frequency (Hz)');
ylabel('Coherence');
title('Ripple-Triggered HPC-Striatum Coherence');


%%%%%%%%%


% === PARAMETERS ===
sampling_rate = 130; % sampling rate in Hz
win_sec = 5; % window around each SWR (±5 sec)
win_samp = win_sec * sampling_rate;
max_lag_sec = 2; % max lag to compute cross-correlation
max_lag_samp = max_lag_sec * sampling_rate;

% === INPUT SIGNALS ===
hpc_da_signal = hpc_sleep_sync.grabDA_z; % 1D array
str_da_signal = striatum_sleep_sync.grabDA_z; % 1D array
swr_times = good_ripples; % array in seconds

swr_samples = swr_times - sleep_start;
swr_samples = round(swr_samples * sampling_rate);

% === STORAGE ===
all_xcorrs = [];

for i = 1:length(swr_samples)
    idx = swr_samples(i);
   
    % Skip if SWR window would go out of bounds
    if idx - win_samp < 1 || idx + win_samp > length(hpc_da_signal)
        continue;
    end

    hpc_win = hpc_da_signal(idx - win_samp : idx + win_samp);
    str_win = str_da_signal(idx - win_samp : idx + win_samp);

    % Z-score each segment (optional but often helpful)
    hpc_win = (hpc_win - mean(hpc_win)) / std(hpc_win);
    str_win = (str_win - mean(str_win)) / std(str_win);

    % Compute cross-correlation (normalized)
    [xc, lags] = xcorr(hpc_win, str_win, max_lag_samp, 'coeff');
    all_xcorrs = [all_xcorrs; xc'];
end


samples = win_sec*sampling_rate;

for j = 1:length(good_ripples)
    [~, ripple_idx] = min(abs(adjusted_ts_hpc - good_ripples(j)));
    start_idx = ripple_idx - samples;
    end_idx = ripple_idx + samples;
    if start_idx >= 1 && end_idx <= height(adjusted_ts_hpc)
        hpc_win = hpc_sleep_sync.grabDA_z(start_idx:end_idx);
        baseline_period = hpc_win(130:440); % using -4s to -1.6s before ripple as baseline period
        baseline_mean = mean(baseline_period);
        baseline_std = std(baseline_period);

        % Z-score the entire trial using the baseline stats
        hpc_win = (hpc_win - baseline_mean) / baseline_std;
    end

    [~, ripple_idx] = min(abs(adjusted_ts_str - good_ripples(j)));
    start_idx = ripple_idx - samples;
    end_idx = ripple_idx + samples;
    if start_idx >= 1 && end_idx <= height(adjusted_ts_str)
        str_win = striatum_sleep_sync.grabDA_z(start_idx:end_idx);
        baseline_period = str_win(130:440); % using -4s to -1.6s before ripple as baseline period
        baseline_mean = mean(baseline_period);
        baseline_std = std(baseline_period);
        
        % Z-score the entire trial using the baseline stats
        str_win = (str_win - baseline_mean) / baseline_std;
    end

    % Compute cross-correlation (normalized)
    [xc, lags] = xcorr(hpc_win, str_win, max_lag_samp, 'coeff');
    all_xcorrs = [all_xcorrs; xc'];
end

% === AVERAGE & SEM ===
mean_xc = mean(all_xcorrs, 1);
sem_xc = std(all_xcorrs, 0, 1) / sqrt(size(all_xcorrs, 1));
lags_sec = lags / sampling_rate;

% === PLOT ===
figure('color','white');
hold on;
fill([lags_sec, fliplr(lags_sec)], ...
     [mean_xc + sem_xc, fliplr(mean_xc - sem_xc)], ...
     [0.9 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(lags_sec, mean_xc, 'm', 'LineWidth', 2);
xlabel('Lag (s)');
ylabel('Cross-correlation coefficient');
title('Cross-correlogram: HPC DA vs STR DA around SWRs');
xline(0, '--k'); % zero lag line
hold off;

%%%%%%%%%%%
% correlations

%% Load data
load('rippleBandPower.mat');       % Contains ripple_power.normed_power_trace, ripple_power.timestamps
load('photometry_hpc.mat');        % Contains photometry_hpc.grabDA_z, photometry_hpc.timestamps

%% Extract signals
ripple_trace = ripple_power.normed_power_trace;
ripple_t = ripple_power.timestamps;

photometry_hpc = hpc_sleep_sync;
photo_signal = photometry_hpc.grabDA_z;
photo_t = photometry_hpc.timestamps;

%% Interpolate ripple power to photometry timebase (if needed)
if abs(median(diff(ripple_t)) - median(diff(photo_t))) > 1e-4
    ripple_interp = interp1(ripple_t, ripple_trace, photo_t, 'linear', 'extrap');
else
    ripple_interp = ripple_trace;
end

%% Find peaks in both signals
minPeakDist = 1; % seconds between peaks

[~, ripple_peak_locs] = findpeaks(ripple_interp, photo_t, 'MinPeakDistance', minPeakDist);
[~, photo_peak_locs] = findpeaks(photo_signal, photo_t, 'MinPeakDistance', minPeakDist);

%% Convert peak times into logical vectors for correlation
binSize = 0.1;  % bin size in seconds for time-alignment

% Define common time vector
common_start = max([photo_t(1), ripple_t(1)]);
common_end = min([photo_t(end), ripple_t(end)]);
time_bins = common_start:binSize:common_end;

% Bin peak occurrences
ripple_peak_binned = histcounts(ripple_peak_locs, time_bins);
photo_peak_binned = histcounts(photo_peak_locs, time_bins);

%% Correlate peak occurrence patterns
[r, p] = corr(ripple_peak_binned', photo_peak_binned');

%% Plot
figure;
subplot(2,1,1)
plot(photo_t, photo_signal, 'g'); hold on;
plot(photo_peak_locs, photo_signal(ismember(photo_t, photo_peak_locs)), 'ko');
xlabel('Time (s)'); ylabel('Photometry (Z)');
title('Photometry Signal with Peaks');

subplot(2,1,2)
plot(photo_t, ripple_interp, 'b'); hold on;
plot(ripple_peak_locs, ripple_interp(ismember(photo_t, ripple_peak_locs)), 'ro');
xlabel('Time (s)'); ylabel('Ripple Power (Z)');
title('Ripple Power Trace with Peaks');

sgtitle(sprintf('Peak Time Correlation: r = %.2f, p = %.3f', r, p));



%% Load data
load('rippleBandPower.mat');       % Contains ripple_power.normed_power_trace, ripple_power.timestamps
load('photometry_hpc.mat');        % Contains photometry_hpc.grabDA_z, photometry_hpc.timestamps

%% Extract signals
ripple_trace = ripple_power.normed_power_trace;
ripple_t = ripple_power.timestamps;

photo_signal = photometry_hpc.grabDA_z;
photo_t = photometry_hpc.timestamps;

%% Resample ripple signal to match photometry timebase
% Use interpolation to align timestamps
ripple_interp = interp1(ripple_t, ripple_trace, photo_t, 'linear', 'extrap');

%% Remove NaNs (if any)
valid_idx = ~isnan(ripple_interp) & ~isnan(photo_signal);
ripple_interp = ripple_interp(valid_idx);
photo_signal = photo_signal(valid_idx);
photo_t = photo_t(valid_idx);

%% Z-score both signals (optional, for normalization)
ripple_z = zscore(ripple_interp);
photo_z = zscore(photo_signal);

%% Compute cross-correlation
[max_lag_sec] = 10;  % maximum lag to test in seconds
sampling_rate = 1 / median(diff(photo_t));  % photometry sampling rate (should be ~50Hz or 100Hz)
max_lag_samples = round(max_lag_sec * sampling_rate);

[xcorr_vals, lags] = xcorr(photo_z, ripple_z, max_lag_samples, 'coeff');  % normalized cross-corr
lag_times = lags / sampling_rate;  % convert lag from samples to seconds

%% Plot cross-correlation
figure;
plot(lag_times, xcorr_vals, 'k', 'LineWidth', 2);
xlabel('Lag (s)');
ylabel('Cross-correlation');
title('Time-lagged Correlation Between Ripple Power and Photometry');
grid on;

%% Find peak correlation and lag
[max_corr, max_idx] = max(xcorr_vals);
best_lag = lag_times(max_idx);

% Annotate plot
hold on;
plot(best_lag, max_corr, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
text(best_lag, max_corr, sprintf('  Max r = %.2f at %.2f s lag', max_corr, best_lag), 'FontSize', 12);




