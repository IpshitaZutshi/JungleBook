%% Plotting Correlations btwn HPC DA and striatum DA 
%
% USAGE
%   Correlation btwn HPC DA and striatum DA across the sleep states
%
% INPUTS 
%    
%
%    =========================================================================

%% PARAMETERS

basepath = pwd;
save_path = '\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\patchTask\figures';
saveLoc = strcat(save_path,'\sleep_phot');

session_info_file = dir(fullfile(basepath, '*sessionInfo.mat'));
load(session_info_file.name);
tokens = regexp(sessionInfo.FileName, 'sess(\d+)$', 'tokens');
session_str = tokens{1}{1}; % string of the number session

sampling_rate = 130; % sampling rate
window = 5;
samples = window * sampling_rate;
max_lag_sec = 2;
max_lag_samp = max_lag_sec * sampling_rate;

states = {'NREM', 'REM', 'Wake'};
colors = {[0.7412    0.8784    0.9176], [0.2431    0.7451    0.8706], [0.1922    0.4000    0.5804]}; % colors for plotting

% Convert SleepStateEpisodes to intervals
% Should be a struct or table with 'state', 'start', and 'stop' fields

%% Loop over states

merge_file = dir(fullfile(basepath,'*.MergePoints.events.mat'));
load(merge_file.name);
sleep_state_file = dir(fullfile(basepath,'*.SleepState.states.mat'));
load(sleep_state_file.name);
ripple_file = dir(fullfile(basepath, '*ripples.events.mat'));
load(ripple_file.name)

% Aggregate all sleep data across sessions and pre/post
agg_all = struct('NREM', [], 'REM', [], 'WAKE', []);
high_corr_times = [];
high_corr_vals = [];

for ii = 2:size(MergePoints.foldernames,2)
    cd(basepath);
    if isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep 'top*']))
        cd([basepath filesep MergePoints.foldernames{ii}]); 
        % set pre_post to be 1 for pre behavior sleep, 2 for post behavior sleep  
        if ii == 1
            pre_post = 1;
        else
            pre_post = 2;
        end
        
        curr_folder = pwd;
        if ~isempty (dir(fullfile(curr_folder, '*HPCSleepPhotomSynced.mat'))) 
            load(dir(fullfile(curr_folder,'*HPCSleepPhotomSynced.mat')).name);
        end

        if ~isempty (dir(fullfile(curr_folder, '*StriatumSleepPhotomSynced.mat')))
            load(dir(fullfile(curr_folder,'*StriatumSleepPhotomSynced.mat')).name);
        end

        if pre_post == 1 
            sleep_start = MergePoints.timestamps(1,1); % start of sleep session
            sleep_end = MergePoints.timestamps(1,2); % time stamp of end of sleep session 
        elseif pre_post == 2
            sleep_start = MergePoints.timestamps(3,1); 
            sleep_end = MergePoints.timestamps(3,2);  
        end
        
        % INSIDE HERE
        figure('color','white'); 
        hold on;
        legend_entries = {};
        
        state_names = fieldnames(SleepState.ints);  % {'WAKEstate', 'NREMstate', 'REMstate'}
        
        for s = 1:length(state_names)
            curr_state_name = state_names{s};
            curr_state_epochs = SleepState.ints.(curr_state_name);  % [N×2] matrix of [start, end] times

            valid_epochs = curr_state_epochs(:,1) > MergePoints.timestamps(ii, 1) & ...
               curr_state_epochs(:,2) < MergePoints.timestamps(ii, 2);
            curr_state_epochs = curr_state_epochs(valid_epochs, :);
        
            if isempty(curr_state_epochs)
                continue;
            end
        
            all_xcorrs = [];
            adjusted_ts_hpc = hpc_sleep_sync.timestamps + sleep_start;
        
            for i = 1:size(curr_state_epochs,1)
                    start_time = curr_state_epochs(i,1);
                    end_time   = curr_state_epochs(i,2);

                    if end_time > adjusted_ts_hpc(end)
                        continue;
                    end

                    [~, start_idx] = min(abs(adjusted_ts_hpc - start_time));
                    [~, end_idx] = min(abs(adjusted_ts_hpc - end_time));

                    % Skip if the window is too short
                    if end_idx - start_idx < 2 * max_lag_samp
                        continue;
                    end
                    if end_idx > length(adjusted_ts_hpc)
                        continue;
                    end
            
                    % Get dopamine signal windows
                    hpc_win = hpc_sleep_sync.grabDA_z(start_idx:end_idx);
                    str_win = striatum_sleep_sync.grabDA_z(start_idx:end_idx);
            
                    % Z-score each segment
                    hpc_win = (hpc_win - mean(hpc_win)) / std(hpc_win);
                    str_win = (str_win - mean(str_win)) / std(str_win);
            
                    % Cross-correlation
                    [xc, lags] = xcorr(hpc_win, str_win, max_lag_samp, 'coeff');
                    all_xcorrs = [all_xcorrs; xc'];
                    
                    state_clean = upper(curr_state_name(1:end-5));  % 'NREM', 'REM', 'WAKE'
                    
                    agg_all.(state_clean) = [agg_all.(state_clean); xc'];

                    % track times of max correlation
                    [~, max_idx] = max(abs(xc));  % or use xc(max_lag_samp+1) for lag=0
                    peak_corr = xc(max_idx);
                    
                    % Save center time of this window
                    center_time = (start_time + end_time) / 2;
                    
                    % Store correlation and timestamp
                    high_corr_times = [high_corr_times; center_time];
                    high_corr_vals = [high_corr_vals; peak_corr];

            end
        
            if ~isempty(all_xcorrs)
                mean_xc = mean(all_xcorrs, 1);
                sem_xc = std(all_xcorrs, 0, 1) / sqrt(size(all_xcorrs,1));
                lags_sec = lags / sampling_rate;
        
                % Plot
                fill([lags_sec, fliplr(lags_sec)], ...
                     [mean_xc + sem_xc, fliplr(mean_xc - sem_xc)], ...
                     colors{s}, 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'HandleVisibility', 'off');
                plot(lags_sec, mean_xc, 'Color', colors{s} * 0.6, 'LineWidth', 2);
                legend_entries{end+1} = curr_state_name(1:end-5);  % e.g., 'WAKE'
            end
        end
        if pre_post == 1
            pp = 'Pre Sleep';
        else
            pp = 'Post Sleep';
        end

        xline(0, '--r');
        xlabel('Lag (s)');
        ylabel('Cross-correlation coefficient');
        title('HPC DA vs STR DA during sleep states - ', pp);
        legend(legend_entries);

        %% Cross-correlation around sharp wave ripples
        % Load ripple events
        
        
            win_sec = 5;  % window around each ripple (±5 seconds)
            win_samples = round(win_sec * sampling_rate);
            lags = -max_lag_samp:max_lag_samp;
            lags_sec = lags / sampling_rate;
            
            all_ripple_xcorrs = [];
        
            for r = 1:size(ripples.peaks, 1)
                center_time = ripples.peaks(r);  % start of ripple
        
                % Align to ripple time
                [~, center_idx] = min(abs(hpc_sleep_sync.timestamps + sleep_start - center_time));
        
                idx_start = center_idx - win_samples;
                idx_end   = center_idx + win_samples;
        
                if idx_start < 1 || idx_end > length(hpc_sleep_sync.timestamps)
                    continue;  % skip if window exceeds bounds
                end
        
                hpc_win = hpc_sleep_sync.grabDA_z(idx_start:idx_end);
                str_win = striatum_sleep_sync.grabDA_z(idx_start:idx_end);
        
                % Z-score each segment
                hpc_win = (hpc_win - mean(hpc_win)) / std(hpc_win);
                str_win = (str_win - mean(str_win)) / std(str_win);
        
                [xc, ~] = xcorr(hpc_win, str_win, max_lag_samp, 'coeff');
                all_ripple_xcorrs = [all_ripple_xcorrs; xc'];
            end
        
            % Plot ripple-aligned xcorrs
            figure('color','white');
            hold on;
        
            if ~isempty(all_ripple_xcorrs)
                mean_xc = mean(all_ripple_xcorrs, 1);
                sem_xc = std(all_ripple_xcorrs, 0, 1) / sqrt(size(all_ripple_xcorrs, 1));
        
                fill([lags_sec, fliplr(lags_sec)], ...
                     [mean_xc + sem_xc, fliplr(mean_xc - sem_xc)], ...
                     [0.2706    0.1020    0.9412], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        
                plot(lags_sec, mean_xc, 'Color', [0.2706    0.1020    0.9412], 'LineWidth', 2);
                xline(0, '--r');
                xlabel('Lag (s)');
                ylabel('Cross-correlation coefficient');
                title('HPC–STR DA Correlation around Ripples');
                grid on;
                set(gcf,'position',[500,200,1120,840])
                %saveas(gcf,[saveLoc,filesep ,'HPC_str_ripple_corr_', session_str, '.png'],'png');
            else
                disp('No valid ripple-aligned segments found.');
            end

    end
end

%% high correlation
% === After session loop: find high correlation times ===
thresh = prctile(high_corr_vals, 95);  % top 5% most correlated
top_corr_times = high_corr_times(high_corr_vals >= thresh);

ripple_win = 0.5; % seconds
ripple_match = false(size(top_corr_times));

for i = 1:length(top_corr_times)
    if any(abs(ripples.peaks - top_corr_times(i)) < ripple_win)
        ripple_match(i) = true;
    end
end

matched_times = top_corr_times(ripple_match);


% === One Combined Aggregate Plot ===
% pre sleep and post sleep session together

colors = {[0.2706    0.1020    0.9412], [0.5451    0.9412    0.9686], [0.2431    0.7020    0.9294]}; % colors for plotting
lags_sec = lags / sampling_rate;

figure('color','white'); hold on;
legend_entries = {};

for s = 1:length(states)
    state = upper(states{s});
    xcs = agg_all.(state);

    if isempty(xcs)
        continue;
    end

    mean_xc = mean(xcs, 1);
    sem_xc = std(xcs, 0, 1) / sqrt(size(xcs, 1));

    fill([lags_sec, fliplr(lags_sec)], ...
         [mean_xc + sem_xc, fliplr(mean_xc - sem_xc)], ...
         colors{s}, 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'HandleVisibility', 'off');

    plot(lags_sec, mean_xc, 'Color', colors{s} * 0.6, 'LineWidth', 2);
    legend_entries{end+1} = state;
end

xline(0, '--r');
xlabel('Lag (s)');
ylabel('Cross-correlation coefficient');
title('HPC ↔ STR Cross-Correlation During Sleep');
legend(legend_entries);
set(gcf,'position',[500,200,1120,840])
grid on;
%saveas(gcf,[saveLoc,filesep ,'HPC_str_sleep_corr_', session_str, '.png'],'png');


%% NEW ABOVE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%{
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
window = 10;                      % ±10 seconds window
win_samples = window * sampling_rate * 2;   % total window in samples
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
good_ripples = ripple_peaks(ripple_peaks > t_common(1) + window & ripple_peaks < t_common(end) - window);

n_ripples = length(good_ripples);
all_coh = [];

for i = 1:n_ripples
    ripple_time = good_ripples(i);
   
    % Get indices around ripple
    idx_center = find(abs(t_common - ripple_time) == min(abs(t_common - ripple_time)), 1);
    idx_start = idx_center - window*sampling_rate;
    idx_end   = idx_center + window*sampling_rate - 1;
   
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
window = 5; % window around each SWR (±5 sec)
win_samp = window * sampling_rate;
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


samples = window*sampling_rate;

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
figure('color','white');
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
figure('color','white');
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



% ripples durin gdiff states

        % figure('color','white'); 
        % hold on;
        % legend_entries = {};
        % 
        % % Get all state names: WAKEstate, NREMstate, REMstate
        % state_names = fieldnames(SleepState.ints);
        % 
        % for s = 1:length(state_names)
        %     curr_state_name = state_names{s};  % e.g., 'WAKEstate'
        %     curr_state_epochs = SleepState.ints.(curr_state_name);  % [N×2] double
        % 
        %     % Skip if there are no epochs in this state
        %     if isempty(curr_state_epochs)
        %         continue;
        %     end
        % 
        %     % === Find SWRs that fall within any epoch of this state ===
        %     swr_in_state = false(size(swr_times));
        %     for i = 1:size(curr_state_epochs,1)
        %         start_time = curr_state_epochs(i,1);
        %         stop_time  = curr_state_epochs(i,2);
        %         swr_in_state = swr_in_state | ...
        %             (swr_times >= start_time & swr_times <= stop_time);
        %     end
        % 
        %     swr_subset = swr_times(swr_in_state);
        %     swr_samples = round(swr_subset * sampling_rate);
        % 
        %     all_xcorrs = [];
        % 
        %     for i = 1:length(swr_samples)
        %         idx = swr_samples(i);
        % 
        %         if idx - samples < 1 || idx + samples > length(hpc_da_signal)
        %             continue;
        %         end
        % 
        %         hpc_win = hpc_da_signal(idx - samples : idx + samples);
        %         str_win = str_da_signal(idx - samples : idx + samples);
        % 
        %         % Z-score
        %         hpc_win = (hpc_win - mean(hpc_win)) / std(hpc_win);
        %         str_win = (str_win - mean(str_win)) / std(str_win);
        % 
        %         [xc, lags] = xcorr(hpc_win, str_win, max_lag_samp, 'coeff');
        %         all_xcorrs = [all_xcorrs; xc'];
        %     end
        % 
        %     % === Plot mean and SEM cross-correlation for this state ===
        %     if ~isempty(all_xcorrs)
        %         mean_xc = mean(all_xcorrs, 1);
        %         sem_xc = std(all_xcorrs, 0, 1) / sqrt(size(all_xcorrs,1));
        %         lags_sec = lags / sampling_rate;
        % 
        %         fill([lags_sec, fliplr(lags_sec)], ...
        %              [mean_xc + sem_xc, fliplr(mean_xc - sem_xc)], ...
        %              colors{s}, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        %         plot(lags_sec, mean_xc, 'Color', colors{s} * 0.6, 'LineWidth', 2);
        % 
        %         % Add to legend (strip off 'state' suffix)
        %         legend_entries{end+1} = curr_state_name(1:end-5);
        %     end
        % end
        % 
        % xline(0, '--k');
        % xlabel('Lag (s)');
        % ylabel('Cross-correlation coefficient');
        % title('Cross-correlogram: HPC DA vs STR DA during SWRs (by sleep state)');
        % legend(legend_entries);


%}
