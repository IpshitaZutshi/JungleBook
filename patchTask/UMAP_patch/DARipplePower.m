

basepath = pwd;
[~, currentFolderName] = fileparts(basepath);
color = 0;

% load ripple file
ripple_file = dir(fullfile(basepath, '*ripples.events.mat'));
load(ripple_file.name)

%% average ripple power around DA peaks - behav
% load hpc ohtometry
HPC_photometry_file = dir(fullfile(basepath, '*PhotometryBehavHPC.mat'));
load(HPC_photometry_file.name);

% load striatum photometry
striatum_photometry_file = dir(fullfile(basepath, '*PhotometryBehavStriatum.mat'));
load(striatum_photometry_file.name);

% load ripple power
ripple_file = dir(fullfile(basepath, '*rippleBandPower.mat'));
load(ripple_file.name);
%ripple_period = ripple_power.normed_power_trace(ripple_power.timestamps <= behav_end & ripple_power.timestamps >= behav_start);

[peakVals, peakLocs] = findpeaks(photometry_hpc.grabDA_z, 'MinPeakHeight', 0.05);
peakTimes = photometry_hpc.timestamps(peakLocs);

% Combine into a two-column matrix: [peak value, timestamp]
da_peaks = [peakVals, peakTimes];

sampling_rate = 1250;
% power around peak
window = 5; % window of time around ripple to average
samples = window*sampling_rate;

% initialize matrix
ripple_matrix = nan(length(da_peaks), (samples*2)+1); 

% average power within a specified time window around peaks
for j = 1:length(da_peaks)
    [~, peak_idx] = min(abs(ripple_power.timestamps - da_peaks(j, 2)));
    start_idx = peak_idx - samples;
    end_idx = peak_idx + samples;
    if start_idx >= 1 && end_idx <= length(ripple_power.timestamps)
        ripple_matrix(j, :) = ripple_power.normed_power_trace(start_idx:end_idx);
        
        % z score
        baseline = ripple_matrix(j, :);

        % Calculate baseline mean and std
        baseline_mean = mean(baseline);
        baseline_std = std(baseline);

        % Z-score the entire trial using the baseline stats
        ripple_matrix(j, :) = (ripple_matrix(j, :) - baseline_mean) / baseline_std;
    end
end

ripple_matrix_2 = ripple_matrix;
ripple_matrix_2(any(isnan(ripple_matrix_2), 2), :) = [];
median_ripple = median(ripple_matrix, 1); % median at each timepoint
time = linspace(-window, window, ((samples*2)+1));

% calculate confidence intervals
N = height(ripple_matrix);                          % Number of Experiments In Data Set
avg_ripple = mean(ripple_matrix, 1);              % Mean Of All Experiments At Each Value Of ,x 
ripple_SEM = std(ripple_matrix, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
ripple_CI95 = bsxfun(@times, ripple_SEM, CI95(:));  % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of exw
smooth_CI95 = smoothdata(ripple_CI95, 2);

% statistical significance of peak before lick
ripple_baseline = ripple_matrix(:,[1:120]);
avg_ripple_baseline = mean(ripple_baseline, 1);
mn = mean(avg_ripple_baseline);
st_d = std(avg_ripple_baseline);
sample_mn = min(median_ripple);
deg_free = length(da_peaks)-1;
t = (sample_mn - mn)/(st_d/(sqrt(deg_free)));

%% Plot 

if color == 0
    % avg_color = [ 0.2392    0.2863    0.9608];
    % conf_color = 'b';
    avg_color = [0.960784313725490, 0.152941176470588, 0.905882352941176];
    conf_color = [0.960784313725490, 0.152941176470588, 0.905882352941176];
else
    avg_color = 'g';
    conf_color = [0.7176    0.9412    0.1020];
end

figure('color','white');
plot(time, median_ripple, 'g', 'LineWidth', 2);
hold on
fill([time,fliplr(time)], [(ripple_CI95(1,:)+median_ripple),fliplr((ripple_CI95(2,:)+median_ripple))], 'b', 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title('Ripple power around DA peaks');
grid on;
hold off


smoothed = smoothdata(median_ripple);
figure('color','white');
plot(time, smoothed, 'g', 'LineWidth', 2);
hold on
fill([time,fliplr(time)], [(smooth_CI95(1,:)+smoothed),fliplr((smooth_CI95(2,:)+smoothed))], 'b', 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score Around Ripples (smoothed)');
grid on;
hold off




%% look at sleep DA
cd(basepath);

% load sleep states
sleep_state_file = dir(fullfile(basepath, '*SleepStateEpisodes.states.mat'));
load(ripple_file.name)

sleep_merge = []; % will hold the merge points (start and stop times) of just sleep sessions
[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);
if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
    load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
    k = 1;
    for ii = 1:size(MergePoints.foldernames,2)
        if isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep 'top*']))
            cd([basepath filesep MergePoints.foldernames{ii}]); 
            curr_folder = pwd;
            % load hpc photometry
            HPC_photometry_file = dir(fullfile(curr_folder, '*HPC*PhotomSynced.mat'));
            load(HPC_photometry_file.name);
            hpc_sleep{k} = hpc_sleep_sync;
            hpc_sleep{k}.timestamps = hpc_sleep{k}.timestamps + MergePoints.timestamps(ii, 1);

            % load striatum photometry
            striatum_photometry_file = dir(fullfile(curr_folder, '*striatum*PhotomSynced.mat'));
            load(striatum_photometry_file.name);
            striatum_sleep{k} = striatum_sleep_sync;
            striatum_sleep{k}.timestamps = striatum_sleep{k}.timestamps + MergePoints.timestamps(ii, 1);

            sleep_merge(k, 1) = MergePoints.timestamps(ii,1); 
            sleep_merge(k, 2) = MergePoints.timestamps(ii,2); 
            k = k+1;
        end
    end
end

cd(basepath);

% combine pre and post sleep
hpc_full_sleep.timestamps = [hpc_sleep{1,1}.timestamps; hpc_sleep{1,2}.timestamps];
hpc_full_sleep.grabDA_z = [hpc_sleep{1,1}.grabDA_z; hpc_sleep{1,2}.grabDA_z];
striatum_full_sleep.timestamps = [striatum_sleep{1,1}.timestamps; striatum_sleep{1,2}.timestamps];
striatum_full_sleep.grabDA_z = [striatum_sleep{1,1}.grabDA_z; striatum_sleep{1,2}.grabDA_z];

% include only NREM
hpc_nrem.timestamps = [];
hpc_nrem.grabDA_z = [];
for i = 1:length(SleepStateEpisodes.ints.NREMepisode)
    hpc_nrem.timestamps = [hpc_nrem.timestamps;...
        hpc_full_sleep.timestamps(hpc_full_sleep.timestamps <= SleepStateEpisodes.ints.NREMepisode(i, 2) & ...
        hpc_full_sleep.timestamps >= SleepStateEpisodes.ints.NREMepisode(i, 1))];
    hpc_nrem.grabDA_z = [hpc_nrem.grabDA_z;...
        hpc_full_sleep.grabDA_z(hpc_full_sleep.timestamps <= SleepStateEpisodes.ints.NREMepisode(i, 2) & ...
        hpc_full_sleep.timestamps >= SleepStateEpisodes.ints.NREMepisode(i, 1))]; 
end

[peakVals, peakLocs] = findpeaks(hpc_nrem.grabDA_z, 'MinPeakHeight', 0.05);
peakTimes = hpc_nrem.timestamps(peakLocs);

% Combine into a two-column matrix: [peak value, timestamp]
da_peaks = [peakVals, peakTimes];

sampling_rate = 1250;
% power around peak
window = 5; % window of time around ripple to average
samples = window*sampling_rate;

% initialize matrix
ripple_matrix = nan(length(da_peaks), (samples*2)+1); 

% average power within a specified time window around peaks
for j = 1:length(da_peaks)
    [~, peak_idx] = min(abs(ripple_power.timestamps - da_peaks(j, 2)));
    start_idx = peak_idx - samples;
    end_idx = peak_idx + samples;
    if start_idx >= 1 && end_idx <= length(ripple_power.timestamps) % issue where if a peak is right at the edge of nrem, this could include samples from a different nrem episode
        ripple_matrix(j, :) = ripple_power.normed_power_trace(start_idx:end_idx);
        
        % z score
        baseline = ripple_matrix(j, :);

        % Calculate baseline mean and std
        baseline_mean = mean(baseline);
        baseline_std = std(baseline);

        % Z-score the entire trial using the baseline stats
        ripple_matrix(j, :) = (ripple_matrix(j, :) - baseline_mean) / baseline_std;
    end
end

%
baseline_start = 5; % seconds before peak (start of baseline window)
baseline_end = 1;   % seconds before peak (end of baseline window)

baseline_samples_start = round(baseline_start * sampling_rate);
baseline_samples_end = round(baseline_end * sampling_rate);

for j = 1:length(da_peaks)
    [~, peak_idx] = min(abs(ripple_power.timestamps - da_peaks(j, 2)));
    start_idx = peak_idx - samples;
    end_idx = peak_idx + samples;
    
    if start_idx >= 1 && end_idx <= length(ripple_power.timestamps)
        segment = ripple_power.normed_power_trace(start_idx:end_idx);
        
        % Define baseline indices relative to segment
        baseline_indices = (samples + baseline_samples_start):(samples + baseline_samples_end);
        
        baseline_data = segment(baseline_indices);
        
        % Compute mean and std of baseline period only
        baseline_mean = mean(baseline_data);
        baseline_std = std(baseline_data);
        
        % Z-score entire segment using baseline mean and std
        ripple_matrix(j, :) = (segment - baseline_mean) / baseline_std;
    else
        ripple_matrix(j, :) = nan(1, (samples*2)+1);
    end
end

ripple_matrix(any(isnan(ripple_matrix), 2), :) = [];
median_ripple = median(ripple_matrix, 1); % median at each timepoint
time = linspace(-window, window, ((samples*2)+1));

% calculate confidence intervals
N = height(ripple_matrix);                          % Number of Experiments In Data Set
avg_ripple = mean(ripple_matrix, 1);              % Mean Of All Experiments At Each Value Of ,x 
ripple_SEM = std(ripple_matrix, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
ripple_CI95 = bsxfun(@times, ripple_SEM, CI95(:));  % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of exw
smooth_CI95 = smoothdata(ripple_CI95, 2);

% statistical significance 
ripple_baseline = ripple_matrix(:,[1:120]);
avg_ripple_baseline = mean(ripple_baseline, 1);
mn = mean(avg_ripple_baseline);
st_d = std(avg_ripple_baseline);
sample_mn = min(median_ripple);
deg_free = length(da_peaks)-1;
t = (sample_mn - mn)/(st_d/(sqrt(deg_free)));

%% Plot 

if color == 0
    % avg_color = [ 0.2392    0.2863    0.9608];
    % conf_color = 'b';
    avg_color = [0.960784313725490, 0.152941176470588, 0.905882352941176];
    conf_color = [0.960784313725490, 0.152941176470588, 0.905882352941176];
else
    avg_color = 'g';
    conf_color = [0.7176    0.9412    0.1020];
end

figure('color','white');
plot(time, median_ripple, 'g', 'LineWidth', 2);
hold on
fill([time,fliplr(time)], [(ripple_CI95(1,:)+median_ripple),fliplr((ripple_CI95(2,:)+median_ripple))], 'b', 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title('Ripple power around DA peaks - NREM');
grid on;
hold off


smoothed = smoothdata(median_ripple);
figure('color','white');
plot(time, smoothed, 'g', 'LineWidth', 2);
hold on
fill([time,fliplr(time)], [(smooth_CI95(1,:)+smoothed),fliplr((smooth_CI95(2,:)+smoothed))], 'b', 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title('Ripple power around DA peaks - NREM (smoothed)');
grid on;
hold off









