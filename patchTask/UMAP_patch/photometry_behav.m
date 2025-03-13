%% Analyze the behavioral data with the photometry data
% 
% Takes in photometry data that has already been synced to intan. Plots
% the times of rewarded and non-rewarded licks on top of the photometry
% data for the entire behavior session.

basepath = pwd;
behaviorFile = dir(fullfile(basepath, '*.TrialBehavior.mat'));
%{
mode = 1; % 0 for downsampled, 1 for regular
%this allows me to easily siwtch btwn down sampled and regular
% N7 session 12 looks good downsamples but not regular
if mode == 0
    sampling_rate = 30;
    trackFile = dir(fullfile(basepath, '*TrackingDS.mat'));
    photometry_file = dir(fullfile(basepath, '*PhotometryBehavDS.mat'));
    load(photometry_file.name);
    load(trackFile.name);
    photom_var = photometry_ds;
    track_var = tracking_ds;
else
    sampling_rate = 130; % sampling rate of photometry set up - 130
    trackFile = dir(fullfile(basepath, '*Tracking.Behavior.mat'));
    photometry_file = dir(fullfile(basepath, '*PhotometryBehav.mat'));
    load(photometry_file.name);
    load(trackFile.name);
    photom_var = photometry;
    track_var = tracking.position;
end

load(behaviorFile.name);

% used to be patch_behav

rewarded_times = [];
nonrewarded_times = [];
% divide rewarded times vs nonrewarded times
% rewarded_times  and nonrewarded_times has first column with time, second column with trial
% number/idx
for i = 1:length(behavTrials.timestamps) 
    if behavTrials.reward_outcome(i) == 0
        nonrewarded_times = [nonrewarded_times; behavTrials.timestamps(i), i];
    else
        rewarded_times = [rewarded_times; behavTrials.timestamps(i), i];
    end
end


%% Plot photometry with behavior
%{
darkerGreen = [0.1098    0.6000    0.2392];
hold on
figure
plot(photom_var.timestamps, photom_var.grabDA_z, 'Color', darkerGreen, 'LineWidth', 2); 
ylabel('z score');
xlabel('time');
xlim([photom_var.timestamps(1), photom_var.timestamps(length(photom_var.timestamps(~isnan(photom_var.timestamps))))]);
for j = 1:length(rewarded_times)
    xline(rewarded_times(j, 1), '-b', 'LineWidth', 2);
end
for k = 1:length(nonrewarded_times)
    xline(nonrewarded_times(k, 1), '-r', 'LineWidth', 2);
end
hold off

%}



%% Average photometry across licks

window = 3.1; % window of time around lick to average
samples = window*sampling_rate;

zscore_matrix = nan(length(rewarded_times), (samples*2)+1); % rewarded
zscore_matrix_non = nan(length(nonrewarded_times), (samples*2)+1); % not rewarded


% average photometry data within a specified time window around rewarded
% licks
for j = 1:length(rewarded_times)
    curr_reward_time = rewarded_times(j, 1);
    [~, reward_idx] = min(abs(photom_var.timestamps - rewarded_times(j, 1)));
    start_idx = reward_idx - samples;
    end_idx = reward_idx + samples;
    if start_idx >= 1 && end_idx <= length(photom_var.timestamps)
        zscore_matrix(j, :) = photom_var.grabDA_z(start_idx:end_idx);

    end
end

zscore_matrix(any(isnan(zscore_matrix), 2), :) = [];
med_z_reward = median(zscore_matrix, 1); % median at each timepoint
time = linspace(-window, window, ((samples*2)+1));

% calculate confidence intervals
N = height(zscore_matrix);                          % Number of eExperimentsn In Data Set
avg_z_reward = mean(zscore_matrix, 1);              % Mean Of All Experiments At Each Value Of ,x 
reward_SEM = std(zscore_matrix, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
reward_CI95 = bsxfun(@times, reward_SEM, CI95(:));  % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of exw


% average photometry data within a specified time window around nonrewarded
% licks
for k = 1:length(nonrewarded_times)
    curr_reward_time = nonrewarded_times(k, 1);
    [~, reward_idx] = min(abs(photom_var.timestamps - nonrewarded_times(k, 1)));
    start_idx = reward_idx - samples;
    end_idx = reward_idx + samples;
    if start_idx >= 1 && end_idx <= length(photom_var.timestamps)
        zscore_matrix_non(k, :) = photom_var.grabDA_z(start_idx:end_idx);

    end
end

zscore_matrix_non(any(isnan(zscore_matrix_non), 2), :) = [];
med_z_no_reward = median(zscore_matrix_non, 1);

% calculate confidence intervals
N = height(zscore_matrix_non);                          % Number of experimentsn In Data Set
avg_z_no_reward = mean(zscore_matrix_non, 1);              % Mean Of All Experiments At Each Value Of ,x 
nonreward_SEM = std(zscore_matrix_non, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
nonreward_CI95 = bsxfun(@times, nonreward_SEM, CI95(:));  % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of exw


%% Plot  average photometry level around licks
figure(2);
subplot(2,1,1)
hold on
plot(time, med_z_reward, 'g', 'LineWidth', 2);
fill([time,fliplr(time)], [(reward_CI95(1,:)+avg_z_reward),fliplr((reward_CI95(2,:)+avg_z_reward))], 'b', 'EdgeColor','none', 'FaceAlpha',0.25)
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score Around Rewards');
grid on;
hold off

subplot(2,1,2)
plot(time, med_z_no_reward, 'g', 'LineWidth', 2);
hold on
fill([time,fliplr(time)], [(nonreward_CI95(1,:)+avg_z_no_reward),fliplr((nonreward_CI95(2,:)+avg_z_no_reward))], 'b', 'EdgeColor','none', 'FaceAlpha',0.25)
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score Around Non-rewarded Licks');
grid on;
hold off

% plot together 
cusGreen = [0.7176    0.9412    0.1020];
cusBlue = [ 0.2392    0.2863    0.9608];
figure;
hold on
plot(time, med_z_reward, 'g', 'LineWidth', 2);
fill([time,fliplr(time)], [(reward_CI95(1,:)+avg_z_reward),fliplr((reward_CI95(2,:)+avg_z_reward))], cusGreen, 'EdgeColor','none', 'FaceAlpha',0.25)
grid on;

plot(time, med_z_no_reward, 'Color', cusBlue, 'LineWidth', 2);
fill([time,fliplr(time)], [(nonreward_CI95(1,:)+avg_z_no_reward),fliplr((nonreward_CI95(2,:)+avg_z_no_reward))], 'b', 'EdgeColor','none', 'FaceAlpha',0.25)
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score Around Licks');
hold off


%% Plot dopamine on track

xq = linspace(min(track_var.x), max(track_var.x), 100); 
yq = linspace(min(track_var.y), max(track_var.y), 100); 
[Xq, Yq] = meshgrid(xq, yq); 

% Interpolate dopamine values onto grid 
Dq = griddata(track_var.x, track_var.y, photom_var.grabDA_z, Xq, Yq, 'cubic'); % use 'cubic' or 'linear' interpolation 

% p values
grabDA_p = arrayfun(@(x) 2*(1-normcdf(abs(x))), photom_var.grabDA_z);
Pq = griddata(track_var.x, track_var.y, grabDA_p, Xq, Yq, 'cubic');

threshold = 0.1;
Pq(Pq > threshold) = 1;
Pq(Pq <= threshold & ~isnan(Pq)) = -1;

%c = uisetcolor

cmap = [0.980392156862745, 0.360784313725490, 0.137254901960784;
    0.709803921568627, 0.909803921568627, 0.960784313725490];
    %0.8, 0.8, 0.8];

figure;
subplot(1, 2, 1)
imagesc([min(track_var.x), max(track_var.x)], [min(track_var.y), max(track_var.y)], Dq); 
set(gca, 'YDir', 'normal'); 
colormap(subplot(1, 2, 1), "parula")
colorbar; 
xlabel('X position'); 
ylabel('Y position'); 
title('grabDA z score');

subplot(1, 2, 2)
imagesc([min(track_var.x), max(track_var.x)], [min(track_var.y), max(track_var.y)], Pq); 
set(gca, 'YDir', 'normal'); 
colormap(subplot(1, 2, 2), cmap);
colorbar('Ticks', [-1, 1], 'TickLabels', {'p <= 0.1', 'p > 0.1'}); 
xlabel('X position'); 
ylabel('Y position'); 
title('p value');

figure;
plot(track_var.speed, photom_var.grabDA_z);
xlabel('Speed'); 
ylabel('grabDA z score'); 


%% Plot dopamine at unexpected reward (reward when mouse is in low patch)

patch_lick = []; % what patch did they lick in
for i = 1:behavTrials.num_trials
    if behavTrials.port(i) == 1 
        patch_lick(i) = 0;
    elseif behavTrials.port(i) == 2
        patch_lick(i) = 0;
    elseif behavTrials.port(i) == 3
        patch_lick(i) = 0;
    elseif behavTrials.port(i) == 4
        patch_lick(i) = 2;
    else
        patch_lick(i) = 1;
    end
end
patch_lick = patch_lick';


lick_low = []; % trial when the mouse licked in the lower probability patch
for i = 1:behavTrials.num_trials
    if behavTrials.patch_number(i) ~= patch_lick(i)
        lick_low(i) = 1; % TRUE - mouse licked in low prob patch
    else
        lick_low(i) = 0; % FALSE - mouse licked in high prob patch
    end
end
lick_low = lick_low';

unexpected_times = [];
for i = 1:behavTrials.num_trials 
    if and((behavTrials.reward_outcome(i) == 1), (lick_low(i) == 1))
        unexpected_times = [unexpected_times; behavTrials.timestamps(i)];
    else
        continue;
    end
end

window = 3.1; % window of time around lick to average
samples = window*sampling_rate;

zscore_matrix_unexpected = nan(length(unexpected_times), (samples*2)+1); % rewarded

% average photometry data within a specified time window around rewarded
% licks
for j = 1:length(unexpected_times)
    curr_reward_time = unexpected_times(j);
    [~, reward_idx] = min(abs(photom_var.timestamps - unexpected_times(j)));
    start_idx = reward_idx - samples;
    end_idx = reward_idx + samples;
    if start_idx >= 1 && end_idx <= length(photom_var.timestamps)
        zscore_matrix_unexpected(j, :) = photom_var.grabDA_z(start_idx:end_idx);

    end
end

zscore_matrix_unexpected(any(isnan(zscore_matrix_unexpected), 2), :) = [];
med_z_unexpected = median(zscore_matrix_unexpected, 1); % median at each timepoint

time = linspace(-window, window, ((samples*2)+1));

% calculate confidence intervals
N = height(zscore_matrix_unexpected);                          % Number of eExperimentsn In Data Set
avg_z_unexpected = mean(zscore_matrix_unexpected, 1);              % Mean Of All Experiments At Each Value Of ,x 
reward_SEM = std(zscore_matrix_unexpected, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
unexpected_CI95 = bsxfun(@times, reward_SEM, CI95(:));  % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of exw

%% Plot  average photometry level around unexpected rewarded licks
figure;
hold on
plot(time, med_z_unexpected, 'g', 'LineWidth', 2);
fill([time,fliplr(time)], [(unexpected_CI95(1,:)+avg_z_unexpected),fliplr((unexpected_CI95(2,:)+avg_z_unexpected))], 'b', 'EdgeColor','none', 'FaceAlpha',0.25)
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score Around Unexpected Rewards');
grid on;
hold off


% find average duration of trial

durs = zeros(behavTrials.num_trials, 1);
for m = 2:behavTrials.num_trials
    durs(m) = behavTrials.timestamps(m) - behavTrials.timestamps(m-1);
end
avg_dur = median(durs(2:end));

durs_reward = zeros(length(rewarded_times), 1);
for m = 2:length(rewarded_times)
    idx = rewarded_times(m, 2);
    durs_reward(m) = behavTrials.timestamps(idx) - behavTrials.timestamps(idx-1);
end
avg_dur_reward = median(durs_reward(2:end));

durs_nonreward = zeros(length(nonrewarded_times), 1);
for m = 2:length(nonrewarded_times)
    idx = nonrewarded_times(m, 2);
    durs_nonreward(m) = behavTrials.timestamps(idx) - behavTrials.timestamps(idx-1);
end
avg_dur_nonreward = median(durs_nonreward(2:end));


x=1;
%}


%% Analyze sleep photometry

% synchronize time stamps for sleep session
%sleep_sync = getSyncPhotometry(photometryData);

% if saveMat
%     C = strsplit(pwd,'\');
%     save([basepath filesep C{end} '.SleepPhotomSynced.mat'],'sleep_sync');
% end

% make this a fxn
%function [avg] = avgAcrossEvents(event)

window = 1; % window of time around lick to average
samples = window*sampling_rate;

% need to find end timestamp of first sleep session. will write code for
% this, but for right now can just use MergePoints

sleep_start = 8335; % start of sleep session
sleep_end = 8445; % time stamp of end of sleep session 
ripple_period = ripples.peaks(ripples.peaks <= sleep_end & ripples.peaks >= sleep_start);
ripple_matrix = nan(length(ripple_period), (samples*2)+1); % ripples

% average photometry data within a specified time window around ripples
adjusted_ts = sleep_sync.timestamps + sleep_start;
for j = 1:length(ripple_period)
    [~, ripple_idx] = min(abs(adjusted_ts - ripple_period(j)));
    start_idx = ripple_idx - samples;
    end_idx = ripple_idx + samples;
    if start_idx >= 1 && end_idx <= height(adjusted_ts)
        ripple_matrix(j, :) = sleep_sync.grabDA_z(start_idx:end_idx);
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


%% Plot photometry against ripples

figure;
plot(time, median_ripple, 'g', 'LineWidth', 2);
hold on
fill([time,fliplr(time)], [(ripple_CI95(1,:)+median_ripple),fliplr((ripple_CI95(2,:)+median_ripple))], 'b', 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score Around Ripples');
grid on;
hold off

smoothed = smoothdata(median_ripple);
plot(time, smoothed, 'g', 'LineWidth', 2);
fill([time,fliplr(time)], [(smooth_CI95(1,:)+smoothed),fliplr((smooth_CI95(2,:)+smoothed))], 'b', 'EdgeColor','none', 'FaceAlpha',0.25)



figure;
plot(sleep_sync(:,1), sleep_sync(:,2), 'g', 'LineWidth', 2);


ripple_durs = zeros(length(ripple_period), 1);
for k = 2:length(ripple_period)
    ripple_durs(k) = ripple_period(k)-ripple_period(k-1);
end



%photom_pre_sleep = photometryData;
% photometry across sleep and behav
figure(4);
subplot(2,1,1)
plot(photom_pre_sleep.timestamps, photom_pre_sleep.grabDA_z)
%xlim([0, photom_pre_sleep.timestamps(size(photom_pre_sleep.timestamps))]);
ylabel('z score');
title('sleep DA');

subplot(2,1,2)
plot(photometryData.timestamps, photometryData.grabDA_z)
%xlim([0, photometryData.timestamps(size(photometryData.timestamps))]);
ylabel('z score');
title('behavior DA');








