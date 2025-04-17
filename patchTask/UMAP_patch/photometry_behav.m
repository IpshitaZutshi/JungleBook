%% Analyze the behavioral data with the photometry data
% 
% Takes in photometry data that has already been synced to intan. Plots
% the times of rewarded and non-rewarded licks on top of the photometry
% data for the entire behavior session.

%% Defaults and Parms

basepath = pwd;
behaviorFile = dir(fullfile(basepath, '*.TrialBehavior.mat'));
session = 18;
<<<<<<< HEAD
=======
%color_scheme = spring;
>>>>>>> 2f4d3a86f04df0da4d5e289cb15f173d8d44909b

color = 1; % 0/blue for striatum, 1/green for HPC

mode = 1; % 0 for downsampled, 1 for regular
%this allows me to easily siwtch btwn down sampled and regular
% N7 session 12 looks good downsamples but not regular
if mode == 0
    sampling_rate = 7; 
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


%% SEPARATIONS
% separate rewarded vs non rewarded
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

% separate high patch vs low patch
high_patch = []; % trials when mouse licks in high patch
low_patch = [];
% high_patch and low_patch has first column with time, second column with trial number/idx
for i = 1:behavTrials.num_trials 
    if behavTrials.patch_number(i) == 0
        if behavTrials.port(i) == 1 || behavTrials.port(i) == 2 || behavTrials.port(i) == 3
            high_patch = [high_patch; behavTrials.timestamps(i), i];
        else
            low_patch = [low_patch; behavTrials.timestamps(i), i];
        end
    elseif behavTrials.patch_number(i) == 1
        if behavTrials.port(i) == 5 || behavTrials.port(i) == 6 || behavTrials.port(i) == 7
            high_patch = [high_patch; behavTrials.timestamps(i), i];
        else
            low_patch = [low_patch; behavTrials.timestamps(i), i];
        end
    end
end

% separate stay vs switch
stay = [];
sswitch = [];
for i = 1:behavTrials.num_trials 
    if behavTrials.stay_switch(i) == 0
        stay = [stay; behavTrials.timestamps(i), i];
    else
        sswitch = [sswitch; behavTrials.timestamps(i), i];
    end
end

%{
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
lick_high = [];
for i = 1:behavTrials.num_trials
    if behavTrials.patch_number(i) ~= patch_lick(i)
        lick_low(i) = 1; % TRUE - mouse licked in low prob patch
        lick_high(i) = 0;
    else
        lick_low(i) = 0; % FALSE - mouse licked in high prob patch
        lick_high(i) = 1;
    end
end
lick_low = lick_low';
lick_high = lick_high';
%}

%% Plot photometry with behavior
%{
darkerGreen = [0.1098    0.6000    0.2392];
hold on
figure('color','white')
ax = gca;
ax.FontSize = 15;
plot(photom_var.timestamps, photom_var.grabDA_z, 'Color', avg_color, 'LineWidth', 2); 
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

window = 5; % window of time around lick to average
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

% statistical significance of peak before lick
reward_baseline = zscore_matrix(:,[1:390]);%135
avg_reward_baseline = mean(reward_baseline, 1);
mn = mean(avg_reward_baseline);
st_d = std(avg_reward_baseline);
sample_mn = max(med_z_reward);
deg_free = length(rewarded_times)-1;
t = (sample_mn - mn)/(st_d/(sqrt(deg_free)));

if color == 0
    % avg_color = [ 0.2392    0.2863    0.9608];
    % conf_color = 'b';
    avg_color = [0.960784313725490, 0.152941176470588, 0.905882352941176];
    conf_color = [0.960784313725490, 0.152941176470588, 0.905882352941176];
else
    avg_color = 'g';
    conf_color = [0.7176    0.9412    0.1020];
end

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


% statistical significance of peak before lick
non_reward_baseline = zscore_matrix_non(:,[1:390]);%135
avg_nonreward_baseline = mean(non_reward_baseline, 1);
mn = mean(avg_nonreward_baseline);
st_d = std(avg_nonreward_baseline);
sample_mn = max(med_z_no_reward);
deg_free = length(nonrewarded_times)-1;
t = (sample_mn - mn)/(st_d/(sqrt(deg_free)));

%% Plot  average photometry level around licks

<<<<<<< HEAD
total_min = min(min(med_z_reward), min(med_z_no_reward));
total_max = max(max(med_z_reward), max(med_z_no_reward));
=======
if color == 0
    % avg_color = [ 0.2392    0.2863    0.9608];
    % conf_color = 'b';
    avg_color = [0.960784313725490, 0.152941176470588, 0.905882352941176];
    conf_color = [0.960784313725490, 0.152941176470588, 0.905882352941176];
else
    avg_color = 'g';
    conf_color = [0.7176    0.9412    0.1020];
end
>>>>>>> 2f4d3a86f04df0da4d5e289cb15f173d8d44909b

figure('color','white');
subplot(2,1,1)
hold on
<<<<<<< HEAD
plot(time, med_z_reward, 'g', 'LineWidth', 2);
fill([time,fliplr(time)], [(reward_CI95(1,:)+med_z_reward),fliplr((reward_CI95(2,:)+med_z_reward))], 'b', 'EdgeColor','none', 'FaceAlpha',0.25)
=======
ax = gca;
ax.FontSize = 15;
plot(time, avg_z_reward, 'color', avg_color, 'LineWidth', 2);
fill([time,fliplr(time)], [(reward_CI95(1,:)+avg_z_reward),fliplr((reward_CI95(2,:)+avg_z_reward))], conf_color, 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
>>>>>>> 2f4d3a86f04df0da4d5e289cb15f173d8d44909b
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score Around Rewards');
ylim([total_min-0.1, total_max+0.1]);
grid on;
hold off

subplot(2,1,2)
plot(time, avg_z_no_reward, 'color', avg_color, 'LineWidth', 2);
hold on
<<<<<<< HEAD
fill([time,fliplr(time)], [(nonreward_CI95(1,:)+med_z_no_reward),fliplr((nonreward_CI95(2,:)+med_z_no_reward))], 'b', 'EdgeColor','none', 'FaceAlpha',0.25)
=======
ax = gca;
ax.FontSize = 15;
fill([time,fliplr(time)], [(nonreward_CI95(1,:)+avg_z_no_reward),fliplr((nonreward_CI95(2,:)+avg_z_no_reward))], conf_color, 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
>>>>>>> 2f4d3a86f04df0da4d5e289cb15f173d8d44909b
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score Around Non-rewarded Licks');
ylim([total_min-0.1, total_max+0.1]);
grid on;
hold off
%{
% plot together 
cusGreen = [0.7176    0.9412    0.1020];
cusBlue = [ 0.2392    0.2863    0.9608];
figure('color','white');
hold on
plot(time, avg_z_reward, 'g', 'LineWidth', 2);
fill([time,fliplr(time)], [(reward_CI95(1,:)+avg_z_reward),fliplr((reward_CI95(2,:)+avg_z_reward))], cusGreen, 'EdgeColor','none', 'FaceAlpha',0.25)
grid on;

plot(time, avg_z_no_reward, 'Color', cusBlue, 'LineWidth', 2);
fill([time,fliplr(time)], [(nonreward_CI95(1,:)+avg_z_no_reward),fliplr((nonreward_CI95(2,:)+avg_z_no_reward))], 'b', 'EdgeColor','none', 'FaceAlpha',0.25)
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score Around Licks');
hold off
%}


%% DA - licks in high patch vs low patch

window = 5; % window of time around lick to average
samples = window*sampling_rate;

zscore_matrix_high = nan(length(high_patch), (samples*2)+1); % rewarded
zscore_matrix_low = nan(length(low_patch), (samples*2)+1); % not rewarded


% average photometry data within a specified time window around rewarded
% licks
for j = 1:length(high_patch)
    if behavTrials.patch_trials(high_patch(j, 2)) >= 30
        curr_high_time = high_patch(j, 1);
        [~, high_idx] = min(abs(photom_var.timestamps - high_patch(j, 1)));
        start_idx = high_idx - samples;
        end_idx = high_idx + samples;
        if start_idx >= 1 && end_idx <= length(photom_var.timestamps)
            zscore_matrix_high(j, :) = photom_var.grabDA_z(start_idx:end_idx);
            baseline = zscore_matrix_high(j, :);
            % Calculate baseline mean and std
            baseline_mean = mean(baseline);
            baseline_std = std(baseline);
           
            % Z-score the entire trial using the baseline stats
            zscore_matrix_high(j, :) = (baseline - baseline_mean) / baseline_std;
        end
    else
        continue
    end
end

zscore_matrix_high(any(isnan(zscore_matrix_high), 2), :) = [];
med_z_high = median(zscore_matrix_high, 1); % median at each timepoint
time = linspace(-window, window, ((samples*2)+1));

% calculate confidence intervals
N = height(zscore_matrix_high);                          % Number of eExperimentsn In Data Set
avg_z_high = mean(zscore_matrix_high, 1);              % Mean Of All Experiments At Each Value Of ,x 
high_SEM = std(zscore_matrix_high, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
high_CI95 = bsxfun(@times, high_SEM, CI95(:));  % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of exw

% statistical significance of peak before lick
% reward_baseline = zscore_matrix_high(:,[1:390]);%135
% avg_reward_baseline = mean(reward_baseline, 1);
% mn = mean(avg_reward_baseline);
% st_d = std(avg_reward_baseline);
% sample_mn = max(med_z_reward);
% deg_free = length(rewarded_times)-1;
% t = (sample_mn - mn)/(st_d/(sqrt(deg_free)));

% average photometry data within a specified time window around nonrewarded
% licks
for k = 1:length(low_patch)
    if behavTrials.patch_trials(low_patch(k, 2)) >= 30
        curr_low_time = low_patch(k, 1);
        [~, low_idx] = min(abs(photom_var.timestamps - low_patch(k, 1)));
        start_idx = low_idx - samples;
        end_idx = low_idx + samples;
        if start_idx >= 1 && end_idx <= length(photom_var.timestamps)
            zscore_matrix_low(k, :) = photom_var.grabDA_z(start_idx:end_idx);
            baseline = zscore_matrix_low(k, :);
            % Calculate baseline mean and std
            baseline_mean = mean(baseline);
            baseline_std = std(baseline);
           
            % Z-score the entire trial using the baseline stats
            zscore_matrix_low(k, :) = (baseline - baseline_mean) / baseline_std;
        end
    else
        continue
    end
end

zscore_matrix_low(any(isnan(zscore_matrix_low), 2), :) = [];
med_z_low = median(zscore_matrix_low, 1);

% calculate confidence intervals
N = height(zscore_matrix_low);                          % Number of experimentsn In Data Set
avg_z_low = mean(zscore_matrix_low, 1);              % Mean Of All Experiments At Each Value Of ,x 
low_SEM = std(zscore_matrix_low, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
low_CI95 = bsxfun(@times, low_SEM, CI95(:));  % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of exw


% statistical significance of peak before lick
% non_reward_baseline = zscore_matrix_low(:,[1:390]);%135
% avg_nonreward_baseline = mean(non_reward_baseline, 1);
% mn = mean(avg_nonreward_baseline);
% st_d = std(avg_nonreward_baseline);
% sample_mn = max(med_z_no_reward);
% deg_free = length(nonrewarded_times)-1;
% t = (sample_mn - mn)/(st_d/(sqrt(deg_free)));

%% Plot photometry level around licks - high vs low

figure('color','white');
subplot(2,1,1)
ax = gca;
ax.FontSize = 15;
hold on
plot(time, med_z_high, 'color', avg_color, 'LineWidth', 2);
fill([time,fliplr(time)], [(high_CI95(1,:)+med_z_high),fliplr((high_CI95(2,:)+med_z_high))], conf_color, 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score in High Patch');
grid on;
hold off

subplot(2,1,2)
plot(time, med_z_low, 'color', avg_color, 'LineWidth', 2);
hold on
ax = gca;
ax.FontSize = 15;
fill([time,fliplr(time)], [(low_CI95(1,:)+med_z_low),fliplr((low_CI95(2,:)+med_z_low))], conf_color, 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score in Low Patch');
grid on;
hold off


%% DA - early vs late in patch

save_path = 'Z:\Buzsakilabspace\LabShare\ZutshiI\patchTask\N11';
saveLoc = strcat(save_path,'\DAChangeDuringPatch');
% if ~isfolder('Maps')
%     mkdir('Maps')
% end    


% Parameters
window = 5; % Time window around event in seconds
samples = window * sampling_rate;
baseline_samples = window * sampling_rate * (4/5);

% Initialize variables
trial_count = 0;
patch_num = 0;

% Loop through all trials
for i = 1:behavTrials.num_trials
   
    % Identify patch switches
    if behavTrials.patch_trials(i) == 1
        % If patch number is not zero, plot the previous patch
        if patch_num > 0
            % Plot the average dopamine signal in increments of 10
            figure('color','white', 'WindowState', 'maximized');
            ax = gca;
            ax.FontSize = 15;
            hold on
            quartile_size = ceil(trial_count / 4);
            cmap = spring(4); 
            for q = 1:4
                % Determine the trial range for each quartile
                start_q = (q - 1) * quartile_size + 1;
                end_q = min(q * quartile_size, trial_count);
               
                % Average the signal within the quartile
                avg_signal = mean(patch_zscore_matrix(start_q:end_q, :), 1);
                smoothed = smoothdata(avg_signal);
                baseline_mean = mean(smoothed);
                baseline_std = std(smoothed);
                normalized_signal = (smoothed - baseline_mean) / baseline_std;
               
                % Plot each quartile with different color
                plot(time, normalized_signal, 'LineWidth', 2, 'Color', cmap(q, :));
            end
           
            % Labeling
            xlabel('Time (s)');
            ylabel('Avg Z-score');
            title(['Patch ', num2str(patch_num), ': Dopamine Signal by Quartiles']);
            grid on;
           
            % Add colorbar with proper trial range labeling
            colormap(spring)
            colorbar_handle = colorbar;
            colorbar_handle.Ticks = linspace(0, 1, 4);
            quartile_labels = arrayfun(@(q) sprintf('Q%d (%d-%d)', q, ...
                (q - 1) * quartile_size + 1, min(q * quartile_size, trial_count)), 1:4, 'UniformOutput', false);
            colorbar_handle.TickLabels = quartile_labels;
            colorbar_handle.Label.String = 'Quartile Range';
           
            hold off
            %saveas(gcf,[saveLoc,filesep ,'Sess_', num2str(session), '_Patch_', num2str(patch_num),'.png'],'png');

        end
       
        % Reset for the new patch
        patch_num = patch_num + 1;
        trial_count = 0;
        patch_zscore_matrix = [];
    end

    % Get the timestamp for the current trial
    curr_time = behavTrials.timestamps(i);
    [~, trial_idx] = min(abs(photom_var.timestamps - curr_time));
   
    % Extract window of data
    start_idx = trial_idx - samples;
    end_idx = trial_idx + samples;
   
    if start_idx >= 1 && end_idx <= length(photom_var.timestamps)
        trial_count = trial_count + 1;
        patch_zscore_matrix(trial_count, :) = photom_var.grabDA_z(start_idx:end_idx);
    end
end

% Plot the final patch
if patch_num > 0 && trial_count > 0
    figure('color','white', 'WindowState', 'maximized');
    ax = gca;
    ax.FontSize = 15;
    hold on
    quartile_size = ceil(trial_count / 4);
    cmap = spring(4); 

    for q = 1:4
        start_q = (q - 1) * quartile_size + 1;
        end_q = min(q * quartile_size, trial_count);
       
        % Average signal within the quartile
        avg_signal = mean(patch_zscore_matrix(start_q:end_q, :), 1);
        smoothed = smoothdata(avg_signal);
        baseline_mean = mean(smoothed);
        baseline_std = std(smoothed);
        normalized_signal = (smoothed - baseline_mean) / baseline_std;
       
        % Plot each quartile with different color
        plot(time, normalized_signal, 'LineWidth', 2, 'Color', cmap(q, :));
    end
   
    % Labeling
    xlabel('Time (s)');
    ylabel('Avg Z-score');
    title(['Patch ', num2str(patch_num), ': Dopamine Signal by Quartiles']);
    grid on;
   
    % Add colorbar with proper trial range labeling
    colormap(spring)
    colorbar_handle = colorbar;
    colorbar_handle.Ticks = linspace(0, 1, 4);
    quartile_labels = arrayfun(@(q) sprintf('Q%d (%d-%d)', q, ...
        (q - 1) * quartile_size + 1, min(q * quartile_size, trial_count)), 1:4, 'UniformOutput', false);
    colorbar_handle.TickLabels = quartile_labels;
    colorbar_handle.Label.String = 'Quartile Range';
   
    hold off
    %saveas(gcf,[saveLoc,filesep ,'Sess_', num2str(session), '_Patch_', num2str(patch_num),'.png'],'png');
end

<<<<<<< HEAD
=======

>>>>>>> 2f4d3a86f04df0da4d5e289cb15f173d8d44909b
%% Stay vs switch

window = 5; % window of time around lick to average
samples = window*sampling_rate;

zscore_matrix_stay = nan(length(high_patch), (samples*2)+1); % rewarded
zscore_matrix_switch = nan(length(low_patch), (samples*2)+1); % not rewarded


% average photometry data within a specified time window around stay
for j = 1:length(stay)
    [~, stay_idx] = min(abs(photom_var.timestamps - stay(j, 1)));
    start_idx = stay_idx - samples;
    end_idx = stay_idx + samples;
    if start_idx >= 1 && end_idx <= length(photom_var.timestamps)
        zscore_matrix_stay(j, :) = photom_var.grabDA_z(start_idx:end_idx);
        baseline = zscore_matrix_stay(j, :);
        baseline_mean = mean(baseline);
        baseline_std = std(baseline);
       
        zscore_matrix_stay(j, :) = (baseline - baseline_mean) / baseline_std;
    end
end

zscore_matrix_stay(any(isnan(zscore_matrix_stay), 2), :) = [];
med_z_stay = median(zscore_matrix_stay, 1); % median at each timepoint
time = linspace(-window, window, ((samples*2)+1));

% calculate confidence intervals
N = height(zscore_matrix_stay);                          % Number of eExperimentsn In Data Set
avg_z_stay = mean(zscore_matrix_stay, 1);              % Mean Of All Experiments At Each Value Of ,x 
SEM = std(zscore_matrix_stay, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
stay_CI95 = bsxfun(@times, SEM, CI95(:));  % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of exw

% average photometry data within a specified time window around nonrewarded
% licks
for k = 1:length(sswitch)
    [~, switch_idx] = min(abs(photom_var.timestamps - sswitch(k, 1)));
    start_idx = switch_idx - samples;
    end_idx = switch_idx + samples;
    if start_idx >= 1 && end_idx <= length(photom_var.timestamps)
        zscore_matrix_switch(k, :) = photom_var.grabDA_z(start_idx:end_idx);
        baseline = zscore_matrix_switch(k, :);
        baseline_mean = mean(baseline);
        baseline_std = std(baseline);
       
        zscore_matrix_switch(k, :) = (baseline - baseline_mean) / baseline_std;
    end
end

zscore_matrix_switch(any(isnan(zscore_matrix_switch), 2), :) = [];
med_z_switch = median(zscore_matrix_switch, 1);

% calculate confidence intervals
N = height(zscore_matrix_switch);                          % Number of experimentsn In Data Set
avg_z_switch = mean(zscore_matrix_switch, 1);              % Mean Of All Experiments At Each Value Of ,x 
SEM = std(zscore_matrix_switch, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
switch_CI95 = bsxfun(@times, SEM, CI95(:));  % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of exw


%% Plot photometry during run stay vs switch

figure('color','white');
subplot(2,1,1)
ax = gca;
ax.FontSize = 15;
hold on
plot(time, med_z_stay, 'color', avg_color, 'LineWidth', 2);
fill([time,fliplr(time)], [(stay_CI95(1,:)+med_z_stay),fliplr((stay_CI95(2,:)+med_z_stay))], conf_color, 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score when mouse stays in patch');
grid on;
hold off

subplot(2,1,2)
plot(time, med_z_switch, 'color', avg_color, 'LineWidth', 2);
hold on
ax = gca;
ax.FontSize = 15;
fill([time,fliplr(time)], [(switch_CI95(1,:)+med_z_switch),fliplr((switch_CI95(2,:)+med_z_switch))], conf_color, 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score when mouse switches patch');
grid on;
hold off

peaks_stay = max(zscore_matrix_stay, [], 2);
peaks_switch = max(zscore_matrix_switch, [], 2);

[h, p, ci, stats] = ttest2(peaks_stay, peaks_switch);

<<<<<<< HEAD
=======

>>>>>>> 2f4d3a86f04df0da4d5e289cb15f173d8d44909b
%% Plot dopamine at unexpected reward (reward when mouse is in low patch)

unexpected_reward = []; % when the mouse is in low patch but gets rewarded
unexpected_non = []; % when the mouse is in high patch and doesn't get rewarded
expected_reward = []; 
expected_non = []; 
for i = 1:length(low_patch) 
    if behavTrials.reward_outcome(low_patch(i, 2)) == 1
        unexpected_reward = [unexpected_reward; low_patch(i, 1), low_patch(i, 2)];
    else
        expected_non = [expected_non; low_patch(i, 1), low_patch(i, 2)];
    end
end

for i = 1:length(high_patch) 
    if behavTrials.reward_outcome(high_patch(i, 2)) == 0
        unexpected_non = [unexpected_non; high_patch(i, 1), high_patch(i, 2)];
    else
        expected_reward = [expected_reward; high_patch(i, 1), high_patch(i, 2)];
    end
end    

window = 3.1; % window of time around lick to average
samples = window*sampling_rate;

zscore_matrix_unexpected = nan(length(unexpected_reward), (samples*2)+1); % rewarded
zscore_matrix_unexpected_non = nan(length(unexpected_non), (samples*2)+1); % unrewarded

% average photometry data within a specified time window around unexpected rewarded
% licks
for j = 1:length(unexpected_reward)
    curr_reward_time = unexpected_reward(j, 1);
    [~, reward_idx] = min(abs(photom_var.timestamps - unexpected_reward(j, 1)));
    start_idx = reward_idx - samples;
    end_idx = reward_idx + samples;
    if start_idx >= 1 && end_idx <= length(photom_var.timestamps)
        zscore_matrix_unexpected(j, :) = photom_var.grabDA_z(start_idx:end_idx);
        baseline = zscore_matrix_unexpected(j, :);
        % Calculate baseline mean and std
        baseline_mean = mean(baseline);
        baseline_std = std(baseline);
       
        % Z-score the entire trial using the baseline stats
        zscore_matrix_unexpected(j, :) = (baseline - baseline_mean) / baseline_std;
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

% unexpected unrewarded
for j = 1:length(unexpected_non)
    curr_reward_time = unexpected_non(j, 1);
    [~, reward_idx] = min(abs(photom_var.timestamps - unexpected_non(j, 1)));
    start_idx = reward_idx - samples;
    end_idx = reward_idx + samples;
    if start_idx >= 1 && end_idx <= length(photom_var.timestamps)
        zscore_matrix_unexpected_non(j, :) = photom_var.grabDA_z(start_idx:end_idx);
        baseline = zscore_matrix_unexpected_non(j, :);
        % Calculate baseline mean and std
        baseline_mean = mean(baseline);
        baseline_std = std(baseline);
       
        % Z-score the entire trial using the baseline stats
        zscore_matrix_unexpected_non(j, :) = (baseline - baseline_mean) / baseline_std;
    end
end

zscore_matrix_unexpected_non(any(isnan(zscore_matrix_unexpected_non), 2), :) = [];
med_z_unexpected_non = median(zscore_matrix_unexpected_non, 1); % median at each timepoint

% calculate confidence intervals
N = height(zscore_matrix_unexpected_non);                          % Number of eExperimentsn In Data Set
avg_z_unexpected_non = mean(zscore_matrix_unexpected_non, 1);              % Mean Of All Experiments At Each Value Of ,x 
reward_SEM = std(zscore_matrix_unexpected_non, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
unexpected_CI95_non = bsxfun(@times, reward_SEM, CI95(:));  % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of exw


%% expected

zscore_expected = nan(length(expected_reward), (samples*2)+1); % rewarded
zscore_expected_non = nan(length(expected_non), (samples*2)+1); % unrewarded

% average photometry data within a specified time window around unexpected rewarded
% licks
for j = 1:length(expected_reward)
    curr_reward_time = expected_reward(j, 1);
    [~, reward_idx] = min(abs(photom_var.timestamps - expected_reward(j, 1)));
    start_idx = reward_idx - samples;
    end_idx = reward_idx + samples;
    if start_idx >= 1 && end_idx <= length(photom_var.timestamps)
        zscore_expected(j, :) = photom_var.grabDA_z(start_idx:end_idx);
        baseline = zscore_expected(j, :);
        % Calculate baseline mean and std
        baseline_mean = mean(baseline);
        baseline_std = std(baseline);
       
        % Z-score the entire trial using the baseline stats
        zscore_expected(j, :) = (baseline - baseline_mean) / baseline_std;
    end
end

zscore_expected(any(isnan(zscore_expected), 2), :) = [];
med_z_expected = median(zscore_expected, 1); % median at each timepoint

time = linspace(-window, window, ((samples*2)+1));

% calculate confidence intervals
N = height(zscore_expected);                          % Number of eExperimentsn In Data Set
avg_z_expected = mean(zscore_expected, 1);              % Mean Of All Experiments At Each Value Of ,x 
reward_SEM = std(zscore_expected, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
expected_CI95 = bsxfun(@times, reward_SEM, CI95(:));  % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of exw

% expected unrewarded
for j = 1:length(expected_non)
    curr_reward_time = expected_non(j, 1);
    [~, reward_idx] = min(abs(photom_var.timestamps - expected_non(j, 1)));
    start_idx = reward_idx - samples;
    end_idx = reward_idx + samples;
    if start_idx >= 1 && end_idx <= length(photom_var.timestamps)
        zscore_expected_non(j, :) = photom_var.grabDA_z(start_idx:end_idx);
        baseline = zscore_expected_non(j, :);
        % Calculate baseline mean and std
        baseline_mean = mean(baseline);
        baseline_std = std(baseline);
       
        % Z-score the entire trial using the baseline stats
        zscore_expected_non(j, :) = (baseline - baseline_mean) / baseline_std;
    end
end

zscore_expected_non(any(isnan(zscore_expected_non), 2), :) = [];
med_z_expected_non = median(zscore_expected_non, 1); % median at each timepoint

% calculate confidence intervals
N = height(zscore_expected_non);                          % Number of eExperimentsn In Data Set
avg_z_expected_non = mean(zscore_expected_non, 1);              % Mean Of All Experiments At Each Value Of ,x 
reward_SEM = std(zscore_expected_non, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
expected_CI95_non = bsxfun(@times, reward_SEM, CI95(:));  % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of exw


%% Plot average photometry level around unexpected rewarded licks
%{
figure('color','white');
subplot(2,1,1)
hold on
ax = gca;
ax.FontSize = 15;
plot(time, med_z_unexpected, 'color', avg_color, 'LineWidth', 2);
fill([time,fliplr(time)], [(unexpected_CI95(1,:)+med_z_unexpected),fliplr((unexpected_CI95(2,:)+med_z_unexpected))], conf_color, 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title(['Z-score Around Rewards in Low Patch'], [num2str(length(unexpected_reward)), ' events']);
grid on;
hold off

subplot(2,1,2)
plot(time, med_z_unexpected_non, 'color', avg_color, 'LineWidth', 2);
hold on
ax = gca;
ax.FontSize = 15;
fill([time,fliplr(time)], [(unexpected_CI95_non(1,:)+med_z_unexpected_non),fliplr((unexpected_CI95_non(2,:)+med_z_unexpected_non))], conf_color, 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title(['Z-score Around Non-rewarded Licks in High Patch'], [num2str(length(unexpected_non)), ' events']);
grid on;
hold off
%}

% all four together

figure('color','white');
subplot(2,2,1)
hold on
ax = gca;
ax.FontSize = 15;
plot(time, med_z_unexpected, 'color', avg_color, 'LineWidth', 2);
fill([time,fliplr(time)], [(unexpected_CI95(1,:)+med_z_unexpected),fliplr((unexpected_CI95(2,:)+med_z_unexpected))], conf_color, 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title(['Z-score Around Rewards in Low Patch'], [num2str(length(unexpected_reward)), ' events']);
grid on;
hold off

subplot(2,2,2)
plot(time, med_z_unexpected_non, 'color', avg_color, 'LineWidth', 2);
hold on
ax = gca;
ax.FontSize = 15;
fill([time,fliplr(time)], [(unexpected_CI95_non(1,:)+med_z_unexpected_non),fliplr((unexpected_CI95_non(2,:)+med_z_unexpected_non))], conf_color, 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title(['Z-score Around Non-rewarded Licks in High Patch'], [num2str(length(unexpected_non)), ' events']);
grid on;
hold off

subplot(2,2,3)
hold on
ax = gca;
ax.FontSize = 15;
plot(time, med_z_expected, 'color', avg_color, 'LineWidth', 2);
fill([time,fliplr(time)], [(expected_CI95(1,:)+med_z_expected),fliplr((expected_CI95(2,:)+med_z_expected))], conf_color, 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title(['Z-score Around Rewards in High Patch'], [num2str(length(expected_reward)), ' events']);
grid on;
hold off

subplot(2,2,4)
plot(time, med_z_expected_non, 'color', avg_color, 'LineWidth', 2);
hold on
ax = gca;
ax.FontSize = 15;
fill([time,fliplr(time)], [(expected_CI95_non(1,:)+med_z_expected_non),fliplr((expected_CI95_non(2,:)+med_z_expected_non))], conf_color, 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title(['Z-score Around Non-rewarded Licks in Low Patch'], [num2str(length(expected_non)), ' events']);
grid on;
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

figure('color','white');
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

figure('color','white');
plot(track_var.speed, photom_var.grabDA_z);
xlabel('Speed'); 
ylabel('grabDA z score'); 


%% find average duration of trial

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


%{
%% DA - early vs late in patch

<<<<<<< HEAD
%% Analyze awake ripples

% Load variables
basepath = pwd;
sampling_rate = 130; % sampling rate of photometry set up - 130

photometry_file = dir(fullfile(basepath, '*PhotometryBehav.mat'));
load(photometry_file.name);

ripple_file = dir(fullfile(basepath, '*ripples.events.mat'));
load(ripple_file.name)

merge_pt_file = dir(fullfile(basepath, '*MergePoints.events.mat'));
load(merge_pt_file.name)

sleep_start = MergePoints.timestamps(2,1); % start of sleep session
sleep_end = MergePoints.timestamps(2,2); % time stamp of end of sleep session 


ripple_period = ripples.peaks(ripples.peaks <= sleep_end & ripples.peaks >= sleep_start);
=======

% Parameters
window = 5; % Time window around event in seconds
samples = window * sampling_rate;
baseline_samples = window * sampling_rate * (4/5);

% Initialize variables
trial_count = 0;
patch_num = 0;

% Loop through all trials
for i = 1:length(behavTrials.patch_trials)
   
    % Identify patch switches
    if behavTrials.patch_trials(i) == 1
        % If patch number is not zero, plot the previous patch
        if patch_num > 0
            % Plot the average dopamine signal in increments of 10
            figure('color','white');
            hold on
            cmap = jet(ceil(trial_count/10)); 
            for b = 1:ceil(trial_count/10)
                avg_signal = mean(patch_zscore_matrix(((b-1)*10 + 1):min(b*10, trial_count), :), 1);
                smoothed = smoothdata(avg_signal);
>>>>>>> 2f4d3a86f04df0da4d5e289cb15f173d8d44909b

                baseline = mean(smoothed(1:baseline_samples));
                normalized_signal = smoothed - baseline;
                plot(time, normalized_signal, 'LineWidth', 2, 'Color', cmap(b, :));
            end
           
            % Labeling
            xlabel('Time (s)');
            ylabel('Avg Z-score');
            title(['Patch ', num2str(patch_num), ': Dopamine Signal in Increments of 10 Trials']);
            grid on;
           
            % Add colorbar with proper trial range labeling
            colormap(jet)
            colorbar_handle = colorbar;
            tick_positions = linspace(0, 1, ceil(trial_count/10));
            trial_ranges = arrayfun(@(x) sprintf('%d-%d', (x-1)*10 + 1, min(x*10, trial_count)), 1:ceil(trial_count/10), 'UniformOutput', false);
            colorbar_handle.Ticks = tick_positions;
            colorbar_handle.TickLabels = trial_ranges;
            colorbar_handle.Label.String = 'Trial Range';
           
            hold off
        end
       
        % Reset for the new patch
        patch_num = patch_num + 1;
        trial_count = 0;
        patch_zscore_matrix = [];
    end

<<<<<<< HEAD
window = 5; % window of time around ripple to average
samples = window*sampling_rate;

% initialize matrix
ripple_matrix = nan(length(ripple_period), (samples*2)+1); % ripples


% average photometry data within a specified time window around ripples
adjusted_ts = sleep_sync.timestamps + sleep_start;
for j = 1:length(ripple_period)
    [~, ripple_idx] = min(abs(adjusted_ts - ripple_period(j)));
    start_idx = ripple_idx - samples;
    end_idx = ripple_idx + samples;
    if start_idx >= 1 && end_idx <= height(adjusted_ts)
        ripple_matrix(j, :) = sleep_sync.grabDA_z(start_idx:end_idx);
        
        % z score
        baseline = ripple_matrix(j, :);

        % Calculate baseline mean and std
        baseline_mean = mean(baseline);
        baseline_std = std(baseline);

        % Z-score the entire trial using the baseline stats
        ripple_matrix(j, :) = (ripple_matrix(j, :) - baseline_mean) / baseline_std;
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

% statistical significance of peak before lick
ripple_baseline = ripple_matrix(:,[1:120]);
avg_ripple_baseline = mean(ripple_baseline, 1);
mn = mean(avg_ripple_baseline);
st_d = std(avg_ripple_baseline);
sample_mn = min(median_ripple);
deg_free = length(ripple_period)-1;
t = (sample_mn - mn)/(st_d/(sqrt(deg_free)));


%% Plot photometry against ripples

figure('color','white');
plot(time, median_ripple, 'g', 'LineWidth', 2);
hold on
fill([time,fliplr(time)], [(ripple_CI95(1,:)+median_ripple),fliplr((ripple_CI95(2,:)+median_ripple))], 'b', 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score Around Ripples - N11 sess 17 Striatum');% (Pre-task Sleep)');
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
title('Average Z-score Around Ripples (Pre-task Sleep)');
grid on;
hold off



=======
    % Get the timestamp for the current trial
    curr_time = behavTrials.timestamps(i);
    [~, trial_idx] = min(abs(photom_var.timestamps - curr_time));
   
    % Extract window of data
    start_idx = trial_idx - samples;
    end_idx = trial_idx + samples;
   
    if start_idx >= 1 && end_idx <= length(photom_var.timestamps)
        trial_count = trial_count + 1;
        patch_zscore_matrix(trial_count, :) = photom_var.grabDA_z(start_idx:end_idx);
    end
end

% Plot the final patch
if patch_num > 0 && trial_count > 0
    figure('color','white');
    hold on
    cmap = jet(ceil(trial_count/10));
    for b = 1:ceil(trial_count/10)
        avg_signal = mean(patch_zscore_matrix(((b-1)*10 + 1):min(b*10, trial_count), :), 1);
        smoothed = smoothdata(avg_signal);

        baseline = mean(smoothed(1:baseline_samples));
        normalized_signal = smoothed - baseline;
        plot(time, normalized_signal, 'LineWidth', 2, 'Color', cmap(b, :));
    end
   
    % Labeling
    xlabel('Time (s)');
    ylabel('Avg Z-score');
    title(['Patch ', num2str(patch_num), ': Dopamine Signal in Increments of 10 Trials']);
    grid on;
   
    % Add colorbar with proper trial range labeling
    colormap(jet)
    colorbar_handle = colorbar;
    tick_positions = linspace(0, 1, ceil(trial_count/10));
    trial_ranges = arrayfun(@(x) sprintf('%d-%d', (x-1)*10 + 1, min(x*10, trial_count)), 1:ceil(trial_count/10), 'UniformOutput', false);
    colorbar_handle.Ticks = tick_positions;
    colorbar_handle.TickLabels = trial_ranges;
    colorbar_handle.Label.String = 'Trial Range';
   
    hold off
end
%}
>>>>>>> 2f4d3a86f04df0da4d5e289cb15f173d8d44909b




