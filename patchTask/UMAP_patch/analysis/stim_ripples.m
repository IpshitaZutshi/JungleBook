%function [sleep_photometry] = dualSleepPhotometry(varargin)

%% Stim photometry
%
% USAGE
%   Looks at ripples around stimulation in post sleep
%
% INPUTS 
%    showfig    true/false - show a summary figure of the results
%         (default:true)
%
%    =========================================================================


%% Load variables

% p = inputParser;
% addParameter(p,'showfig',true,@islogical)
% addParameter(p,'saveMat',true,@islogical)
% 
% parse(p,varargin{:})
% showfig = p.Results.showfig;
% saveMat = p.Results.saveMat;

direc = '\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\patchTask';

save_path = '\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\patchTask\figures';
saveLoc = strcat(save_path,'\behavior');

fold = 3; % set fold to 3 if you want to look at post sleep, 1 if you want to look at pre

sessions_N17 = {'N17\N17_250430_sess9', ...
    'N17\N17_250501_sess10', ...
    'N17\N17_250509_sess15', ...
    'N17\N17_250511_sess17', ...
    'N17\N17_250512_sess18', ...
    'N17\N17_250513_sess19', ...
    'N17\N17_250519_sess21', ...
    'N17\N17_250520_sess22'};

sessions = sessions_N17;


% Initialize matrix to store ripple counts for each session
all_ripple_counts = [];
session_labels = {};


for s = 1:length(sessions)
    sessionpath = [direc filesep sessions{s}];
    cd(sessionpath);
    disp(['Collecting data from ' sessions{s}]);


    basepath = pwd;
    [~, currentFolderName] = fileparts(basepath);
    
    sampling_rate = 130; % sampling rate of photometry set up - 130
    
    ripple_file = dir(fullfile(basepath, '*ripples.events.mat'));
    load(ripple_file.name)
    
    [sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);
    if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
        load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
        
        tokens = regexp(sessionInfo.FileName, 'sess(\d+)$', 'tokens');
        sessionStr = tokens{1}{1};
        sleep_photometry.session = sessionStr; 

        if isempty(dir([basepath filesep MergePoints.foldernames{fold} filesep 'top*']))
            cd([basepath filesep MergePoints.foldernames{fold}]); 
            % set pre_post to be 1 for pre behavior sleep, 2 for post behavior sleep  
    
            pre_post = 2;
            
            curr_folder = pwd;
            % get HPC sleep photometry data
            if ~isempty (dir(fullfile(curr_folder, '*HPC*photometry.mat'))) 
                photometry_file = dir(fullfile(curr_folder, '*HPCSleepPhotomSynced.mat'));
                load(photometry_file.name);
                HPC_photometry_file = dir(fullfile(curr_folder, '*HPC*_photometry.mat'));
                load(HPC_photometry_file.name);
                hpc_stim = photometryData;
            end
    
            % get striatum sleep photometry data
            if ~isempty (dir(fullfile(curr_folder, '*striatum*photometry.mat')))
                photometry_file = dir(fullfile(curr_folder, '*StriatumSleepPhotomSynced.mat'));
                load(photometry_file.name);
                striatum_photometry_file = dir(fullfile(curr_folder, '*striatum*_photometry.mat'));
                load(striatum_photometry_file.name);
                striatum_stim = photometryData;
            end
                
            sleep_start = MergePoints.timestamps(fold,1); 
            sleep_end = MergePoints.timestamps(fold,2); 

            if ~isempty (dir(fullfile(curr_folder, '*striatum*photometry.mat')))
                striatum_photometry_file = dir(fullfile(curr_folder, '*striatum*_photometry.mat'));
                load(striatum_photometry_file.name);
                
                timestamps = striatum_sleep_sync.timestamps; % vector of timestamps
                stim_vector = striatum_stim.stim_time;      % 0/1 vector
            else
                HPC_photometry_file = dir(fullfile(curr_folder, '*HPC*_photometry.mat'));
                load(HPC_photometry_file.name);
                
                timestamps = hpc_sleep_sync.timestamps; % vector of timestamps
                stim_vector = hpc_stim.stim_time;      % 0/1 vector
            end
    
         
            
            ripple_times = ripples.peaks;            % ripple timestamps
            
            % Initialize logical array for ripples during stimulation
            ripples_during_stim = false(size(ripple_times));
            
            % For each ripple, find nearest photometry timestamp index
            for i = 1:length(ripple_times)
                % Find index of closest timestamp in photometry data to ripple time
                [~, idx] = min(abs(timestamps - ripple_times(i)));
            
                % Check if stimulation was ON at this timestamp
                ripples_during_stim(i) = stim_vector(idx) == 1;
            end
            
            % Display results
            num_ripples = length(ripple_times);
            num_ripples_stim = sum(ripples_during_stim);
            fprintf('Total ripples: %d\nRipples during stimulation: %d\n', num_ripples, num_ripples_stim);
            
            % Optionally, get the ripple times during stimulation
            ripple_times_stim = ripple_times(ripples_during_stim);
    
            counts = [sum(ripples_during_stim), sum(~ripples_during_stim)];
    
            figure('Color','w');
            bar(categorical({'Stim ON','Stim OFF'}), counts, 0.5);
            ylabel('Ripple Count');
            title('Number of Ripples During vs. Outside Stimulation');
    
            %% DA during stimulation
            % Find indices where stim_time switches from 0 to 1 (rising edges)
            stim_diff = diff([0; striatum_stim.stim_time(:)]);  % prepend 0 for correct diff length
            stim_onset_indices = find(stim_diff == 1);
            
            % Select timestamps at those indices
            stim_times = striatum_sleep_sync.timestamps(stim_onset_indices);
    
            window = 6; % window of time around ripple to average
            samples = window*sampling_rate;
    
            % initialize matrix
            stim_matrix = nan(length(stim_times), (samples*2)+1); 
    
    
            % average photometry data within a specified time window around stim
            adjusted_ts = striatum_sleep_sync.timestamps + sleep_start;
            stim_times = stim_times + sleep_start;
            for j = 1:length(stim_times)
                [~, stim_idx] = min(abs(adjusted_ts - stim_times(j)));
                start_idx = stim_idx - samples;
                end_idx = stim_idx + samples;
                if start_idx >= 1 && end_idx <= height(adjusted_ts)
                    stim_matrix(j, :) = striatum_sleep_sync.grabDA_z(start_idx:end_idx);
                    
                    % z score
                    baseline = stim_matrix(j, :);
                    % using -4s to -1.6s before ripple as baseline
                    % period
                    baseline_period = stim_matrix(j, 1:650); % Define baseline window before ripple event
                    baseline = baseline_period;
    
                    % Calculate baseline mean and std
                    baseline_mean = mean(baseline);
                    baseline_std = std(baseline);
    
                    % Z-score the entire trial using the baseline stats
                    stim_matrix(j, :) = (stim_matrix(j, :) - baseline_mean) / baseline_std;
                end
            end
    
            stim_matrix(any(isnan(stim_matrix), 2), :) = [];   
    
            if pre_post == 1 
                pp = 'pre_sleep';
            elseif pre_post == 2
                pp = 'post_sleep';  
            end
    
            tokens = regexp(sessionInfo.FileName, 'sess(\d+)$', 'tokens');
            sessionStr = tokens{1}{1};
    
            median_ripple = median(stim_matrix, 1); % median at each timepoint
            time = linspace(-window, window, ((samples*2)+1));
    
            % calculate confidence intervals
            N = height(stim_matrix);                          % Number of Experiments In Data Set
            avg_ripple = mean(stim_matrix, 1);              % Mean Of All Experiments At Each Value Of ,x 
            ripple_SEM = std(stim_matrix, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
            CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
            ripple_CI95 = bsxfun(@times, ripple_SEM, CI95(:));  % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of exw
            smooth_CI95 = smoothdata(ripple_CI95, 2);
    
    
            %% Plot photometry during stim
    
            striatum_color = [0.909803921568627   0.290196078431373   0.454901960784314];
        
            figure('color','white');
            plot(time, median_ripple, 'Color', striatum_color, 'LineWidth', 2);
            hold on
            fill([time,fliplr(time)], [(ripple_CI95(1,:)+median_ripple),fliplr((ripple_CI95(2,:)+median_ripple))], striatum_color, 'EdgeColor','none', 'FaceAlpha',0.25)
            xline(0, '--r', 'LineWidth', 1)
            xlabel('time (s)');
            ylabel('avg z-score');
            title(['Average Z-score Around stim - striatum'], [num2str(size(stim_matrix, 1)), ' stims']);
            grid on;
            hold off
    
            %% how many ripples around stim
    
            % Define time windows (seconds)
            pre_window = [-3 0];
            post_window = [0 3];
            late_post_window = [3 6];
            
            % Initialize counts
            count_pre = 0;
            count_post = 0;
            count_late_post = 0;
            
            % Loop over all stim times and count ripples in each window
            for i = 1:length(stim_times)
                t = stim_times(i);
                
                % Count ripples in pre-window
                count_pre = count_pre + sum(ripple_times >= (t + pre_window(1)) & ripple_times < (t + pre_window(2)));
                
                % Count ripples in post-window
                count_post = count_post + sum(ripple_times >= (t + post_window(1)) & ripple_times < (t + post_window(2)));
                
                % Count ripples in late post-window
                count_late_post = count_late_post + sum(ripple_times >= (t + late_post_window(1)) & ripple_times < (t + late_post_window(2)));
            end
            
            % Prepare data for bar plot
            ripple_counts = [count_pre, count_post, count_late_post];
            
            % Plot bar graph
            figure('Color','w');
            bar(categorical({'-3 to 0 s', '0 to 3 s', '3 to 6 s'}), ripple_counts, 0.5);
            ylabel('Number of Ripples');
            title('Ripple Counts Around Stimulation Onset');
            grid on;

            all_ripple_counts(end+1, :) = [count_pre, count_post, count_late_post];
            session_labels{end+1} = sessions{s};
        end
    end
end

%% dot plot
% Convert to array
counts_array = all_ripple_counts;

% X positions for time windows
x = [1, 2, 3];
labels = {'-3 to 0 s', '0 to 3 s', '3 to 6 s'};

figure('Color','w'); hold on;
colors = lines(size(counts_array,1));
%colors = jet(size(counts_array,1));

% Plot dots for each session
for sess = 1:size(counts_array,1)
    plot(x, counts_array(sess,:), '-', 'Color', colors(sess,:), 'LineWidth', 1.5, 'HandleVisibility','off');
    scatter(x, counts_array(sess,:), 60, 'filled', 'MarkerFaceColor', colors(sess,:), 'DisplayName', session_labels{sess});
end


xticks(x);
xticklabels(labels);
ylabel('Ripple Count');
title('Ripples Around Stim Onset per Session');
legend('Location','bestoutside');
grid on;



% Define analysis window
window = 5;  % seconds before and after stimulation onset
samples = window * sampling_rate;  % samples per side

% Find stimulation onsets (0 → 1 transitions)
stim_diff = diff([0; stim_vector(:)]);
stim_onsets_idx = find(stim_diff == 1);
stim_onsets_time = timestamps(stim_onsets_idx);  % time values of stim onsets

% Initialize matrix to hold aligned DA traces
stim_aligned_DA = nan(length(stim_onsets_time), samples*2 + 1);

% Extract DA traces around each stim onset
for i = 1:length(stim_onsets_idx)
    idx = stim_onsets_idx(i);
    start_idx = idx - samples;
    end_idx = idx + samples;
    
    if start_idx > 0 && end_idx <= length(hpc_sleep_sync.grabDA_z)
        stim_aligned_DA(i, :) = hpc_sleep_sync.grabDA_z(start_idx:end_idx);
    end
end

% Remove any rows with NaNs (e.g., boundary issues)
stim_aligned_DA(any(isnan(stim_aligned_DA), 2), :) = [];

% Average across stimulations
avg_stim_trace = mean(stim_aligned_DA, 1);
sem_stim_trace = std(stim_aligned_DA, 0, 1) ./ sqrt(size(stim_aligned_DA, 1));
t = linspace(-window, window, size(stim_aligned_DA, 2));

% Plot
figure('Color','w');
hold on;
plot(t, avg_stim_trace, 'k', 'LineWidth', 2);
fill([t fliplr(t)], ...
     [avg_stim_trace + sem_stim_trace, fliplr(avg_stim_trace - sem_stim_trace)], ...
     [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
xline(0, '--r', 'LineWidth', 1.5);
xlabel('Time (s) relative to stim onset');
ylabel('Avg. DA signal (z-score)');
title('Average Photometry Trace Around Stimulation');
grid on;


x=1;

%end


