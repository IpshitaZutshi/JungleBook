function [sleep_photometry] = dualSleepPhotometry(varargin)

%% Analyze the sleep photometry data for mice with two simultaneous recording sites
%
% USAGE
%   Takes in photometry data from the pre- and post- sleep sessions that has 
%   not already been synced to intan. Plots the averaged photometry data around 
%   ripple events for the entire session.
%
% INPUTS 
%    showfig    true/false - show a summary figure of the results
%         (default:true)
%
%    =========================================================================


%% Load variables

p = inputParser;
addParameter(p,'showfig',true,@islogical)
addParameter(p,'saveMat',true,@islogical)

parse(p,varargin{:})
showfig = p.Results.showfig;
saveMat = p.Results.saveMat;

basepath = pwd;
[~, currentFolderName] = fileparts(basepath);
save_mat = 0;
sleep_photometry = [];

sampling_rate = 130; % sampling rate of photometry set up - 130

ripple_file = dir(fullfile(basepath, '*ripples.events.mat'));
load(ripple_file.name)

[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);
if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
    load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
    
    tokens = regexp(sessionInfo.FileName, 'sess(\d+)$', 'tokens');
    sessionStr = tokens{1}{1};
    sleep_photometry.session = sessionStr; 
    
    for ii = 1:size(MergePoints.foldernames,2)
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
            % sync HPC sleep photometry data
            if ~isempty (dir(fullfile(curr_folder, '*HPC*photometry.mat'))) 
                if isempty(dir(fullfile(curr_folder, '*HPCSleepPhotomSynced.mat'))) 
                    HPC_photometry_file = dir(fullfile(curr_folder, '*HPC*_photometry.mat'));
                    load(HPC_photometry_file.name);
                
                    % synchronize time stamps for sleep session
                    disp('Synchronizing HPC photometry');
                    hpc_sleep_sync = getSyncPhotometry(photometryData);
                
                    C = strsplit(pwd,'\');
                    save([curr_folder filesep C{end} '.HPCSleepPhotomSynced.mat'],'hpc_sleep_sync');
                else
                    photometry_file = dir(fullfile(curr_folder, '*HPCSleepPhotomSynced.mat'));
                    load(photometry_file.name);
                    save_mat = 0;
                end
            end
    
            % sync striatum sleep photometry data
            if ~isempty (dir(fullfile(curr_folder, '*striatum*photometry.mat')))
                if isempty(dir(fullfile(curr_folder, '*StriatumSleepPhotomSynced.mat')))
                    striatum_photometry_file = dir(fullfile(curr_folder, '*striatum*_photometry.mat'));
                    load(striatum_photometry_file.name);
                
                    % synchronize time stamps for sleep session
                    disp('Synchronizing striatum photometry');
                    striatum_sleep_sync = getSyncPhotometry(photometryData);
                
                    C = strsplit(pwd,'\');
                    save([curr_folder filesep C{end} '.StriatumSleepPhotomSynced.mat'],'striatum_sleep_sync');
                else
                    photometry_file = dir(fullfile(curr_folder, '*StriatumSleepPhotomSynced.mat'));
                    load(photometry_file.name);
                    save_mat = 0;
                end
            end

            if pre_post == 1 
                sleep_start = MergePoints.timestamps(1,1); % start of sleep session
                sleep_end = MergePoints.timestamps(1,2); % time stamp of end of sleep session 
            elseif pre_post == 2
                sleep_start = MergePoints.timestamps(3,1); 
                sleep_end = MergePoints.timestamps(3,2);  
            end

%{
    sleep_start = MergePoints.timestamps(2,1); 
    sleep_end = MergePoints.timestamps(2,2); 
%}
            % average photometry around ripples - HPC
            if ~isempty(dir(fullfile(pwd, '*HPCSleepPhotomSynced.mat')))
                ripple_period = ripples.peaks(ripples.peaks <= sleep_end & ripples.peaks >= sleep_start);

                % find average time between ripple events
                ripple_durs = zeros(length(ripple_period), 1);
                for k = 2:length(ripple_period)
                    ripple_durs(k) = ripple_period(k)-ripple_period(k-1);
                end
                median(ripple_durs)
                
                window = 5; % window of time around ripple to average
                samples = window*sampling_rate;
    
                % initialize matrix
                ripple_matrix = nan(length(ripple_period), (samples*2)+1); 
    
    
                % average photometry data within a specified time window around ripples
                adjusted_ts = hpc_sleep_sync.timestamps + sleep_start;
                for j = 1:length(ripple_period)
                    [~, ripple_idx] = min(abs(adjusted_ts - ripple_period(j)));
                    start_idx = ripple_idx - samples;
                    end_idx = ripple_idx + samples;
                    if start_idx >= 1 && end_idx <= height(adjusted_ts)
                        ripple_matrix(j, :) = hpc_sleep_sync.grabDA_z(start_idx:end_idx);
                        
                        % z score
                        baseline = ripple_matrix(j, :);
                        % using -4s to -1.6s before ripple as baseline
                        % period
                        baseline_period = ripple_matrix(j, 130:440); % Define baseline window before ripple event
                        baseline = baseline_period;
    
                        % Calculate baseline mean and std
                        baseline_mean = mean(baseline);
                        baseline_std = std(baseline);
    
                        % Z-score the entire trial using the baseline stats
                        ripple_matrix(j, :) = (ripple_matrix(j, :) - baseline_mean) / baseline_std;
                    end
                end
    
                ripple_matrix(any(isnan(ripple_matrix), 2), :) = [];

                % save data
                if pre_post == 1 
                    pp = 'pre_sleep';
                elseif pre_post == 2
                    pp = 'post_sleep';  
                end
    
                if save_mat == 1
                    save([basepath filesep 'HPC_' pp '_' sessionStr],'ripple_matrix');
                    %save([basepath filesep 'HPC_' pp '_sleep_' sessionStr],'ripple_matrix');
                end

                sleep_photometry.hpc.(pp) = ripple_matrix;
    
                % prepare plot
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
                ripple_baseline = ripple_matrix(:,1:120);
                avg_ripple_baseline = mean(ripple_baseline, 1);
                mn = mean(avg_ripple_baseline);
                st_d = std(avg_ripple_baseline);
                sample_mn = min(median_ripple);
                deg_free = length(ripple_period)-1;
                t = (sample_mn - mn)/(st_d/(sqrt(deg_free)));
    
                %% Plot photometry against ripples
                
                if showfig
                    hpc_color = [0.588235294117647   0.800000000000000   0.345098039215686];
    
                    figure('color','white');
                    plot(time, median_ripple, 'Color', hpc_color, 'LineWidth', 2);
                    hold on
                    fill([time,fliplr(time)], [(ripple_CI95(1,:)+median_ripple),fliplr((ripple_CI95(2,:)+median_ripple))], hpc_color, 'EdgeColor','none', 'FaceAlpha',0.25)
                    xline(0, '--r', 'LineWidth', 1)
                    xlabel('time (s)');
                    ylabel('avg z-score');
                    title(['Average Z-score Around Ripples - HPC'], [num2str(size(ripple_matrix, 1)), ' ripples']);
                    grid on;
                    hold off
                    
                    
                    smoothed = smoothdata(median_ripple);
                    figure('color','white');
                    plot(time, smoothed, 'Color', hpc_color, 'LineWidth', 2);
                    hold on
                    fill([time,fliplr(time)], [(smooth_CI95(1,:)+smoothed),fliplr((smooth_CI95(2,:)+smoothed))], hpc_color, 'EdgeColor','none', 'FaceAlpha',0.25)
                    xline(0, '--r', 'LineWidth', 1)
                    xlabel('time (s)');
                    ylabel('avg z-score');
                    title(['Average Z-score Around Ripples (smoothed) - HPC'], [num2str(size(ripple_matrix, 1)), ' ripples']);
                    grid on;
                    hold off
                end
            end


            % average photometry around ripples - striatum
            if ~isempty(dir(fullfile(pwd, '*StriatumSleepPhotomSynced.mat')))
                ripple_period = ripples.peaks(ripples.peaks <= sleep_end & ripples.peaks >= sleep_start);

                % find average time between ripple events
                ripple_durs = zeros(length(ripple_period), 1);
                for k = 2:length(ripple_period)
                    ripple_durs(k) = ripple_period(k)-ripple_period(k-1);
                end
                median(ripple_durs)
                
                window = 5; % window of time around ripple to average
                samples = window*sampling_rate;
    
                % initialize matrix
                ripple_matrix = nan(length(ripple_period), (samples*2)+1); 
    
    
                % average photometry data within a specified time window around ripples
                adjusted_ts = striatum_sleep_sync.timestamps + sleep_start;
                for j = 1:length(ripple_period)
                    [~, ripple_idx] = min(abs(adjusted_ts - ripple_period(j)));
                    start_idx = ripple_idx - samples;
                    end_idx = ripple_idx + samples;
                    if start_idx >= 1 && end_idx <= height(adjusted_ts)
                        ripple_matrix(j, :) = striatum_sleep_sync.grabDA_z(start_idx:end_idx);
                        
                        % z score
                        baseline = ripple_matrix(j, :);
                        % using -4s to -1.6s before ripple as baseline
                        % period
                        baseline_period = ripple_matrix(j, 130:440); % Define baseline window before ripple event
                        baseline = baseline_period;

                        % Calculate baseline mean and std
                        baseline_mean = mean(baseline);
                        baseline_std = std(baseline);

                        % Z-score the entire trial using the baseline stats
                        ripple_matrix(j, :) = (ripple_matrix(j, :) - baseline_mean) / baseline_std;
                    end
                end
    
                ripple_matrix(any(isnan(ripple_matrix), 2), :) = [];   
    
                if pre_post == 1 
                    pp = 'pre_sleep';
                elseif pre_post == 2
                    pp = 'post_sleep';  
                end
    
                tokens = regexp(sessionInfo.FileName, 'sess(\d+)$', 'tokens');
                sessionStr = tokens{1}{1};
    
                if save_mat == 1
                    save([basepath filesep 'striatum_' pp '_' sessionStr],'ripple_matrix');
                end

                sleep_photometry.striatum.(pp) = ripple_matrix;
    
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
                ripple_baseline = ripple_matrix(:,1:120);
                avg_ripple_baseline = mean(ripple_baseline, 1);
                mn = mean(avg_ripple_baseline);
                st_d = std(avg_ripple_baseline);
                sample_mn = min(median_ripple);
                deg_free = length(ripple_period)-1;
                t = (sample_mn - mn)/(st_d/(sqrt(deg_free)));
    
                %% Plot photometry against ripples

                if showfig
                    striatum_color = [0.909803921568627   0.290196078431373   0.454901960784314];
                
                    figure('color','white');
                    plot(time, median_ripple, 'Color', striatum_color, 'LineWidth', 2);
                    hold on
                    fill([time,fliplr(time)], [(ripple_CI95(1,:)+median_ripple),fliplr((ripple_CI95(2,:)+median_ripple))], striatum_color, 'EdgeColor','none', 'FaceAlpha',0.25)
                    xline(0, '--r', 'LineWidth', 1)
                    xlabel('time (s)');
                    ylabel('avg z-score');
                    title(['Average Z-score Around Ripples - striatum'], [num2str(size(ripple_matrix, 1)), ' ripples']);
                    grid on;
                    hold off
                    
                    smoothed = smoothdata(median_ripple);
                    figure('color','white');
                    plot(time, smoothed, 'Color', striatum_color, 'LineWidth', 2);
                    hold on
                    fill([time,fliplr(time)], [(smooth_CI95(1,:)+smoothed),fliplr((smooth_CI95(2,:)+smoothed))], striatum_color, 'EdgeColor','none', 'FaceAlpha',0.25)
                    xline(0, '--r', 'LineWidth', 1)
                    xlabel('time (s)');
                    ylabel('avg z-score');
                    title(['Average Z-score Around Ripples (smoothed) - striatum'], [num2str(size(ripple_matrix, 1)), ' ripples']);
                    grid on;
                    hold off
                end
            end
        end
    end
end

% Combine and plot pre and post sleep averages for HPC and Striatum

if showfig && isfield(sleep_photometry, 'hpc') && isfield(sleep_photometry, 'striatum')
    figure('color', 'white');
    hold on;

    time = linspace(-window, window, size(sleep_photometry.hpc.pre_sleep, 2));

    % Average pre and post for HPC
    avg_hpc = mean([sleep_photometry.hpc.pre_sleep; sleep_photometry.hpc.post_sleep], 1);
    hpc_color = [0.588 0.8 0.345];
    plot(time, avg_hpc, 'Color', hpc_color, 'LineWidth', 2, 'DisplayName', 'HPC Avg Pre/Post Sleep');

    xline(0, '--r', 'LineWidth', 1)
    xlabel('Time around ripple (s)');
    ylabel('Average z-scored DA signal');
    title('Average Photometry Around Ripples');
    % legend('Location', 'best');
    grid on;
    hold off;

    figure('color', 'white');
    hold on;

    time = linspace(-window, window, size(sleep_photometry.hpc.pre_sleep, 2));

    % Average pre and post for Striatum
    avg_striatum = mean([sleep_photometry.striatum.pre_sleep; sleep_photometry.striatum.post_sleep], 1);
    striatum_color = [0.91 0.29 0.45];
    plot(time, avg_striatum, 'Color', striatum_color, 'LineWidth', 2, 'DisplayName', 'Striatum Avg Pre/Post Sleep');

    xline(0, '--r', 'LineWidth', 1)
    xlabel('Time around ripple (s)');
    ylabel('Average z-scored DA signal');
    title('Average Photometry Around Ripples');
    % legend('Location', 'best');
    grid on;
    hold off;
end

% % Assuming HPC photometry data is loaded and synced, i.e.
% % hpc_sleep_sync.timestamps - vector of time stamps for dopamine data
% % hpc_sleep_sync.grabDA_z - dopamine signal (z-scored)
% % ripples.peaks - vector of ripple peak timestamps
% 
% % 1) Create binary vector marking ripple events aligned to dopamine timestamps
% dt = mean(diff(hpc_sleep_sync.timestamps)); % sampling interval in seconds
% ripple_binary = zeros(size(hpc_sleep_sync.timestamps));
% 
% for i = 1:length(ripples.peaks)
%     % Find nearest dopamine timestamp index for each ripple peak
%     [~, idx] = min(abs(hpc_sleep_sync.timestamps - ripples.peaks(i)));
%     ripple_binary(idx) = 1;
% end
% 
% % 2) Dopamine signal vector
% dopamine_signal = hpc_sleep_sync.grabDA_z;
% 
% % 3) Calculate Pearson correlation (zero lag)
% [rho, pval] = corr(dopamine_signal, ripple_binary);
% 
% fprintf('Correlation between dopamine and ripple events: r = %.3f, p = %.3g\n', rho, pval);
% 
% % 4) Cross-correlation for temporal lag analysis (±5 seconds window)
% max_lag_sec = 5;
% max_lag_samples = round(max_lag_sec/dt);
% [xcorr_vals, lags] = xcorr(dopamine_signal - mean(dopamine_signal), ripple_binary - mean(ripple_binary), max_lag_samples, 'coeff');
% lags_sec = lags * dt;
% 
% % 5) Plot cross-correlation
% figure('colo', 'white)');
% plot(lags_sec, xcorr_vals, 'LineWidth', 2);
% xlabel('Lag (seconds)');
% ylabel('Cross-correlation coefficient');
% title('Cross-correlation between dopamine signal and ripple events');
% grid on;
% xline(0, '--r');


x=1;

end

%% full trace

% BLUE = NREM
% RED = REM
% YELLOW = WAKE
%{
if exist('full_photometry')
    sleep_state_file = dir(fullfile(basepath, '..', '*SleepState.states.mat'));
    load(sleep_state_file.name)
    adjusted_ts = sleep_sync.timestamps + sleep_start;
    
    darkerGreen = [0.1098    0.6000    0.2392];
    striatum_pink = [0.960784313725490, 0.152941176470588, 0.905882352941176];
    figure('color','white')
    ax = gca;
    ax.FontSize = 15;
    hold on
    plot(adjusted_ts, sleep_sync.grabDA_z, 'Color', darkerGreen, 'LineWidth', 2);
    for ii = 1:size(SleepState.ints.NREMstate, 1)
        y = [-4, -4, 6, 6];
        x = [SleepState.ints.NREMstate(ii, 1), SleepState.ints.NREMstate(ii, 2), ...
            SleepState.ints.NREMstate(ii, 2), SleepState.ints.NREMstate(ii, 1)];
        patch(x, y, [0.301960784313725   0.745098039215686   0.933333333333333], 'FaceAlpha', 0.3, 'EdgeAlpha', 0)
    end
    
    for ii = 1:size(SleepState.ints.REMstate, 1)
        y = [-4, -4, 6, 6];
        x = [SleepState.ints.REMstate(ii, 1), SleepState.ints.REMstate(ii, 2), ...
            SleepState.ints.REMstate(ii, 2), SleepState.ints.REMstate(ii, 1)];
        patch(x, y, [1.000000000000000   0.266666666666667                   0], 'FaceAlpha', 0.3, 'EdgeAlpha', 0)
    end
    
    for ii = 1:size(SleepState.ints.WAKEstate, 1)
        y = [-4, -4, 6, 6];
        x = [SleepState.ints.WAKEstate(ii, 1), SleepState.ints.WAKEstate(ii, 2), ...
            SleepState.ints.WAKEstate(ii, 2), SleepState.ints.WAKEstate(ii, 1)];
        patch(x, y, [1.000000000000000   1.000000000000000   0.066666666666667], 'FaceAlpha', 0.3, 'EdgeAlpha', 0)
    end
    xlim([adjusted_ts(1), adjusted_ts(end)]);
    ylim([min(y), max(y)])
    hold off

    sleep_state_file = dir(fullfile(basepath, '..', '*SleepState.states.mat'));
    load(sleep_state_file.name)
    adjusted_ts = sleep_sync.timestamps + sleep_start;
    
    darkerGreen = [0.1098    0.6000    0.2392];
    striatum_pink = [0.960784313725490, 0.152941176470588, 0.905882352941176];
    figure('color','white')
    ax = gca;
    ax.FontSize = 15;
    hold on
    plot(adjusted_ts, sleep_sync.grabDA_z, 'Color', darkerGreen, 'LineWidth', 2);
    for ii = 1:size(SleepState.ints.NREMstate, 1)
        y = [-4, -4, 6, 6];
        x = [SleepState.ints.NREMstate(ii, 1), SleepState.ints.NREMstate(ii, 2), ...
            SleepState.ints.NREMstate(ii, 2), SleepState.ints.NREMstate(ii, 1)];
        patch(x, y, [0.301960784313725   0.745098039215686   0.933333333333333], 'FaceAlpha', 0.3, 'EdgeAlpha', 0)
    end
    
    for ii = 1:size(SleepState.ints.REMstate, 1)
        y = [-4, -4, 6, 6];
        x = [SleepState.ints.REMstate(ii, 1), SleepState.ints.REMstate(ii, 2), ...
            SleepState.ints.REMstate(ii, 2), SleepState.ints.REMstate(ii, 1)];
        patch(x, y, [1.000000000000000   0.266666666666667                   0], 'FaceAlpha', 0.3, 'EdgeAlpha', 0)
    end
    
    for ii = 1:size(SleepState.ints.WAKEstate, 1)
        y = [-4, -4, 6, 6];
        x = [SleepState.ints.WAKEstate(ii, 1), SleepState.ints.WAKEstate(ii, 2), ...
            SleepState.ints.WAKEstate(ii, 2), SleepState.ints.WAKEstate(ii, 1)];
        patch(x, y, [1.000000000000000   1.000000000000000   0.066666666666667], 'FaceAlpha', 0.3, 'EdgeAlpha', 0)
    end
    xlim([adjusted_ts(1), adjusted_ts(end)]);
    ylim([min(y), max(y)])
    hold off
end


if color == 0
                    % avg_color = [ 0.2392    0.2863    0.9608];
                    % conf_color = 'b';
                    avg_color = [0.960784313725490, 0.152941176470588, 0.905882352941176];
                    conf_color = [0.960784313725490, 0.152941176470588, 0.905882352941176];
                else
                    avg_color = 'g';
                    conf_color = [0.7176    0.9412    0.1020];
                end

%}
