% z score across behavior sessions
%
% USAGE
%    z score photometry data from across behavior sessions and plot the 
%    dopamine traces around rewarded or nonrewarded licks, in the HPC and striatum
% 
%
% INPUTS 
%    
%
%    =========================================================================

direc = '\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\patchTask';

save_path = '\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\patchTask\figures';
saveLoc = strcat(save_path,'\behavior');

sessions_N11 = {'N11\Final\N11_250312_sess16', ...
    'N11\Final\N11_250313_sess17', ...
    'N11\Final\N11_250314_sess18', ...
    'N11\Final\N11_250318_sess19' ...
    'N11\Final\N11_250319_sess20', ...
    'N11\Final\N11_250321_sess22', ...
    'N11\Final\N11_250331_sess24', ...
    'N11\Final\N11_250401_sess25', ...
    'N11\Final\N11_250403_sess26', ...
    'N11\Final\N11_250407_sess27', ...
    'N11\Final\N11_250408_sess28', ...
    'N11\Final\N11_250410_sess30', ...
    'N11\Final\N11_250411_sess31'};

sessions_N14 = {'N14\N14_250429_sess8', ...
    'N14\N14_250501_sess10', ...
    'N14\N14_250507_sess13', ...
    'N14\N14_250509_sess15', ...
    'N14\N14_250510_sess16', ...
    'N14\N14_250511_sess17', ...
    'N14\N14_250512_sess18', ...
    'N14\N14_250519_sess20', ...
    'N14\N14_250520_sess21', ...
    'N14\N14_250521_sess22'};

sessions_N17 = {'N17\N17_250430_sess9', ...
    'N17\N17_250501_sess10', ...
    'N17\N17_250509_sess15', ...
    'N17\N17_250510_sess16', ...
    'N17\N17_250511_sess17', ...
    'N17\N17_250512_sess18', ...
    'N17\N17_250513_sess19', ...
    'N17\N17_250519_sess21'};

sessions = cat(2, sessions_N11, sessions_N14, sessions_N17);
sessions = cat(2, sessions_N11, sessions_N17);

% Initialize arrays to accumulate data
all_rewarded_hpc = {};
all_nonrewarded_hpc = {};
all_expected_reward_hpc = {};
all_unexpected_reward_hpc = {};
all_expected_nonreward_hpc = {};
all_unexpected_nonreward_hpc = {};
all_awake_ripples_hpc = {};

all_rewarded_striatum = {};
all_nonrewarded_striatum = {};
all_expected_reward_striatum = {};
all_unexpected_reward_striatum = {};
all_expected_nonreward_striatum = {};
all_unexpected_nonreward_striatum = {};
all_awake_ripples_striatum = {};

all_ripples_near_rewarded = {};
all_ripples_near_nonrewarded = {};
all_ripples_near_expected_reward = {};
all_ripples_near_unexpected_reward = {};


for s = 1:length(sessions)
    sessionpath = [direc filesep sessions{s}];
    cd(sessionpath);
    disp(['Collecting data from ' sessions{s}]);

    %fprintf('Running session: %s\n', sessions(s).name);
    session_data = photometry_behav('showfig', false); 

    % Accumulate HPC
    if isfield(session_data, 'hpc')
        if isfield(session_data.hpc, 'rewarded')
            all_rewarded_hpc{end+1} = session_data.hpc.rewarded;
        end
        if isfield(session_data.hpc, 'nonrewarded')
            all_nonrewarded_hpc{end+1} = session_data.hpc.nonrewarded;
        end
        if isfield(session_data.hpc, 'awake_ripples')
            all_awake_ripples_hpc{end+1} = session_data.hpc.awake_ripples;
        end
        if isfield(session_data.hpc, 'expected_reward')
            all_expected_reward_hpc{end+1} = session_data.hpc.expected_reward;
        end
        if isfield(session_data.hpc, 'expected_no_reward')
            all_expected_nonreward_hpc{end+1} = session_data.hpc.expected_no_reward;
        end
        if isfield(session_data.hpc, 'unexpected_reward')
            all_unexpected_reward_hpc{end+1} = session_data.hpc.unexpected_reward;
        end
        if isfield(session_data.hpc, 'unexpected_no_reward')
            all_unexpected_nonreward_hpc{end+1} = session_data.hpc.unexpected_no_reward;
        end
    end
    
    % Accumulate STRIATUM
    if isfield(session_data, 'striatum')
        if isfield(session_data.striatum, 'rewarded')
            all_rewarded_striatum{end+1} = session_data.striatum.rewarded;
        end
        if isfield(session_data.striatum, 'nonrewarded')
            all_nonrewarded_striatum{end+1} = session_data.striatum.nonrewarded;
        end
        if isfield(session_data.striatum, 'awake_ripples')
            all_awake_ripples_striatum{end+1} = session_data.striatum.awake_ripples;
        end
        if isfield(session_data.striatum, 'expected_reward')
            all_expected_reward_striatum{end+1} = session_data.striatum.expected_reward;
        end
        if isfield(session_data.striatum, 'expected_no_reward')
            all_expected_nonreward_striatum{end+1} = session_data.striatum.expected_no_reward;
        end
        if isfield(session_data.striatum, 'unexpected_reward')
            all_unexpected_reward_striatum{end+1} = session_data.striatum.unexpected_reward;
        end
        if isfield(session_data.striatum, 'unexpected_no_reward')
            all_unexpected_nonreward_striatum{end+1} = session_data.striatum.unexpected_no_reward;
        end
    end

    
        %% RIpples around events
    behaviorFile = dir(fullfile(sessionpath, '*.TrialBehavior.mat'));
    load(behaviorFile.name);

    ripple_file = dir(fullfile(sessionpath, '*.ripples.events.mat'));
    load(ripple_file.name);
    
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
    
    ripple_times = ripples.peaks; % all ripple event timestamps in this session
    
    % Define a time window around events to check ripples
    window_sec = 2;
    
    % Function to count ripples near events
    countRipplesNearEvents = @(event_times) arrayfun(@(t) ...
        sum(ripple_times >= t - window_sec & ripple_times <= t + window_sec), ...
        event_times);
    
    % Count ripples near each event type
    ripples_near_expected_reward = countRipplesNearEvents(expected_reward(:, 1));
    ripples_near_unexpected_reward = countRipplesNearEvents(unexpected_reward(:, 1));
    ripples_near_expected_non = countRipplesNearEvents(expected_non(:, 1));
    ripples_near_unexpected_non = countRipplesNearEvents(unexpected_non(:, 1));
    
    % Store these counts per session for later analysis
    all_ripples_near_expected_reward{s} = ripples_near_expected_reward;
    all_ripples_near_unexpected_reward{s} = ripples_near_unexpected_reward;
    all_ripples_near_expected_non{s} = ripples_near_expected_non;
    all_ripples_near_unexpected_non{s} = ripples_near_unexpected_non;

end

%% average and plot across sessions

%% Rewarded vs non-rewarded

% Concatenate and average for HPC
rewarded_hpc_all = cat(1, all_rewarded_hpc{:});
nonrewarded_hpc_all = cat(1, all_nonrewarded_hpc{:});

        % rewarded
avg_rewarded_hpc = mean(rewarded_hpc_all, 1);
med_rewarded_hpc = median(rewarded_hpc_all, 1); % median at each timepoint

            % confidence intervals
N = height(rewarded_hpc_all);                          % # of experiments
a_SEM = std(rewarded_hpc_all, 1)/sqrt(N);         % Standard Error Of The Mean 
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
a_CI95 = bsxfun(@times, a_SEM, CI95(:));  % Calculate 95% Confidence Intervals 

        % nonrewarded
avg_nonrewarded_hpc = mean(nonrewarded_hpc_all, 1);
med_nonrewarded_hpc = median(nonrewarded_hpc_all, 1);

            % confidence intervals
N = height(nonrewarded_hpc_all);                          % Number of eExperimentsn In Data Set
b_SEM = std(nonrewarded_hpc_all, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
b_CI95 = bsxfun(@times, b_SEM, CI95(:));  % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of exw


%% PLOT reward/nonrewarded hpc

% Create time vector
sampling_rate = 130;
window = 5;
samples = window * sampling_rate;
time = linspace(-window, window, 2*samples+1);

plot_color = [0.588235294117647   0.800000000000000   0.345098039215686]; % green
all_axes = gobjects(2,1);  % Preallocate axes handles
y_min = inf;
y_max = -inf;

figure('color','white');
subplot(2,1,1)
hold on
plot(time, med_rewarded_hpc, 'color', plot_color, 'LineWidth', 2);
fill([time,fliplr(time)], [(a_CI95(1,:)+med_rewarded_hpc),fliplr((a_CI95(2,:)+med_rewarded_hpc))], plot_color, 'EdgeColor','none', 'FaceAlpha',0.25)
ax = gca;
ax.FontSize = 15;
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score Around Rewards');
grid on;
% Store axes handle
all_axes(1) = gca;  % for subplot 1
% Update Y-axis limits
current_ylim = ylim(gca);
y_min = min(y_min, current_ylim(1));
y_max = max(y_max, current_ylim(2));
hold off

subplot(2,1,2)
hold on
plot(time, med_nonrewarded_hpc, 'color', plot_color, 'LineWidth', 2);
fill([time,fliplr(time)], [(b_CI95(1,:)+med_nonrewarded_hpc),fliplr((b_CI95(2,:)+med_nonrewarded_hpc))], plot_color, 'EdgeColor','none', 'FaceAlpha',0.25)
ax = gca;
ax.FontSize = 15;
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score Around Non-rewarded Licks');
grid on;
all_axes(2) = gca; 
% Update Y-axis limits
current_ylim = ylim(gca);
y_min = min(y_min, current_ylim(1));
y_max = max(y_max, current_ylim(2));
hold off

% Set all to same y-limits
for i = 1:2
    ylim(all_axes(i), [y_min y_max]);
end

% Link all axes
linkaxes(all_axes, 'xy');



% Concatenate and average for STRIATUM
rewarded_striatum_all = cat(1, all_rewarded_striatum{:});
nonrewarded_striatum_all = cat(1, all_nonrewarded_striatum{:});

        % rewarded
avg_rewarded_striatum = mean(rewarded_striatum_all, 1);
med_rewarded_striatum = median(rewarded_striatum_all, 1);

            % confidence intervals
N = height(rewarded_striatum_all);                          % # of experiments
a_SEM = std(rewarded_striatum_all, 1)/sqrt(N);         % Standard Error Of The Mean 
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
a_CI95 = bsxfun(@times, a_SEM, CI95(:));  % Calculate 95% Confidence Intervals 

        % nonrewarded
avg_nonrewarded_striatum = mean(nonrewarded_striatum_all, 1);
med_nonrewarded_striatum = median(nonrewarded_striatum_all, 1);

            % confidence intervals
N = height(nonrewarded_striatum_all);                          % Number of eExperimentsn In Data Set
b_SEM = std(nonrewarded_striatum_all, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
b_CI95 = bsxfun(@times, b_SEM, CI95(:));  % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of exw


%% PLOT Reward/nonrewarded str

% Create time vector
sampling_rate = 130;
window = 5;
samples = window * sampling_rate;
time = linspace(-window, window, 2*samples+1);

plot_color = [0.909803921568627   0.290196078431373   0.454901960784314]; % pink
all_axes = gobjects(2,1);  % Preallocate axes handles
y_min = inf;
y_max = -inf;

figure('color','white');
subplot(2,1,1)
hold on
plot(time, med_rewarded_striatum, 'color', plot_color, 'LineWidth', 2);
fill([time,fliplr(time)], [(a_CI95(1,:)+med_rewarded_striatum),fliplr((a_CI95(2,:)+med_rewarded_striatum))], plot_color, 'EdgeColor','none', 'FaceAlpha',0.25)
ax = gca;
ax.FontSize = 15;
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score Around Rewards');
grid on;
% Store axes handle
all_axes(1) = gca;  % for subplot 1
% Update Y-axis limits
current_ylim = ylim(gca);
y_min = min(y_min, current_ylim(1));
y_max = max(y_max, current_ylim(2));
hold off

subplot(2,1,2)
hold on
plot(time, med_nonrewarded_striatum, 'color', plot_color, 'LineWidth', 2);
fill([time,fliplr(time)], [(b_CI95(1,:)+med_nonrewarded_striatum),fliplr((b_CI95(2,:)+med_nonrewarded_striatum))], plot_color, 'EdgeColor','none', 'FaceAlpha',0.25)
ax = gca;
ax.FontSize = 15;
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score Around Non-rewarded Licks');
grid on;
all_axes(2) = gca; 
% Update Y-axis limits
current_ylim = ylim(gca);
y_min = min(y_min, current_ylim(1));
y_max = max(y_max, current_ylim(2));
hold off

% Set all to same y-limits
for i = 1:2
    ylim(all_axes(i), [y_min y_max]);
end

% Link all axes
linkaxes(all_axes, 'xy');

%% plotted on same axes
customGreen = [152, 194, 9] / 255;
customRed = [238, 75, 43] / 255;

figure('color','white');
subplot(2,1,1)
plot(time, avg_rewarded_hpc, 'Color', customGreen, 'LineWidth', 2);
hold on
plot(time, avg_nonrewarded_hpc, 'Color', customRed, 'LineWidth', 2);
title('HPC: Rewarded vs Non-rewarded');
xlabel('Time (s)');
ylabel('Avg z-score');
legend('Rewarded', 'Non-rewarded');
grid on;

subplot(2,1,2)
plot(time, avg_rewarded_striatum, 'Color', customGreen, 'LineWidth', 2);
hold on
plot(time, avg_nonrewarded_striatum, 'Color', customRed, 'LineWidth', 2);
title('Striatum: Rewarded vs Non-rewarded');
xlabel('Time (s)');
ylabel('Avg z-score');
legend('Rewarded', 'Non-rewarded');
grid on;
set(gcf,'position',[500,200,1120,840])
saveas(gcf,[saveLoc,filesep ,'reward_overview_N11_N17','.png'],'png');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Awake ripples

% Concatenate and average for HPC
awake_ripples_hpc_all = cat(1, all_awake_ripples_hpc{:});

        % rewarded
avg_ripples_hpc = mean(awake_ripples_hpc_all, 1);
med_ripples_hpc = median(awake_ripples_hpc_all, 1); % median at each timepoint

            % confidence intervals
N = height(awake_ripples_hpc_all);                          % # of experiments
a_SEM = std(awake_ripples_hpc_all, 1)/sqrt(N);         % Standard Error Of The Mean 
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
a_CI95 = bsxfun(@times, a_SEM, CI95(:));  % Calculate 95% Confidence Intervals 
smooth_CI95 = smoothdata(a_CI95, 2);
       
% PLOT

% Create time vector
sampling_rate = 130;
window = 5;
samples = window * sampling_rate;
time = linspace(-window, window, 2*samples+1);

plot_color = [0.588235294117647   0.800000000000000   0.345098039215686]; % green

figure('color','white');
plot(time, med_ripples_hpc, 'Color', plot_color, 'LineWidth', 2);
hold on
fill([time,fliplr(time)], [(a_CI95(1,:)+med_ripples_hpc),fliplr((a_CI95(2,:)+med_ripples_hpc))], plot_color, 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title(['Average Z-score Around Ripples'], [num2str(size(awake_ripples_hpc_all, 1)), ' ripples']);
grid on;
set(gcf,'position',[500,200,1120,840])
saveas(gcf,[saveLoc,filesep ,'hpc_awake_ripples_N11_N17','.png'],'png');
hold off

smoothed = smoothdata(med_ripples_hpc);
figure('color','white');
plot(time, smoothed, 'Color', plot_color, 'LineWidth', 2);
hold on
fill([time,fliplr(time)], [(smooth_CI95(1,:)+smoothed),fliplr((smooth_CI95(2,:)+smoothed))], plot_color, 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title(['Average Z-score Around Ripples'], [num2str(size(awake_ripples_hpc_all, 1)), ' ripples']);
grid on;
set(gcf,'position',[500,200,1120,840])
saveas(gcf,[saveLoc,filesep ,'hpc_awake_ripples_smooth_N11_N17','.png'],'png');
hold off


% Concatenate and average for STRIATUM
awake_ripples_striatum_all = cat(1, all_awake_ripples_striatum{:});

        % rewarded
avg_ripples_striatum = mean(awake_ripples_striatum_all, 1);
med_ripples_striatum = median(awake_ripples_striatum_all, 1); % median at each timepoint

            % confidence intervals
N = height(awake_ripples_striatum_all);                          % # of experiments
a_SEM = std(awake_ripples_striatum_all, 1)/sqrt(N);         % Standard Error Of The Mean 
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
b_CI95 = bsxfun(@times, a_SEM, CI95(:));  % Calculate 95% Confidence Intervals 
smooth_CI95 = smoothdata(b_CI95, 2);
       
% PLOT

% Create time vector
sampling_rate = 130;
window = 5;
samples = window * sampling_rate;
time = linspace(-window, window, 2*samples+1);

plot_color = [0.909803921568627   0.290196078431373   0.454901960784314]; % pink

figure('color','white');
plot(time, med_ripples_striatum, 'Color', plot_color, 'LineWidth', 2);
hold on
fill([time,fliplr(time)], [(b_CI95(1,:)+med_ripples_striatum),fliplr((b_CI95(2,:)+med_ripples_striatum))], plot_color, 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title(['Average Z-score Around Ripples'], [num2str(size(awake_ripples_striatum_all, 1)), ' ripples']);
grid on;
set(gcf,'position',[500,200,1120,840])
saveas(gcf,[saveLoc,filesep ,'str_awake_ripples_N11_N17','.png'],'png');
hold off

smoothed = smoothdata(med_ripples_striatum);
figure('color','white');
plot(time, smoothed, 'Color', plot_color, 'LineWidth', 2);
hold on
fill([time,fliplr(time)], [(smooth_CI95(1,:)+smoothed),fliplr((smooth_CI95(2,:)+smoothed))], plot_color, 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title(['Average Z-score Around Ripples'], [num2str(size(awake_ripples_striatum_all, 1)), ' ripples']);
grid on;
set(gcf,'position',[500,200,1120,840])
saveas(gcf,[saveLoc,filesep ,'str_awake_ripples_smooth_N11_N17','.png'],'png');
hold off



%% EXPECTATIONS

%% Plot average z-scored dopamine signal: Expected/Unexpected Reward & Non-Reward (HPC)
% Create time vector
sampling_rate = 130;
window = 3.1;
samples = window * sampling_rate;
time = linspace(-window, window, 2*samples+1);

% Concatenate trials
expected_reward = cat(1, all_expected_reward_hpc{:});
unexpected_reward = cat(1, all_unexpected_reward_hpc{:});
expected_non = cat(1, all_expected_nonreward_hpc{:});
unexpected_non = cat(1, all_unexpected_nonreward_hpc{:});

% Compute medians
med_expected_reward = median(expected_reward, 1);
med_unexpected_reward = median(unexpected_reward, 1);
med_expected_non = median(expected_non, 1);
med_unexpected_non = median(unexpected_non, 1);

% Confidence intervals
CI95 = tinv([0.025 0.975], size(expected_reward, 1)-1);

expected_CI95 = bsxfun(@times, std(expected_reward, 1)/sqrt(size(expected_reward,1)), CI95(:));
unexpected_CI95 = bsxfun(@times, std(unexpected_reward, 1)/sqrt(size(unexpected_reward,1)), CI95(:));
expected_CI95_non = bsxfun(@times, std(expected_non, 1)/sqrt(size(expected_non,1)), CI95(:));
unexpected_CI95_non = bsxfun(@times, std(unexpected_non, 1)/sqrt(size(unexpected_non,1)), CI95(:));

% Plot setup
plot_color = [0.588235294117647   0.800000000000000   0.345098039215686]; % green
all_axes = gobjects(4,1);
y_min = inf;
y_max = -inf;

figure('color','white');

subplot(2,2,1)
hold on
ax = gca; ax.FontSize = 15;
plot(time, med_unexpected_reward, 'color', plot_color, 'LineWidth', 2);
fill([time,fliplr(time)], ...
     [(unexpected_CI95(1,:) + med_unexpected_reward), fliplr((unexpected_CI95(2,:) + med_unexpected_reward))], ...
     plot_color, 'EdgeColor','none', 'FaceAlpha',0.25);
xline(0, '--r', 'LineWidth', 1)
xlim([-3.1, 3.1]);
xlabel('Time (s)');
ylabel('Avg z-score');
title(sprintf('Unexpected Reward (%d trials)', size(unexpected_reward,1)));
grid on;
all_axes(1) = ax;
current_ylim = ylim; y_min = min(y_min, current_ylim(1)); y_max = max(y_max, current_ylim(2));
hold off

subplot(2,2,2)
hold on
ax = gca; ax.FontSize = 15;
plot(time, med_unexpected_non, 'color', plot_color, 'LineWidth', 2);
fill([time,fliplr(time)], ...
     [(unexpected_CI95_non(1,:) + med_unexpected_non), fliplr((unexpected_CI95_non(2,:) + med_unexpected_non))], ...
     plot_color, 'EdgeColor','none', 'FaceAlpha',0.25);
xline(0, '--r', 'LineWidth', 1)
xlim([-3.1, 3.1]);
xlabel('Time (s)');
ylabel('Avg z-score');
title(sprintf('Unexpected Non-reward (%d trials)', size(unexpected_non,1)));
grid on;
all_axes(2) = ax;
current_ylim = ylim; y_min = min(y_min, current_ylim(1)); y_max = max(y_max, current_ylim(2));
hold off

subplot(2,2,3)
hold on
ax = gca; ax.FontSize = 15;
plot(time, med_expected_reward, 'color', plot_color, 'LineWidth', 2);
fill([time,fliplr(time)], ...
     [(expected_CI95(1,:) + med_expected_reward), fliplr((expected_CI95(2,:) + med_expected_reward))], ...
     plot_color, 'EdgeColor','none', 'FaceAlpha',0.25);
xline(0, '--r', 'LineWidth', 1)
xlim([-3.1, 3.1]);
xlabel('Time (s)');
ylabel('Avg z-score');
title(sprintf('Expected Reward (%d trials)', size(expected_reward,1)));
grid on;
all_axes(3) = ax;
current_ylim = ylim; y_min = min(y_min, current_ylim(1)); y_max = max(y_max, current_ylim(2));
hold off

subplot(2,2,4)
hold on
ax = gca; ax.FontSize = 15;
plot(time, med_expected_non, 'color', plot_color, 'LineWidth', 2);
fill([time,fliplr(time)], ...
     [(expected_CI95_non(1,:) + med_expected_non), fliplr((expected_CI95_non(2,:) + med_expected_non))], ...
     plot_color, 'EdgeColor','none', 'FaceAlpha',0.25);
xline(0, '--r', 'LineWidth', 1)
xlim([-3.1, 3.1]);
xlabel('Time (s)');
ylabel('Avg z-score');
title(sprintf('Expected Non-reward (%d trials)', size(expected_non,1)));
grid on;
all_axes(4) = ax;
current_ylim = ylim; y_min = min(y_min, current_ylim(1)); y_max = max(y_max, current_ylim(2));
hold off
set(gcf,'position',[200,200,1600,840])

% Set all y-limits the same
for i = 1:4
    ylim(all_axes(i), [y_min y_max]);
end
linkaxes(all_axes, 'xy');

% Optional save
saveas(gcf, [saveLoc, filesep, 'expected_unexpected_hpc_N11_N17','.png'], 'png');

%% striatum
% Concatenate trials
expected_reward = cat(1, all_expected_reward_striatum{:});
unexpected_reward = cat(1, all_unexpected_reward_striatum{:});
expected_non = cat(1, all_expected_nonreward_striatum{:});
unexpected_non = cat(1, all_unexpected_nonreward_striatum{:});

% Compute medians
med_expected_reward = median(expected_reward, 1);
med_unexpected_reward = median(unexpected_reward, 1);
med_expected_non = median(expected_non, 1);
med_unexpected_non = median(unexpected_non, 1);

% Confidence intervals
CI95 = tinv([0.025 0.975], size(expected_reward, 1)-1);

expected_CI95 = bsxfun(@times, std(expected_reward, 1)/sqrt(size(expected_reward,1)), CI95(:));
unexpected_CI95 = bsxfun(@times, std(unexpected_reward, 1)/sqrt(size(unexpected_reward,1)), CI95(:));
expected_CI95_non = bsxfun(@times, std(expected_non, 1)/sqrt(size(expected_non,1)), CI95(:));
unexpected_CI95_non = bsxfun(@times, std(unexpected_non, 1)/sqrt(size(unexpected_non,1)), CI95(:));

% Plot setup
plot_color = [0.909803921568627   0.290196078431373   0.454901960784314]; % pink
all_axes = gobjects(4,1);
y_min = inf;
y_max = -inf;

figure('color','white');

subplot(2,2,1)
hold on
ax = gca; ax.FontSize = 15;
plot(time, med_unexpected_reward, 'color', plot_color, 'LineWidth', 2);
fill([time,fliplr(time)], ...
     [(unexpected_CI95(1,:) + med_unexpected_reward), fliplr((unexpected_CI95(2,:) + med_unexpected_reward))], ...
     plot_color, 'EdgeColor','none', 'FaceAlpha',0.25);
xline(0, '--r', 'LineWidth', 1)
xlim([-3.1, 3.1]);
xlabel('Time (s)');
ylabel('Avg z-score');
title(sprintf('Unexpected Reward (%d trials)', size(unexpected_reward,1)));
grid on;
all_axes(1) = ax;
current_ylim = ylim; y_min = min(y_min, current_ylim(1)); y_max = max(y_max, current_ylim(2));
hold off

subplot(2,2,2)
hold on
ax = gca; ax.FontSize = 15;
plot(time, med_unexpected_non, 'color', plot_color, 'LineWidth', 2);
fill([time,fliplr(time)], ...
     [(unexpected_CI95_non(1,:) + med_unexpected_non), fliplr((unexpected_CI95_non(2,:) + med_unexpected_non))], ...
     plot_color, 'EdgeColor','none', 'FaceAlpha',0.25);
xline(0, '--r', 'LineWidth', 1)
xlim([-3.1, 3.1]);
xlabel('Time (s)');
ylabel('Avg z-score');
title(sprintf('Unexpected Non-reward (%d trials)', size(unexpected_non,1)));
grid on;
all_axes(2) = ax;
current_ylim = ylim; y_min = min(y_min, current_ylim(1)); y_max = max(y_max, current_ylim(2));
hold off

subplot(2,2,3)
hold on
ax = gca; ax.FontSize = 15;
plot(time, med_expected_reward, 'color', plot_color, 'LineWidth', 2);
fill([time,fliplr(time)], ...
     [(expected_CI95(1,:) + med_expected_reward), fliplr((expected_CI95(2,:) + med_expected_reward))], ...
     plot_color, 'EdgeColor','none', 'FaceAlpha',0.25);
xline(0, '--r', 'LineWidth', 1)
xlim([-3.1, 3.1]);
xlabel('Time (s)');
ylabel('Avg z-score');
title(sprintf('Expected Reward (%d trials)', size(expected_reward,1)));
grid on;
all_axes(3) = ax;
current_ylim = ylim; y_min = min(y_min, current_ylim(1)); y_max = max(y_max, current_ylim(2));
hold off

subplot(2,2,4)
hold on
ax = gca; ax.FontSize = 15;
plot(time, med_expected_non, 'color', plot_color, 'LineWidth', 2);
fill([time,fliplr(time)], ...
     [(expected_CI95_non(1,:) + med_expected_non), fliplr((expected_CI95_non(2,:) + med_expected_non))], ...
     plot_color, 'EdgeColor','none', 'FaceAlpha',0.25);
xline(0, '--r', 'LineWidth', 1)
xlim([-3.1, 3.1]);
xlabel('Time (s)');
ylabel('Avg z-score');
title(sprintf('Expected Non-reward (%d trials)', size(expected_non,1)));
grid on;
all_axes(4) = ax;
current_ylim = ylim; y_min = min(y_min, current_ylim(1)); y_max = max(y_max, current_ylim(2));
hold off
set(gcf,'position',[200,200,1600,840])

% Set all y-limits the same
for i = 1:4
    ylim(all_axes(i), [y_min y_max]);
end
linkaxes(all_axes, 'xy');

% Optional save
saveas(gcf, [saveLoc, filesep, 'expected_unexpected_striatum_N11_N17','.png'], 'png');


%% ripples around events
% Parameters
window_sec = 2; % time window before and after event
bin_size = 0.1; % bin size for histogram (100 ms)
bins = -window_sec:bin_size:window_sec;

% Select event types for plotting
eventTypes = {'expected_reward', 'unexpected_reward', 'expected_non', 'unexpected_non'};
eventLabels = {'Expected Reward', 'Unexpected Reward', 'Expected Non-reward', 'Unexpected Non-reward'};
eventData = {expected_reward(:,1), unexpected_reward(:,1), expected_non(:,1), unexpected_non(:,1)};

figure;
for i = 1:length(eventTypes)
    event_times = eventData{i};
    
    % Collect ripple times relative to each event for raster plot
    relRippleTimes = [];
    for ev = 1:length(event_times)
        % Ripples within window around this event
        ripplesAround = ripple_times(ripple_times >= event_times(ev)-window_sec & ripple_times <= event_times(ev)+window_sec);
        % Relative ripple times
        relTimes = ripplesAround - event_times(ev);
        relRippleTimes = [relRippleTimes; [relTimes, ev*ones(size(relTimes))]];
    end
    
    % Plot raster of ripple events
    subplot(length(eventTypes), 2, (i-1)*2 + 1);
    if ~isempty(relRippleTimes)
        plot(relRippleTimes(:,1), relRippleTimes(:,2), 'k.', 'MarkerSize', 10);
    end
    xlim([-window_sec window_sec]);
    ylim([0 length(event_times)+1]);
    xlabel('Time from event (s)');
    ylabel('Trial');
    title([eventLabels{i} ' Ripple Raster']);
    hold on;
    plot([0 0], ylim, 'r--'); % event time line
    
    % Calculate histogram (ripple counts per bin, normalized by number of events)
    ripple_counts = histcounts(relRippleTimes(:,1), bins);
    ripple_rate = ripple_counts / length(event_times) / bin_size; % ripples/sec
    
    % Plot histogram of ripple rate
    subplot(length(eventTypes), 2, (i-1)*2 + 2);
    bar(bins(1:end-1) + bin_size/2, ripple_rate, 'k');
    xlim([-window_sec window_sec]);
    xlabel('Time from event (s)');
    ylabel('Ripple rate (Hz)');
    title([eventLabels{i} ' Ripple Rate']);
    hold on;
    plot([0 0], ylim, 'r--'); % event time line
end

figure;
hold on;

% Plot ripple times
scatter(ripple_times, ones(size(ripple_times)), 10, 'k', 'filled'); % black dots for ripples

% Overlay event times as vertical lines with different colors
y_limits = [0.8 1.2];

% Plot each event type
eventColorMap = containers.Map(...
    {'expected_reward', 'unexpected_reward', 'expected_non', 'unexpected_non'}, ...
    {[0 0.6 0], [1 0 1], [0 0 1], [1 0 0]}); % green, magenta, blue, red

eventTypeVars = {expected_reward, unexpected_reward, expected_non, unexpected_non};
eventNames = {'expected_reward', 'unexpected_reward', 'expected_non', 'unexpected_non'};

for i = 1:length(eventNames)
    ev_times = eventTypeVars{i}(:,1);
    for t = 1:length(ev_times)
        plot([ev_times(t), ev_times(t)], y_limits, '-', 'Color', eventColorMap(eventNames{i}), 'LineWidth', 1);
    end
end

% Formatting
ylim(y_limits);
xlabel('Time (s)');
yticks([]); % remove y-axis labels
title('Ripple Events and Behavioral Event Timeline');

legend({'Ripples', 'Expected Reward', 'Unexpected Reward', 'Expected Non-reward', 'Unexpected Non-reward'}, ...
    'Location', 'northoutside', 'Orientation', 'horizontal');
