% z score across sleep sessions
%
% USAGE
%    z score photometry data from across sleep sessions and plot the 
%    dopamine traces around ripples, in the HPC and striatum
% 
%
% INPUTS 
%    need to load in the z score matrices around ripples, and then update code to reflect proper session
%    numbers
%
%    =========================================================================

%% HPC

%% combine sessions

wndow = 5; % seconds
sampling_rate = 130;
samples = wndow * sampling_rate;
baseline_idx = 1:520;

time = linspace(-wndow, wndow, ((samples*2)+1));

% Preallocate matrices for normalized data
norm_session1 = zeros(size(ripple_matrix_sleep_18));
norm_session2 = zeros(size(ripple_matrix_sleep_19));
norm_session3 = zeros(size(ripple_matrix_sleep_22));
norm_session4 = zeros(size(ripple_matrix_sleep_HPC_23));
%norm_session5 = zeros(size(ripple_matrix_sleep_24));
norm_session6 = zeros(size(ripple_matrix_sleep_26));
norm_session7 = zeros(size(ripple_matrix_sleep_28));
norm_session8 = zeros(size(ripple_matrix_sleep_30));
% N17
norm_session9 = zeros(size(ripple_matrix_pre_HPC_15));
norm_session10 = zeros(size(ripple_matrix_post_HPC_15));
norm_session11 = zeros(size(ripple_matrix_pre_HPC_21));
norm_session12 = zeros(size(ripple_matrix_post_HPC_21));


% Normalize each trial in each session
sessions = {ripple_matrix_sleep_18, ripple_matrix_sleep_19, ripple_matrix_sleep_22, ... 
    ripple_matrix_sleep_HPC_23, ripple_matrix_sleep_26, ... 
    ripple_matrix_sleep_28, ripple_matrix_sleep_30, ripple_matrix_pre_HPC_15, ...
    ripple_matrix_post_HPC_15, ripple_matrix_pre_HPC_21, ripple_matrix_post_HPC_21};
norm_sessions = {norm_session1, norm_session2, norm_session3, norm_session4,... 
    norm_session6, norm_session7, norm_session8, norm_session9, norm_session10,...
    norm_session11, norm_session12};


for s = 1:length(sessions)
    current_session = sessions{s};
   
    for trial = 1:size(current_session,1)
        % Extract baseline for current trial
        baseline = current_session(trial, :);
       
        % Calculate baseline mean and std
        baseline_mean = mean(baseline);
        baseline_std = std(baseline);
       
        % Z-score the entire trial using the baseline stats
        norm_sessions{s}(trial, :) = (current_session(trial, :) - baseline_mean) / baseline_std;
    end
end



%% PLOT

figure('color','white');
hold on
ax = gca;
ax.FontSize = 15;

for ii=1:length(norm_sessions)
    current_norm = norm_sessions{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_ripple = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    ripple_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    ripple_CI95 = bsxfun(@times, ripple_SEM, CI95(:)); 
    median_ripple = median(current_norm, 1); % median at each timepoint

    cmap = summer(length(sessions)+1);
    plot(time, median_ripple, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(ripple_CI95(1,:)+median_ripple),fliplr((ripple_CI95(2,:)+median_ripple))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title(['Average Z-score Around Ripples'], [num2str(size(norm_sessions{1}, 1)), ...
        ' ripples, ', num2str(size(norm_sessions{2}, 1)), ' ripples, ', ...
        num2str(size(norm_sessions{3}, 1)), ' ripples, ', num2str(size(norm_sessions{4}, 1)), ' ripples, ',... 
        num2str(size(norm_sessions{5}, 1)), ' ripples, ', num2str(size(norm_sessions{6}, 1)), ' ripples, ', ...
        num2str(size(norm_sessions{7}, 1)), ' ripples, ']);%, num2str(size(norm_sessions{8}, 1)), ' ripples, ']);
end
grid on;
hold off

% smoothed
figure('color','white');
hold on
ax = gca;
ax.FontSize = 15;

for ii=1:length(norm_sessions)
    current_norm = norm_sessions{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_ripple = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    ripple_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    ripple_CI95 = bsxfun(@times, ripple_SEM, CI95(:)); 
    median_ripple = median(current_norm, 1); % median at each timepoint

    cmap = summer(length(sessions)+1);
    smoothed = smoothdata(median_ripple);
    plot(time, smoothed, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(ripple_CI95(1,:)+smoothed),fliplr((ripple_CI95(2,:)+smoothed))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title(['Average Z-score Around Ripples'], [num2str(size(norm_sessions{1}, 1)), ...
        ' ripples, ', num2str(size(norm_sessions{2}, 1)), ' ripples, ', ...
        num2str(size(norm_sessions{3}, 1)), ' ripples, ', num2str(size(norm_sessions{4}, 1)), ' ripples, ',... 
        num2str(size(norm_sessions{5}, 1)), ' ripples, ', num2str(size(norm_sessions{6}, 1)), ' ripples, ', ...
        num2str(size(norm_sessions{7}, 1)), ' ripples, ']);%, num2str(size(norm_sessions{8}, 1)), ' ripples, ']);
end
grid on;
hold off

figure('color','white');
hold on
ax = gca;
ax.FontSize = 15;

for ii=1:length(norm_sessions)
    current_norm = norm_sessions{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_ripple = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    ripple_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    ripple_CI95 = bsxfun(@times, ripple_SEM, CI95(:)); 
    median_ripple = median(current_norm, 1); % median at each timepoint

    cmap = summer(length(sessions)+1);
    smoothed = smoothdata(median_ripple);
    plot(time, smoothed, 'color', cmap(ii, :), 'LineWidth', 2);
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title(['Average Z-score Around Ripples'], [num2str(size(norm_sessions{1}, 1)), ...
        ' ripples, ', num2str(size(norm_sessions{2}, 1)), ' ripples, ', ...
        num2str(size(norm_sessions{3}, 1)), ' ripples, ', num2str(size(norm_sessions{4}, 1)), ' ripples, ',... 
        num2str(size(norm_sessions{5}, 1)), ' ripples, ', num2str(size(norm_sessions{6}, 1)), ' ripples, ', ...
        num2str(size(norm_sessions{7}, 1)), ' ripples, ']);%, num2str(size(norm_sessions{8}, 1)), ' ripples, ']);
end
grid on;
hold off

%% striatum


% Preallocate matrices for normalized data
norm_session1 = zeros(size(ripple_matrix_sleep_16));
norm_session2 = zeros(size(ripple_matrix_sleep_17));
norm_session3 = zeros(size(ripple_matrix_sleep_20));
norm_session4 = zeros(size(ripple_matrix_sleep_striatum_23));
% N17
norm_session5 = zeros(size(ripple_matrix_pre_striatum_15));
norm_session6 = zeros(size(ripple_matrix_post_striatum_15));
norm_session7 = zeros(size(ripple_matrix_pre_striatum_21));
norm_session8 = zeros(size(ripple_matrix_post_striatum_21));


% Normalize each trial in each session
% sessions = {ripple_matrix_sleep_16, ripple_matrix_sleep_17, ripple_matrix_sleep_20, ripple_matrix_sleep_striatum_23};
% norm_sessions = {norm_session1, norm_session2, norm_session3, norm_session4};
sessions = {ripple_matrix_sleep_16, ripple_matrix_sleep_17, ripple_matrix_sleep_striatum_23, ...
    ripple_matrix_pre_striatum_15, ripple_matrix_post_striatum_15, ...
    ripple_matrix_pre_striatum_21, ripple_matrix_post_striatum_21};
norm_sessions = {norm_session1, norm_session2, norm_session4, norm_session5, norm_session6, ...
    norm_session7, norm_session8};

for s = 1:length(sessions)
    current_session = sessions{s};

    for trial = 1:size(current_session,1)
        % Extract baseline for current trial
        baseline = current_session(trial, :);

        % Calculate baseline mean and std
        baseline_mean = mean(baseline);
        baseline_std = std(baseline);

        % Z-score the entire trial using the baseline stats
        norm_sessions{s}(trial, :) = (current_session(trial, :) - baseline_mean) / baseline_std;

        %  % Extract baseline for current trial
        % baseline = current_session(trial, :);
        % 
        % % Z-score the entire trial using the baseline stats
        % norm_sessions{s}(trial, :) = baseline;
    end
end



%% PLOT

figure('color','white');
hold on
ax = gca;
ax.FontSize = 15;

for ii=1:length(norm_sessions)
    current_norm = norm_sessions{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_ripple = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    ripple_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    ripple_CI95 = bsxfun(@times, ripple_SEM, CI95(:));
    median_ripple = median(current_norm, 1); % median at each timepoint

    cmap = spring(length(sessions)+1);
    smoothed = smoothdata(median_ripple);
    plot(time, median_ripple, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(ripple_CI95(1,:)+median_ripple),fliplr((ripple_CI95(2,:)+median_ripple))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title(['Average Z-score Around Ripples'], [num2str(size(norm_sessions{1}, 1)), ...
        ' ripples, ', num2str(size(norm_sessions{2}, 1)), ' ripples, ', num2str(size(norm_sessions{3}, 1)), ' ripples, ']);

end
grid on;
hold off

% smoothed
figure('color','white');
hold on
ax = gca;
ax.FontSize = 15;

for ii=1:length(norm_sessions)
    current_norm = norm_sessions{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_ripple = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    ripple_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    ripple_CI95 = bsxfun(@times, ripple_SEM, CI95(:));
    median_ripple = median(current_norm, 1); % median at each timepoint

    cmap = spring(length(sessions)+1);
    plot(time, median_ripple, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(ripple_CI95(1,:)+median_ripple),fliplr((ripple_CI95(2,:)+median_ripple))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title(['Average Z-score Around Ripples'], [num2str(size(norm_sessions{1}, 1)), ...
        ' ripples, ', num2str(size(norm_sessions{2}, 1)), ' ripples, ', num2str(size(norm_sessions{3}, 1)), ' ripples, ']);

end
grid on;
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% awake ripples
%% combine sessions

wndow = 5; % seconds
sampling_rate = 130;
samples = wndow * sampling_rate;

time = linspace(-wndow, wndow, ((samples*2)+1));

% Preallocate matrices for normalized data
norm_session1 = zeros(size(ripple_matrix_18));
norm_session2 = zeros(size(ripple_matrix_19));
norm_session3 = zeros(size(ripple_matrix_22));

% Normalize each trial in each session
sessions = {ripple_matrix_18, ripple_matrix_19, ripple_matrix_22};
norm_sessions_awake_ripple = {norm_session1, norm_session2, norm_session3};

for s = 1:length(sessions)
    current_session = sessions{s};
   
    for trial = 1:size(current_session,1)
        % Extract baseline for current trial
        baseline = current_session(trial, :);
       
        % Calculate baseline mean and std
        baseline_mean = mean(baseline);
        baseline_std = std(baseline);
       
        % Z-score the entire trial using the baseline stats
        norm_sessions_awake_ripple{s}(trial, :) = (current_session(trial, :) - baseline_mean) / baseline_std;
    end
end



%% PLOT

figure('color','white');
hold on
ax = gca;
ax.FontSize = 15;

for ii=1:length(norm_sessions_awake_ripple)
    current_norm = norm_sessions_awake_ripple{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_ripple = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    ripple_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    ripple_CI95 = bsxfun(@times, ripple_SEM, CI95(:)); 
    median_ripple = median(current_norm, 1); % median at each timepoint

    cmap = summer(5);
    plot(time, median_ripple, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(ripple_CI95(1,:)+median_ripple),fliplr((ripple_CI95(2,:)+median_ripple))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title(['Average Z-score Around Ripples'], [num2str(size(norm_sessions_awake_ripple{1}, 1)), ...
        ' ripples, ', num2str(size(norm_sessions_awake_ripple{2}, 1)), ' ripples, ', ...
        num2str(size(norm_sessions_awake_ripple{3}, 1)), ' ripples ']);
end
grid on;
hold off


%% striatum


% Preallocate matrices for normalized data
norm_session1 = zeros(size(ripple_matrix_16));
norm_session2 = zeros(size(ripple_matrix_17));
norm_session3 = zeros(size(ripple_matrix_20));

% Normalize each trial in each session
sessions = {ripple_matrix_16, ripple_matrix_17, ripple_matrix_20};
norm_sessions_awake_ripple = {norm_session1, norm_session2, norm_session3};

for s = 1:length(sessions)
    current_session = sessions{s};
   
    for trial = 1:size(current_session,1)
        % Extract baseline for current trial
        baseline = current_session(trial, :);
       
        % Calculate baseline mean and std
        baseline_mean = mean(baseline);
        baseline_std = std(baseline);
       
        % Z-score the entire trial using the baseline stats
        norm_sessions_awake_ripple{s}(trial, :) = (current_session(trial, :) - baseline_mean) / baseline_std;
    end
end



%% PLOT

figure('color','white');
hold on
ax = gca;
ax.FontSize = 15;

for ii=1:length(norm_sessions_awake_ripple)
    current_norm = norm_sessions_awake_ripple{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_ripple = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    ripple_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    ripple_CI95 = bsxfun(@times, ripple_SEM, CI95(:));
    median_ripple = median(current_norm, 1); % median at each timepoint

    cmap = spring(4);
    plot(time, median_ripple, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(ripple_CI95(1,:)+median_ripple),fliplr((ripple_CI95(2,:)+median_ripple))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title(['Average Z-score Around Ripples'], [num2str(size(norm_sessions_awake_ripple{1}, 1)), ...
        ' ripples, ', num2str(size(norm_sessions_awake_ripple{2}, 1)), ' ripples, ', num2str(size(norm_sessions_awake_ripple{3}, 1)), ' ripples ']);

end
grid on;
hold off
