
%% unexpected vs expected rewards - HPC
%% combine sessions

window = 3.1; % seconds
sampling_rate = 130;
samples = window * sampling_rate;

time = linspace(-window, window, ((samples*2)+1));

% unexpected reward
norm_session1 = zscore_matrix_unexpected_18;
norm_session2 = zscore_matrix_unexpected_19;
norm_session3 = zscore_matrix_unexpected_22;

norm_sessions_unexpected = {norm_session1, norm_session2, norm_session3};


% unexpected no reward
norm_session1 = zscore_matrix_unexpected_non_18;
norm_session2 = zscore_matrix_unexpected_non_19;
norm_session3 = zscore_matrix_unexpected_non_22;

norm_sessions_unexpected_non = {norm_session1, norm_session2, norm_session3};

% expected reward
norm_session1 = zscore_expected_18;
norm_session2 = zscore_expected_19;
norm_session3 = zscore_expected_22;

norm_sessions_expected = {norm_session1, norm_session2, norm_session3};

% expected no reward
norm_session1 = zscore_expected_non_18;
norm_session2 = zscore_expected_non_19;
norm_session3 = zscore_expected_non_22;

norm_sessions_expected_non = {norm_session1, norm_session2, norm_session3};

%% PLOT 

figure('color','white');
subplot(2, 2, 1)
hold on
ax = gca;
ax.FontSize = 15;
 
% unexpected
for ii=1:length(norm_sessions_unexpected)
    current_norm = norm_sessions_unexpected{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 
    med_z_reward = median(current_norm, 1); % median at each timepoint

    cmap = summer(4);
    plot(time, med_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+med_z_reward),fliplr((reward_CI95(2,:)+med_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title(['Z-score Around Rewards in Low Patch'], [num2str(size(norm_sessions_unexpected{1}, 1)), ...
        ' events, ', num2str(size(norm_sessions_unexpected{2}, 1)), ' events, ', ...
        num2str(size(norm_sessions_unexpected{3}, 1)), ' events ']);
end
ylim([-1.5, 2]);
xlim([-3.1,3.1]);
grid on;
hold off

% unexpected nonreward
subplot(2, 2, 2)
hold on
ax = gca;
ax.FontSize = 15;

for ii=1:length(norm_sessions_unexpected_non)
    current_norm = norm_sessions_unexpected_non{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 
    med_z_reward = median(current_norm, 1); % median at each timepoint

    cmap = summer(4);
    plot(time, med_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+med_z_reward),fliplr((reward_CI95(2,:)+med_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title(['Z-score Around Non-Rewarded Licks in High Patch'], [num2str(size(norm_sessions_unexpected_non{1}, 1)), ...
        ' events, ', num2str(size(norm_sessions_unexpected_non{2}, 1)), ' events, ', ...
        num2str(size(norm_sessions_unexpected_non{3}, 1)), ' events ']);
end
ylim([-1.5, 2]);
xlim([-3.1,3.1]);
grid on;
hold off

subplot(2, 2, 3)
hold on
ax = gca;
ax.FontSize = 15;
 
% unexpected reward
for ii=1:length(norm_sessions_expected)
    current_norm = norm_sessions_expected{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 
    med_z_reward = median(current_norm, 1); % median at each timepoint

    cmap = summer(4);
    plot(time, med_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+med_z_reward),fliplr((reward_CI95(2,:)+med_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title(['Z-score Around Rewards in High Patch'], [num2str(size(norm_sessions_expected{1}, 1)), ...
        ' events, ', num2str(size(norm_sessions_expected{2}, 1)), ' events, ', ...
        num2str(size(norm_sessions_expected{3}, 1)), ' events ']);
end
ylim([-1.5, 2]);
xlim([-3.1,3.1]);
grid on;
hold off

% expected nonreward
subplot(2, 2, 4)
hold on
ax = gca;
ax.FontSize = 15;

for ii=1:length(norm_sessions_expected_non)
    current_norm = norm_sessions_expected_non{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 
    med_z_reward = median(current_norm, 1); % median at each timepoint

    cmap = summer(4);
    plot(time, med_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+med_z_reward),fliplr((reward_CI95(2,:)+med_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title(['Z-score Around Non-Rewarded Licks in Low Patch'], [num2str(size(norm_sessions_expected_non{1}, 1)), ...
        ' events, ', num2str(size(norm_sessions_expected_non{2}, 1)), ' events, ', ...
        num2str(size(norm_sessions_expected_non{3}, 1)), ' events ']);
end
ylim([-1.5, 2]);
xlim([-3.1,3.1]);
grid on;
hold off




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% unexpected vs expected rewards - striatum
%% combine sessions

window = 3.1; % seconds
sampling_rate = 130;
samples = window * sampling_rate;

time = linspace(-window, window, ((samples*2)+1));

% unexpected reward
norm_session1 = zscore_matrix_unexpected_16;
norm_session2 = zscore_matrix_unexpected_17;
norm_session3 = zscore_matrix_unexpected_20;

norm_sessions_unexpected = {norm_session1, norm_session2, norm_session3};


% unexpected no reward
norm_session1 = zscore_matrix_unexpected_non_16;
norm_session2 = zscore_matrix_unexpected_non_17;
norm_session3 = zscore_matrix_unexpected_non_20;

norm_sessions_unexpected_non = {norm_session1, norm_session2, norm_session3};

% expected reward
norm_session1 = zscore_expected_16;
norm_session2 = zscore_expected_17;
norm_session3 = zscore_expected_20;

norm_sessions_expected = {norm_session1, norm_session2, norm_session3};

% expected no reward
norm_session1 = zscore_expected_non_16;
norm_session2 = zscore_expected_non_17;
norm_session3 = zscore_expected_non_20;

norm_sessions_expected_non = {norm_session1, norm_session2, norm_session3};

%% PLOT 

figure('color','white');
subplot(2, 2, 1)
hold on
ax = gca;
ax.FontSize = 15;
 
% unexpected
for ii=1:length(norm_sessions_unexpected)
    current_norm = norm_sessions_unexpected{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 
    med_z_reward = median(current_norm, 1); % median at each timepoint

    cmap = spring(4);
    plot(time, med_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+med_z_reward),fliplr((reward_CI95(2,:)+med_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title(['Z-score Around Rewards in Low Patch'], [num2str(size(norm_sessions_unexpected{1}, 1)), ...
        ' events, ', num2str(size(norm_sessions_unexpected{2}, 1)), ' events, ', ...
        num2str(size(norm_sessions_unexpected{3}, 1)), ' events ']);
end
ylim([-1.5, 3.2]);
xlim([-3.1,3.1]);
grid on;
hold off

% unexpected nonreward
subplot(2, 2, 2)
hold on
ax = gca;
ax.FontSize = 15;

for ii=1:length(norm_sessions_unexpected_non)
    current_norm = norm_sessions_unexpected_non{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 
    med_z_reward = median(current_norm, 1); % median at each timepoint

    cmap = spring(4);
    plot(time, med_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+med_z_reward),fliplr((reward_CI95(2,:)+med_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title(['Z-score Around Non-Rewarded Licks in High Patch'], [num2str(size(norm_sessions_unexpected_non{1}, 1)), ...
        ' events, ', num2str(size(norm_sessions_unexpected_non{2}, 1)), ' events, ', ...
        num2str(size(norm_sessions_unexpected_non{3}, 1)), ' events ']);
end
ylim([-1.5, 3.2]);
xlim([-3.1,3.1]);
grid on;
hold off

subplot(2, 2, 3)
hold on
ax = gca;
ax.FontSize = 15;
 
% unexpected reward
for ii=1:length(norm_sessions_expected)
    current_norm = norm_sessions_expected{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 
    med_z_reward = median(current_norm, 1); % median at each timepoint

    cmap = spring(4);
    plot(time, med_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+med_z_reward),fliplr((reward_CI95(2,:)+med_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title(['Z-score Around Rewards in High Patch'], [num2str(size(norm_sessions_expected{1}, 1)), ...
        ' events, ', num2str(size(norm_sessions_expected{2}, 1)), ' events, ', ...
        num2str(size(norm_sessions_expected{3}, 1)), ' events ']);
end
ylim([-1.5, 3.2]);
xlim([-3.1,3.1]);
grid on;
hold off

% expected nonreward
subplot(2, 2, 4)
hold on
ax = gca;
ax.FontSize = 15;

for ii=1:length(norm_sessions_expected_non)
    current_norm = norm_sessions_expected_non{ii};
    N = height(current_norm);                          % Number of eExperimentsn In Data Set
    avg_z_reward = mean(current_norm, 1);              % Mean Of All Experiments At Each Value Of ,x 
    reward_SEM = std(current_norm, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    reward_CI95 = bsxfun(@times, reward_SEM, CI95(:)); 
    med_z_reward = median(current_norm, 1); % median at each timepoint

    cmap = spring(4);
    plot(time, med_z_reward, 'color', cmap(ii, :), 'LineWidth', 2);
    fill([time,fliplr(time)], [(reward_CI95(1,:)+med_z_reward),fliplr((reward_CI95(2,:)+med_z_reward))], cmap(ii, :), 'EdgeColor','none', 'FaceAlpha',0.25)
    xline(0, '--r', 'LineWidth', 1)
    xlabel('time (s)');
    ylabel('avg z-score');
    title(['Z-score Around Non-Rewarded Licks in Low Patch'], [num2str(size(norm_sessions_expected_non{1}, 1)), ...
        ' events, ', num2str(size(norm_sessions_expected_non{2}, 1)), ' events, ', ...
        num2str(size(norm_sessions_expected_non{3}, 1)), ' events ']);
end
ylim([-1.5, 3.2]);
xlim([-3.1,3.1]);
grid on;
hold off