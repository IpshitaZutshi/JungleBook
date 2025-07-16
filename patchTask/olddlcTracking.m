
function [DLC_tracking] = dlcTracking(varargin)

p = inputParser;
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'plotfig',true,@islogical);
addParameter(p,'forceRun',true,@islogical);
%addParameter(p,'updatedIntan',true,@islogical);

parse(p,varargin{:});
saveMat = p.Results.saveMat;
plotfig = p.Results.plotfig;
forceRun = p.Results.forceRun;
%updatedIntan = p.Results.updatedIntan;

basepath = pwd;
[~, currentFolderName] = fileparts(basepath);


tracking_file = dir(fullfile(basepath, '*Tracking.Behavior.mat'));
load(tracking_file.name);

behavior_file = dir(fullfile(basepath, '*TrialBehavior.Events.mat'));
load(behavior_file.name);

j=2;%76570
while j < size(tracking.position.x, 1)
    if tracking.position.x(j) >= tracking.position.x(j-1)+0.9 || tracking.position.y(j) >= tracking.position.y(j-1)+0.9 || ...
            tracking.position.x(j) <= tracking.position.x(j-1)-0.9 || tracking.position.y(j) <= tracking.position.y(j-1)-0.9
        start_cut = j;
        
        k=j+1;
        while k < size(tracking.position.x, 1)
            % middle section of track
            if tracking.position.y(k) <= 110 && tracking.position.x(k) > 92 ... 
                    && tracking.position.x(k) < 98 && tracking.position.y(k) >= 35 % these are rough numbers -> and only for this session
                break;  % End of bad block (jump back)

            % T section of track
            elseif tracking.position.y(k) > 110 && tracking.position.y(k) < 123 ...
                    && tracking.position.x(k) <= 108 && tracking.position.x(k) >= 79
                break;  % End of bad block (jump back)

            % home
            elseif tracking.position.y(k) < 35 && tracking.position.x(k) > 90 ... 
                    && tracking.position.x(k) < 97
                break;
            else
                k = k + 1;
            end
        end

        end_cut = k - 1;

        % Delete the bad block
        tracking.position.x(start_cut:end_cut, :) = [];
        tracking.position.y(start_cut:end_cut, :) = [];
       % cut_rows = [cut_rows; (start_cut:end_cut)'];

        % Adjust j after deletion
        j = start_cut;
        j = j + 1;
    else
        j = j + 1;
    end
end


% manual deletion

figure('color','white');
colormap(jet)
scatter(tracking.position.x, tracking.position.y, 20, tracking.position.x, 'filled') % 20 is point size
colorbar




%% Look at trajectories

forward_tracking.x = [];
forward_tracking.y = [];
forward_tracking.direction = nan(size(tracking.position.x, 1), 1);
forward_tracking.trial = nan(size(tracking.position.x, 1), 1);

num_trials = length(behavTrials.choiceTS);

for i = 1:num_trials
    curr_trial_tracking_x = tracking.position.x(tracking.timestamps >= behavTrials.timestamps(i, 1) ... 
        & tracking.timestamps <= behavTrials.choiceTS(i));
    curr_trial_tracking_y = tracking.position.y(tracking.timestamps >= behavTrials.timestamps(i, 1) ... 
        & tracking.timestamps <= behavTrials.choiceTS(i));

    forward_tracking.x = [forward_tracking.x; curr_trial_tracking_x];
    forward_tracking.y = [forward_tracking.y; curr_trial_tracking_y];

    forward_tracking.direction(tracking.timestamps >= behavTrials.timestamps(i, 1) ... 
        & tracking.timestamps <= behavTrials.choiceTS(i)) = behavTrials.choice(i);
    forward_tracking.trial(tracking.timestamps >= behavTrials.timestamps(i, 1) ... 
        & tracking.timestamps <= behavTrials.choiceTS(i)) = i;

end

forward_tracking.direction = forward_tracking.direction(~isnan(forward_tracking.direction));
forward_tracking.trial = forward_tracking.trial(~isnan(forward_tracking.trial));

figure('color','white');
colormap(cool)
scatter(forward_tracking.x, forward_tracking.y, 20, forward_tracking.direction, 'filled') % 20 is point size
colorbar


%% filter correct vs incorrect

corr_forward_tracking.x = [];
corr_forward_tracking.y = [];
corr_forward_tracking.direction = nan(size(tracking.position.x, 1), 1);
corr_forward_tracking.trial = nan(size(tracking.position.x, 1), 1);
trial_lens = [];
left_trials = {};
right_trials = {};

incorr_forward_tracking.x = [];
incorr_forward_tracking.y = [];
incorr_forward_tracking.direction = nan(size(tracking.position.x, 1), 1);
incorr_forward_tracking.trial = nan(size(tracking.position.x, 1), 1);

num_trials = length(behavTrials.choiceTS);

for i = 1:num_trials
    if behavTrials.correct(i) == 1
        curr_trial_tracking_x = tracking.position.x(tracking.timestamps >= behavTrials.timestamps(i, 1) ... 
            & tracking.timestamps <= behavTrials.choiceTS(i));
        curr_trial_tracking_y = tracking.position.y(tracking.timestamps >= behavTrials.timestamps(i, 1) ... 
            & tracking.timestamps <= behavTrials.choiceTS(i));
    
        corr_forward_tracking.x = [corr_forward_tracking.x; curr_trial_tracking_x];
        corr_forward_tracking.y = [corr_forward_tracking.y; curr_trial_tracking_y];
    
        corr_forward_tracking.direction(tracking.timestamps >= behavTrials.timestamps(i, 1) ... 
            & tracking.timestamps <= behavTrials.choiceTS(i)) = behavTrials.choice(i);
        corr_forward_tracking.trial(tracking.timestamps >= behavTrials.timestamps(i, 1) ... 
            & tracking.timestamps <= behavTrials.choiceTS(i)) = i;

        trial_lens = [trial_lens; length(curr_trial_tracking_x)];

        if behavTrials.choice(i) == 1
            left_trials{end+1} = [curr_trial_tracking_x(:), curr_trial_tracking_y(:)];
        else
            right_trials{end+1} = [curr_trial_tracking_x(:), curr_trial_tracking_y(:)];
        end
    else
        curr_trial_tracking_x = tracking.position.x(tracking.timestamps >= behavTrials.timestamps(i, 1) ... 
            & tracking.timestamps <= behavTrials.choiceTS(i));
        curr_trial_tracking_y = tracking.position.y(tracking.timestamps >= behavTrials.timestamps(i, 1) ... 
            & tracking.timestamps <= behavTrials.choiceTS(i));
    
        incorr_forward_tracking.x = [incorr_forward_tracking.x; curr_trial_tracking_x];
        incorr_forward_tracking.y = [incorr_forward_tracking.y; curr_trial_tracking_y];
    
        incorr_forward_tracking.direction(tracking.timestamps >= behavTrials.timestamps(i, 1) ... 
            & tracking.timestamps <= behavTrials.choiceTS(i)) = behavTrials.choice(i);
        incorr_forward_tracking.trial(tracking.timestamps >= behavTrials.timestamps(i, 1) ... 
            & tracking.timestamps <= behavTrials.choiceTS(i)) = i;
    end

end

corr_forward_tracking.direction = corr_forward_tracking.direction(~isnan(corr_forward_tracking.direction));
corr_forward_tracking.trial = corr_forward_tracking.trial(~isnan(corr_forward_tracking.trial));

incorr_forward_tracking.direction = incorr_forward_tracking.direction(~isnan(incorr_forward_tracking.direction));
incorr_forward_tracking.trial = incorr_forward_tracking.trial(~isnan(incorr_forward_tracking.trial));

figure('color','white');
subplot(1, 2, 1)
colormap(jet)
scatter(corr_forward_tracking.x, corr_forward_tracking.y, 20, corr_forward_tracking.direction, 'filled') % 20 is point size
title('Correct trials');
colorbar

subplot(1, 2, 2)
colormap(jet)
scatter(incorr_forward_tracking.x, incorr_forward_tracking.y, 20, incorr_forward_tracking.direction, 'filled') % 20 is point size
title('Incorrect trials');
colorbar


%% separate left vs right (all)
% left_trials = [];
% right_trials = [];
% 
% for kk = 1:num_trials
%     if behavTrials.choice(kk) == 0 
%         right_trials = [right_trials; kk];
%     elseif behavTrials.choice(kk) == 1
%         left_trials = [left_trials; kk];
%     end
% end


%% separate left vs right (correct only)


%% Find average trajectories

num_pts = 482; % number of points for interpolation

interp_trials = @(trials) cellfun(@(traj) ...
    [interp1(1:size(traj,1), traj(:,1), linspace(1, size(traj,1), num_pts), 'linear', 'extrap')', ...
     interp1(1:size(traj,1), traj(:,2), linspace(1, size(traj,1), num_pts), 'linear', 'extrap')'], ...
    trials, 'UniformOutput', false);

left_interp = interp_trials(left_trials);
right_interp = interp_trials(right_trials);

left_matrix = cat(3, left_interp{:});
right_matrix = cat(3, right_interp{:});

left_mean = mean(left_matrix, 3, 'omitnan');
right_mean = mean(right_matrix, 3, 'omitnan');

figure('color','white'); 
hold on;
plot(left_mean(:,1), left_mean(:,2), 'r-', 'LineWidth', 2); % left = red
plot(right_mean(:,1), right_mean(:,2), 'b-', 'LineWidth', 2); % right = blue
%legend('Left Trials', 'Right Trials');
title('Average Trajectories');


% calculate confidence intervals
N = size(left_matrix, 3);                          % Number of Experiments In Data Set
left_SEM = std(left_matrix, 3, 'omitnan')/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
left_CI95 = bsxfun(@times, left_SEM, CI95(:));  % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of exw
%smooth_CI95 = smoothdata(ripple_CI95, 2);

%% second approach
% Parameters
bin_size = 0.1;  % Y-axis bin size

left_avg_x = [];
left_avg_y = [];

right_avg_x = [];
right_avg_y = [];

num_pts = 482;

interp_trials = @(trials) cellfun(@(traj) ...
    [interp1(1:size(traj,1), traj(:,1), linspace(1, size(traj,1), num_pts), 'linear', 'extrap')', ...
     interp1(1:size(traj,1), traj(:,2), linspace(1, size(traj,1), num_pts), 'linear', 'extrap')'], ...
    trials, 'UniformOutput', false);

left_interp = interp_trials(left_trials);
right_interp = interp_trials(right_trials);

% Combine all y-values to find the global min and max for binning
all_left_y = cell2mat(cellfun(@(traj) traj(:, 2), left_interp, 'UniformOutput', false));
all_right_y = cell2mat(cellfun(@(traj) traj(:, 2), right_interp, 'UniformOutput', false));

min_y = min(min(min(all_left_y)), min(min(all_right_y)));
max_y = max(max(max(all_left_y)), max(max(all_right_y)));

% Define bin edges and initialize cell arrays to hold x-values
y_bin_edges = min_y:bin_size:max_y;
num_bins = length(y_bin_edges) - 1;
left_bins = cell(num_bins, 1);
right_bins = cell(num_bins, 1);


% Fill bins for LEFT trials
for t = 1:length(left_trials)
    x = left_trials{t}(:,1);
    y = left_trials{t}(:,2);
    for i = 1:num_bins
        mask = y >= y_bin_edges(i) & y < y_bin_edges(i+1);
        left_bins{i} = [left_bins{i}; x(mask)];
    end
end

% Fill bins for RIGHT trials
for t = 1:length(right_trials)
    x = right_trials{t}(:,1);
    y = right_trials{t}(:,2);
    for i = 1:num_bins
        mask = y >= y_bin_edges(i) & y < y_bin_edges(i+1);
        right_bins{i} = [right_bins{i}; x(mask)];
    end
end

% Compute mean x-values per bin
left_avg_x = cellfun(@(x) mean(x, 'omitnan'), left_bins);
right_avg_x = cellfun(@(x) mean(x, 'omitnan'), right_bins);
smoothed_left = smoothdata(left_avg_x);
smoothed_right = smoothdata(right_avg_x);

% Compute y-bin centers
y_bin_centers = (y_bin_edges(1:end-1) + y_bin_edges(2:end)) / 2;

% Plot
figure('Color', 'white'); 
hold on;
plot(smoothed_left, y_bin_centers, 'r-', 'LineWidth', 2);
plot(smoothed_right, y_bin_centers, 'b-', 'LineWidth', 2);
xlabel('X Position');
ylabel('Y Position');
legend({'Left', 'Right'});
title('Average Trajectories by Binned Y');

% Compute 95% confidence intervals for left and right
left_std = cellfun(@(x) std(x, 'omitnan'), left_bins);
right_std = cellfun(@(x) std(x, 'omitnan'), right_bins);

left_n = cellfun(@(x) sum(~isnan(x)), left_bins);
right_n = cellfun(@(x) sum(~isnan(x)), right_bins);

left_se = left_std ./ sqrt(left_n);
right_se = right_std ./ sqrt(right_n);

% 95% CI = mean ± 1.96 * SE
left_ci_upper = left_avg_x + 1.96 * left_se;
left_ci_lower = left_avg_x - 1.96 * left_se;

right_ci_upper = right_avg_x + 1.96 * right_se;
right_ci_lower = right_avg_x - 1.96 * right_se;

y_bin_centers = y_bin_centers';


% Plot
figure('Color', 'white'); 
hold on;

% Mask out NaNs for plotting left CI
valid_left = ~isnan(left_ci_upper) & ~isnan(left_ci_lower) & ~isnan(y_bin_centers);
left_ci_upper_plot = left_ci_upper(valid_left);
left_ci_lower_plot = left_ci_lower(valid_left);
left_y_plot = y_bin_centers(valid_left);

% Mask out NaNs for plotting right CI
valid_right = ~isnan(right_ci_upper) & ~isnan(right_ci_lower) & ~isnan(y_bin_centers);
right_ci_upper_plot = right_ci_upper(valid_right);
right_ci_lower_plot = right_ci_lower(valid_right);
right_y_plot = y_bin_centers(valid_right);

% Shaded confidence interval areas
valid_left = ~isnan(left_ci_upper) & ~isnan(left_ci_lower) & ~isnan(left_avg_x);
left_y = y_bin_centers(valid_left);
left_upper = left_ci_upper(valid_left);
left_lower = left_ci_lower(valid_left);
left_mean = left_avg_x(valid_left);
left_lower_smooth = smoothdata(left_lower);
left_upper_smooth = smoothdata(left_upper);

% RIGHT
valid_right = ~isnan(right_ci_upper) & ~isnan(right_ci_lower) & ~isnan(right_avg_x);
right_y = y_bin_centers(valid_right);
right_upper = right_ci_upper(valid_right);
right_lower = right_ci_lower(valid_right);
right_mean = right_avg_x(valid_right);
right_lower_smooth = smoothdata(right_lower);
right_upper_smooth = smoothdata(right_upper);

% Plot shaded confidence bands
fill([left_upper_smooth; flipud(left_lower_smooth)], [left_y; flipud(left_y)], ...
     'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([right_upper_smooth; flipud(right_lower_smooth)], [right_y; flipud(right_y)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% Mean lines
plot(smoothed_left, y_bin_centers, 'r-', 'LineWidth', 2);
plot(smoothed_right, y_bin_centers, 'b-', 'LineWidth', 2);

xlabel('X Position');
ylabel('Y Position');
legend({'Left CI', 'Right CI', 'Left Mean', 'Right Mean'}, 'Location', 'Best');
title('Average Trajectories with 95% Confidence Intervals');

%% split end portion
%% Bin size and threshold
bin_size = 0.1;
y_threshold = 112;

all_left_x = cell2mat(cellfun(@(traj) traj(:, 1), left_interp, 'UniformOutput', false));
all_right_x = cell2mat(cellfun(@(traj) traj(:, 1), right_interp, 'UniformOutput', false));
min_x = min(min(min(all_left_x)), min(min(all_right_x)));
max_x = max(max(max(all_left_x)), max(max(all_right_x)));

% Initialize bins
y_bin_edges = floor(min_y):bin_size:y_threshold;
x_bin_edges = floor(min_x):bin_size:ceil(max_x);  

y_bin_centers = (y_bin_edges(1:end-1) + y_bin_edges(2:end)) / 2;
x_bin_centers = (x_bin_edges(1:end-1) + x_bin_edges(2:end)) / 2;

% Initialize left bins
left_x_by_y = cell(length(y_bin_centers),1);
left_y_by_x = cell(length(x_bin_centers),1);

% Populate bins for left trials
for i = 1:length(left_trials)
    traj = left_trials{i};
    for j = 1:size(traj,1)
        x = traj(j,1);
        y = traj(j,2);

        if y < y_threshold
            bin_idx = find(y >= y_bin_edges(1:end-1) & y < y_bin_edges(2:end), 1);
            if ~isempty(bin_idx)
                left_x_by_y{bin_idx} = [left_x_by_y{bin_idx}; x];
            end
        else
            bin_idx = find(x >= x_bin_edges(1:end-1) & x < x_bin_edges(2:end), 1);
            if ~isempty(bin_idx)
                left_y_by_x{bin_idx} = [left_y_by_x{bin_idx}; y];
            end
        end
    end
end

% Compute means and CIs
left_mean_x_by_y = cellfun(@(v) mean(v, 'omitnan'), left_x_by_y);
left_ci_x_upper = cellfun(@(v) prctile(v, 97.5), left_x_by_y);
left_ci_x_lower = cellfun(@(v) prctile(v, 2.5), left_x_by_y);

left_mean_y_by_x = cellfun(@(v) mean(v, 'omitnan'), left_y_by_x);
left_ci_y_upper = cellfun(@(v) prctile(v, 97.5), left_y_by_x);
left_ci_y_lower = cellfun(@(v) prctile(v, 2.5), left_y_by_x);

% Repeat for right trials
right_x_by_y = cell(length(y_bin_centers),1);
right_y_by_x = cell(length(x_bin_centers),1);

for i = 1:length(right_trials)
    traj = right_trials{i};
    for j = 1:size(traj,1)
        x = traj(j,1);
        y = traj(j,2);

        if y < y_threshold
            bin_idx = find(y >= y_bin_edges(1:end-1) & y < y_bin_edges(2:end), 1);
            if ~isempty(bin_idx)
                right_x_by_y{bin_idx} = [right_x_by_y{bin_idx}; x];
            end
        else
            bin_idx = find(x >= x_bin_edges(1:end-1) & x < x_bin_edges(2:end), 1);
            if ~isempty(bin_idx)
                right_y_by_x{bin_idx} = [right_y_by_x{bin_idx}; y];
            end
        end
    end
end

right_mean_x_by_y = cellfun(@(v) mean(v, 'omitnan'), right_x_by_y);
right_ci_x_upper = cellfun(@(v) prctile(v, 97.5), right_x_by_y);
right_ci_x_lower = cellfun(@(v) prctile(v, 2.5), right_x_by_y);

right_mean_y_by_x = cellfun(@(v) mean(v, 'omitnan'), right_y_by_x);
right_ci_y_upper = cellfun(@(v) prctile(v, 97.5), right_y_by_x);
right_ci_y_lower = cellfun(@(v) prctile(v, 2.5), right_y_by_x);

%% Plotting
y_bin_centers = y_bin_centers';
x_bin_centers = x_bin_centers';

figure('Color','white'); hold on;

% Left: x vs y before turn
fill([left_ci_x_upper fliplr(left_ci_x_lower)], [y_bin_centers fliplr(y_bin_centers)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(left_mean_x_by_y, y_bin_centers, 'r-', 'LineWidth', 2);

% Left: x vs y after turn (x bins)
fill([x_bin_centers fliplr(x_bin_centers)], [left_ci_y_upper fliplr(left_ci_y_lower)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(x_bin_centers, left_mean_y_by_x, 'r-', 'LineWidth', 2);

% Right: x vs y before turn
fill([right_ci_x_upper fliplr(right_ci_x_lower)], [y_bin_centers fliplr(y_bin_centers)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(right_mean_x_by_y, y_bin_centers, 'b-', 'LineWidth', 2);

% Right: x vs y after turn (x bins)
fill([x_bin_centers fliplr(x_bin_centers)], [right_ci_y_upper fliplr(right_ci_y_lower)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(x_bin_centers, right_mean_y_by_x, 'b-', 'LineWidth', 2);

xlabel('X'); ylabel('Y');
title('Average Trajectories with Turn Handling and 95% CI');
legend('Left CI', 'Left Mean', 'Right CI', 'Right Mean');

%% Plotting
y_bin_centers = y_bin_centers(:);
x_bin_centers = x_bin_centers(:);

figure('Color','white'); 
hold on;


% --- Plot Left Trials (Y < 112): x vs y ---
valid_idx = ~isnan(left_ci_x_lower) & ~isnan(left_ci_x_upper) & ~isnan(left_mean_x_by_y);
fill([left_ci_x_upper(valid_idx); flipud(left_ci_x_lower(valid_idx))], ...
     [y_bin_centers(valid_idx); flipud(y_bin_centers(valid_idx))], ...
     'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
smooth = smoothdata(left_mean_x_by_y);
plot(smooth(valid_idx), y_bin_centers(valid_idx), 'r-', 'LineWidth', 2);

% --- Plot Left Trials (Y > 112): y vs x ---
valid_idx = ~isnan(left_ci_y_lower) & ~isnan(left_ci_y_upper) & ~isnan(left_mean_y_by_x);
fill([x_bin_centers(valid_idx); flipud(x_bin_centers(valid_idx))], ...
     [left_ci_y_upper(valid_idx); flipud(left_ci_y_lower(valid_idx))], ...
     'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
smooth = smoothdata(x_bin_centers);
plot(smooth(valid_idx), left_mean_y_by_x(valid_idx), 'r-', 'LineWidth', 2);

% --- Plot Right Trials (Y < 112): x vs y ---
valid_idx = ~isnan(right_ci_x_lower) & ~isnan(right_ci_x_upper) & ~isnan(right_mean_x_by_y);
fill([right_ci_x_upper(valid_idx); flipud(right_ci_x_lower(valid_idx))], ...
     [y_bin_centers(valid_idx); flipud(y_bin_centers(valid_idx))], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
smooth = smoothdata(right_mean_x_by_y);
plot(smooth(valid_idx), y_bin_centers(valid_idx), 'b-', 'LineWidth', 2);

% --- Plot Right Trials (Y > 112): y vs x ---
valid_idx = ~isnan(right_ci_y_lower) & ~isnan(right_ci_y_upper) & ~isnan(right_mean_y_by_x);
fill([x_bin_centers(valid_idx); flipud(x_bin_centers(valid_idx))], ...
     [right_ci_y_upper(valid_idx); flipud(right_ci_y_lower(valid_idx))], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
smooth = smoothdata(x_bin_centers);
plot(smooth(valid_idx), right_mean_y_by_x(valid_idx), 'b-', 'LineWidth', 2);

xlabel('X');
ylabel('Y');
title('Average Trajectories with Turn Handling and 95% CI');
legend('Left CI (before)', 'Left Mean (before)', ...
       'Left CI (after)', 'Left Mean (after)', ...
       'Right CI (before)', 'Right Mean (before)', ...
       'Right CI (after)', 'Right Mean (after)');

%% Smoothed
y_bin_centers = y_bin_centers(:);
x_bin_centers = x_bin_centers(:);

figure('Color','white'); 
hold on;

% --- Plot Left Trials (Y < 112): x vs y ---
valid_idx = ~isnan(left_ci_x_lower) & ~isnan(left_ci_x_upper) & ~isnan(left_mean_x_by_y);
smooth_mean = smoothdata(left_mean_x_by_y);
smooth_ci_upper = smoothdata(left_ci_x_upper);
smooth_ci_lower = smoothdata(left_ci_x_lower);

fill([smooth_ci_upper(valid_idx); flipud(smooth_ci_lower(valid_idx))], ...
     [y_bin_centers(valid_idx); flipud(y_bin_centers(valid_idx))], ...
     'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(smooth_mean(valid_idx), y_bin_centers(valid_idx), 'r-', 'LineWidth', 2);

% --- Plot Left Trials (Y > 112): y vs x ---
valid_idx = ~isnan(left_ci_y_lower) & ~isnan(left_ci_y_upper) & ~isnan(left_mean_y_by_x);
smooth_mean = smoothdata(left_mean_y_by_x);
smooth_ci_upper = smoothdata(left_ci_y_upper);
smooth_ci_lower = smoothdata(left_ci_y_lower);

fill([x_bin_centers(valid_idx); flipud(x_bin_centers(valid_idx))], ...
     [smooth_ci_upper(valid_idx); flipud(smooth_ci_lower(valid_idx))], ...
     'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(x_bin_centers(valid_idx), smooth_mean(valid_idx), 'r-', 'LineWidth', 2);

% --- Plot Right Trials (Y < 112): x vs y ---
valid_idx = ~isnan(right_ci_x_lower) & ~isnan(right_ci_x_upper) & ~isnan(right_mean_x_by_y);
smooth_mean = smoothdata(right_mean_x_by_y);
smooth_ci_upper = smoothdata(right_ci_x_upper);
smooth_ci_lower = smoothdata(right_ci_x_lower);

fill([smooth_ci_upper(valid_idx); flipud(smooth_ci_lower(valid_idx))], ...
     [y_bin_centers(valid_idx); flipud(y_bin_centers(valid_idx))], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(smooth_mean(valid_idx), y_bin_centers(valid_idx), 'b-', 'LineWidth', 2);

% --- Plot Right Trials (Y > 112): y vs x ---
valid_idx = ~isnan(right_ci_y_lower) & ~isnan(right_ci_y_upper) & ~isnan(right_mean_y_by_x);
smooth_mean = smoothdata(right_mean_y_by_x);
smooth_ci_upper = smoothdata(right_ci_y_upper);
smooth_ci_lower = smoothdata(right_ci_y_lower);

fill([x_bin_centers(valid_idx); flipud(x_bin_centers(valid_idx))], ...
     [smooth_ci_upper(valid_idx); flipud(smooth_ci_lower(valid_idx))], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(x_bin_centers(valid_idx), smooth_mean(valid_idx), 'b-', 'LineWidth', 2);

xlabel('X');
ylabel('Y');
title('Average Trajectories with Turn Handling and Smoothed 95% CI');
legend('Left CI (before)', 'Left Mean (before)', ...
       'Left CI (after)', 'Left Mean (after)', ...
       'Right CI (before)', 'Right Mean (before)', ...
       'Right CI (after)', 'Right Mean (after)');


%% at each time point
%{
% Compute p-values at each time point (x-dimension)
p_vals_x = nan(num_pts, 1);
p_vals_y = nan(num_pts, 1);

for t = 1:num_pts
    % Get x values for left and right at time t
    left_x = squeeze(left_matrix(t, 1, :)); % x at time t for all left trials 
    right_x = squeeze(right_matrix(t, 1, :)); % x at time t for all right trials

    left_y = squeeze(left_matrix(t, 2, :)); % y at time t for all left trials
    right_y = squeeze(right_matrix(t, 2, :)); % y at time t for all right trials

    % t-test
    [~, p_vals_x(t)] = ttest2(left_x, right_x, 'Vartype', 'unequal');
    [~, p_vals_y(t)] = ttest2(left_y, right_y, 'Vartype', 'unequal');
end

% plot p-values 
figure('color','white');
subplot(2,1,1);
plot(p_vals_x, 'r'); hold on;
yline(0.05, '--k');
title('p-values for x-dimension');
xlabel('Time point');
ylabel('p-value');

subplot(2,1,2);
plot(p_vals_y, 'b'); hold on;
yline(0.05, '--k');
title('p-values for y-dimension');
xlabel('Time point');
ylabel('p-value');


left_x = squeeze(left_matrix(:, 1, :)); % [time, trials]
left_y = squeeze(left_matrix(:, 2, :));
right_x = squeeze(right_matrix(:, 1, :));
right_y = squeeze(right_matrix(:, 2, :));

left_mean_x = mean(left_x, 2, 'omitnan');
%left_sem_x = std(left_x, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(left_x), 2));
left_sem_x = 1.96 * std(left_x, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(left_x), 2));

left_mean_y = mean(left_y, 2, 'omitnan');
%left_sem_y = std(left_y, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(left_y), 2));
left_sem_y = 1.96 * std(left_y, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(left_y), 2));

right_mean_x = mean(right_x, 2, 'omitnan');
%right_sem_x = std(right_x, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(right_x), 2));
right_sem_x = 1.96 * std(right_x, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(right_x), 2));

right_mean_y = mean(right_y, 2, 'omitnan');
%right_sem_y = std(right_y, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(right_y), 2));
right_sem_y = 1.96 * std(right_y, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(right_y), 2));

% Plot with confidence intervals
figure('color','white');
hold on;
left_x_shade = [left_mean_x - left_sem_x; flipud(left_mean_x + left_sem_x)];
left_y_shade = [left_mean_y - left_sem_y; flipud(left_mean_y + left_sem_y)];

% Plot shaded region
fill([left_mean_x; flipud(left_mean_x)], left_y_shade, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(left_mean_x, left_mean_y, 'r-', 'LineWidth', 2);

% Plot Right with shaded error (blue)
right_x_shade = [right_mean_x - right_sem_x; flipud(right_mean_x + right_sem_x)];
right_y_shade = [right_mean_y - right_sem_y; flipud(right_mean_y + right_sem_y)];

% Plot shaded region
fill([right_mean_x; flipud(right_mean_x)], right_y_shade, 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(right_mean_x, right_mean_y, 'b-', 'LineWidth', 2);

title('Average Trajectories with 95% CI');
xlabel('X Position');
ylabel('Y Position');
% axis equal;

%}

num_pts = 482; % number of path-based sample points
n_left = length(left_trials);
n_right = length(right_trials);

left_matrix = nan(num_pts, 2, n_left);
right_matrix = nan(num_pts, 2, n_right);

for i = 1:n_left
    traj = left_trials{i}; % traj is [N x 2] matrix of [x, y] coordinates

    % Compute cumulative path length
    d = sqrt(sum(diff(traj).^2, 2));        % distances between points
    cumdist = [0; cumsum(d)];               % cumulative distance
    total_length = cumdist(end);

    if total_length == 0 || length(cumdist) < 2
        resampled = repmat(traj(1,:), num_pts, 1); % handle degenerate cases
    else
        % Remove duplicate cumdist values
        [cumdist_unique, ia] = unique(cumdist, 'stable');
        traj_unique = traj(ia, :);

        % Resample at equal distances along the path
        target_dists = linspace(0, total_length, num_pts);
        resampled_x = interp1(cumdist_unique, traj_unique(:,1), target_dists, 'linear');
        resampled_y = interp1(cumdist_unique, traj_unique(:,2), target_dists, 'linear');
        resampled = [resampled_x', resampled_y'];
    end

    left_matrix(:,:,i) = resampled;
end

for i = 1:n_right
    traj = right_trials{i};

    d = sqrt(sum(diff(traj).^2, 2));
    cumdist = [0; cumsum(d)];
    total_length = cumdist(end);

    if total_length == 0 || length(cumdist) < 2
        resampled = repmat(traj(1,:), num_pts, 1);
    else
        % Remove duplicate cumdist values
        [cumdist_unique, ia] = unique(cumdist, 'stable');
        traj_unique = traj(ia, :);

        target_dists = linspace(0, total_length, num_pts);
        resampled_x = interp1(cumdist_unique, traj_unique(:,1), target_dists, 'linear');
        resampled_y = interp1(cumdist_unique, traj_unique(:,2), target_dists, 'linear');
        resampled = [resampled_x', resampled_y'];
    end

    right_matrix(:,:,i) = resampled;
end


left_x = squeeze(left_matrix(:, 1, :)); % [time, trials]
left_y = squeeze(left_matrix(:, 2, :));
right_x = squeeze(right_matrix(:, 1, :));
right_y = squeeze(right_matrix(:, 2, :));

left_mean_x = mean(left_x, 2, 'omitnan');
%left_sem_x = std(left_x, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(left_x), 2));
left_sem_x = 1.96 * std(left_x, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(left_x), 2));

left_mean_y = mean(left_y, 2, 'omitnan');
%left_sem_y = std(left_y, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(left_y), 2));
left_sem_y = 1.96 * std(left_y, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(left_y), 2));

right_mean_x = mean(right_x, 2, 'omitnan');
%right_sem_x = std(right_x, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(right_x), 2));
right_sem_x = 1.96 * std(right_x, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(right_x), 2));

right_mean_y = mean(right_y, 2, 'omitnan');
%right_sem_y = std(right_y, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(right_y), 2));
right_sem_y = 1.96 * std(right_y, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(right_y), 2));

% Plot with confidence intervals
figure('color','white');
hold on;
left_x_shade = [left_mean_x - left_sem_x; flipud(left_mean_x + left_sem_x)];
left_y_shade = [left_mean_y - left_sem_y; flipud(left_mean_y + left_sem_y)];

% Plot shaded region
fill([left_mean_x; flipud(left_mean_x)], left_y_shade, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(left_mean_x, left_mean_y, 'r-', 'LineWidth', 2);

% Plot Right with shaded error (blue)
right_x_shade = [right_mean_x - right_sem_x; flipud(right_mean_x + right_sem_x)];
right_y_shade = [right_mean_y - right_sem_y; flipud(right_mean_y + right_sem_y)];

% Plot shaded region
fill([right_mean_x; flipud(right_mean_x)], right_y_shade, 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(right_mean_x, right_mean_y, 'b-', 'LineWidth', 2);

title('Average Trajectories with 95% CI');
xlabel('X Position');
ylabel('Y Position');
% axis equal;

figure('color','white');
hold on;
% Compute tangent vectors
dx = gradient(left_mean_x);
dy = gradient(left_mean_y);
tangents = [dx, dy];
tangent_norms = sqrt(sum(tangents.^2, 2));

% Normalize tangents
tangents_unit = tangents ./ tangent_norms;

% Get normals by rotating tangents 90 degrees counterclockwise
normals = [-tangents_unit(:,2), tangents_unit(:,1)];

% Compute upper and lower bounds using SEM in the normal direction
left_upper = [left_mean_x, left_mean_y] + normals .* left_sem_y;
left_lower = [left_mean_x, left_mean_y] - normals .* left_sem_y;

% Plot ribbon with patch
x_fill = [left_upper(:,1); flipud(left_lower(:,1))];
y_fill = [left_upper(:,2); flipud(left_lower(:,2))];
patch(x_fill, y_fill, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot the mean line
plot(left_mean_x, left_mean_y, 'r-', 'LineWidth', 2);


% Repeat for right trajectory
dx = gradient(right_mean_x);
dy = gradient(right_mean_y);
tangents = [dx, dy];
tangent_norms = sqrt(sum(tangents.^2, 2));
tangents_unit = tangents ./ tangent_norms;
normals = [-tangents_unit(:,2), tangents_unit(:,1)];

right_upper = [right_mean_x, right_mean_y] + normals .* right_sem_y;
right_lower = [right_mean_x, right_mean_y] - normals .* right_sem_y;

x_fill = [right_upper(:,1); flipud(right_lower(:,1))];
y_fill = [right_upper(:,2); flipud(right_lower(:,2))];
patch(x_fill, y_fill, 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(right_mean_x, right_mean_y, 'b-', 'LineWidth', 2);

title('Mean Trajectories with 95% CI (Orthogonal)');
xlabel('X Position'); ylabel('Y Position'); 
axis equal;



%{
% dlc_file = dir('*T22'); 
% tracking_data = csvread(dlc_file, 2, 1);
tracking_data = csvread('\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\Cue_task-Lucy-2025-02-24\T22_41022_150631 - T22_41022_150631DLC_resnet50_Cue_taskFeb24shuffle1_100000.csv', 2, 0);

for i = 1:size(tracking_data, 1)
    if tracking_data(i, 4) >= 0.98 && tracking_data(i, 7) >= 0.98 && tracking_data(i, 10) >= 0.98 ...
            && tracking_data(i, 13) >= 0.98 && tracking_data(i, 16) >= 0.98
        cutoff = i;
        break
    else
        continue
    end
end

tracking_data = tracking_data(cutoff:end, :);

tracking_data_og=tracking_data;

sz = size(tracking_data, 1);
cut_rows = [];

% delete uncertainty
j = 1;
while j <= sz
    if tracking_data(j, 13) < 0.99
        tracking_data(j, :) = [];
        sz = sz - 1;
    else
        j = j + 1;
    end
end

tracking_data_og_2=tracking_data;

j=2;
while j < size(tracking_data, 1)
    if tracking_data(j, 11) >= tracking_data(j-1, 11)+10 || tracking_data(j, 12) >= tracking_data(j-1, 12)+10 || ...
            tracking_data(j, 11) <= tracking_data(j-1, 11)-10 || tracking_data(j, 12) <= tracking_data(j-1, 12)-10
        start_cut = j;
        cut_pt = tracking_data(j, 11);
        
        k=j+1;
        % while k < size(tracking_data, 1)
        %     if tracking_data(k, 11) >= tracking_data(k-1, 11)+10 || tracking_data(k, 12) >= tracking_data(k-1, 12)+10 || ...
        %             tracking_data(k, 11) <= tracking_data(k-1, 11)-10 || tracking_data(k, 12) <= tracking_data(k-1, 12)-10
        %         break;  % End of bad block (jump back)
        %     else
        %         k = k + 1;
        %     end
        % end

        while k < size(tracking_data, 1)
            if tracking_data(k, 12) <= 850 && tracking_data(k, 11) > 680 ... 
                    && tracking_data(k, 11) < 730 % these are rough numbers -> and only for this session
                break;  % End of bad block (jump back)
            elseif tracking_data(k, 12) > 850 && tracking_data(k, 12) < 890 
                break;  % End of bad block (jump back)
            else
                k = k + 1;
            end
        end

        end_cut = k - 1;

        % Delete the bad block
        tracking_data(start_cut:end_cut, :) = [];
        cut_rows = [cut_rows; (start_cut:end_cut)'];

        % Adjust j after deletion
        j = start_cut;
        j = j + 1;
    else
        j = j + 1;
    end
end


figure('color','white');
colormap(jet)
scatter(tracking_data(:,11), tracking_data(:,12), 20, tracking_data(:,11), 'filled') % 20 is point size
colorbar

%}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




end








