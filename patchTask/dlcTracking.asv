
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
        tracking.timestamps(start_cut:end_cut, :) = [];
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


%% Find average trajectories - interpolating by time - BAD WAY OF DOING IT

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
title('Average Trajectories - interpolating across time');
axis equal;




%% interpolate to distance

% COLORS 
custom_green = [0.3294    0.7216    0.4353];
custom_blue = [0.2235    0.4353    0.6784];

num_pts = 482; % number of path-based sample points
n_left = length(left_trials);
n_right = length(right_trials);


for i = 1:n_left
    traj = left_trials{i}; % traj is [N x 2] matrix of [x, y] coordinates
    % Remove rows with NaNs or Infs
    traj = traj(all(isfinite(traj), 2), :);

    % Compute cumulative path length
    d = sqrt(sum(diff(traj).^2, 2));        % distances between points
    cumdist = [0; cumsum(d)];               % cumulative distance
    total_length = cumdist(end);

    norm_cumdist = cumdist / cumdist(end); % normalized to [0, 1]

    if total_length == 0 || length(cumdist) < 2
        resampled = repmat(traj(1,:), num_pts, 1); % handle degenerate cases
    else
        % Remove duplicate cumdist values
        [cumdist_unique, ia] = unique(norm_cumdist, 'stable');
        traj_unique = traj(ia, :);

        % Resample at equal distances along the path
        % target_dists = linspace(0, total_length, num_pts);
        target_dists = linspace(0, 1, num_pts); % uniform sampling from 0 to 1
        resampled_x = interp1(cumdist_unique, traj_unique(:,1), target_dists, 'linear');
        resampled_y = interp1(cumdist_unique, traj_unique(:,2), target_dists, 'linear');
        resampled = [resampled_x', resampled_y'];
    end

    left_matrix(:,:,i) = resampled;
end

for i = 1:n_right
    traj = right_trials{i};
    % Remove rows with NaNs or Infs
    traj = traj(all(isfinite(traj), 2), :);
    % figure;
    % plot(right_trials{1,i}(:,1),right_trials{1,i}(:,2))

    d = sqrt(sum(diff(traj).^2, 2));
    cumdist = [0; cumsum(d)];
    total_length = cumdist(end);

    norm_cumdist = cumdist / cumdist(end); % normalized to [0, 1]

    if total_length == 0 || length(cumdist) < 2
        resampled = repmat(traj(1,:), num_pts, 1);
    else
        % Remove duplicate cumdist values
        [cumdist_unique, ia] = unique(norm_cumdist, 'stable');
        traj_unique = traj(ia, :);

        % target_dists = linspace(0, total_length, num_pts);
        target_dists = linspace(0, 1, num_pts); % uniform sampling from 0 to 1
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

left_mean_x = median(left_x, 2, 'omitnan');
left_mean_y = median(left_y, 2, 'omitnan');
left_sem_x = 1.96 * std(left_x, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(left_x), 2));
left_sem_y = 1.96 * std(left_y, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(left_y), 2));

right_mean_x = median(right_x, 2, 'omitnan');
right_sem_x = 1.96 * std(right_x, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(right_x), 2));
right_mean_y = median(right_y, 2, 'omitnan');
right_sem_y = 1.96 * std(right_y, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(right_y), 2));

%left_sem_x = std(left_x, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(left_x), 2));
%left_sem_y = std(left_y, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(left_y), 2));
%right_sem_x = std(right_x, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(right_x), 2));
%right_sem_y = std(right_y, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(right_y), 2));

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

sem_vector = [left_sem_x, left_sem_y];

% Project SEM onto the normal
%sem_along_normal = sum(sem_vector .* normals, 2); % Dot product per frame
sem_along_normal = dot(sem_vector, normals, 2);  % this gives signed magnitude

% Then apply as magnitude in normal direction
left_upper = [left_mean_x, left_mean_y] + normals .* sem_along_normal;
left_lower = [left_mean_x, left_mean_y] - normals .* sem_along_normal;


% Compute upper and lower bounds using SEM in the normal direction
% left_upper = [left_mean_x, left_mean_y] + normals .* left_sem_y;
% left_lower = [left_mean_x, left_mean_y] - normals .* left_sem_y;

% Plot ribbon with patch
x_fill = [left_upper(:,1); flipud(left_lower(:,1))];
y_fill = [left_upper(:,2); flipud(left_lower(:,2))];
patch(x_fill, y_fill, custom_green, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot the mean line
plot(left_mean_x, left_mean_y, 'color', custom_green, 'LineWidth', 2);

% Repeat for right trajectory
dx = gradient(right_mean_x);
dy = gradient(right_mean_y);
tangents = [dx, dy];
tangent_norms = sqrt(sum(tangents.^2, 2));
tangents_unit = tangents ./ tangent_norms;
normals = [-tangents_unit(:,2), tangents_unit(:,1)];

sem_vector = [right_sem_x, right_sem_y];

% Project SEM onto the normal
sem_along_normal = sum(sem_vector .* normals, 2); % Dot product per frame

% Then apply as magnitude in normal direction
% right_upper = [right_mean_x, right_mean_y] + normals .* sem_along_normal;
% right_lower = [right_mean_x, right_mean_y] - normals .* sem_along_normal;

right_upper = [right_mean_x, right_mean_y] + normals .* right_sem_y;
right_lower = [right_mean_x, right_mean_y] - normals .* right_sem_y;

x_fill = [right_upper(:,1); flipud(right_lower(:,1))];
y_fill = [right_upper(:,2); flipud(right_lower(:,2))];
patch(x_fill, y_fill, custom_blue, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(right_mean_x, right_mean_y, 'color', custom_blue, 'LineWidth', 2);

% track lines
plot([79 107.35], [119.5 119.5], 'k-', 'LineWidth', 2);  % Black line, width 2
plot([79 79], [119.5 112.5], 'k-', 'LineWidth', 2);  % Black line, width 2
plot([107.35 107.35], [119.5 112.5], 'k-', 'LineWidth', 2);  % Black line, width 2
plot([79 89.5], [112.5 112.5], 'k-', 'LineWidth', 2);  % Black line, width 2
plot([96.5 107.35], [112.5 112.5], 'k-', 'LineWidth', 2);  % Black line, width 2
plot([89.5 89.5], [112.5 20.5], 'k-', 'LineWidth', 2);  % Black line, width 2
plot([96.5 96.5], [112.5 20.5], 'k-', 'LineWidth', 2);  % Black line, width 2
plot([89.5 96.5], [20.5 20.5], 'k-', 'LineWidth', 2);  % Black line, width 2

% plot([79 119.5], [107 119.5], 'k-', 'LineWidth', 2);  % Black line, width 2
% plot([79 119.5], [79 112.5], 'k-', 'LineWidth', 2);  % Black line, width 2
% plot([107 119.5], [107 112.5], 'k-', 'LineWidth', 2);  % Black line, width 2
% plot([79 112.5], [89.5 112.5], 'k-', 'LineWidth', 2);  % Black line, width 2
% plot([96.5 112.5], [107 112.5], 'k-', 'LineWidth', 2);  % Black line, width 2
% plot([89.5 112.5], [89.5 12.5], 'k-', 'LineWidth', 2);  % Black line, width 2
% plot([96.5 112.5], [96.5 12.5], 'k-', 'LineWidth', 2);  % Black line, width 2
% plot([89.5 12.5], [96.5 12.5], 'k-', 'LineWidth', 2);  % Black line, width 2

xlabel('X Position'); ylabel('Y Position'); 
axis equal;
xlim([70, 115]); 
ylim([12, 125]); 
hold off;

%% example trajectory
f = figure('color','white');
f.Position(1:4) = [50, 50, 862.5, 750];
%f.Position(1:4) = [50, 50, 1150, 1000];
hold on;

plot(left_x(:,1), left_y(:,1), 'color', custom_green, 'LineWidth', 4);
plot(right_x(:,3), right_y(:,3), 'color', custom_blue, 'LineWidth', 4);

xlabel('X Position'); ylabel('Y Position'); 
axis equal;
xlim([75, 110]); 
ylim([30, 118]); 
hold off;


%% p value

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



m=1;
end








