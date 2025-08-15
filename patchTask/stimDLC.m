% analyzing stim vs no stim trajectories


function stimDLC

direc = '\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\Cue task';

stim_sessions = {'T20\Final\T20_241101_171908', ...
    'T20\Final\T20_241102_163911', ...
    'T22\Final\T22_241022_150631', ...
    'T22\Final\T22_241029_145741'};%,... 
    % 'T22\Final\T22_241114_143845', ...
    % 'T22\Final\T22_241119_144043', ...
    % 'T22\Final\T22_241209_173821', ...
    % 'T22\Final\T22_241209_182341', ... 
    % 'T22\Final\T22_241213_160938'};

left_trials_all = {};
right_trials_all = {};
left_trials_all_stim = {};
right_trials_all_stim = {};

for s = 1:length(stim_sessions)
    basepath = [direc filesep stim_sessions{s}];
    cd(basepath);
   
    % Load tracking and behavior data  
    tracking_file = dir(fullfile(basepath, '*Tracking.Behavior.mat'));
    load(tracking_file.name);

    behavior_file = dir(fullfile(basepath, '*TrialBehavior.Events.mat'));
    load(behavior_file.name);

    % clean session
    j=2;
    while j < size(tracking.position.x, 1)
        if tracking.position.x(j) >= tracking.position.x(j-1)+0.85 || tracking.position.y(j) >= tracking.position.y(j-1)+0.85 || ...
                tracking.position.x(j) <= tracking.position.x(j-1)-0.85 || tracking.position.y(j) <= tracking.position.y(j-1)-0.85
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

    figure('color','white');
    colormap(jet)
    scatter(tracking.position.x, tracking.position.y, 20, tracking.position.x, 'filled') % 20 is point size
    colorbar

    % Initialize session-specific trial storage
    left_trials = {};
    right_trials = {};
    left_trials_stim = {};
    right_trials_stim = {};

    % Loop through trials
    num_trials = length(behavTrials.choiceTS);
    for i = 1:num_trials
        if behavTrials.stim(i) == 0
            if behavTrials.correct(i) == 1
                curr_x = tracking.position.x(tracking.timestamps >= behavTrials.timestamps(i, 1) ...
                    & tracking.timestamps <= behavTrials.choiceTS(i));
                curr_y = tracking.position.y(tracking.timestamps >= behavTrials.timestamps(i, 1) ...
                    & tracking.timestamps <= behavTrials.choiceTS(i));
                if behavTrials.choice(i) == 1
                    left_trials{end+1} = [curr_x(:), curr_y(:)];
                else
                    right_trials{end+1} = [curr_x(:), curr_y(:)];
                end   
            end
        elseif behavTrials.stim(i) == 1
            if behavTrials.correct(i) == 1
                curr_x = tracking.position.x(tracking.timestamps >= behavTrials.timestamps(i, 1) ...
                    & tracking.timestamps <= behavTrials.choiceTS(i));
                curr_y = tracking.position.y(tracking.timestamps >= behavTrials.timestamps(i, 1) ...
                    & tracking.timestamps <= behavTrials.choiceTS(i));
                if behavTrials.choice(i) == 1
                    left_trials_stim{end+1} = [curr_x(:), curr_y(:)];
                else
                    right_trials_stim{end+1} = [curr_x(:), curr_y(:)];
                end   
            end
    end

    % Append to master trial lists
    left_trials_all = [left_trials_all, left_trials];
    right_trials_all = [right_trials_all, right_trials];
    left_trials_all_stim = [left_trials_all_stim, left_trials_stim];
    right_trials_all_stim = [right_trials_all_stim, right_trials_stim];
end


%% interpolate to distance - no stim

num_pts = 482; % number of path-based sample points
n_left = length(left_trials_all);
n_right = length(right_trials_all);


for i = 1:n_left
    traj = left_trials_all{i}; % traj is [N x 2] matrix of [x, y] coordinates
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
    traj = right_trials_all{i};
    % Remove rows with NaNs or Infs
    traj = traj(all(isfinite(traj), 2), :);

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
right_mean_y = median(right_y, 2, 'omitnan');
right_sem_x = 1.96 * std(right_x, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(right_x), 2));
right_sem_y = 1.96 * std(right_y, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(right_y), 2));

% COLORS
no_stim_color = [0.8902    0.8235    0.4353];
%no_stim_color = [1.0000    0.9020    0.4275];
%no_stim_color = [1.0000    0.5412    0.3569];
%no_stim_color = [0.9176    0.3216    0.4353];
stim_color = [0.4431    0.5373    1.0000];
custom_green = [0.3294    0.7216    0.4353];
custom_blue = [0.2235    0.4353    0.6784];

f = figure('color','white');
f.Position(1:4) = [50, 50, 900, 450];
hold on;
subplot(2, 2, [1,3])
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
patch(x_fill, y_fill, custom_green, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot the mean line
hold on;
plot(left_mean_x, left_mean_y, 'color', custom_green, 'LineWidth', 2.5);

% Repeat for right trajectory
dx = gradient(right_mean_x);
dy = gradient(right_mean_y);
tangents = [dx, dy];
tangent_norms = sqrt(sum(tangents.^2, 2));
tangents_unit = tangents ./ tangent_norms;
normals = [-tangents_unit(:,2), tangents_unit(:,1)];

right_upper = [right_mean_x, right_mean_y] + normals .* right_sem_y;
right_lower = [right_mean_x, right_mean_y] - normals .* right_sem_y;
hold on;
x_fill = [right_upper(:,1); flipud(right_lower(:,1))];
y_fill = [right_upper(:,2); flipud(right_lower(:,2))];
patch(x_fill, y_fill, custom_blue, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(right_mean_x, right_mean_y, 'color', custom_blue, 'LineWidth', 2.5);
%plot(93.2679, 100.341, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); %
%median nose ALL
%plot(93.2447, 106.85, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); % median body ALL

plot(93.2406, 98.6168, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); %
%median nose CORRECT
%plot(93.264, 107.249, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); % median body CORRECT

title('Baseline');
% xlabel('X Position'); 
% ylabel('Y Position'); 
xlim([78, 107]); 
ylim([30, 118]);
%axis equal;
hold off;

%% interpolate to distance - WITH STIM 

num_pts = 482; % number of path-based sample points
n_left = length(left_trials_all_stim);
n_right = length(right_trials_all_stim);


for i = 1:n_left
    traj = left_trials_all_stim{i}; % traj is [N x 2] matrix of [x, y] coordinates
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
    traj = right_trials_all_stim{i};
    % Remove rows with NaNs or Infs
    traj = traj(all(isfinite(traj), 2), :);

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
right_mean_y = median(right_y, 2, 'omitnan');
right_sem_x = 1.96 * std(right_x, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(right_x), 2));
right_sem_y = 1.96 * std(right_y, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(right_y), 2));

% figure('color','white');
subplot(2,2,[2,4])
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
patch(x_fill, y_fill, custom_green, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot the mean line
plot(left_mean_x, left_mean_y,  'color', custom_green, 'LineWidth', 2.5);

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
patch(x_fill, y_fill, custom_blue, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(right_mean_x, right_mean_y, 'color', custom_blue, 'LineWidth', 2.5);
%plot(93.2024, 96.2067, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); %
%median nose ALL
%plot(93.1851, 105.512, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); % median body ALL

plot(93.2029, 92.8714, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); %
%median nose CORRECT
%plot(93.1861, 105.203, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); % median body CORRECT

title('Stim');
% xlabel('X Position'); 
% ylabel('Y Position'); 
xlim([78, 107]); 
ylim([30, 118]); 
%axis equal;

% nose all
baseline_pt_divergence = (100.341 - 34.001) / (116.5164 - 34.001);
stim_pt_divergence = (96.2067 - 34.0047) / (116.7728 - 34.0047);

% nose correct
baseline_pt_divergence = (98.6168 - 34.001) / (116.5164 - 34.001);
stim_pt_divergence = (92.8714 - 34.0047) / (116.7728 - 34.0047);

% body all
baseline_pt_divergence = (106.85 - 34.001) / (116.5164 - 34.001);
stim_pt_divergence = (105.512 - 34.0047) / (116.7728 - 34.0047);

% body correct
baseline_pt_divergence = (107.249 - 34.001) / (116.5164 - 34.001);
stim_pt_divergence = (105.203 - 34.0047) / (116.7728 - 34.0047);

m = 1;

end

