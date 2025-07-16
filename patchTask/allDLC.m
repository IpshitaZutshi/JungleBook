

function allDLC

direc = '\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\Cue task';

sessions = {'T20\Final\T20_241030_153043', ...
    'T20\Final\T20_241031_173852', ...
    'T20\Final\T20_241101_171908', ...
    'T20\Final\T20_241102_163911', ...
    'T20\Final\T20_241103_131201', ...
    'T20\Final\T20_241107_145457', ...
    'T20\Final\T20_241111_154144', ...
    'T20\Final\T20_241114_133417', ...
    'T20\Final\T20_241115_172023', ...
    'T20\Final\T20_241117_132754', ...
    'T20\Final\T20_241118_144118', ...
    'T20\Final\T20_241119_124533', ...
    'T20\Final\T20_241213_152020', ...
    'T22\Final\T22_241022_150631', ... 
    'T22\Final\T22_241023_175232', ...
    'T22\Final\T22_241024_153843', ...
    'T22\Final\T22_241029_145741', ...
    'T22\Final\T22_241030_144753', ...
    'T22\Final\T22_241102_173932', ...
    'T22\Final\T22_241103_123921', ...
    'T22\Final\T22_241104_150209', ...
    'T22\Final\T22_241113_150043', ...
    'T22\Final\T22_241114_143845', ...
    'T22\Final\T22_241119_144043', ...
    'T22\Final\T22_241209_173821', ...
    'T22\Final\T22_241209_182341', ... 
    'T22\Final\T22_241209_190750', ...
    'T22\Final\T22_241212_154105', ...
    'T22\Final\T22_241213_160938'};

% sessions = {'T20\Final\T20_241101_171908', ...
%     'T20\Final\T20_241102_163911', ...
%     'T22\Final\T22_241022_150631', ... 
%     'T22\Final\T22_241029_145741'};

% these sessions have screwed up tracking
% 'T22\Final\T22_241101_164241'
% 'T22\Final\T22_241031_162840', ...
% 'T22\Final\T22_241216_120243'

left_trials_all = {};
right_trials_all = {};

for s = 1:length(sessions)
    basepath = [direc filesep sessions{s}];
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

    % forward_tracking.x = [];
    % forward_tracking.y = [];
    % forward_tracking.direction = nan(size(tracking.position.x, 1), 1);
    % forward_tracking.trial = nan(size(tracking.position.x, 1), 1);
    % 
    % num_trials = length(behavTrials.choiceTS);
    % 
    % for i = 1:num_trials
    %     curr_trial_tracking_x = tracking.position.x(tracking.timestamps >= behavTrials.timestamps(i, 1) ... 
    %         & tracking.timestamps <= behavTrials.choiceTS(i));
    %     curr_trial_tracking_y = tracking.position.y(tracking.timestamps >= behavTrials.timestamps(i, 1) ... 
    %         & tracking.timestamps <= behavTrials.choiceTS(i));
    % 
    %     forward_tracking.x = [forward_tracking.x; curr_trial_tracking_x];
    %     forward_tracking.y = [forward_tracking.y; curr_trial_tracking_y];
    % 
    %     forward_tracking.direction(tracking.timestamps >= behavTrials.timestamps(i, 1) ... 
    %         & tracking.timestamps <= behavTrials.choiceTS(i)) = behavTrials.choice(i);
    %     forward_tracking.trial(tracking.timestamps >= behavTrials.timestamps(i, 1) ... 
    %         & tracking.timestamps <= behavTrials.choiceTS(i)) = i;
    % 
    % end
    % 
    % forward_tracking.direction = forward_tracking.direction(~isnan(forward_tracking.direction));
    % forward_tracking.trial = forward_tracking.trial(~isnan(forward_tracking.trial));

    % Initialize session-specific trial storage
    left_trials = {};
    right_trials = {};

    % Loop through trials
    num_trials = length(behavTrials.choiceTS);
    for i = 1:num_trials
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
    end

    % Append to master trial lists
    left_trials_all = [left_trials_all, left_trials];
    right_trials_all = [right_trials_all, right_trials];
end


%% interpolate to distance

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
custom_green = [0.3294    0.7216    0.4353];
custom_blue = [0.2235    0.4353    0.6784];

f = figure('color','white');
f.Position(1:4) = [50, 50, 862.5, 750];
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
plot(left_mean_x, left_mean_y,'color', custom_green, 'LineWidth', 5);

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

plot(right_mean_x, right_mean_y, 'color', custom_blue, 'LineWidth', 5);

% title('Mean Trajectories with 95% CI');
% xlabel('X Position'); 
% ylabel('Y Position'); 
% axis equal;

% % Create polyshapes from the confidence interval ribbons
% left_ci_poly = polyshape(x_fill, y_fill);  % x_fill, y_fill from left patch
% right_ci_poly = polyshape([right_upper(:,1); flipud(right_lower(:,1))], ...
%                           [right_upper(:,2); flipud(right_lower(:,2))]);
% 
% % Find the intersection of the two CI polygons
% ci_intersection = intersect(left_ci_poly, right_ci_poly);
% 
% % If intersection is not empty, find the last point
% if ~isempty(ci_intersection.Vertices)
%     % Option 1: Use the point in the intersection polygon with the highest index in the trajectory
%     % Find the closest point in the intersection to the end of the left/right trajectory
%     end_point = [left_mean_x(end), left_mean_y(end)];
%     distances = vecnorm(ci_intersection.Vertices - end_point, 2, 2);
%     [~, idx] = min(distances);
%     last_point = ci_intersection.Vertices(idx, :);
% 
%     % Plot the red dot
%     plot(last_point(1), last_point(2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
% end

%plot(93.4247, 106.687, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'g'); % mean
plot(93.4221, 108.158, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); %
% median nose
%plot(93.4234, 108.917, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % median body

xlim([78, 107]); 
ylim([30, 118]); 
axis equal;

point_of_divergence = (108.158 - 33.9509) / (116.5964 - 33.9509);

m=1;

end

