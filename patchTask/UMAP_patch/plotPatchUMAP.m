function plotPatchUMAP(varargin)

p = inputParser;
addParameter(p,'min_trial',0,@isnumeric);
addParameter(p,'max_trial',1000,@isnumeric);
parse(p,varargin{:});
min_trial = p.Results.min_trial; % confused about this
max_trial = p.Results.max_trial;

basepath = pwd;
umapFile = dir(fullfile(basepath, '*speed.csv'));
behavFile = dir(fullfile(basepath, '*.position_behavior_speed.mat'));
behaviorFile = dir(fullfile(basepath, '*.TrialBehavior.mat'));
trackFile = dir(fullfile(basepath, '*Tracking.Behavior.mat'));

currentFolder = pwd;
[~, currentFolderName] = fileparts(currentFolder);

umap_path = umapFile.folder;
umap_name = umapFile.name;

%% SET PARAMETERS
dim1 = 1;
dim2 = 2;
dim3 = 3;
dotSize = 5;
A = -59.8364; E = -4.0112; % azimuth and elevation for viewing
%% load Umap result
Umap_results = readtable(fullfile(umap_path, umap_name));
Umap_results = table2array(Umap_results);

%% load position direction and behavior information
load(behavFile.name);
load(behaviorFile.name);
load(trackFile.name);
timestamp_beh = timestamp_beh';

pos_plot = position_y_all;
pos_plot(isnan(pos_plot))=0; % deal with nan

lick_plot = licked_port_ds;
lick_plot(isnan(lick_plot))=0; % deal with nan

trial_num_plot = trial_num_ds;
trial_num_plot(isnan(trial_num_plot))=0; % deal with nan

patch_num_plot = patch_num_ds;
patch_num_plot(isnan(patch_num_plot))=0; % deal with nan

trial_patch_num_plot = patch_trial_num_ds;
trial_patch_num_plot(isnan(trial_patch_num_plot))=0; % deal with nan

speed_plot = speed_all;
speed_plot(isnan(speed_plot))=0; % deal with nan

speed_dir_plot = speed_dir;
speed_dir_plot(isnan(speed_dir_plot))=0; % deal with nan

direction_plot = direction;
direction_plot(isnan(direction_plot))=0;

traject_dir_plot = trajectory_dir_ds;
traject_dir_plot(isnan(traject_dir_plot))=0;

%outcome_plot = outcome_ds;
outcome_ds(isnan(outcome_ds))=0; %plot

patch_type_plot = patch_type_ds;
patch_type_plot(isnan(patch_type_plot))=0;

%outcome_lick = outcome_plot + 1; % delete
outcome_lick_time = zeros(size(outcome_ds));%plot
outcome_lick_time(lick_plot ~= 0) = outcome_ds(lick_plot ~= 0);%lick

outcome_lick_high_patch = zeros(size(outcome_ds));%plot
outcome_lick_high_patch(patch_type_plot==1) = outcome_lick_time(patch_type_plot==1);



%% Filter trials if necessary
filtered_idx = (trial_num_plot > min_trial) & (trial_num_plot < max_trial);

% Filter variables
pos_plot_filtered = pos_plot(filtered_idx);
lick_plot_filtered = lick_plot(filtered_idx);
trial_num_plot_filtered = trial_num_plot(filtered_idx);
direction_plot_filtered = direction_plot(filtered_idx);
speed_plot_filtered = speed_plot(filtered_idx);
Umap_results_filtered = Umap_results(filtered_idx, :);
event_number_filtered = 1:length(pos_plot_filtered);
outcome_filtered = outcome_ds(filtered_idx);%plot
outcome_lick_time_filtered = outcome_lick_time(filtered_idx);
outcome_lick_high_patch_filtered = outcome_lick_high_patch(filtered_idx);

%% Create multiplot

figure(1);
t = tiledlayout(2, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
title(t, currentFolderName);

ax1 = nexttile;
ax2 = nexttile;
ax3 = nexttile;
ax4 = nexttile;
ax5 = nexttile;
ax6 = nexttile;

Link = linkprop([ax1, ax2, ax3, ax4, ax5, ax6], {'CameraPosition', 'CameraTarget', 'CameraUpVector', 'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', Link);

color = 'jet';
dotSize = 5;
% Plot1: colored by position
hold on
scatter3(ax1, Umap_results_filtered(:,dim1), Umap_results_filtered(:,dim2), Umap_results_filtered(:,dim3), dotSize, pos_plot_filtered, 'filled');
colormap(ax1, color);
view(ax1, A, E);
grid(ax1, 'off');
axis(ax1, 'tight');
axis(ax1, 'off');
colorbar(ax1);
title(ax1, 'Position');
hold off

% Plot 2: Colored by chosen port
scatter3(ax2, Umap_results_filtered(:,dim1), Umap_results_filtered(:,dim2), Umap_results_filtered(:,dim3), dotSize, lick_plot_filtered, 'filled');
colormap(ax2, color);
view(ax2, A, E);
grid(ax2, 'off');
axis(ax2, 'tight');
axis(ax2, 'off');
colorbar(ax2);
title(ax2, 'Chosen port');

% Plot 3: Colored by trial number
scatter3(ax3, Umap_results_filtered(:,dim1), Umap_results_filtered(:,dim2), Umap_results_filtered(:,dim3), dotSize, trial_num_plot_filtered, 'filled');
colormap(ax3, color);
view(ax3, A, E);
grid(ax3, 'off');
axis(ax3, 'tight');
axis(ax3, 'off');
colorbar(ax3);
title(ax3, 'Trial number');

% Plot 4: Colored by direction
scatter3(ax4, Umap_results_filtered(:,dim1), Umap_results_filtered(:,dim2), Umap_results_filtered(:,dim3), dotSize, direction_plot_filtered, 'filled');
colormap(ax4, color);
view(ax4, A, E);
grid(ax4, 'off');
axis(ax4, 'tight');
axis(ax4, 'off');
colorbar(ax4);
title(ax4, 'Direction');

%Plot 5: Colored by speed
scatter3(ax5, Umap_results_filtered(:,dim1), Umap_results_filtered(:,dim2), Umap_results_filtered(:,dim3), dotSize, speed_plot_filtered, 'filled');
colormap(ax5, color);
view(ax5, A, E);
grid(ax5, 'off');
axis(ax5, 'tight');
axis(ax5, 'off');
colorbar(ax5);
% caxis(ax5,[min(speed_plot_filtered) 30]);
title(ax5, 'Speed');


%Plot 6: Only show times of licks
% find timestamps when trial switches
trial_lengths = zeros(behavTrials.num_trials, 2);
trial = 1;
len = 1;
for i = 2:length(lick_plot_filtered)
    if i == length(lick_plot_filtered)
        len = len + 1;
        trial_lengths(trial, 1) = len;
        trial_lengths(trial, 2) = i;
    elseif lick_plot_filtered(i) == lick_plot_filtered(i-1)
        len = len + 1;
        continue
    else
        trial_lengths(trial, 1) = len;
        trial_lengths(trial, 2) = i-1;
        trial = trial + 1;
        len = 1;
    end
end

umap_minus = Umap_results_filtered;
umap_minus(trial_lengths(:,2), :) = [];
scatter3(ax6, umap_minus(:,dim1), umap_minus(:,dim2), umap_minus(:,dim3), ...
    dotSize, [0.85 0.85 0.85], 'filled', 'MarkerFaceAlpha', 0.4);
view(ax6, A, E);

cmap = jet(7);
colormap(ax6, cmap);
% caxis(ax2, [1 7]);  
grid(ax6, 'off');
axis(ax6, 'tight');
axis(ax6, 'off');
title(ax6, 'Licked Ports');
hold(ax6, 'on');
colorbar(ax6, 'Ticks', [0, 0.2, 0.35, 0.5, 0.65, 0.8, 1], 'TickLabels', {'1', '2', '3', '4','5', '6', '7'});

for i = 1:7
    idx_1 = behavTrials.port == i;
    idx_2 = trial_lengths(idx_1, 2);
    scatter3(ax6, Umap_results_filtered(idx_2, dim1), Umap_results_filtered(idx_2, dim2), Umap_results_filtered(idx_2, dim3), ...
        dotSize, cmap(i, :), 'filled', 'MarkerFaceAlpha', 0.8);
end

% % Save figure
% figFile = fullfile('/Users/lauraribalta/Desktop', [currentFolderName '_UMAP_Multiplot.fig']);
% savefig(figFile);

%% PLOT OUTCOME
figure(2);
t = tiledlayout(2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact'); % 1 fila, 3 columnas
title(t, 'Multiplot');

% Subplot 1: Outcome
ax1 = nexttile;
view(ax1, A, E);
grid(ax1, 'off');
axis(ax1, 'tight');
axis(ax1, 'off');
title(ax1, 'Outcome');
hold(ax1, 'on');

customGreen = [152, 194, 9] / 255;
customRed = [238, 75, 43] / 255;
nonreward_idx = outcome_lick_time_filtered == 1; % idx of all timepoints in the trial
reward_idx = outcome_lick_time_filtered == 2;

scatter3(ax1, Umap_results_filtered(nonreward_idx, dim1), Umap_results_filtered(nonreward_idx, dim2), Umap_results_filtered(nonreward_idx, dim3), ...
    dotSize, customRed, 'filled', 'MarkerFaceAlpha', 0.8, 'DisplayName', 'Non-Rewarded');

scatter3(ax1, Umap_results_filtered(reward_idx, dim1), Umap_results_filtered(reward_idx, dim2), Umap_results_filtered(reward_idx, dim3), ...
    dotSize, customGreen, 'filled', 'MarkerFaceAlpha', 0.8, 'DisplayName', 'Rewarded');

legend(ax1, 'show', 'Location', 'bestoutside');

% Subplot 2: Licked ports
% find timestamps when trial switches
trial_lengths = zeros(behavTrials.num_trials, 2);
trial = 1;
len = 1;
for i = 2:length(lick_plot_filtered)
    if i == length(lick_plot_filtered)
        len = len + 1;
        trial_lengths(trial, 1) = len;
        trial_lengths(trial, 2) = i;
    elseif lick_plot_filtered(i) == lick_plot_filtered(i-1)
        len = len + 1;
        continue
    else
        trial_lengths(trial, 1) = len;
        trial_lengths(trial, 2) = i-1;
        trial = trial + 1;
        len = 1;
    end
end

ax2 = nexttile;
umap_minus = Umap_results_filtered;
umap_minus(trial_lengths(:,2), :) = [];
scatter3(ax2, umap_minus(:,dim1), umap_minus(:,dim2), umap_minus(:,dim3), ...
    dotSize, [0.85 0.85 0.85], 'filled', 'MarkerFaceAlpha', 0.4);
view(ax2, A, E);

cmap = jet(7);
colormap(ax2, cmap);
% caxis(ax2, [1 7]);  
grid(ax2, 'off');
axis(ax2, 'tight');
axis(ax2, 'off');
title(ax2, 'Licked Ports');
hold(ax2, 'on');
colorbar(ax2, 'Ticks', [0, 0.2, 0.35, 0.5, 0.65, 0.8, 1], 'TickLabels', {'1', '2', '3', '4','5', '6', '7'});

tracker = [];

for i = 1:7
    idx_1 = behavTrials.port == i;
    idx_2 = trial_lengths(idx_1, 2);
    scatter3(ax2, Umap_results_filtered(idx_2, dim1), Umap_results_filtered(idx_2, dim2), Umap_results_filtered(idx_2, dim3), ...
        dotSize, cmap(i, :), 'filled', 'MarkerFaceAlpha', 0.8);
    tracker = [tracker; idx_2];
end


% Subplot 3: Outcome at lick
ax3 = nexttile;
scatter3(ax3, umap_minus(:,dim1), umap_minus(:,dim2), umap_minus(:,dim3), ...
    dotSize, [0.85 0.85 0.85], 'filled', 'MarkerFaceAlpha', 0.4);
hold(ax3, 'on');

nonreward_idx_2 = [];
reward_idx_2 = [];
for i = 1:behavTrials.num_trials
    if reward_idx(trial_lengths(i, 2)) == 1
        reward_idx_2 = [reward_idx_2; trial_lengths(i, 2)]; % idx of just the lick
    else
        nonreward_idx_2 = [nonreward_idx_2; trial_lengths(i, 2)];
    end
end


scatter3(ax3, Umap_results_filtered(nonreward_idx_2, dim1), Umap_results_filtered(nonreward_idx_2, dim2), Umap_results_filtered(nonreward_idx_2, dim3), ...
    dotSize, customRed, 'filled', 'MarkerFaceAlpha', 0.8, 'DisplayName', 'Non-Rewarded');

scatter3(ax3, Umap_results_filtered(reward_idx_2, dim1), Umap_results_filtered(reward_idx_2, dim2), Umap_results_filtered(reward_idx_2, dim3), ...
    dotSize, customGreen, 'filled', 'MarkerFaceAlpha', 0.8, 'DisplayName', 'Rewarded');

view(ax3, A, E);
grid(ax3, 'off');
axis(ax3, 'tight');
axis(ax3, 'off');
title(ax3, 'Outcome at Lick');
hold(ax3, 'off');

% Subplot 4: Only trials in high probability patch
ax4 = nexttile;
scatter3(ax4, Umap_results_filtered(:,dim1), Umap_results_filtered(:,dim2), Umap_results_filtered(:,dim3), ...
    dotSize, [0.85 0.85 0.85], 'filled', 'MarkerFaceAlpha', 0.4);
view(ax4, A, E);
grid(ax4, 'off');
axis(ax4, 'tight');
axis(ax4, 'off');
title(ax4, 'Outcome (High probability patch)');
hold(ax4, 'on');


scatter3(ax4, Umap_results_filtered(nonreward_idx, dim1), Umap_results_filtered(nonreward_idx, dim2), Umap_results_filtered(nonreward_idx, dim3), ...
    dotSize, customRed, 'filled', 'MarkerFaceAlpha', 0.8, 'DisplayName', 'Non-Rewarded');

scatter3(ax4, Umap_results_filtered(reward_idx, dim1), Umap_results_filtered(reward_idx, dim2), Umap_results_filtered(reward_idx, dim3), ...
    dotSize, customGreen, 'filled', 'MarkerFaceAlpha', 0.8, 'DisplayName', 'Rewarded');

hold(ax4, 'off');

Link = linkprop([ax1, ax2, ax3, ax4], {'CameraPosition', 'CameraTarget', 'CameraUpVector', 'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', Link);

%% PORT LOCATIONS
% positions along track where ports are located. determined from
% pos_plot_filtered
port_5_location = 91.3;
port_6_location = 111.7;
port_7_location = 118.7;


 %% PLOT 3 PORT TRAJECTORIES
 ports = [1, 2, 3];
 permutations = perms(ports); % Generate all permutations
 permutations = [permutations; 1, 2, NaN; 1, 3, NaN; 2, 1, NaN; 2, 3, NaN; 3, 1, NaN; 3, 2, NaN];
 %permutations = [permutations; 5, 6, NaN; 5, 7, NaN; 6, 5, NaN; 6, 7, NaN; 7, 6, NaN; 7, 5, NaN];

 % DOES THIS START AT 5, go to 7, or is it starting with the run TO 5, then
 % to 7
 % maybe plot in  a red dot where th elick is
 
  % sliding window

 mode = 0; % 1 is when you just want ONE run. 0 is when you want ALL runs 
 
 for i = 1:size(permutations, 1)
     perm_found = zeros(1, behavTrials.num_trials); % a logical vector saying if the permutation occurs in those trials (1) or does not (0)
     selected_trial_vector = [];
     find_first_position = NaN;
     find_last_position = NaN;
     permutation = permutations(i, :);
     permutation = permutation(~isnan(permutation));
     window_size = length(permutation); 
     match_found = false;  
 
     for start_idx = 1:(behavTrials.num_trials - window_size + 1)
         current_window = behavTrials.port(start_idx:start_idx + window_size - 1);
         current_window = current_window';
         if isequal(current_window, permutation)
             perm_found(start_idx:start_idx + window_size - 1) = 1; % should it be window size -1
             match_found = true;
             if mode == 1
                 break; % Exit inner loop after finding the first match 
             else
                 continue; %  -> use BREAK for one trial, CONTINUE for all trials
             end
         end
     end
     if ~match_found
         disp(['No match found for combination ', mat2str(permutation)]); 
         continue;
     end
     selected_trial_vector = zeros(1, length(lick_plot_filtered));
     idx_found_1 = 0 ;
     new_trials = [1; (trial_lengths(1:(length(trial_lengths)-1), 2))+1];  % Get the positions where a new trial starts
 
     for ii = 1:behavTrials.num_trials
         if perm_found(ii) == 1
             idx_found_1 = idx_found_1 + 1;  
             if idx_found_1 == 1 
                 find_first_position = new_trials(ii);
            
             elseif idx_found_1 == length(permutation)
                 if ii ~= behavTrials.num_trials
                    find_last_position = new_trials(ii+1);
                 end
                 if mode == 0
                    selected_trial_vector(find_first_position:find_last_position-1) = 1;
                 end
                 idx_found_1 = 0;
             end
         end
     end
     
     if mode == 1
        selected_trial_vector(find_first_position:find_last_position-1) = 1;
        % use when you just want one run for each permutation to be shown
     end

     %Plot the entire manifold in gray
     selected_trials = find(selected_trial_vector == 1);
     selected_data = Umap_results_filtered(selected_trials, :);
     selected_ports = lick_plot_filtered(selected_trials);
     selected_events = event_number_filtered(selected_trials); 
     selected_location = pos_plot_filtered(selected_trials);

     % for labelling colorbar
    % first_port_idx = selected_trials(1);
    % second_port_idx = selected_trials(round((length(selected_trials))/2));
    % third_port_idx = selected_trials(length(selected_trials));

     filt = Umap_results_filtered;
     filt(selected_trials, :) = [];
     figure;
     clf;
     scatter3(filt(:,dim1), filt(:,dim2), filt(:,dim3), dotSize, [0.85 0.85 0.85],'filled','MarkerFaceAlpha',0.4);
     hold on;
     view(A, E);
     scatter3(selected_data(:, dim1), selected_data(:, dim2), selected_data(:, dim3), dotSize * 4, selected_location , 'filled');
     colormap('parula');
     grid('off');
     axis('tight');
     axis('off');
     hold('on');
     colorbar;
     %colorbar('Ticks', [first_port_idx, second_port_idx, third_port_idx], 'TickLabels', {string(permutation(1)), string(permutation(2)), string(permutation(3))});
     permutation_str = strjoin(string(permutation), ' ');
     title(permutation_str);
 
 end
 
 
%% PLOT DELIBERATE vs EXPLORE
% Look at behavior using getPatchBehavior, and choose trials to inspect

% while in high patch, look at trials where action seems to be very
% deliberate (such as going back and forth between two highly rewarded
% ports
delib = [265, 296]; % indexes in behavTrials where mouse is being deliberate
delib = [delib; 124, 151]; % add as many instances as you want

delib_times = zeros(size(delib));

selected_trial_vector = zeros(1, length(lick_plot_filtered));

for j = 1:size(delib, 1)
    t_start = delib(j, 1);
    t_end = delib(j, 2);

    % find idx of closest timestamp
    [~, x_idx_1] = min(abs(behavTrials.timestamps(t_start)-timestamp_beh));
    [~, x_idx_2] = min(abs(behavTrials.timestamps(t_end)-timestamp_beh));

    delib_times(j, 1) = x_idx_1;
    delib_times(j, 2) = x_idx_2;

    selected_trial_vector(x_idx_1:x_idx_2) = 1;
    together = selected_trial_vector;
end


selected_trials = find(selected_trial_vector == 1);
selected_data = Umap_results_filtered(selected_trials, :);
selected_delib = Umap_results_filtered(selected_trials, :);
selected_ports = lick_plot_filtered(selected_trials);
selected_location = pos_plot_filtered(selected_trials);

filt = Umap_results_filtered; 
filt(selected_trials, :) = [];
figure;
clf;
scatter3(filt(:,dim1), filt(:,dim2), filt(:,dim3), dotSize, [0.85 0.85 0.85],'filled','MarkerFaceAlpha',0.4);
hold on;
view(A, E);
scatter3(selected_data(:, dim1), selected_data(:, dim2), selected_data(:, dim3), dotSize * 4, selected_location , 'filled');
colormap('parula');
grid('off');
axis('tight');
axis('off');
hold('on');
colorbar;
%colorbar('Ticks', [min(selected_location), max(selected_location)], 'TickLabels', {'Port 6', 'Port 7'}); % fix locations
title('Deliberate running between ports 5, 6, and 7');


% EXPLORE
% look for those same ports in sequence when the mouse is clearly exploring
% more and running between ports randomly or in sequence
%explore = [5, 7]; % indexes where mouse is exploring
%explore = [explore; 17, 23];
%explore = [explore; 29, 32];
%explore = [explore; 38, 41];
explore = [420, 423]; % indexes where mouse is exploring
%explore = [explore; 420, 423];
explore = [explore; 429, 433];
explore = [explore; 439, 443];
explore = [explore; 449, 452];
explore = [explore; 462, 465];

explore_times = zeros(size(explore));

selected_trial_vector = zeros(1, length(lick_plot_filtered));

for j = 1:size(explore, 1)
    t_start = explore(j, 1);
    t_end = explore(j, 2);

    % find idx of closest timestamp
    [~, x_idx_1] = min(abs(behavTrials.timestamps(t_start)-timestamp_beh));
    [~, x_idx_2] = min(abs(behavTrials.timestamps(t_end)-timestamp_beh));

    explore_times(j, 1) = x_idx_1;
    explore_times(j, 2) = x_idx_2;

    selected_trial_vector(x_idx_1:x_idx_2) = 1;
    together(x_idx_1:x_idx_2) = 2;
end

selected_trials = find(selected_trial_vector == 1);
selected_data = Umap_results_filtered(selected_trials, :);
selected_explore = Umap_results_filtered(selected_trials, :);
selected_ports = lick_plot_filtered(selected_trials);
selected_location = pos_plot_filtered(selected_trials);
selected_events = event_number_filtered(selected_trials); 


filt = Umap_results_filtered; 
filt(selected_trials, :) = [];
figure;
clf;
scatter3(filt(:,dim1), filt(:,dim2), filt(:,dim3), dotSize, [0.85 0.85 0.85],'filled','MarkerFaceAlpha',0.4);
hold on;
view(A, E);
scatter3(selected_data(:, dim1), selected_data(:, dim2), selected_data(:, dim3), dotSize * 4, selected_location , 'filled');
colormap('parula');
grid('off');
axis('tight');
axis('off');
hold('on');
colorbar('Ticks', [port_5_location, port_6_location, port_7_location], 'TickLabels', {'Port 5', 'Port 6', 'Port 7'}); % fix locations
title('Exploratory running between ports 5, 6, and 7 (by location)');

figure(9);
clf;
scatter3(filt(:,dim1), filt(:,dim2), filt(:,dim3), dotSize, [0.85 0.85 0.85],'filled','MarkerFaceAlpha',0.4);
hold on;
view(A, E);
scatter3(selected_data(:, dim1), selected_data(:, dim2), selected_data(:, dim3), dotSize * 4, selected_events , 'filled');
colormap('parula');
grid('off');
axis('tight');
axis('off');
hold('on');
colorbar;
title('Exploratory running between ports 5, 6, and 7 (by time)');


% deliberate and exploratory plotted together
selected_trials = find(together ~= 0);
selected_data = Umap_results_filtered(selected_trials, :);
selected_ports = lick_plot_filtered(selected_trials);
selected_location = pos_plot_filtered(selected_trials);
selected_events = event_number_filtered(selected_trials); 

customBlue = [0.2235, 0.2235, 0.8588];
customYellow = [0.9412, 0.9412, 0.3843];

filt = Umap_results_filtered; 
filt(selected_trials, :) = [];
figure;
clf;
scatter3(filt(:,dim1), filt(:,dim2), filt(:,dim3), dotSize, [0.85 0.85 0.85],'filled','MarkerFaceAlpha',0.4);
hold on;
view(A, E);
scatter3(selected_delib(:, dim1), selected_delib(:, dim2), selected_delib(:, dim3), dotSize * 4, customBlue, 'filled', 'DisplayName', 'Deliberate');
scatter3(selected_explore(:, dim1), selected_explore(:, dim2), selected_explore(:, dim3), dotSize * 4, customYellow, 'filled', 'DisplayName', 'Exploring');
%colormap('parula');
grid('off');
axis('tight');
axis('off');
hold('on');
title('Deliberate vs. exploratory running');
legend('show', 'Location', 'bestoutside');



%% Plot tracking 
figure
hold on
plot(timestamp_beh, position_y_all); 
scatter(timestamp_beh(reward_idx_2), position_y_all(reward_idx_2), 36, customGreen, 'filled');
scatter(timestamp_beh(nonreward_idx_2), position_y_all(nonreward_idx_2), 36, customRed, 'filled');

ylabel('Location');
xlabel('Time');
hold off


%% Create plot with just one direction

direction_filtered = Umap_results_filtered(direction_plot_filtered==0, :);

% Filter variables
pos_plot_filtered = pos_plot_filtered';
lick_plot_filtered = lick_plot_filtered';
trial_num_plot_filtered = trial_num_plot_filtered';
%direction_plot_filtered = direction_plot_filtered';
%speed_plot_filtered = speed_plot_filtered';
outcome_filtered = outcome_filtered';
outcome_lick_time_filtered = outcome_lick_time_filtered';
outcome_lick_high_patch_filtered = outcome_lick_high_patch_filtered';

direc_pos_plot_filtered = pos_plot_filtered(direction_plot_filtered==0, :);
direc_lick_plot_filtered = lick_plot_filtered(direction_plot_filtered==0, :);
direc_trial_num_plot_filtered = trial_num_plot_filtered(direction_plot_filtered==0, :);
direc_direction_plot_filtered = direction_plot_filtered(direction_plot_filtered==0, :);
direc_speed_plot_filtered = speed_plot_filtered(direction_plot_filtered==0, :);
direc_event_number_filtered = 1:length(direc_pos_plot_filtered);
direc_outcome_filtered = outcome_filtered(direction_plot_filtered==0, :);%plot
direc_outcome_lick_time_filtered = outcome_lick_time_filtered(direction_plot_filtered==0, :);
direc_outcome_lick_high_patch_filtered = outcome_lick_high_patch_filtered(direction_plot_filtered==0, :);


figure;
t = tiledlayout(2, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
title('Single Direction');

ax1 = nexttile;
ax2 = nexttile;
ax3 = nexttile;
ax4 = nexttile;
ax5 = nexttile;
ax6 = nexttile;

color = 'jet';
dotSize = 5;
% Plot1: colored by position
hold on
scatter3(ax1, direction_filtered(:,dim1), direction_filtered(:,dim2), direction_filtered(:,dim3), dotSize, direc_pos_plot_filtered, 'filled');
colormap(ax1, color);
view(ax1, A, E);
grid(ax1, 'off');
axis(ax1, 'tight');
axis(ax1, 'off');
colorbar(ax1);
title(ax1, 'Position');
hold off

% Plot 2: Colored by chosen port
scatter3(ax2, direction_filtered(:,dim1), direction_filtered(:,dim2), direction_filtered(:,dim3), dotSize, direc_lick_plot_filtered, 'filled');
colormap(ax2, color);
view(ax2, A, E);
grid(ax2, 'off');
axis(ax2, 'tight');
axis(ax2, 'off');
colorbar(ax2);
title(ax2, 'Chosen port');

% Plot 3: Colored by trial number
scatter3(ax3, direction_filtered(:,dim1), direction_filtered(:,dim2), direction_filtered(:,dim3), dotSize, direc_trial_num_plot_filtered, 'filled');
colormap(ax3, color);
view(ax3, A, E);
grid(ax3, 'off');
axis(ax3, 'tight');
axis(ax3, 'off');
colorbar(ax3);
title(ax3, 'Trial number');

% Plot 4: Colored by direction
scatter3(ax4, direction_filtered(:,dim1), direction_filtered(:,dim2), direction_filtered(:,dim3), dotSize, direc_direction_plot_filtered, 'filled');
colormap(ax4, color);
view(ax4, A, E);
grid(ax4, 'off');
axis(ax4, 'tight');
axis(ax4, 'off');
colorbar(ax4);
title(ax4, 'Direction');

%Plot 5: Colored by speed
scatter3(ax5, direction_filtered(:,dim1), direction_filtered(:,dim2), direction_filtered(:,dim3), dotSize, direc_speed_plot_filtered, 'filled');
colormap(ax5, color);
view(ax5, A, E);
grid(ax5, 'off');
axis(ax5, 'tight');
axis(ax5, 'off');
colorbar(ax5);
% caxis(ax5,[min(speed_plot_filtered) 30]);
title(ax5, 'Speed');


%Plot 6: Only show times of licks
umap_minus = direction_filtered;
umap_minus(trial_lengths(:,2), :) = [];
scatter3(ax6, umap_minus(:,dim1), umap_minus(:,dim2), umap_minus(:,dim3), ...
    dotSize, [0.85 0.85 0.85], 'filled', 'MarkerFaceAlpha', 0.4);
view(ax6, A, E);

cmap = jet(7);
colormap(ax6, cmap);
% caxis(ax2, [1 7]);  
grid(ax6, 'off');
axis(ax6, 'tight');
axis(ax6, 'off');
title(ax6, 'Licked Ports');
hold(ax6, 'on');
colorbar(ax6, 'Ticks', [0, 0.2, 0.35, 0.5, 0.65, 0.8, 1], 'TickLabels', {'1', '2', '3', '4','5', '6', '7'});

for i = 1:7
    idx_1 = behavTrials.port == i;
    idx_2 = trial_lengths(idx_1, 2);
    scatter3(ax6, direction_filtered(idx_2, dim1), direction_filtered(idx_2, dim2), direction_filtered(idx_2, dim3), ...
        dotSize, cmap(i, :), 'filled', 'MarkerFaceAlpha', 0.8);
end


x = 3;









 % set a better default view

