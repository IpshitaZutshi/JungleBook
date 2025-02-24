function plotPatchUMAP(varargin)

p = inputParser;
addParameter(p,'min_trial',0,@isnumeric);
addParameter(p,'max_trial',1000,@isnumeric);
parse(p,varargin{:});
min_trial = p.Results.min_trial; % confused about this
max_trial = p.Results.max_trial;

basepath = pwd;
umapFile = dir(fullfile(basepath, '*.csv'));
behavFile = dir(fullfile(basepath, '*.position_behavior_speed1.mat'));
behavior = dir(fullfile(basepath, '*.TrialBehavior.mat'));

currentFolder = pwd;
[~, currentFolderName] = fileparts(currentFolder);

umap_path = umapFile.folder;
behav_file = behavFile.name;
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
load(behav_file);
load(behavior.name);

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

licks_plot = licks_ds; 
licks_plot(isnan(licks_plot))=0;

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
licks_plot_filtered = licks_plot(filtered_idx); % REDUNDANT
Umap_results_filtered = Umap_results(filtered_idx, :);
event_number_filered = 1:length(pos_plot_filtered);
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


%Plot the entire manifold in gray
scatter3(ax6, Umap_results_filtered(:,dim1), Umap_results_filtered(:,dim2), Umap_results_filtered(:,dim3), dotSize, [0.85 0.85 0.85],'filled','MarkerFaceAlpha',0.4);
view(ax6, A, E);

cmap = jet(7);
colormap(ax6, cmap);
%clim(ax6, [1 7]);  % was caxis
grid(ax6, 'off');
axis(ax6, 'tight');
axis(ax6, 'off');
title(ax6, 'Licked Ports');
hold(ax6, 'on');
colorbar(ax6);
    for i = 1:7
        idx = lick_plot_filtered == i;
        scatter3(ax6, Umap_results_filtered(idx, dim1), Umap_results_filtered(idx, dim2), Umap_results_filtered(idx, dim3), ...
                 dotSize, cmap(i, :), 'filled', 'MarkerFaceAlpha', 0.8);
    end
%hold(ax6, 'off');

% % Save figure
% figFile = fullfile('/Users/lauraribalta/Desktop', [currentFolderName '_UMAP_Multiplot.fig']);
% savefig(figFile);

%% PLOT OUTCOME
figure(2);
t = tiledlayout(2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact'); % 1 fila, 3 columnas
title(t, 'Multiplot');

% Subplot 1: Outcome
ax1 = nexttile;
%scatter3(ax1, Umap_results_filtered(:,dim1), Umap_results_filtered(:,dim2), Umap_results_filtered(:,dim3), ...
   % dotSize, [0.85 0.85 0.85], 'filled', 'MarkerFaceAlpha', 0.4);
view(ax1, A, E);
grid(ax1, 'off');
axis(ax1, 'tight');
axis(ax1, 'off');
title(ax1, 'Outcome');
hold(ax1, 'on');


customGreen = [152, 194, 9] / 255;
customRed = [238, 75, 43] / 255;
idx1 = outcome_lick_time_filtered == 1; 
idx2 = outcome_lick_time_filtered == 2;

scatter3(ax1, Umap_results_filtered(idx1, dim1), Umap_results_filtered(idx1, dim2), Umap_results_filtered(idx1, dim3), ...
    dotSize, customRed, 'filled', 'MarkerFaceAlpha', 0.8, 'DisplayName', 'Non-Rewarded');

scatter3(ax1, Umap_results_filtered(idx2, dim1), Umap_results_filtered(idx2, dim2), Umap_results_filtered(idx2, dim3), ...
    dotSize, customGreen, 'filled', 'MarkerFaceAlpha', 0.8, 'DisplayName', 'Rewarded');

legend(ax1, 'show', 'Location', 'bestoutside');

% Subplot 2: Licked ports
ax2 = nexttile;
scatter3(ax2, Umap_results_filtered(:,dim1), Umap_results_filtered(:,dim2), Umap_results_filtered(:,dim3), ...
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
colorbar(ax2);

for i = 1:7
    idx = lick_plot_filtered == i;
    scatter3(ax2, Umap_results_filtered(idx, dim1), Umap_results_filtered(idx, dim2), Umap_results_filtered(idx, dim3), ...
        dotSize, cmap(i, :), 'filled', 'MarkerFaceAlpha', 0.8);
end


% Subplot 3: Trial Number
ax3 = nexttile;
scatter3(ax3, Umap_results_filtered(:,dim1), Umap_results_filtered(:,dim2), Umap_results_filtered(:,dim3), ...
    dotSize, trial_num_plot_filtered, 'filled');
colormap(ax3, 'jet'); % Puedes cambiar el colormap si prefieres otro
view(ax3, A, E);
grid(ax3, 'off');
axis(ax3, 'tight');
axis(ax3, 'off');
colorbar(ax3);
title(ax3, 'Trial Number');

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

customGreen = [152, 194, 9] / 255;
customRed = [238, 75, 43] / 255;
idx1 = outcome_lick_high_patch_filtered == 1; 
idx2 = outcome_lick_high_patch_filtered == 2;

scatter3(ax4, Umap_results_filtered(idx1, dim1), Umap_results_filtered(idx1, dim2), Umap_results_filtered(idx1, dim3), ...
    dotSize, customRed, 'filled', 'MarkerFaceAlpha', 0.8, 'DisplayName', 'Non-Rewarded');

scatter3(ax4, Umap_results_filtered(idx2, dim1), Umap_results_filtered(idx2, dim2), Umap_results_filtered(idx2, dim3), ...
    dotSize, customGreen, 'filled', 'MarkerFaceAlpha', 0.8, 'DisplayName', 'Rewarded');

hold(ax4, 'off');

Link = linkprop([ax1, ax2, ax3, ax4], {'CameraPosition', 'CameraTarget', 'CameraUpVector', 'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', Link);


 %% PLOT TRAJECTORIES
 ports = [5, 6, 7];
 combinations = perms(ports); % Generate all permutations
 
 change_indices = [true, diff(lick_plot_filtered) ~= 0]; % start of new trial i believe
 individual_numbers = lick_plot_filtered(change_indices); % lick_plot_filtered is cutting out first trial
 % individual_numbers is just the port number with one cut out -> not sure why this variable
 % is necessary
 logical_vector = zeros(1, length(individual_numbers));
 
 
  % sliding window
 
 for i = 1:size(combinations, 1)
     matched_combination = [];
     selected_trial_vector = [];
     find_first_position = NaN;
     find_last_position = NaN;
     combination = combinations(i, :);
     % combination = [5, 7];
     window_size = length(combination);
     match_found = false;  
 
     for start_idx = 1:(length(individual_numbers) - window_size + 1)
         current_window = individual_numbers(start_idx:start_idx + window_size - 1);
         if isequal(current_window, combination)
             logical_vector(start_idx:start_idx + window_size) = 1; % should it be window size -1
             matched_combination = combination;
             match_found = true;
             break; % Exit inner loop after finding the first match
         end
     end
     if ~match_found
         disp(['No match found for combination ', mat2str(combination)]); 
         continue;
     end
     selected_trial_vector = zeros(1, length(lick_plot_filtered));
     idx_found_1 = 0 ;
     positions_of_ones = find(change_indices);  % Get the positions where change_indices is 1 -> a new trial starts
 
     for ii = 1:length(logical_vector)
         if logical_vector(ii) == 1
             idx_found_1 = idx_found_1 + 1;  
             if idx_found_1 == 1 
                 find_first_position = positions_of_ones(ii);
             end
 
             if idx_found_1 == length(combination)
                 find_last_position = positions_of_ones(ii+1);
             end
         end
     end
     selected_trial_vector(find_first_position:find_last_position-1) = 1;
 
     %Plot the entire manifold in gray
     selected_trials = find(selected_trial_vector == 1);
     selected_data = Umap_results_filtered(selected_trials, :);
     selected_ports = lick_plot_filtered(selected_trials);
     selected_events = event_number_filered(selected_trials);
 
     figure;
     clf;
     scatter3(Umap_results_filtered(:,dim1), Umap_results_filtered(:,dim2), Umap_results_filtered(:,dim3), dotSize, [0.85 0.85 0.85],'filled','MarkerFaceAlpha',0.4);
     hold on;
     view(A, E);
     scatter3(selected_data(:, dim1), selected_data(:, dim2), selected_data(:, dim3), dotSize * 4, selected_events , 'filled');
     colormap('parula');
     grid('off');
     axis('tight');
     axis('off');
     title('Licked Ports');
     hold('on');
     colorbar;
     combination_str = strjoin(string(combination), ' ');
     title(combination_str);
 
 end
 
 
 % plot random trials, not only the first one
 % set a better default view