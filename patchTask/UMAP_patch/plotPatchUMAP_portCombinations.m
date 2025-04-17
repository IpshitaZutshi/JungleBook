function plotPatchUMAP_portCombinations(varargin)

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
A = 0; E = 90; % azimuth and elevation for viewing

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

pos_x_filtered = position_x_all(filtered_idx);
pos_y_filtered = position_y_all(filtered_idx);

time_filtered = timestamp_beh(filtered_idx);
[~,x]=min(abs(time_filtered - behavTrials.timestamps(5)));

% Create multiplot

figure;
set(gcf,'Color','w')
set(gcf,'Position',[22 65 1650 940])
set(gcf,'Renderer','painters')
t = tiledlayout(2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
title(t, currentFolderName);

ax1 = nexttile;
ax2 = nexttile;
ax3 = nexttile;
ax4 = nexttile;

Link = linkprop([ax1, ax2, ax3, ax4], {'CameraPosition', 'CameraTarget', 'CameraUpVector', 'XLim', 'YLim', 'ZLim'});
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

% Plot 3: Colored by direction
scatter3(ax3, Umap_results_filtered(:,dim1), Umap_results_filtered(:,dim2), Umap_results_filtered(:,dim3), dotSize, direction_plot_filtered, 'filled');
colormap(ax3, color);
view(ax3, A, E);
grid(ax3, 'off');
axis(ax3, 'tight');
axis(ax3, 'off');
colorbar(ax3);
title(ax3, 'Direction');

%Plot 4: Only show times of licks
scatter3(ax4, Umap_results_filtered(:,dim1), Umap_results_filtered(:,dim2), Umap_results_filtered(:,dim3), ...
    dotSize, [0.85 0.85 0.85], 'filled', 'MarkerFaceAlpha', 0.4);
view(ax4, A, E);

cmap = jet(7);
colormap(ax4, cmap); 
grid(ax4, 'off');
axis(ax4, 'tight');
axis(ax4, 'off');
title(ax4, 'Licked Ports');
hold(ax4, 'on');
colorbar(ax4, 'Ticks', [0, 0.2, 0.35, 0.5, 0.65, 0.8, 1], 'TickLabels', {'1', '2', '3', '4','5', '6', '7'});

for i = 1:7
    ts = behavTrials.timestamps(behavTrials.port == i);
    for tt = 1:length(ts)
        [~, tsIdx] = min(abs(time_filtered-ts(tt)));
        scatter3(ax4, Umap_results_filtered(tsIdx, dim1), Umap_results_filtered(tsIdx, dim2), Umap_results_filtered(tsIdx, dim3), ...
            20, cmap(i, :), 'filled');
    end
end

% expPath = 'C:\Users\Ipshita\NYU Langone Health Dropbox\Ipshita Zutshi\Patch Task\Figures\';
% saveas(gcf,strcat(expPath,'UMAP_N7_sess23.png'));
% saveas(gcf,strcat(expPath,'UMAP_N7_sess23.eps'),'epsc');
% saveas(gcf,strcat(expPath,'UMAP_N7_sess23.fig'));

%% PLOT specific port trajectories
% Define permutations
permutations = {[7,6,5], [7,5], [7,6,7]};

%results = cell(size(permutations));  % Preallocate results

% for p = 1:length(permutations)
%     perm = permutations{p};
%     len = length(perm);
%     matches = [];
%     
%     % Loop through behavTrials.port and check for occurrences
%     for i = 1:(length(behavTrials.port) - len + 1)
%         if isequal(behavTrials.port(i:i+len-1)', perm)
%             matches = [matches; i, i+len-1];  % Store start and end indices
%         end
%     end
%     
%     results{p} = matches;  % Store results in a cell array
% end

results{1} = [112 114; 258 259; 272 273];
results{2} = [344 346; 337 338];


trial1 = [7 4 48];
col_gray = [0.3 0.3 0.3; 0.5 0.5 0.5; 0.8 0.8 0.8];
col_red = [107/255 0/255 0/255; 186/255 0/255 0/255; 254/255 77/255 77/255];

for ii = 1:length(results)

    %trial=trial1(ii);

    figHandle = figure;
    set(gcf,'Color','w')
    set(gcf,'Position',[53 260 1820 655])
    set(gcf,'Renderer','painters')

    for kk = 1:size(results{ii},1)

        %Plot the trajectory of the mouse
        tsWin = behavTrials.timestamps(results{ii}(kk,:));
        [~,startIdx] = min(abs(tracking.timestamps-tsWin(1)));
        [~,endIdx] = min(abs(tracking.timestamps-tsWin(2)));
        
        subplot(5,size(results{ii},1),[1:1:size(results{ii},1)])
        hold on
        plot(tracking.position.y(startIdx-2:endIdx+2),tracking.position.x(startIdx-2:endIdx+2),'Color',col_gray(kk,:),'LineWidth',1.5)
        for pp = results{ii}(kk,1):1:results{ii}(kk,2)
            [~,lickIdx] = min(abs(tracking.timestamps-behavTrials.timestamps(pp)));
            scatter(tracking.position.y(lickIdx),tracking.position.x(lickIdx),40,col_red(kk,:),"filled");
        end
        ylim([0 7])
        xlim([30 123])       
        
        ax1 = subplot(5,size(results{ii},1),size(results{ii},1)+kk);
        hold on
        yPos = tracking.position.y(startIdx-2:endIdx+2);
        tsAxis = linspace(tsWin(1)-2,tsWin(2)+2,length(yPos));
        scatter(tracking.position.y(startIdx-2:endIdx+2),tracking.position.x(startIdx-2:endIdx+2),10,tsAxis',"filled");
        for pp = results{ii}(kk,1):1:results{ii}(kk,2)
            [~,lickIdx] = min(abs(tracking.timestamps-behavTrials.timestamps(pp)));
            scatter(tracking.position.y(lickIdx),tracking.position.x(lickIdx),40,col_red(kk,:),"filled");
        end
        ylim([0 7])
        xlim([30 123])
        colormap(ax1,'jet');
        
        %Plot the UMAP in gray
        ax1 = subplot(5,size(results{ii},1),[size(results{ii},1)*2+kk size(results{ii},1)*3+kk size(results{ii},1)*4+kk],'Parent', figHandle);
        hold on
        scatter3(Umap_results_filtered(:,dim1),Umap_results_filtered(:,dim2),Umap_results_filtered(:,dim3),2,[0.85 0.85 0.85],'filled','MarkerFaceAlpha',0.4);
        view(A,E)
        grid off;
        axis off;
        axis tight

        %Now plot the trial block of interest
        [~,startidx] = min(abs(time_filtered-tsWin(1)));
        [~,endidx] = min(abs(time_filtered-tsWin(2)));
        plot_ind_pos = zeros(1,length(pos_plot_filtered));
        plot_ind_pos(startidx:1:endidx) = 1;
        plot_ind_pos = logical(plot_ind_pos);
        tsAxis = linspace(tsWin(1),tsWin(2),sum(plot_ind_pos));
   
        scatter3(Umap_results_filtered(plot_ind_pos,dim1),Umap_results_filtered(plot_ind_pos,dim2),Umap_results_filtered(plot_ind_pos,dim3),40,tsAxis','filled'); 
        for pp = results{ii}(kk,1):1:results{ii}(kk,2)
            [~,lickIdx] = min(abs(time_filtered-behavTrials.timestamps(pp)));
            scatter3(Umap_results_filtered(lickIdx,dim1),Umap_results_filtered(lickIdx,dim2),Umap_results_filtered(lickIdx,dim3),60,[186/255 0/255 0/255],'filled'); 
        end
        colormap(ax1,'jet');
        colorbar    
        %title(num2str(results{ii}(kk,1):1:results{ii}(kk,2)))

    end
        
    % 
    expPath = 'C:\Users\Ipshita\NYU Langone Health Dropbox\Ipshita Zutshi\Patch Task\Figures';
    saveas(gcf,strcat(expPath,'UMAP_N7_sess23_singletrials_',num2str(ii),'.png'));
    saveas(gcf,strcat(expPath,'UMAP_N7_sess23_singletrials_',num2str(ii),'.eps'),'epsc');
    saveas(gcf,strcat(expPath,'UMAP_N7_sess23_singletrials_',num2str(ii),'.fig'));
end


