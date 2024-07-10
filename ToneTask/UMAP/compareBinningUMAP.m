
fig2 = figure;
TRIAL_TYPE = [0 1 2 3 4 5];
%col = 'jet';
col = [238/255 67/255 69/255;...
    241/255 114/255 42/255;...
    247/255 149/255 33/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];

umap_name = 'behavior_speed_0_smooth_5_bin_0.1_correct_only';
file = dir('*.position_behavior_speed_0_smooth_5_bin_0.1_correct_only.mat');
behav_file = file.name;

manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,...
    'numrow',1,'numcol',2,'rowloc',1,'colloc',1,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'error',1)

umap_name = 'behavior_speed_1_smooth_1_bin_0.1';
file = dir('*.position_behavior_speed_1_smooth_1_bin_0.1.mat');
behav_file = file.name;

manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,...
    'numrow',3,'numcol',2,'rowloc',2,'colloc',1,'col',col,'TRIAL_TYPE', TRIAL_TYPE)

umap_name = 'behavior_speed_1_smooth_5_bin_0.02';
file = dir('*.position_behavior_speed_1_smooth_5_bin_0.02.mat');
behav_file = file.name;
% 
manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,...
    'numrow',3,'numcol',2,'rowloc',3,'colloc',1,'col',col,'TRIAL_TYPE', TRIAL_TYPE)

% manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,...
%     'numrow',3,'numcol',2,'rowloc',3,'colloc',1,'col',col,'TRIAL_TYPE', TRIAL_TYPE)
% 
% ax1 = subplot(3,2,1,'Parent',fig2);
% ax2 = subplot(3,2,2,'Parent',fig2);
% Link = linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
% setappdata(gcf, 'StoreTheLink', Link);
% 
% 
% ax1 = subplot(3,2,3,'Parent',fig2);
% ax2 = subplot(3,2,4,'Parent',fig2);
% Link = linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
% setappdata(gcf, 'StoreTheLink', Link);
% 
% 
% ax1 = subplot(3,2,5,'Parent',fig2);
% ax2 = subplot(3,2,6,'Parent',fig2);
% Link = linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
% setappdata(gcf, 'StoreTheLink', Link);

% trial1 = 68;
% trial2 = 66;
% umap_name = 'behavior_speed_1_smooth_5_bin_0.02';
% file = dir('*.position_behavior_speed_1_smooth_5_bin_0.02.mat');
% behav_file = file.name;
% A = -271.4603;
% E = -90;
% manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,...
%     'numrow',2,'numcol',3,'rowloc',1,'colloc',1,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'singleTrial',true,'tsWin',behavTrials.timestamps(trial1,:),'A',A,'E',E)
% 
% manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,...
%     'numrow',2,'numcol',3,'rowloc',2,'colloc',1,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'singleTrial',true,'tsWin',behavTrials.timestamps(trial2,:),'A',A,'E',E)
% 
% %51, 26
% umap_name = 'behavior_speed_1_smooth_1_bin_0.1';
% file = dir('*.position_behavior_speed_1_smooth_1_bin_0.1.mat');
% behav_file = file.name;
% A = -179.4570;
% E = -90;
% manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,...
%     'numrow',2,'numcol',3,'rowloc',1,'colloc',2,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'singleTrial',true,'tsWin',behavTrials.timestamps(trial1,:),'A',A,'E',E)
% 
% manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,...
%     'numrow',2,'numcol',3,'rowloc',2,'colloc',2,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'singleTrial',true,'tsWin',behavTrials.timestamps(trial2,:),'A',A,'E',E)
% 
% 
% umap_name = 'behavior_speed_1_smooth_5_bin_0.1';
% file = dir('*.position_behavior_speed_1_smooth_5_bin_0.1.mat');
% behav_file = file.name;
% A = -93.25;
% E = -84.0809;
% manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,...
%     'numrow',2,'numcol',3,'rowloc',1,'colloc',3,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'singleTrial',true,'tsWin',behavTrials.timestamps(trial1,:),'A',A,'E',E)
% 
% manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,...
%     'numrow',2,'numcol',3,'rowloc',2,'colloc',3,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'singleTrial',true,'tsWin',behavTrials.timestamps(trial2,:),'A',A,'E',E)
