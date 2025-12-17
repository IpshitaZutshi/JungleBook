fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')

numrows = 4;
numcol = 2;

%% UMAP of probe trial
TRIAL_TYPE = [0:5];
col = [83/255 0/255 0/255;...
    184/255 15/255 10/255;...
    241/255 114/255 42/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];

umap_path = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ39_220707_sess17\manifold';
sessname  = strsplit(umap_path,'\');
behav_file = strcat(umap_path,'\',sessname{7},'.position_behavior_speed_1_smooth_5_bin_0.1.mat');
load(behav_file)
umap_name = 'behavior_speed_1_smooth_5_bin_0.1';

file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*Tracking.Behavior.mat']);
load(file.name);

% Add stim mask
stimMask = zeros(size(tracking.timestamps));
for ii = 1:(length(behavTrials.timestamps)-1)
    posTrials = tracking.timestamps >= behavTrials.timestamps(ii,1) & ...
                tracking.timestamps < behavTrials.timestamps(ii+1,1); % maybe i should keep all timestamps in behavTrials

    stimMask(posTrials) =  behavTrials.stim(ii);
end
stim_ds = interp1(tracking.timestamps,stimMask,timestamp_beh,'nearest');

A = -45;
E = 2.7;

manifoldPlot_stim('figHandle',fig2,'umap_path',umap_path,'behav_file',behav_file,'umap_name',umap_name,...
    'numrow',numrows,'numcol',numcol,'rowloc',1,'colloc',1,'stim',false,'error',false,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'stim_ds',stim_ds)
manifoldPlot_stim('figHandle',fig2,'umap_path',umap_path,'behav_file',behav_file,'umap_name',umap_name,...
    'numrow',numrows,'numcol',numcol,'rowloc',2,'colloc',1,'stim',false,'error',true,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'stim_ds',stim_ds)
manifoldPlot_stim('figHandle',fig2,'umap_path',umap_path,'behav_file',behav_file,'umap_name',umap_name,...
    'numrow',numrows,'numcol',numcol,'rowloc',3,'colloc',1,'stim',true,'error',false,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'stim_ds',stim_ds)
manifoldPlot_stim('figHandle',fig2,'umap_path',umap_path,'behav_file',behav_file,'umap_name',umap_name,...
    'numrow',numrows,'numcol',numcol,'rowloc',4,'colloc',1,'stim',true,'error',true,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'stim_ds',stim_ds)

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')

numrows = 4;
numcol = 2;
%% UMAP of probe trial
TRIAL_TYPE = [0:5];
col = [83/255 0/255 0/255;...
    184/255 15/255 10/255;...
    241/255 114/255 42/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];


umap_path = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ43_220915_sess13\manifold';
sessname  = strsplit(umap_path,'\');
behav_file = strcat(umap_path,'\',sessname{7},'.position_behavior_speed_1_smooth_5_bin_0.1.mat');
load(behav_file)
umap_name = 'behavior_speed_1_smooth_5_bin_0.1';

file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*Tracking.Behavior.mat']);
load(file.name);

% Add stim mask
stimMask = zeros(size(tracking.timestamps));
for ii = 1:(length(behavTrials.timestamps)-1)
    posTrials = tracking.timestamps >= behavTrials.timestamps(ii,1) & ...
                tracking.timestamps < behavTrials.timestamps(ii+1,1); % maybe i should keep all timestamps in behavTrials

    stimMask(posTrials) =  behavTrials.stim(ii);
end
stim_ds = interp1(tracking.timestamps,stimMask,timestamp_beh,'nearest');

A = -224;
E = 3.6;

manifoldPlot_stim('figHandle',fig2,'umap_path',umap_path,'behav_file',behav_file,'umap_name',umap_name,...
    'numrow',numrows,'numcol',numcol,'rowloc',1,'colloc',1,'stim',false,'error',false,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'stim_ds',stim_ds)
manifoldPlot_stim('figHandle',fig2,'umap_path',umap_path,'behav_file',behav_file,'umap_name',umap_name,...
    'numrow',numrows,'numcol',numcol,'rowloc',2,'colloc',1,'stim',false,'error',true,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'stim_ds',stim_ds)
manifoldPlot_stim('figHandle',fig2,'umap_path',umap_path,'behav_file',behav_file,'umap_name',umap_name,...
    'numrow',numrows,'numcol',numcol,'rowloc',3,'colloc',1,'stim',true,'error',false,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'stim_ds',stim_ds)
manifoldPlot_stim('figHandle',fig2,'umap_path',umap_path,'behav_file',behav_file,'umap_name',umap_name,...
    'numrow',numrows,'numcol',numcol,'rowloc',4,'colloc',1,'stim',true,'error',true,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'stim_ds',stim_ds)

