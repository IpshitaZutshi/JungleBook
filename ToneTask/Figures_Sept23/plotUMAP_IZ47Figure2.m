fig2 = figure;
set(fig2, 'Renderer','painters')
set(fig2,'Position',[21 177 1366 579])

%umap_path = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\finalSessions\Task\IZ47_230710_sess25\manifold';
umap_path = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ47_230626_sess15\manifold';
sessname  = strsplit(umap_path,'\');
behav_file = strcat(umap_path,'\',sessname{7},'.position_behavior_speed_1_smooth_5_bin_0.1.mat');
umap_name = 'behavior_speed_1_smooth_5_bin_0.1';
% A = 0;%-89;%
% E = 90;%7.68;%
A = -1.17; E = -31.94;
% 
% % First, entire manifold
% TRIAL_TYPE = [0 1 2 3 4 5 6 7 8];
% manifoldPlot('figHandle',fig2,'umap_path',umap_path,'behav_file',behav_file,'umap_name',umap_name,...
%     'numrow',1,'numcol',2,'rowloc',1,'colloc',1,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E)

% Next, tone trials
TRIAL_TYPE = [0:5];
%col = 'jet';
fig2 = figure;
set(fig2, 'Renderer','painters')
set(fig2,'Position',[21 177 1366 579])
set(fig2,'Color','k')

col = [117/255 26/255 51/255;...
    184/255 15/255 10/255;...
    241/255 114/255 42/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];

manifoldPlot('figHandle',fig2,'umap_path',umap_path,'behav_file',behav_file,'umap_name',umap_name,'poscutOff',7,'speedThresh',2,...
    'numrow',1,'numcol',3,'rowloc',1,'colloc',1,'col',col,'addFreq',true,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E)
% 
% % Next, tone 6 and no tone trials
% TRIAL_TYPE = [5 6 8];
% col = [0/255 0/255 255/255;...
%     0/255 0/255 0/255;...
%    224/255 163/255 46/255];
% A = 78;%-179
% E = -3.8;%-22
% manifoldPlot('figHandle',fig2,'umap_path',umap_path,'behav_file',behav_file,'umap_name',umap_name,...
%     'numrow',3,'numcol',3,'rowloc',1,'colloc',3,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E)

% ax1 = subplot(3,3,1,'Parent',fig2);
% ax2 = subplot(3,3,2,'Parent',fig2);
% ax3 = subplot(3,3,3,'Parent',fig2);
% ax4 = subplot(3,3,4,'Parent',fig2);
% ax5 = subplot(3,3,5,'Parent',fig2);
% ax6 = subplot(3,3,6,'Parent',fig2);
% ax7 = subplot(3,3,7,'Parent',fig2);
% ax8 = subplot(3,3,8,'Parent',fig2);
% ax9 = subplot(3,3,9,'Parent',fig2);
% 
% Link = linkprop([ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
% setappdata(gcf, 'StoreTheLink', Link);

% view(ax1,-2,90)
% view(ax2,-2,90)
%view(ax5,-2,270)

