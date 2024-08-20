fig2 = figure;
set(fig2, 'Renderer','painters')
set(fig2,'Position',[21 177 1366 579])

umap_path = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ43_220915_sess13\manifold';
%umap_path = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\finalSessions\Task\IZ48_230714_sess28\manifold';
cd(umap_path)
sessname  = strsplit(umap_path,'\');
umap_name = 'behavior_speed_1_smooth_5_bin_0.1';
behav_file = strcat(umap_path,'\',sessname{7},'.position_behavior_speed_1_smooth_5_bin_0.1.mat');
A = 81.8716;%-89;%-179;
E = 3.7413;%7.68;%-22;

% Next, tone non-stim trials
TRIAL_TYPE = [0:5];
%col = 'jet';
col = [238/255 67/255 69/255;...
    241/255 114/255 42/255;...
    247/255 149/255 33/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];
manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,'addPosPlot',true,'poscutOff',0,'speedThresh',1,'addFreq',true,...
    'numrow',1,'numcol',3,'rowloc',1,'colloc',1,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E)

% Next, tone stim trials
% TRIAL_TYPE = [0:5];
% %col = 'jet';
% col = [238/255 67/255 69/255;...
%     241/255 114/255 42/255;...
%     247/255 149/255 33/255;...
%     249/255 197/255 81/255;...
%     143/255 189/255 107/255;...
%     87/255 116/255 144/255];
% manifoldPlot_error('figHandle',fig2,'umap_path',umap_path,'behav_file',behav_file,...
%     'numrow',2,'numcol',3,'rowloc',2,'colloc',1,'probe',true,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E)


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

