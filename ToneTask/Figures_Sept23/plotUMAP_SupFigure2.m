umap_path = {'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\finalSessions\Task\IZ39_220714_sess18\manifold',...
    'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\finalSessions\Task\IZ40_220714_sess18\manifold',...
    'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\finalSessions\Task\IZ43_220828_sess4\manifold',...
    'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\finalSessions\Task\IZ44_220830_sess7\manifold',...
    'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\finalSessions\Task\IZ48_230714_sess28\manifold'};...
%'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\finalSessions\Task\IZ47_230710_sess25\manifold'

warning off

A_all = [0.82,93.4153,-182.76,191,23];
E_all = [-90,90,-1.55,-11,1.55];

A_all2 = [176.3896,176.3896,-154,191,-179];
E_all2 = [17.3750,17.3750,2.33,-11,-22];

for ii = 1:length(umap_path)

    fig2 = figure;
    set(fig2, 'Renderer','painters')
    set(fig2,'Position',[1 41 1920 970])
        
    sessname  = strsplit(umap_path{ii},'\');
    behav_file = strcat(umap_path{ii},'\',sessname{8},'.position_behavior_speed_1_smooth_5.mat');
    umap_name = 'behavior_speed_1_smooth_5';

    if ii == 5
        behav_file = strcat(umap_path{ii},'\',sessname{8},'.position_behavior_speed_1_smooth_10.mat');
        umap_name = 'behavior_speed_1_smooth_10';
    end
    A = A_all(ii);
    E = E_all(ii);
    % 
    % % First, entire manifold
    % TRIAL_TYPE = [0 1 2 3 4 5 6 7 8];
    % manifoldPlot('figHandle',fig2,'umap_path',umap_path,'behav_file',behav_file,...
    %     'numrow',3,'numcol',3,'rowloc',1,'colloc',1,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E)
    
    % Next, tone trials
    TRIAL_TYPE = [0:5];
    %col = 'jet';
    col = [238/255 67/255 69/255;...
        241/255 114/255 42/255;...
        247/255 149/255 33/255;...
        249/255 197/255 81/255;...
        143/255 189/255 107/255;...
        87/255 116/255 144/255];
    manifoldPlot('figHandle',fig2,'umap_path',umap_path{ii},'umap_name',umap_name,'behav_file',behav_file,'addFreq',true,...
         'numrow',2,'numcol',3,'rowloc',1,'colloc',1,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E)


    % Next, tone 6 and no tone trials
    TRIAL_TYPE = [5 6 8];
    col = [0/255 0/255 255/255;...
        0/255 0/255 0/255;...
       224/255 163/255 46/255];
    A = A_all2(ii);
    E = E_all2(ii);
    % manifoldPlot('figHandle',fig2,'umap_path',umap_path{ii},'umap_name',umap_name,'behav_file',behav_file,...
    %     'numrow',length(umap_path),'numcol',4,'rowloc',ii,'colloc',3,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E)
    manifoldPlot('figHandle',fig2,'umap_path',umap_path{ii},'umap_name',umap_name,'behav_file',behav_file,...
    'numrow',2,'numcol',3,'rowloc',2,'colloc',1,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E)
    % ax1 = subplot(1,4,3,'Parent',fig2);
    % ax2 = subplot(1,4,4,'Parent',fig2);
    % Link = linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
    % setappdata(gcf, 'StoreTheLink', Link);

    % view(ax1,-2,90)
    % view(ax2,-2,90)
    %view(ax5,-2,270)

    %% Save figure
    % expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';
    % saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure2_UMAP',num2str(ii),'.png'));
    % saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure2_UMAP',num2str(ii),'.eps'),'epsc');
    % saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure2_UMAP',num2str(ii),'.fig'));
    % close all
end




