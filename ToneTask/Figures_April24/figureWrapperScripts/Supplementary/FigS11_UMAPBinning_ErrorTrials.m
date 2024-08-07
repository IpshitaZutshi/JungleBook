function FigS11_UMAPBinning_ErrorTrials

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[680 65 1225 962]);

numrow = 6;
numcol = 6;

TRIAL_TYPE = [0 1 2 3 4 5];
col = [83/255 0/255 0/255;...
        184/255 15/255 10/255;...
        241/255 114/255 42/255;...
        249/255 197/255 81/255;...
        143/255 189/255 107/255;...
        87/255 116/255 144/255];

%% IZ47
umap_path = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ47_230707_sess24\manifold';
cd(umap_path)
umap_name = 'behavior_speed_1_smooth_1_bin_0.02';
file = dir('*.position_behavior_speed_1_smooth_1_bin_0.02.mat');
behav_file = file.name;
A = -118;
E = 80;
manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,'poscutOff',5,'speedThresh',2,...
    'numrow',numrow,'numcol',numcol,'rowloc',1,'colloc',1,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'dotSize',2)

umap_name = 'behavior_speed_1_smooth_5_bin_0.02';
file = dir('*.position_behavior_speed_1_smooth_5_bin_0.02.mat');
behav_file = file.name;
A = -2.5881;
E = 23.8308;
manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,'poscutOff',5,'speedThresh',2,...
    'numrow',numrow,'numcol',numcol,'rowloc',2,'colloc',1,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'dotSize',2)

umap_name = 'behavior_speed_1_smooth_1_bin_0.1';
file = dir('*.position_behavior_speed_1_smooth_1_bin_0.1.mat');
behav_file = file.name;
A = -185.7439;
E = -77.644;
manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,'poscutOff',5,'speedThresh',2,...
    'numrow',numrow,'numcol',numcol,'rowloc',3,'colloc',1,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'dotSize',2)

umap_name = 'behavior_speed_1_smooth_5_bin_0.1';
file = dir('*.position_behavior_speed_1_smooth_5_bin_0.1.mat');
behav_file = file.name;
A = -74;
E = -90;
manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,'poscutOff',5,'speedThresh',2,...
    'numrow',numrow,'numcol',numcol,'rowloc',4,'colloc',1,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'dotSize',4)

%% IZ48
umap_path = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ48_230628_sess17\manifold';
cd(umap_path)
umap_name = 'behavior_speed_1_smooth_1_bin_0.02';
file = dir('*.position_behavior_speed_1_smooth_1_bin_0.02.mat');
behav_file = file.name;
A = 70;
E = -9;
manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,'poscutOff',5,'speedThresh',2,...
    'numrow',numrow,'numcol',numcol,'rowloc',1,'colloc',3,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'dotSize',2)

umap_name = 'behavior_speed_1_smooth_5_bin_0.02';
file = dir('*.position_behavior_speed_1_smooth_5_bin_0.02.mat');
behav_file = file.name;
A = 11;
E = -90;
manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,'poscutOff',5,'speedThresh',2,...
    'numrow',numrow,'numcol',numcol,'rowloc',2,'colloc',3,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'dotSize',2)


umap_name = 'behavior_speed_1_smooth_1_bin_0.1';
file = dir('*.position_behavior_speed_1_smooth_1_bin_0.1.mat');
behav_file = file.name;
A = -208;
E = 90;
manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,'poscutOff',5,'speedThresh',2,...
    'numrow',numrow,'numcol',numcol,'rowloc',3,'colloc',3,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'dotSize',2)


umap_name = 'behavior_speed_1_smooth_5_bin_0.1';
file = dir('*.position_behavior_speed_1_smooth_5_bin_0.1.mat');
behav_file = file.name;
A = 179;
E = 90;
manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,'poscutOff',5,'speedThresh',2,...
    'numrow',numrow,'numcol',numcol,'rowloc',4,'colloc',3,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'dotSize',2)


%% IZ44
umap_path = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ44_220830_sess7\manifold';
cd(umap_path)
umap_name = 'behavior_speed_1_smooth_1_bin_0.02';
file = dir('*.position_behavior_speed_1_smooth_1_bin_0.02.mat');
behav_file = file.name;
A = 70;
E = -9;
manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,'poscutOff',5,'speedThresh',2,...
    'numrow',numrow,'numcol',numcol,'rowloc',1,'colloc',5,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'dotSize',2)

umap_name = 'behavior_speed_1_smooth_5_bin_0.02';
file = dir('*.position_behavior_speed_1_smooth_5_bin_0.02.mat');
behav_file = file.name;
A = 20;
E = -90;
manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,'poscutOff',5,'speedThresh',2,...
    'numrow',numrow,'numcol',numcol,'rowloc',2,'colloc',5,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'dotSize',2)

umap_name = 'behavior_speed_1_smooth_1_bin_0.1';
file = dir('*.position_behavior_speed_1_smooth_1_bin_0.1.mat');
behav_file = file.name;
A = 24;
E = -90;
manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,'poscutOff',5,'speedThresh',2,...
    'numrow',numrow,'numcol',numcol,'rowloc',3,'colloc',5,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'dotSize',2)

umap_name = 'behavior_speed_1_smooth_5_bin_0.1';
file = dir('*.position_behavior_speed_1_smooth_5_bin_0.1.mat');
behav_file = file.name;
A = -17;
E = -1.8;
manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,'poscutOff',5,'speedThresh',2,...
    'numrow',numrow,'numcol',numcol,'rowloc',4,'colloc',5,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'dotSize',3)

%% IZ39
% umap_path = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ39_220702_sess14\manifold';
% cd(umap_path)
% umap_name = 'behavior_speed_1_smooth_1_bin_0.02';
% file = dir('*.position_behavior_speed_1_smooth_1_bin_0.02.mat');
% behav_file = file.name;
% A = 183;
% E = 0;
% manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,...
%     'numrow',numrow,'numcol',numcol,'rowloc',1,'colloc',7,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'dotSize',2)
% 
% umap_name = 'behavior_speed_1_smooth_5_bin_0.02';
% file = dir('*.position_behavior_speed_1_smooth_5_bin_0.02.mat');
% behav_file = file.name;
% A = -157;
% E = 6;
% manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,...
%     'numrow',numrow,'numcol',numcol,'rowloc',2,'colloc',7,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'dotSize',2)
% 
% umap_name = 'behavior_speed_1_smooth_1_bin_0.1';
% file = dir('*.position_behavior_speed_1_smooth_1_bin_0.1.mat');
% behav_file = file.name;
% A = -2;
% E = -90;
% manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,...
%     'numrow',numrow,'numcol',numcol,'rowloc',3,'colloc',7,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'dotSize',5)
% 
% umap_name = 'behavior_speed_1_smooth_5_bin_0.1';
% file = dir('*.position_behavior_speed_1_smooth_5_bin_0.1.mat');
% behav_file = file.name;
% A = -117;
% E = -72;
% manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,...
%     'numrow',numrow,'numcol',numcol,'rowloc',4,'colloc',7,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'dotSize',5)



%% For the same sessions, add a panel below showing only incorrect trials
umap_path = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ47_230707_sess24\manifold';
cd(umap_path)
umap_name = 'behavior_speed_1_smooth_5_bin_0.1';
file = dir('*.position_behavior_speed_1_smooth_5_bin_0.1.mat');
behav_file = file.name;
A = -74;
E = -90;
manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,'poscutOff',5,'speedThresh',2,...'dotSize',2
    'numrow',numrow,'numcol',numcol,'rowloc',5,'colloc',1,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'dotSize',4,'error',0)

umap_path = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ48_230628_sess17\manifold';
cd(umap_path)
umap_name = 'behavior_speed_1_smooth_5_bin_0.1';
file = dir('*.position_behavior_speed_1_smooth_5_bin_0.1.mat');
behav_file = file.name;
A = 179;
E = 90;
manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,'poscutOff',5,'speedThresh',2,...
    'numrow',numrow,'numcol',numcol,'rowloc',5,'colloc',3,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'dotSize',2,'error',0)

umap_path = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ44_220830_sess7\manifold';
cd(umap_path)
umap_name = 'behavior_speed_1_smooth_5_bin_0.1';
file = dir('*.position_behavior_speed_1_smooth_5_bin_0.1.mat');
behav_file = file.name;
A = -17;
E = -1.8;
manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,'poscutOff',5,'speedThresh',2,...
    'numrow',numrow,'numcol',numcol,'rowloc',5,'colloc',5,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'dotSize',3,'error',0)

%% Finally manifolds generated separately for correct and incorrect trials
umap_path = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ47_230707_sess24\manifold';
cd(umap_path)

umap_name = 'behavior_speed_0_smooth_5_bin_0.1_correct_only';
file = dir('*.position_behavior_speed_0_smooth_5_bin_0.1_correct_only.mat');
behav_file = file.name;
A = -96.3633;
E = 70;
manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,'poscutOff',5,'speedThresh',2,...'dotSize',2
    'numrow',numrow,'numcol',numcol,'rowloc',6,'colloc',1,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'dotSize',4,'error',1)

umap_name = 'behavior_speed_0_smooth_5_bin_0.1_error_only';
file = dir('*.position_behavior_speed_0_smooth_5_bin_0.1_error_only.mat');
behav_file = file.name;
A = -375.87;
E = -1.94;
manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,'poscutOff',5,'speedThresh',2,...'dotSize',2
    'numrow',numrow,'numcol',numcol,'rowloc',6,'colloc',3,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E,'dotSize',4,'error',0)


expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
saveas(gcf,strcat(expPath,'SupFigure11_UMAPbinning.png'));
saveas(gcf,strcat(expPath,'SupFigure11_UMAPbinning.eps'),'epsc');
saveas(gcf,strcat(expPath,'SupFigure11_UMAPbinning.fig'));
end
