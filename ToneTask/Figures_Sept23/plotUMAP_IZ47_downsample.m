fig2 = figure;
set(fig2, 'Renderer','painters')
set(fig2,'Position',[21 177 1366 579])

umap_path = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\finalSessions\Task\IZ47_230710_sess25\manifold';
%umap_path = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\finalSessions\Task\IZ48_230714_sess28\manifold';
cellnum = 108;
seed = 6;

sessname  = strsplit(umap_path,'\');
behav_file = strcat(umap_path,'\',sessname{8},'.position_behavior_downsample.mat');
umap_name = strcat('behavior_downsample_cell_',num2str(cellnum),'_seed_',num2str(seed));
A = 23;%-89;%-179;
E = 1.55;%7.68;%-22;

% Next, tone non-probe trials
TRIAL_TYPE = [0:5];
%col = 'jet';
col = [238/255 67/255 69/255;...
    241/255 114/255 42/255;...
    247/255 149/255 33/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];
manifoldPlot_probe('figHandle',fig2,'umap_path',umap_path,'umap_name',umap_name,'behav_file',behav_file,...
    'numrow',1,'numcol',2,'rowloc',1,'colloc',1,'probe',false,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A,'E',E)

