function FigS12_UMAP_acrossMice

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[680 65 1877 962]);

sess = {'IZ39_220714_sess18','IZ40_220708_sess17','IZ43_220828_sess4','IZ44_220830_sess7','IZ47_230707_sess24','IZ48_230628_sess17'};

col1 = [83/255 0/255 0/255;...
    184/255 15/255 10/255;...
    241/255 114/255 42/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];

col2 = [0/255 0/255 255/255;...
    0/255 0/255 0/255;...
   224/255 163/255 46/255];

A1 = [-181.1478,-100,-182,363.1123,-74,179];
E1 = [90,-31.0830,-2,6.4,90,90];

A2 = [-183.3154, 131.88, 37.5821,-11.1170,-182.7542,51.5817];
E2 = [12.1462, 14.1623, -7.3622,7.5846,-53.9952,90];


numrow = 6;
numcol = 5;

for ss = 1:length(sess)
    umap_path = strcat('Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\',sess{ss},'\manifold');
    cd(umap_path)     
    umap_name = 'behavior_speed_1_smooth_5_bin_0.1';
    file = dir('*.position_behavior_speed_1_smooth_5_bin_0.1.mat');
    behav_file = file.name;
    TRIAL_TYPE = [0 1 2 3 4 5];
    manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,'addPosPlot',true,'addFreq',true,'poscutOff',5,...
        'numrow',numrow,'numcol',numcol,'rowloc',ss,'colloc',1,'col',col1,'TRIAL_TYPE', TRIAL_TYPE,'A',A1(ss),'E',E1(ss),'dotSize',2,'plotcolorbar',true)

    TRIAL_TYPE = [5 6 8];
    manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,'addPosPlot',true,'addFreq',false,...
        'numrow',numrow,'numcol',numcol,'rowloc',ss,'colloc',4,'col',col2,'TRIAL_TYPE', TRIAL_TYPE,'A',A2(ss),'E',E2(ss),'dotSize',2)
    
end


expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
saveas(gcf,strcat(expPath,'SupFigure12_UMAPacrossMice.png'));
saveas(gcf,strcat(expPath,'SupFigure12_UMAPacrossMice.eps'),'epsc');
saveas(gcf,strcat(expPath,'SupFigure12_UMAPacrossMice.fig'));
end
