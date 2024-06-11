function FigS14_UMAP_rotationData

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[680 65 1000 962]);

numrow = 5;
numcol = 3;

[clock,counterclock] = extractRotationChoice;
ax1 =  subplot(numrow,numcol,10);
imagesc(clock);
ylabel('Sess num')
xlabel('Ports')
colormap(ax1,'viridis')
title('Fraction clock-wise')
colorbar

ax1 =  subplot(numrow,numcol,13);
imagesc(counterclock);
ylabel('Sess num')
xlabel('Ports')
colormap(ax1,'viridis')
colorbar
title('Fraction counter clock-wise')

sess = {'IZ39_220702_sess14','IZ43_220913_sess11','IZ43_220828_sess4','IZ47_230707_sess24','IZ48_230705_sess22'};

TRIAL_TYPE = [0 1 2 3 4 5];
col = [238/255 67/255 69/255;...
    241/255 114/255 42/255;...
    247/255 149/255 33/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];

A = [-161,-1.5,-182,-65,-271];
E = [3,-3,-2,90,-90];

for ss = 1:length(sess)
    umap_path = strcat('Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\',sess{ss},'\manifold');
    cd(umap_path)    
    umap_name = 'behavior_speed_1_smooth_5_bin_0.1';
    file = dir('*.position_behavior_speed_1_smooth_5_bin_0.1.mat');
    behav_file = file.name;
    load(behav_file);
    load(strcat('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\rotationData\',sess{ss},'.rotation.mat'))
    Umap_results = readtable([umap_path, '\Umap_',umap_name,'.csv']);
    Umap_results = table2array(Umap_results);
    manifoldPlot('figHandle',fig2,'umap_name',umap_name,'behav_file',behav_file,'addPosPlot',true,'poscutOff',0,'speedThresh',-2,...
        'numrow',numrow,'numcol',numcol,'rowloc',ss,'colloc',2,'col',col,'TRIAL_TYPE', TRIAL_TYPE,'A',A(ss),'E',E(ss),'dotSize',5)     
    projectTurnsonManifold(Umap_results,timestamp_beh,trackRotation.timestamps,trackRotation.direction,fig2,numrow,numcol,ss,2)
end


expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
saveas(gcf,strcat(expPath,'SupFigure14_UMAP_rotation.png'));
saveas(gcf,strcat(expPath,'SupFigure14_UMAP_rotation.eps'),'epsc');
saveas(gcf,strcat(expPath,'SupFigure14_UMAP_rotation.fig'));
end

function projectTurnsonManifold(Umap_results,man_ts,points_ts,direction,fighandle,numrow,numcol,rowloc,colloc)

for ii = 1:length(points_ts)
    [~,idx] = min(abs(man_ts-points_ts(ii)));
    % Just subtract 2 frames to get the point before the turn
    idx = idx-1;
    subplot(numrow,numcol,numcol*(rowloc-1)+colloc,'Parent',fighandle)
    hold on
    if direction(ii) == 1
        scatter3(Umap_results(idx,1),Umap_results(idx,2),Umap_results(idx,3),20,'c','filled');
    elseif direction(ii) == -1
        scatter3(Umap_results(idx,1),Umap_results(idx,2),Umap_results(idx,3),20,'m','filled');
    end

    subplot(numrow,numcol,numcol*(rowloc-1)+colloc+1,'Parent',fighandle)
    hold on
    if direction(ii) == 1
        scatter3(Umap_results(idx,1),Umap_results(idx,2),Umap_results(idx,3),20,'c','filled');
    elseif direction(ii) == -1
        scatter3(Umap_results(idx,1),Umap_results(idx,2),Umap_results(idx,3),20,'m','filled');
    end

end
end

