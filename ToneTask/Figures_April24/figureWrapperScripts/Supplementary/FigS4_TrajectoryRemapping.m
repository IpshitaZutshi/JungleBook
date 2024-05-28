function FigS4_TrajectoryRemapping

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[1000 253 852 673]);
YlGnBu=cbrewer('div', 'RdYlBu', 11);
YlGnBu = YlGnBu(end:-1:1,:);
%rows = 22, columns = 6
numrows = 8;
numcol = 10;


%% Panel F: 2-D maps of the same neurons from A-C
% Cell 1
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220919_sess14';
cd(sessloc)
file = dir(['*.rateMapsAvg2D.cellinfo.mat']);
load(file.name);
cellNum = 7;
maxRate = mean([max(firingMaps.forward.rateMaps{cellNum}{1},[],'all'),max(firingMaps.forward.rateMaps{cellNum}{4},[],'all'),max(firingMaps.forward.rateMaps{cellNum}{5},[],'all')]);

subplot(numrows,numcol,[1 2])
hold off
h = imagesc(firingMaps.forward.rateMaps{cellNum}{1});
colormap(YlGnBu)
set(h, 'AlphaData', ~isnan(firingMaps.forward.rateMaps{cellNum}{1}))
caxis([0 maxRate-5])
set(gca,'YDir','normal')
axis off
title(strcat(num2str(maxRate),' Hz'))

subplot(numrows,numcol,[1+1*numcol 2+1*numcol])
h = imagesc(firingMaps.forward.rateMaps{cellNum}{5});
colormap(YlGnBu)
set(h, 'AlphaData', ~isnan(firingMaps.forward.rateMaps{cellNum}{5}))
caxis([0 maxRate-5])
set(gca,'YDir','normal')
axis off
colorbar

subplot(numrows,numcol,[1+2*numcol 2+2*numcol])
h = imagesc(firingMaps.forward.rateMaps{cellNum}{4});
colormap(YlGnBu)
set(h, 'AlphaData', ~isnan(firingMaps.forward.rateMaps{cellNum}{4}))
caxis([0 maxRate-5])
set(gca,'YDir','normal')
axis off

corrLin = corrcoef(firingMaps.forward.rateMaps{cellNum}{1}(:)',firingMaps.forward.rateMaps{cellNum}{5}(:)','rows','pairwise');
title(num2str(corrLin(1,2)))

% Cell 2
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220919_sess14';
cd(sessloc)
file = dir(['*.rateMapsAvg2D.cellinfo.mat']);
load(file.name);
cellNum = 75;
maxRate = mean([max(firingMaps.forward.rateMaps{cellNum}{1},[],'all'),max(firingMaps.forward.rateMaps{cellNum}{4},[],'all'),max(firingMaps.forward.rateMaps{cellNum}{5},[],'all')]);

subplot(numrows,numcol,[3 4])
hold off
h = imagesc(firingMaps.forward.rateMaps{cellNum}{1});
colormap(YlGnBu)
caxis([0 maxRate-5])
set(h, 'AlphaData', ~isnan(firingMaps.forward.rateMaps{cellNum}{1}))
set(gca,'YDir','normal')
axis off
title(strcat(num2str(maxRate),' Hz'))

subplot(numrows,numcol,[3+1*numcol 4+1*numcol])
h = imagesc(firingMaps.forward.rateMaps{cellNum}{5});
colormap(YlGnBu)
set(h, 'AlphaData', ~isnan(firingMaps.forward.rateMaps{cellNum}{5}))
caxis([0 maxRate-5])
set(gca,'YDir','normal')
axis off

subplot(numrows,numcol,[3+2*numcol 4+2*numcol])
h = imagesc(firingMaps.forward.rateMaps{cellNum}{4});
colormap(YlGnBu)
set(h, 'AlphaData', ~isnan(firingMaps.forward.rateMaps{cellNum}{4}))
caxis([0 maxRate-5])
set(gca,'YDir','normal')
axis off
corrLin = corrcoef(firingMaps.forward.rateMaps{cellNum}{1}(:)',firingMaps.forward.rateMaps{cellNum}{5}(:)','rows','pairwise');
title(num2str(corrLin(1,2)))

% Cell 3
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ44\Final\IZ44_220830_sess7';
cd(sessloc)
file = dir(['*.rateMapsAvg2D.cellinfo.mat']);
load(file.name);
cellNum = 114;
maxRate = mean([max(firingMaps.forward.rateMaps{cellNum}{1},[],'all'),max(firingMaps.forward.rateMaps{cellNum}{4},[],'all'),max(firingMaps.forward.rateMaps{cellNum}{5},[],'all')]);
axis off

subplot(numrows,numcol,[5 6])
hold off
h = imagesc(firingMaps.forward.rateMaps{cellNum}{1});
colormap(YlGnBu)
set(h, 'AlphaData', ~isnan(firingMaps.forward.rateMaps{cellNum}{1}))
set(gca,'YDir','normal')
caxis([0 maxRate-5])
title(strcat(num2str(maxRate),' Hz'))
axis off

subplot(numrows,numcol,[5+1*numcol 6+1*numcol])
h = imagesc(firingMaps.forward.rateMaps{cellNum}{5});
colormap(YlGnBu)
set(h, 'AlphaData', ~isnan(firingMaps.forward.rateMaps{cellNum}{5}))
set(gca,'YDir','normal')
caxis([0 maxRate-5])
axis off

subplot(numrows,numcol,[5+2*numcol 6+2*numcol])
h = imagesc(firingMaps.forward.rateMaps{cellNum}{4});
colormap(YlGnBu)
set(h, 'AlphaData', ~isnan(firingMaps.forward.rateMaps{cellNum}{4}))
set(gca,'YDir','normal')
caxis([0 maxRate-5])
axis off
corrLin = corrcoef(firingMaps.forward.rateMaps{cellNum}{1}(:)',firingMaps.forward.rateMaps{cellNum}{5}(:)','rows','pairwise');
title(num2str(corrLin(1,2)))

%% Correlation of 2D maps
Summary = compileMice2DPlaceFields;
col = [52/243 52/243 52/243;...
    180/243 180/243 180/243;...
    56/243 61/243 150/243;...
    80/243 91/243 166/243;...
    ];
    
subplot(numrows,numcol,[7 8 7+1*numcol 8+1*numcol 7+2*numcol 8+2*numcol])
Stats.Maps2D = groupStats({Summary.linCorr,Summary.linlinEndCorr,Summary.toneNoToneCorr,Summary.tonelinEndCorr},[],...
    'inAxis',true,'plotType','boxplot','color',col,'labelSummary',false);

subplot(numrows,numcol,[9 10 9+1*numcol 10+1*numcol 9+2*numcol 10+2*numcol])
Stats.binNumber = groupStats({Summary.linCorrvalidPairs,Summary.linlinEndCorrvalidPairs,Summary.toneNoToneCorrvalidPairs,Summary.tonelinEndCorrvalidPairs},[],...
    'inAxis',true,'plotType','boxplot','color',col,'labelSummary',false);


expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
saveas(gcf,strcat(expPath,'SupFigure4A_2DTrajectory.png'));
saveas(gcf,strcat(expPath,'SupFigure4A_2DTrajectory.eps'),'epsc');
saveas(gcf,strcat(expPath,'SupFigure4A_2DTrajectory.fig'));
save(strcat(expPath,'SupFigure4A_2DTrajectory.mat'),'Stats'); 

%%Now analyze trajectory
fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[680 42 975 962]);

%rows = 22, columns = 6
numrows = 18;
numcol = 9;
%% Panel A: Example trajectory modulation cell ACgN % 13 rows, 2 columns each

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220919_sess14';
cd(sessloc)
cellNum = 7;
file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file.name);
plotExampleCell(1,1,fig2,numrows,numcol,cellNum)

subplot(numrows,numcol,[1 2])
corrLin = corrcoef(firingMaps.forward.rateMaps{cellNum}{1}',firingMaps.forward.rateMaps{cellNum}{7}','rows','complete');
title(num2str(corrLin(1,2)))
axis off

%% Panel B: Example trajectory modulation cell ACgN % 13 rows, 2 columns each

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220919_sess14';
cd(sessloc)
cellNum = 75;
file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file.name);
plotExampleCell(1,3, fig2,numrows,numcol,cellNum)

subplot(numrows,numcol,[3 4])
corrLin = corrcoef(firingMaps.forward.rateMaps{cellNum}{1}',firingMaps.forward.rateMaps{cellNum}{7}','rows','complete');
title(num2str(corrLin(1,2)))
axis off

%% Panel C: Example trajectory 1
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ44\Final\IZ44_220830_sess7';
cd(sessloc)
cellNum = 114;
file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file.name);
plotExampleCell(1,5,fig2,numrows,numcol,cellNum)

subplot(numrows,numcol,[5 6])
corrLin = corrcoef(firingMaps.forward.rateMaps{cellNum}{1}',firingMaps.forward.rateMaps{cellNum}{7}','rows','complete');
title(num2str(corrLin(1,2)))
axis off

%% Panel D: Trajectory variance
Summary = calculateTrajectoryVariance(2);
col = [52/243 52/243 52/243;...
    56/243 61/243 150/243;...
    122/243 122/243 122/243];

subplot(numrows,numcol,[7 8 7+numcol 8+numcol])
avgVar = nanmean(Summary.stdLin,1);
stderrVar = nanstd(Summary.stdLin,[],1)/sqrt(size(Summary.stdLin,1));
fill([(1:1:length(avgVar))'; (length(avgVar):-1:1)'],[avgVar'-stderrVar';flipud(avgVar'+stderrVar')],col(1,:),'linestyle','none','FaceAlpha',0.2);
hold on
plot(avgVar,'color',col(1,:),'LineWidth',1.5)  

avgVar = nanmean(Summary.stdTone6,1);
stderrVar = nanstd(Summary.stdTone6,[],1)/sqrt(size(Summary.stdTone6,1));
fill([(1:1:length(avgVar))'; (length(avgVar):-1:1)'],[avgVar'-stderrVar';flipud(avgVar'+stderrVar')],col(2,:),'linestyle','none','FaceAlpha',0.2);
hold on
plot(avgVar,'color',col(2,:),'LineWidth',1.5)     

avgVar = nanmean(Summary.stdLinEnd,1);
stderrVar = nanstd(Summary.stdLinEnd,[],1)/sqrt(size(Summary.stdLinEnd,1));
fill([(1:1:length(avgVar))'; (length(avgVar):-1:1)'],[avgVar'-stderrVar';flipud(avgVar'+stderrVar')],col(3,:),'linestyle','none','FaceAlpha',0.2);
hold on
plot(avgVar,'color',col(3,:),'LineWidth',1.5)     
box off
set(gca,'XTickLabels',[])
ylim([0 1.2])
title('ACgN')

subplot(numrows,numcol,[7+2*numcol 8+2*numcol 7+3*numcol 8+3*numcol])
avgVar = nanmean(Summary.posLin,1);
stderrVar = nanstd(Summary.posLin,[],1)/sqrt(size(Summary.posLin,1));
fill([(1:1:length(avgVar))'; (length(avgVar):-1:1)'],[avgVar'-stderrVar';flipud(avgVar'+stderrVar')],col(1,:),'linestyle','none','FaceAlpha',0.2);
hold on
plot(avgVar,'color',col(1,:),'LineWidth',1.5)  

avgVar = nanmean(Summary.posTone6,1);
stderrVar = nanstd(Summary.posTone6,[],1)/sqrt(size(Summary.posTone6,1));
fill([(1:1:length(avgVar))'; (length(avgVar):-1:1)'],[avgVar'-stderrVar';flipud(avgVar'+stderrVar')],col(2,:),'linestyle','none','FaceAlpha',0.2);
hold on
plot(avgVar,'color',col(2,:),'LineWidth',1.5)     

avgVar = nanmean(Summary.posLinEnd,1);
stderrVar = nanstd(Summary.posLinEnd,[],1)/sqrt(size(Summary.posLinEnd,1));
fill([(1:1:length(avgVar))'; (length(avgVar):-1:1)'],[avgVar'-stderrVar';flipud(avgVar'+stderrVar')],col(3,:),'linestyle','none','FaceAlpha',0.2);
hold on
plot(avgVar,'color',col(3,:),'LineWidth',1.5)     
box off
set(gca,'XTickLabels',[])
%ylim([0 1.2])
ylabel('Avg x - position (cm)')

%% Panel F: Clustering trajectories
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ44\Final\IZ44_220830_sess7';
%'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ44\Final\IZ44_220828_sess5';
%'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ39\Final\IZ39_220707_sess17';
cd(sessloc)
cellNum = 114;%32;%38;%44;
file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);
file = dir(['*.rateMapsTrial.cellinfo.mat']);
load(file(1).name);   

[posXtrial,posXLabel,T,positions,lin_begin,lin_end,tone6,posX,posY] = clusterTrajectories('plotfig',false);

subplot(numrows,numcol,[7+4*numcol 8+4*numcol 7+5*numcol 8+5*numcol])
idxLin = posXtrial(posXLabel==1);
idxLinT = T(posXLabel==1);
for ii = 1:length(idxLin)  
    if idxLinT(ii)==1
        plot(positions.forward{idxLin(ii)}(:,2),positions.forward{idxLin(ii)}(:,1),'Color',[0.6 0.6 0.6])
    elseif idxLinT(ii)==2
        plot(positions.forward{idxLin(ii)}(:,2),positions.forward{idxLin(ii)}(:,1),'Color','b')
    end
    
    hold on
    posID = ismember(positions.forward{idxLin(ii)}(:,4),spikeData.posIdx{cellNum});
    plot(positions.forward{idxLin(ii)}(posID,2),positions.forward{idxLin(ii)}(posID,1),'r.')
end

xlim([0 122])
ylim([0 6])
%box off

idxT1 = idxLin(idxLinT==1);
Maps{1} = [];
for ii = 1:length(idxT1)
    Maps{1} = [Maps{1};firingMaps.forward.rateMaps{cellNum}{idxT1(ii)}];
end

idxT2 = idxLin(idxLinT==2);
Maps{2} = [];
for ii = 1:length(idxT2)
    Maps{2} = [Maps{2};firingMaps.forward.rateMaps{cellNum}{idxT2(ii)}];
end

maxFR = 20;%max([nanmean(Maps{1}) nanmean(Maps{2})]);


YlGnBu=cbrewer('seq', 'YlGnBu', 11);

ax1 = subplot(numrows,numcol,[7+6*numcol 8+6*numcol]);
imagesc(nanmean(Maps{1},1))
colormap(ax1,YlGnBu);
caxis([0 maxFR])
axis off


ax1 = subplot(numrows,numcol,[7+7*numcol 8+7*numcol]);
imagesc(nanmean(Maps{2},1))
colormap(ax1,YlGnBu)
caxis([0 maxFR])
box off
set(gca,'XTickLabels',[],'YTickLabels',[])
xlabel(num2str(maxFR))


subplot(numrows,numcol,[7+8*numcol 8+8*numcol 7+9*numcol 8+9*numcol])
idxLin = posXtrial(posXLabel==2);
idxLinT = T(posXLabel==2);
for ii = 1:length(idxLin)  
    if idxLinT(ii)==1
        plot(positions.forward{idxLin(ii)}(:,2),positions.forward{idxLin(ii)}(:,1),'Color',[0.6 0.6 0.6])
    elseif idxLinT(ii)==2
        plot(positions.forward{idxLin(ii)}(:,2),positions.forward{idxLin(ii)}(:,1),'Color','b')
    end
    hold on
    posID = ismember(positions.forward{idxLin(ii)}(:,4),spikeData.posIdx{cellNum});
    plot(positions.forward{idxLin(ii)}(posID,2),positions.forward{idxLin(ii)}(posID,1),'r.')
end
xlim([0 122])
ylim([0 6])
box off

idxT1 = idxLin(idxLinT==1);
Maps{1} = [];
for ii = 1:length(idxT1)
    Maps{1} = [Maps{1};firingMaps.forward.rateMaps{cellNum}{idxT1(ii)}];
end

idxT2 = idxLin(idxLinT==2);
Maps{2} = [];
for ii = 1:length(idxT2)
    Maps{2} = [Maps{2};firingMaps.forward.rateMaps{cellNum}{idxT2(ii)}];
end

maxFR = 20;%max([nanmean(Maps{1}) nanmean(Maps{2})]);

ax1 = subplot(numrows,numcol,[7+10*numcol 8+10*numcol]);
imagesc(nanmean(Maps{1},1))
colormap(ax1,YlGnBu)
caxis([0 maxFR])
axis off

ax1 = subplot(numrows,numcol,[7+11*numcol 8+11*numcol]);
imagesc(nanmean(Maps{2},1))
colormap(ax1,YlGnBu)
caxis([0 maxFR])
box off
set(gca,'XTickLabels',[],'YTickLabels',[])
xlabel(num2str(maxFR))

% No tone 2
subplot(numrows,numcol,[7+12*numcol 8+12*numcol 7+13*numcol 8+13*numcol])
idxLin = posXtrial(posXLabel==3);
idxLinT = T(posXLabel==3);
for ii = 1:length(idxLin)  
    if idxLinT(ii)==1
        plot(positions.forward{idxLin(ii)}(:,2),positions.forward{idxLin(ii)}(:,1),'Color',[0.6 0.6 0.6])
    elseif idxLinT(ii)==2
        plot(positions.forward{idxLin(ii)}(:,2),positions.forward{idxLin(ii)}(:,1),'Color','b')
    end
    hold on
    posID = ismember(positions.forward{idxLin(ii)}(:,4),spikeData.posIdx{cellNum});
    plot(positions.forward{idxLin(ii)}(posID,2),positions.forward{idxLin(ii)}(posID,1),'r.')
end
xlim([0 122])
ylim([0 6])
box off

idxT1 = idxLin(idxLinT==1);
Maps{1} = [];
for ii = 1:length(idxT1)
    Maps{1} = [Maps{1};firingMaps.forward.rateMaps{cellNum}{idxT1(ii)}];
end

idxT2 = idxLin(idxLinT==2);
Maps{2} = [];
for ii = 1:length(idxT2)
    Maps{2} = [Maps{2};firingMaps.forward.rateMaps{cellNum}{idxT2(ii)}];
end

maxFR = 20;%max([nanmean(Maps{1}) nanmean(Maps{2})]);

ax1 = subplot(numrows,numcol,[7+14*numcol 8+14*numcol]);
imagesc(nanmean(Maps{1},1))
colormap(ax1,YlGnBu)
caxis([0 maxFR])
axis off


ax1 = subplot(numrows,numcol,[7+15*numcol 8+15*numcol]);
imagesc(nanmean(Maps{2},1))
colormap(ax1,YlGnBu)
caxis([0 maxFR])
box off
set(gca,'XTickLabels',[],'YTickLabels',[])
xlabel(num2str(maxFR))

%% Finally, summary of trajectory

Summary = sessionRemappingTrajectorynoToneVer2('plotfig',false,'matchTrials',false);

col = [52/243 52/243 52/243;...
    160/243 160/243 160/243;... 
    214/243 126/243 44/243;...
    224/243 163/243 46/243];

subplot(numrows,numcol,[4+15*numcol 5+15*numcol 6+15*numcol 4+16*numcol 5+16*numcol 6+16*numcol 4+17*numcol 5+17*numcol 6+17*numcol])
Stats = groupStats([{Summary.comp1},{Summary.comp2},{Summary.comp3},{Summary.comp4}],[],'color',col,'inAxis',true','plotType','boxplot','labelSummary',false);
xticks('auto')
xticklabels({'','SC,ST','SC,DT','DC,ST','DC,DT'})

Summary = sessionRemappingTrajectorynoToneVer2('plotfig',false,'matchTrials',true);

col = [52/243 52/243 52/243;...
    160/243 160/243 160/243;... 
    214/243 126/243 44/243;...
    224/243 163/243 46/243];

subplot(numrows,numcol,[1+15*numcol 2+15*numcol 3+15*numcol 1+16*numcol 2+16*numcol 3+16*numcol 1+17*numcol 2+17*numcol 3+17*numcol])
Stats = groupStats([{Summary.comp1},{Summary.comp2},{Summary.comp3},{Summary.comp4}],[],'color',col,'inAxis',true','plotType','boxplot','labelSummary',false);
xticks('auto')
xticklabels({'','SC,ST','SC,DT','DC,ST','DC,DT'})

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
saveas(gcf,strcat(expPath,'SupFigure4B_TrajectoryCluster.png'));
saveas(gcf,strcat(expPath,'SupFigure4B_TrajectoryCluster.eps'),'epsc');
saveas(gcf,strcat(expPath,'SupFigure4B_TrajectoryCluster.fig'));
save(strcat(expPath,'SupFigure4B_TrajectoryCluster.mat'),'Stats'); 

end