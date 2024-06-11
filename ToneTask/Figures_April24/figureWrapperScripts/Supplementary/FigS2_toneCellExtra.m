function FigS2_toneCellExtra

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[22 465 1892 488]);

numrows = 1;
numcol = 10;

%% Calculate place cells
%Summary = compileAvgTonePlaceMaps('plotfig',false,'savefig',false);
load('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\compilePlaceFields.mat')

%% Fraction of tone tuned cells
ss=2;
sessIndices = unique(Summary.AllsessID(Summary.AllsessType==(ss-1)));
for ll = 1:length(sessIndices)
    sumActive(ll) = sum(Summary.AllsessType==(ss-1) & Summary.AllcellType == 1 & nanmean(Summary.AllspaceMapAvg,2)>0.2 & Summary.AllsessID==sessIndices(ll)); 
    sumMaps(ll) = sum(Summary.AllsessType==(ss-1) & Summary.AllcellType == 1 & Summary.AlltoneField & Summary.AlltoneCorr>0.1 & Summary.AllsessID==sessIndices(ll)); 
end
sessFractTuned = [sumMaps./sumActive];
sessActive = sumActive;
subplot(numrows,numcol,1)
bar(sessFractTuned)
line([0 37],[median(sessFractTuned) median(sessFractTuned)])
title('Non spatial cells')

subplot(numrows,numcol,2)
bar(sessActive)
line([0 37],[median(sessActive) median(sessActive)])
title('Active cells')


%% Panel B: Spatial tuning of tone cells during no tone trials

Summary = compareActionNoToneTrials('plotfig',false,'savefig',false);

linPos = linspace(1,122,50);
c = linspace(0,1.05,50);

freqExp = log10(22000/1000);
for ii = 1:length(c)
    linTone(ii) = (1000*(10.^(freqExp*c(ii))));
end

normlinMap = (Summary.AllrateMapLin-nanmean(Summary.AllrateMapLin,2))./nanstd(Summary.AllrateMapLin,[],2);
normlinMap1 = (Summary.AllrateMapLin1-nanmean(Summary.AllrateMapLin1,2))./nanstd(Summary.AllrateMapLin1,[],2);
normlinMap2 = (Summary.AllrateMapLin2-nanmean(Summary.AllrateMapLin2,2))./nanstd(Summary.AllrateMapLin2,[],2);
normlinMapEnd = (Summary.AllrateMapLinEnd-nanmean(Summary.AllrateMapLinEnd,2))./nanstd(Summary.AllrateMapLinEnd,[],2);
normspaceMap = (Summary.AllrateMapSpace-nanmean(Summary.AllrateMapSpace,2))./nanstd(Summary.AllrateMapSpace,[],2);
normtoneMap = (Summary.AllrateMaptone-nanmean(Summary.AllrateMaptone,2))./nanstd(Summary.AllrateMaptone,[],2);

[maxLin,idxLin] = max(Summary.AllrateMapLin,[],2);
[~,sortidx] = sort(idxLin,'ascend');
    
BuPu=cbrewer('seq', 'BuPu', 11);
YlGnBu=cbrewer('seq', 'YlGnBu', 11);

ax1 = subplot(numrows,numcol,3);
imagesc(linPos, 1:length(sortidx),normlinMap(sortidx,:))
colormap(ax1, YlGnBu)
ylabel(strcat('Cell ID ',num2str(length(sortidx))))
caxis([-1 3])
title('No tone I')

ax1 = subplot(numrows,numcol,4);
imagesc(linPos, 1:length(sortidx),normlinMap1(sortidx,:))
colormap(ax1, YlGnBu)
caxis([-1 3])
title('No tone I first half')

ax1 = subplot(numrows,numcol,5);
imagesc(linPos, 1:length(sortidx),normlinMap2(sortidx,:))
colormap(ax1, YlGnBu)
caxis([-1 3])
title('No tone I second half')

ax1 = subplot(numrows,numcol,6);
imagesc(linPos, 1:length(sortidx),normlinMapEnd(sortidx,:))
colormap(ax1, YlGnBu)
caxis([-1 3])
title('No tone II')

ax1 = subplot(numrows,numcol,7);
imagesc(linPos, 1:length(sortidx),normspaceMap(sortidx,:))
colormap(ax1, YlGnBu)
ylabel(strcat('Cell ID ',num2str(length(sortidx))))
caxis([-1 3])
title('Space')

ax1 = subplot(numrows,numcol,8);
imagesc(linPos, 1:length(sortidx),normtoneMap(sortidx,:))
colormap(ax1, BuPu)
caxis([-1 3])
title('Tone')

ax1 = subplot(numrows,numcol,9);
bar([sum(Summary.AlllinField)./length(Summary.AlllinField) sum(Summary.AlllinFieldEnd)./length(Summary.AlllinFieldEnd)])
box off
%ylim([0 0.5])

subplot(numrows,numcol,10)
Summary = plotToneBehaviorStatistics;
[r,p] = corr([sessFractTuned' Summary.performance'],'Type','Pearson','rows','complete');
scatter(sessFractTuned,Summary.performance,5,'filled')
hold on
lsline
box off
xlabel('Fraction non-spatial')
ylabel('Behavior performance')
title(strcat(num2str(r(1,2)),'|',num2str(p(1,2))))


%% Save figure
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
saveas(gcf,strcat(expPath,'SupFigure2_toneCellLinTrial.png'));
saveas(gcf,strcat(expPath,'SupFigure2_toneCellLinTrial.eps'),'epsc');
saveas(gcf,strcat(expPath,'SupFigure2_toneCellLinTrial.fig'));


%% Next, make a separate figure to plot what these cells do in error trials

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[680 42 975 962]);

numrows = 21;
numcol = 18;

%% Panel A1: Example tone tuned cell 1

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220828_sess4';
cd(sessloc)
cellNum = 28;
file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*.rateMapsAvgnotLog.cellinfo.mat']);
load(file.name);
file = dir(['*.rateMapsAvgError.cellinfo.mat']);
Maps = load(file.name);
errorMaps = Maps.firingMaps;
plotExampleCellErrorTrials(1,1, spikeData, tracking, behavTrials, firingMaps,errorMaps,fig2,numrows,numcol,cellNum,0)
plotExampleCellErrorTrials(1,9, spikeData, tracking, behavTrials, firingMaps,errorMaps,fig2,numrows,numcol,cellNum,1)

%% Panel A2: Example tone tuned cell 2

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230707_sess24';
cd(sessloc)
cellNum = 275;
file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*.rateMapsAvgnotLog.cellinfo.mat']);
load(file.name);
file = dir(['*.rateMapsAvgError.cellinfo.mat']);
Maps = load(file.name);
errorMaps = Maps.firingMaps;
plotExampleCellErrorTrials(1,3, spikeData, tracking, behavTrials, firingMaps,errorMaps,fig2,numrows,numcol,cellNum,0)
plotExampleCellErrorTrials(1,11, spikeData, tracking, behavTrials, firingMaps,errorMaps,fig2,numrows,numcol,cellNum,1)


%% Panel A3: Example tone tuned cell 3

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ44\Final\IZ44_220830_sess7';
cd(sessloc)
cellNum = 167;
file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*.rateMapsAvgnotLog.cellinfo.mat']);
load(file.name);
file = dir(['*.rateMapsAvgError.cellinfo.mat']);
Maps = load(file.name);
errorMaps = Maps.firingMaps;
plotExampleCellErrorTrials(1,5, spikeData, tracking, behavTrials, firingMaps,errorMaps,fig2,numrows,numcol,cellNum,0)
plotExampleCellErrorTrials(1,13, spikeData, tracking, behavTrials, firingMaps,errorMaps,fig2,numrows,numcol,cellNum,1)


%% Panel A4: Example tone tuned cell 4

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220901_sess8';
cd(sessloc)
cellNum = 46;
file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*.rateMapsAvgnotLog.cellinfo.mat']);
load(file.name);
file = dir(['*.rateMapsAvgError.cellinfo.mat']);
Maps = load(file.name);
errorMaps = Maps.firingMaps;
plotExampleCellErrorTrials(1,7, spikeData, tracking, behavTrials, firingMaps,errorMaps,fig2,numrows,numcol,cellNum,0)
plotExampleCellErrorTrials(1,15, spikeData, tracking, behavTrials, firingMaps,errorMaps,fig2,numrows,numcol,cellNum,1)


%% Panel B: Heatmaps of error and no error tuning
load('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\compilePlaceFields.mat')
ss = 2;
idxSess = Summary.AllsessType==(ss-1) & Summary.AllcellType == 1;
idxMaps = idxSess & Summary.AlltoneField & Summary.AlltoneCorr>0.1;  

selectedtoneMap = Summary.AlltoneMap(idxMaps,:);
selectedtoneErrorMap1 = Summary.AlltoneMapError1(idxMaps,:);
selectedtoneErrorMap2 = Summary.AlltoneMapError2(idxMaps,:);

[maxTone,idxTone] = max(selectedtoneMap,[],2);

normtoneMap = (selectedtoneMap-nanmean(selectedtoneMap,2))./nanstd(selectedtoneMap,[],2);
normtoneMapError1 = (selectedtoneErrorMap1-nanmean(selectedtoneErrorMap1,2))./nanstd(selectedtoneErrorMap1,[],2);
normtoneMapError2 = (selectedtoneErrorMap2-nanmean(selectedtoneErrorMap2,2))./nanstd(selectedtoneErrorMap2,[],2);

[~,sortidx] = sort(idxTone,'ascend');

BuPu=cbrewer('seq', 'BuPu', 11);

ax1 = subplot(numrows,numcol,[(numcol*16)+1 (numcol*16)+2 (numcol*17)+1 (numcol*17)+2 (numcol*18)+1 (numcol*18)+2 (numcol*19)+1 (numcol*19)+2]);
imagesc(linPos, 1:length(sortidx),normtoneMap(sortidx,:))
colormap(ax1, BuPu)
caxis([-1 3])
title('Tone')

ax1 = subplot(numrows,numcol,[(numcol*16)+3 (numcol*16)+4 (numcol*17)+3 (numcol*17)+4 (numcol*18)+3 (numcol*18)+4 (numcol*19)+3 (numcol*19)+4]);
imagesc(linPos, 1:length(sortidx),normtoneMapError1(sortidx,:))
colormap(ax1, BuPu)
caxis([-1 3])
axis off
title('Tone E1')

ax1 = subplot(numrows,numcol,[(numcol*16)+5 (numcol*16)+6 (numcol*17)+5 (numcol*17)+6 (numcol*18)+5 (numcol*18)+6 (numcol*19)+5 (numcol*19)+6]);
imagesc(linTone, 1:length(sortidx),normtoneMapError2(sortidx,:))
colormap(ax1, BuPu)
caxis([-1 3])
xlim([1000 25000])
xscale log
title('Tone E2')

%% Panel C: Quantification of error and no error tuning

col = [52/243 52/243 52/243;...
    156/243 156/243 156/243];
    
data{1} = Summary.AlltoneCorrError1(idxMaps,:);
data{2} = Summary.AlltoneCorrError2(idxMaps,:);

subplot(numrows,numcol,[7+(numcol*16) 8+(numcol*16) 7+(numcol*17) 8+(numcol*17)  7+(numcol*18) 8+(numcol*18)])
Stats.ErrorCorr= groupStats(data,[],...
    'inAxis',true,'plotType','boxplot','color',col,'labelSummary',false);

%% Save figure

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
saveas(gcf,strcat(expPath,'SupFigure2_toneCellerrorTrials.png'));
saveas(gcf,strcat(expPath,'SupFigure2_toneCellerrorTrials.eps'),'epsc');
saveas(gcf,strcat(expPath,'SupFigure2_toneCellerrorTrials.fig'));
save(strcat(expPath,'SupFigure2_toneCellerrorTrials.mat'),'Stats'); 

end