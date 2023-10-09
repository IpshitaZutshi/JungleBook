function plotToneCells_Figure4

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[680 42 975 962]);
groupStyle = 1; %0 if group by lickloc, otherwise by toneGain

%rows = 22, columns = 6
numrows = 21;
numcol = 18;

%% Panel A1: Example tone tuned cell 1

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ39\Final\IZ39_220624_sess10';
cd(sessloc)
cellNum = 19;%13
file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file.name);
plotExampleCell(1,1, spikeData, tracking, behavTrials, firingMaps,fig2,numrows,numcol,cellNum)

file = dir(['*.rateMapsAvgError.cellinfo.mat']);
Maps = load(file.name);
errorMaps = Maps.firingMaps;
plotExampleCellErrorTrials(1,9, spikeData, tracking, behavTrials, firingMaps,errorMaps, fig2,numrows,numcol,cellNum,groupStyle)

%% Panel A2: Example tone tuned cell 2

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230707_sess24';
%'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220828_sess4';
cd(sessloc)
cellNum = 275;%4;%28
file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file.name);
plotExampleCell(1,3, spikeData, tracking, behavTrials, firingMaps,fig2,numrows,numcol,cellNum)

file = dir(['*.rateMapsAvgError.cellinfo.mat']);
Maps = load(file.name);
errorMaps = Maps.firingMaps;
plotExampleCellErrorTrials(1,11, spikeData, tracking, behavTrials, firingMaps,errorMaps,fig2,numrows,numcol,cellNum,groupStyle)


%% Panel A3: Example tone tuned cell 3

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230710_sess25';
%'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220828_sess4';
%'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230707_sess24';
%
cd(sessloc)
cellNum = 95;%134;%46;
file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file.name);
plotExampleCell(1,5, spikeData, tracking, behavTrials, firingMaps,fig2,numrows,numcol,cellNum)

file = dir(['*.rateMapsAvgError.cellinfo.mat']);
Maps = load(file.name);
errorMaps = Maps.firingMaps;
plotExampleCellErrorTrials(1,13, spikeData, tracking, behavTrials, firingMaps,errorMaps,fig2,numrows,numcol,cellNum,groupStyle)


%% Panel D: Example tone tuned cell 4
% 
 sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220901_sess8';
%'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220828_sess4';
cd(sessloc)
cellNum = 46;%134;
file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file.name);
plotExampleCell(1,7, spikeData, tracking, behavTrials, firingMaps,fig2,numrows,numcol,cellNum)

file = dir(['*.rateMapsAvgError.cellinfo.mat']);
Maps = load(file.name);
errorMaps = Maps.firingMaps;
plotExampleCellErrorTrials(1,15, spikeData, tracking, behavTrials, firingMaps,errorMaps,fig2,numrows,numcol,cellNum,groupStyle)


%% Calculate place cells
Summary = compileAvgTonePlaceMaps('plotfig',false,'savefig',false);

%% Panel B: Heat maps

linPos = linspace(1,122,50);
linTone = linspace(2000,22000,50);
    
ss = 2;
idxSess = Summary.AllsessType==(ss-1) & Summary.AllcellType == 1;
idxActive = Summary.AllsessType==(ss-1) & Summary.AllcellType == 1 & nanmean(Summary.AllspaceMapAvg,2)>0.2;
idxMaps = idxSess & Summary.AlltoneField & Summary.AlltoneCorr>0.1;  

selectedlinMap = Summary.AlllinMapInit(idxMaps,:);
selectedlinMapEnd = Summary.AlllinMapEnd(idxMaps,:);
selectedspaceMap = Summary.AllspaceMapAvg(idxMaps,:);
selectedtoneMap = Summary.AlltoneMap(idxMaps,:);

[maxTone,idxTone] = max(selectedtoneMap,[],2);

normlinMap = (selectedlinMap-nanmean(selectedlinMap,2))./nanstd(selectedlinMap,[],2);
normlinMapEnd = (selectedlinMapEnd-nanmean(selectedlinMapEnd,2))./nanstd(selectedlinMapEnd,[],2);
normtoneMap = (selectedtoneMap-nanmean(selectedtoneMap,2))./nanstd(selectedtoneMap,[],2);
normspaceMap = (selectedspaceMap-nanmean(selectedspaceMap,2))./nanstd(selectedspaceMap,[],2);

[~,sortidx] = sort(idxTone,'ascend');

BuPu=cbrewer('seq', 'BuPu', 11);
YlGnBu=cbrewer('seq', 'YlGnBu', 11);

ax1 = subplot(numrows,numcol,[(numcol*16)+1 (numcol*16)+2 (numcol*17)+1 (numcol*17)+2 (numcol*18)+1 (numcol*18)+2 (numcol*19)+1 (numcol*19)+2]);
imagesc(linPos, 1:length(sortidx),normlinMap(sortidx,:))
colormap(ax1, YlGnBu)
ylabel(strcat('Cell ID ',num2str(length(sortidx))))
caxis([-1 3])
title('No tone I')

ax1 = subplot(numrows,numcol,[(numcol*16)+3 (numcol*16)+4 (numcol*17)+3 (numcol*17)+4 (numcol*18)+3 (numcol*18)+4 (numcol*19)+3 (numcol*19)+4]);
imagesc(linPos, 1:length(sortidx),normspaceMap(sortidx,:))
colormap(ax1, YlGnBu)
caxis([-1 3])
axis off
title('Tone - space')


ax1 = subplot(numrows,numcol,[(numcol*16)+5 (numcol*16)+6 (numcol*17)+5 (numcol*17)+6 (numcol*18)+5 (numcol*18)+6 (numcol*19)+5 (numcol*19)+6]);
imagesc(linTone, 1:length(sortidx),normtoneMap(sortidx,:))
colormap(ax1, BuPu)
colorbar
caxis([-1 3])
axis off
title('Tone - freq')

ax1 = subplot(numrows,numcol,[(numcol*16)+7 (numcol*16)+8 (numcol*17)+7 (numcol*17)+8 (numcol*18)+7 (numcol*18)+8 (numcol*19)+7 (numcol*19)+8]);
imagesc(linPos, 1:length(sortidx),normlinMapEnd(sortidx,:))
colormap(ax1, YlGnBu)
caxis([-1 3])
axis off
title('No tone II')

% %% Panel F: Correlation between space and tone tuning
% 
% ss = 2;
% idxSess = Summary.AllsessType==(ss-1) & Summary.AllcellType == 1;
% 
% idxTone = idxSess & Summary.AlltoneField & Summary.AlltoneCorr>0.2;  
% idxSpace = idxSess & Summary.AllspaceField;  
% 
% diffSpace = Summary.AllspatialCorr(idxSpace)-Summary.AlltoneCorr(idxSpace);
% diffTone = Summary.AlltoneCorr(idxTone)-Summary.AllspatialCorr(idxTone);
% 
% subplot(numrows, numcol, [9 10 11 numcol+9 numcol+10 numcol+11])
% histogram(diffSpace,[-1:0.1:1],'Normalization','probability','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeColor','none')
% line([median(diffSpace) median(diffSpace)],[0 0.15],'Color','k','LineWidth',2)
% hold on
% histogram(diffTone,[-1:0.1:1],'Normalization','probability','FaceColor',[96/243 60/243 108/243],'FaceAlpha',0.5,'EdgeColor','none')
% line([median(diffTone) median(diffTone)],[0 0.15],'Color',[96/243 60/243 108/243],'LineWidth',2)

%% Panel D: Heatmaps of error and no error tuning

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

ax1 = subplot(numrows,numcol,[(numcol*16)+9 (numcol*16)+10 (numcol*17)+9 (numcol*17)+10 (numcol*18)+9 (numcol*18)+10 (numcol*19)+9 (numcol*19)+10]);
imagesc(linPos, 1:length(sortidx),normtoneMap(sortidx,:))
colormap(ax1, BuPu)
caxis([-1 3])
title('Tone')

ax1 = subplot(numrows,numcol,[(numcol*16)+11 (numcol*16)+12 (numcol*17)+11 (numcol*17)+12 (numcol*18)+11 (numcol*18)+12 (numcol*19)+11 (numcol*19)+12]);
imagesc(linPos, 1:length(sortidx),normtoneMapError1(sortidx,:))
colormap(ax1, BuPu)
caxis([-1 3])
axis off
title('Tone E1')

ax1 = subplot(numrows,numcol,[(numcol*16)+13 (numcol*16)+14 (numcol*17)+13 (numcol*17)+14 (numcol*18)+13 (numcol*18)+14 (numcol*19)+13 (numcol*19)+14]);
imagesc(linTone, 1:length(sortidx),normtoneMapError2(sortidx,:))
colormap(ax1, BuPu)
caxis([-1 3])
axis off
title('Tone E2')

%% Panel E: Quantification of error and no error tuning

col = [52/243 52/243 52/243;...
    156/243 156/243 156/243];
    
data{1} = Summary.AlltoneCorrError1(idxMaps,:);
data{2} = Summary.AlltoneCorrError2(idxMaps,:);

subplot(numrows,numcol,[15+(numcol*16) 16+(numcol*16) 15+(numcol*17) 16+(numcol*17)  15+(numcol*18) 16+(numcol*18)])
Fig4Stats.ErrorCorr= groupStats(data,[],...
    'inAxis',true,'plotType','boxplot','color',col,'labelSummary',false);

%% Save figure

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';
saveas(gcf,strcat(expPath,'\Compiled\Figures\Figure4_taskTuning',num2str(groupStyle),'.png'));
saveas(gcf,strcat(expPath,'\Compiled\Figures\Figure4_taskTuning',num2str(groupStyle),'.eps'),'epsc');
saveas(gcf,strcat(expPath,'\Compiled\Figures\Figure4_taskTuning',num2str(groupStyle),'.fig'));
save(strcat(expPath,'\Compiled\Figures\Figure4_taskTuning',num2str(groupStyle),'.mat'),'Fig4Stats'); 

end