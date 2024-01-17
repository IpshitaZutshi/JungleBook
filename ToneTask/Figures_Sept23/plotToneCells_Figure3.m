function plotToneCells_Figure3

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[680 42 975 962]);

numrows = 33;
numcol = 8;

sessID = [29 47 40 31];
cellID = [28 275 167 46];    

%% Panel A1: Example tone tuned cell 1
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220828_sess4';
%sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ39\Final\IZ39_220624_sess10';
cd(sessloc)
cellNum = 28;%19;
file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file.name);
plotExampleCell(1,1, spikeData, tracking, behavTrials, firingMaps,fig2,numrows,numcol,cellNum)


%% Panel A2: Example tone tuned cell 2

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230707_sess24';
%'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220828_sess4';
cd(sessloc)
cellNum = 275;%4;
file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file.name);
plotExampleCell(1,3, spikeData, tracking, behavTrials, firingMaps,fig2,numrows,numcol,cellNum)

%% Panel A3: Example tone tuned cell 3
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ44\Final\IZ44_220830_sess7';
%sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230710_sess25';
%'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220828_sess4';
%'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230707_sess24';

cd(sessloc)
cellNum = 167;%95;%134;%46;
file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file.name);
plotExampleCell(1,5, spikeData, tracking, behavTrials, firingMaps,fig2,numrows,numcol,cellNum)

%% Panel A4: Example tone tuned cell 4

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


%% Calculate place cells
Summary = compileAvgTonePlaceMaps('plotfig',false,'savefig',false);

%% Panel E: Heat maps
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
selectedsessID = Summary.AllsessID(idxMaps);
selectedcellID = Summary.AllcellID(idxMaps);

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
ax1 = subplot(numrows,numcol,[(numcol*21)+1 (numcol*21)+2 (numcol*22)+1 (numcol*22)+2 (numcol*23)+1 (numcol*23)+2]);
plot(linPos, nanmean(normlinMap,1));
ylim([-0.5 1])
box off

ax1 = subplot(numrows,numcol,[(numcol*16)+3 (numcol*16)+4 (numcol*17)+3 (numcol*17)+4 (numcol*18)+3 (numcol*18)+4 (numcol*19)+3 (numcol*19)+4]);
imagesc(linPos, 1:length(sortidx),normspaceMap(sortidx,:))
colormap(ax1, YlGnBu)
caxis([-1 3])
axis off
title('Tone - space')
ax1 = subplot(numrows,numcol,[(numcol*21)+3 (numcol*21)+4 (numcol*22)+3 (numcol*22)+4 (numcol*23)+3 (numcol*23)+4]);
plot(linPos, nanmean(normspaceMap,1));
ylim([-0.5 1])
box off

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
ax1 = subplot(numrows,numcol,[(numcol*21)+7 (numcol*21)+8 (numcol*22)+7 (numcol*22)+8 (numcol*23)+7 (numcol*23)+8]);
plot(linPos, nanmean(normlinMapEnd,1));
ylim([-0.5 1])
box off

% Identify example cells on the heatmaps
sortsessID = selectedsessID(sortidx);
sortcellID = selectedcellID(sortidx);
col = [119/255 221/255 229/255;...
    122/255 201/255 229/255;...
    38/255 169/255 224/255;...
    73/255 136/255 189/255];

for ll = 1:length(sessID)
    cellPos(ll) = find(sortsessID == sessID(ll) & sortcellID == cellID(ll));
    line([0 linPos(end)],[cellPos(ll) cellPos(ll)],'Color',col(ll,:),'LineWidth',2)
    hold on
end


%% Fraction of tone tuned cells
subplot(numrows,numcol,[(numcol*29)+7 (numcol*29)+8 (numcol*30)+7 (numcol*30)+8 (numcol*31)+7 (numcol*31)+8 (numcol*32)+7 (numcol*32)+8]);
ss=2;
sessIndices = unique(Summary.AllsessID(Summary.AllsessType==(ss-1)));
for ll = 1:length(sessIndices)
    sumActive(ll) = sum(Summary.AllsessType==(ss-1) & Summary.AllcellType == 1 & nanmean(Summary.AllspaceMapAvg,2)>0.2 & Summary.AllsessID==sessIndices(ll)); 
    sumMaps(ll) = sum(Summary.AllsessType==(ss-1) & Summary.AlltoneField & Summary.AlltoneCorr>0.1 & Summary.AllsessID==sessIndices(ll)); 
end
sessFractTuned = [sumMaps./sumActive];
sessActive = sumActive;
barh(sessFractTuned)
hold on
line([median(sumMaps./sumActive) median(sumMaps./sumActive)],[0 38],'Color','k','LineWidth',1.5)
box off

%% Panel F: Spatial tuning of tone cells during no tone trials

Summary = compareActionNoToneTrials('plotfig',false,'savefig',false);

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

ax1 = subplot(numrows,numcol,[(numcol*24)+1 (numcol*24)+2 (numcol*25)+1 (numcol*25)+2 (numcol*26)+1 (numcol*26)+2 (numcol*27)+1 (numcol*27)+2]);
imagesc(linPos, 1:length(sortidx),normlinMap(sortidx,:))
colormap(ax1, YlGnBu)
ylabel(strcat('Cell ID ',num2str(length(sortidx))))
caxis([-1 3])
title('No tone I')

ax1 = subplot(numrows,numcol,[(numcol*24)+3 (numcol*24)+4 (numcol*25)+3 (numcol*25)+4 (numcol*26)+3 (numcol*26)+4 (numcol*27)+3 (numcol*27)+4]);
imagesc(linPos, 1:length(sortidx),normlinMap1(sortidx,:))
colormap(ax1, YlGnBu)
caxis([-1 3])
title('No tone I first half')

ax1 = subplot(numrows,numcol,[(numcol*24)+5 (numcol*24)+6 (numcol*25)+5 (numcol*25)+6 (numcol*26)+5 (numcol*26)+6 (numcol*27)+5 (numcol*27)+6]);
imagesc(linPos, 1:length(sortidx),normlinMap2(sortidx,:))
colormap(ax1, YlGnBu)
caxis([-1 3])
title('No tone I second half')

ax1 = subplot(numrows,numcol,[(numcol*24)+7 (numcol*24)+8 (numcol*25)+7 (numcol*25)+8 (numcol*26)+7 (numcol*26)+8 (numcol*27)+7 (numcol*27)+8]);
imagesc(linPos, 1:length(sortidx),normlinMapEnd(sortidx,:))
colormap(ax1, YlGnBu)
caxis([-1 3])
title('No tone II')

ax1 = subplot(numrows,numcol,[(numcol*29)+1 (numcol*29)+2 (numcol*30)+1 (numcol*30)+2 (numcol*31)+1 (numcol*31)+2 (numcol*32)+1 (numcol*32)+2]);
imagesc(linPos, 1:length(sortidx),normspaceMap(sortidx,:))
colormap(ax1, YlGnBu)
ylabel(strcat('Cell ID ',num2str(length(sortidx))))
caxis([-1 3])
title('Space')

ax1 = subplot(numrows,numcol,[(numcol*29)+3 (numcol*29)+4 (numcol*30)+3 (numcol*30)+4 (numcol*31)+3 (numcol*31)+4 (numcol*32)+3 (numcol*32)+4]);
imagesc(linPos, 1:length(sortidx),normtoneMap(sortidx,:))
colormap(ax1, BuPu)
caxis([-1 3])
title('Tone')

subplot(numrows,numcol,[(numcol*29)+5 (numcol*29)+6 (numcol*30)+5 (numcol*30)+6 (numcol*31)+5 (numcol*31)+6 (numcol*32)+5 (numcol*32)+6]);
bar([sum(Summary.AlllinField)./length(Summary.AlllinField) sum(Summary.AlllinFieldEnd)./length(Summary.AlllinFieldEnd)])
box off
%ylim([0 0.5])


%% Save figure
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\Figure3_toneTuning.png'));
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\Figure3_toneTuning.eps'),'epsc');
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\Figure3_toneTuning.fig'));

end