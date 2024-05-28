function Fig2_plotSingleCells

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[680 42 975 962]);

numrows = 28;
numcol = 14;

sessID = [29 47 40 31];
cellID = [28 275 167 46];    

%% Panel A: Examples of single cells
%Tone-like
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220828_sess4';
cd(sessloc)
cellNum = 28;
plotExampleCell(1,1,fig2,numrows,numcol,cellNum)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230707_sess24';
cd(sessloc)
cellNum = 275;
plotExampleCell(1,3,fig2,numrows,numcol,cellNum)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ44\Final\IZ44_220830_sess7';
cd(sessloc)
cellNum = 167;
plotExampleCell(1,5,fig2,numrows,numcol,cellNum)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220901_sess8';
cd(sessloc)
cellNum = 46;
plotExampleCell(1,7,fig2,numrows,numcol,cellNum)

%Place-like
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220830_sess6';
cd(sessloc)
cellNum = 41;
plotExampleCell(1,9,fig2,numrows,numcol,cellNum)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ44\Final\IZ44_220828_sess5';
cd(sessloc)
cellNum = 38;
plotExampleCell(1,11,fig2,numrows,numcol,cellNum)

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230707_sess24';
cd(sessloc)
cellNum = 47;
plotExampleCell(1,13,fig2,numrows,numcol,cellNum)


%% Calculate place/tone cells
%Summary = compileAvgTonePlaceMaps('plotfig',false,'savefig',false,'saveMat',true);
load('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\compilePlaceFields.mat')

%% Panel B: Heat maps
linPos = linspace(1,122,50);
c = linspace(0,1.05,50);

freqExp = log10(22000/1000);
for ii = 1:length(c)
    linTone(ii) = (1000*(10.^(freqExp*c(ii))));
end

    
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

%% Plot tone tuned cells

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
title('Tone - freq')
xscale log

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

subplot(numrows,numcol,[(numcol*16)+1 (numcol*16)+2 (numcol*17)+1 (numcol*17)+2 (numcol*18)+1 (numcol*18)+2 (numcol*19)+1 (numcol*19)+2]);
for ll = 1:length(sessID)
    cellPos(ll) = find(sortsessID == sessID(ll) & sortcellID == cellID(ll));
    line([0 linPos(end)],[cellPos(ll) cellPos(ll)],'Color',col(ll,:),'LineWidth',2)
    hold on
end

%% Plot place tuned cells

ss = 2;
idxSess = Summary.AllsessType==(ss-1) & Summary.AllcellType == 1;
%idxMaps{1} = idxSess & Summary.AlllinField & Summary.AlllinCorr>0.1;
idxMaps = [];
idxMaps{1} = idxSess & Summary.AllspaceField & Summary.AllspatialCorr>0.1;

for ii = 1%:2            
    selectedlinMap = Summary.AlllinMapInit(idxMaps{ii},:);
    selectedlinMapEnd = Summary.AlllinMapEnd(idxMaps{ii},:);
    selectedspaceMap = Summary.AllspaceMap(idxMaps{ii},:);
    
    [maxLin,idxLin] = max(selectedlinMap,[],2);
    [maxSpace,idxSpace] = max(selectedspaceMap,[],2);    
    
    normlinMap = (selectedlinMap-nanmean(selectedlinMap,2))./nanstd(selectedlinMap,[],2);
    normlinMapEnd = (selectedlinMapEnd-nanmean(selectedlinMapEnd,2))./nanstd(selectedlinMapEnd,[],2);
    normspaceMap = (selectedspaceMap-nanmean(selectedspaceMap,2))./nanstd(selectedspaceMap,[],2);
      
    % if ii ==1            
    %     [~,sortidx] = sort(idxLin,'ascend');
    % elseif ii == 2            
        [~,sortidx] = sort(idxSpace,'ascend');
    %end
         
    ax1 = subplot(numrows,numcol,[(numcol*16)+9+(numcol*(ii-1)) (numcol*16)+10+(numcol*(ii-1)) (numcol*17)+9+(numcol*(ii-1)) (numcol*17)+10+(numcol*(ii-1)) (numcol*18)+9+(numcol*(ii-1)) (numcol*18)+10+(numcol*(ii-1)) (numcol*19)+9+(numcol*(ii-1)) (numcol*19)+10+(numcol*(ii-1))]);
    imagesc(linPos, 1:length(sortidx),normlinMap(sortidx,:))
    colormap(ax1,YlGnBu)
    ylabel(strcat('Cell ID ',num2str(length(sortidx))))
    caxis([-1 3.5])
    if ii==1
        title('No tone I')
    end
    
    ax1 = subplot(numrows,numcol,[(numcol*16)+11+(numcol*(ii-1)) (numcol*16)+12+(numcol*(ii-1)) (numcol*17)+11+(numcol*(ii-1)) (numcol*17)+12+(numcol*(ii-1)) (numcol*18)+11+(numcol*(ii-1)) (numcol*18)+12+(numcol*(ii-1)) (numcol*19)+11+(numcol*(ii-1)) (numcol*19)+12+(numcol*(ii-1))]);
    imagesc(linPos, 1:length(sortidx),normspaceMap(sortidx,:))
    colormap(ax1,YlGnBu)
    ylabel(strcat('Cell ID ',num2str(length(sortidx))))
    xlabel('Position')
    axis off
    caxis([-1 3.5])
    if ii==1
        title('Tone')
    end
    
    ax1 = subplot(numrows,numcol,[(numcol*16)+13+(numcol*(ii-1)) (numcol*16)+14+(numcol*(ii-1)) (numcol*17)+13+(numcol*(ii-1)) (numcol*17)+14+(numcol*(ii-1)) (numcol*18)+13+(numcol*(ii-1)) (numcol*18)+14+(numcol*(ii-1)) (numcol*19)+13+(numcol*(ii-1)) (numcol*19)+14+(numcol*(ii-1))]);       
    imagesc(linPos, 1:length(sortidx),normlinMapEnd(sortidx,:))
    colormap(ax1,YlGnBu)
    ylabel('Cell ID')
    xlabel('Position')
    axis off
    caxis([-1 3.5])
    if ii == 1
        title('No tone II')
    end
end   

%% Scatter of place versus tone tuning
subplot(numrows,numcol,[(numcol*25)+1 (numcol*25)+2 (numcol*26)+1 (numcol*26)+2 (numcol*27)+1 (numcol*27)+2]);
idxSess = Summary.AllsessType==1 & Summary.AllcellType == 1 & Summary.AllspacepeakRate > 2 & Summary.AlltonepeakRate > 2;
scatter(Summary.AllspatialCorr(idxSess),Summary.AlltoneCorr(idxSess),2,[0.5 0.5 0.5],'filled','MarkerFaceAlpha',0.4)
hold on
idxSess = Summary.AllsessType==1 & Summary.AllcellType == 1 & (Summary.AllspaceField==1) & (Summary.AllspatialCorr >0.1) ...
    & Summary.AllspacepeakRate > 2 & Summary.AlltonepeakRate > 2;
scatter(Summary.AllspatialCorr(idxSess),Summary.AlltoneCorr(idxSess),2,'k','filled')
idxSess = Summary.AllsessType==1 & Summary.AllcellType == 1 & (Summary.AlltoneField==1) & (Summary.AlltoneCorr >0.1) ...
    & Summary.AllspacepeakRate > 2 & Summary.AlltonepeakRate > 2;
scatter(Summary.AllspatialCorr(idxSess),Summary.AlltoneCorr(idxSess),2,'m','filled')
xlim([-0.2 1])
ylim([-0.2 1])
box off
xlabel('Spatial corr')
ylabel('Tone corr')
refline(1)

%% PGAM mutual information distribution - mixed selectivity
SummaryPGAM = compileTuningProportion('plotfig',false);
idxBoth  =  SummaryPGAM.sigAll(:,1) & SummaryPGAM.sigAll(:,4);
infoExtract = SummaryPGAM.mutInfoAll(idxBoth,:);
subplot(numrows,numcol,[(numcol*25)+4 (numcol*25)+5 (numcol*26)+4 (numcol*26)+5 (numcol*27)+4 (numcol*27)+5])
data = log2((infoExtract(:,1)./infoExtract(:,4)));%./(infoExtract(:,1)+infoExtract(:,4));
histogram(data,-3:0.3:5,'Normalization','probability')
xlabel('Log(2) of mutual info')
ylabel('Proportion')
box off

% Save figure
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\MainFigures\';
saveas(gcf,strcat(expPath,'Figure2A_plotSingleCells.png'));
saveas(gcf,strcat(expPath,'Figure2A_plotSingleCells.eps'),'epsc');
saveas(gcf,strcat(expPath,'Figure2A_plotSingleCells.fig'));


%% Make a new figure with the assembly stuff
load('Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220828_sess4\sessAssemblyPlot.mat')

fig2 = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
%set(fig2,'Position',[1921 41 1920 970]);

cellID = find(keepCells);
toneCell = toneCellLog(keepCells);
placeCell = placeCellLog(keepCells);

[~,maxPSTH]  = max(psth1,[],2);
[~,sortidx] = sort(maxPSTH,'descend');

%Plot a "place assembly" - assembly 4
subplot(2,5,1)
activation = Vectors(:,4);
[~,sortVec] = sort(activation,'ascend');
sortedTone = logical(toneCell(sortVec));
sortedPlace = logical(placeCell(sortVec));
sortedAct = activation(sortVec);
stem(1:length(cellID),sortedAct,'filled','color',[0.5 0.5 0.5],'MarkerSize',2)
hold on
stem(find(sortedTone),sortedAct(sortedTone),'filled','color','m','MarkerSize',2,'LineWidth',1.5)
stem(find(sortedPlace),sortedAct(sortedPlace),'filled','color','k','MarkerSize',2,'LineWidth',1.5)
ylim([-0.3 0.55])
xlim([1 length(cellID)])
camroll(90)
box off
set(gca,'XDir','reverse')

%Plot a "Tone assembly" - assembly 5
subplot(2,5,2)
activation = Vectors(:,5);
[~,sortVec] = sort(activation,'ascend');
sortedTone = logical(toneCell(sortVec));
sortedPlace = logical(placeCell(sortVec));
sortedAct = activation(sortVec);
stem(1:length(cellID),sortedAct,'filled','color',[0.5 0.5 0.5],'MarkerSize',2)
hold on
stem(find(sortedTone),sortedAct(sortedTone),'filled','color','m','MarkerSize',2,'LineWidth',1.5)
stem(find(sortedPlace),sortedAct(sortedPlace),'filled','color','k','MarkerSize',2,'LineWidth',1.5)
ylim([-0.3 0.55])
xlim([1 length(cellID)])
camroll(90)
box off
set(gca,'XDir','reverse')

%Show PSTH
t = linspace(-3, 3, 61);
ax = subplot(2,5,[3 4]);
imagesc(t,1:size(psth1,1),zscore(psth1(sortidx,:),[],2))
set(gca,'YDir','normal')
RdPu=cbrewer('seq', 'RdPu', 11);
colormap(ax,RdPu)
caxis([-1 5])
colorbar
hold on
line([0 0],[1 size(psth1,1)],'Color','k','LineWidth',1.5)
line([t(1) t(end)],[7 7])
line([t(1) t(end)],[17 17])

% Fraction of tone/place
subplot(2,5,5);
for aa = 1:length(assembliesID)    
    fractTonePlace(aa,1) = sum(Vectors(logical(toneCell),aa));
    fractTonePlace(aa,2) = sum(Vectors(logical(placeCell),aa));   
end
plot(fractTonePlace(sortedtoKeep,1),1:length(sortedtoKeep),'Color','m','LineWidth',1.5)
hold on
plot(fractTonePlace(sortedtoKeep,2),1:length(sortedtoKeep),'Color','k','LineWidth',1.5)
ylim([1 length(sortedtoKeep)])

% Plot average
%[AllPlaceWt, AllToneWt, tBins] = calcSessionCellTypeAssemblies;
load('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Assemblycontribution.mat')
zero_rows = all(AllPlaceWt == 0, 2);
AllPlaceWt = AllPlaceWt(~zero_rows, :);
AllToneWt = AllToneWt(~zero_rows, :);

subplot(2,5,[6 7])
meanpsth = nanmean(AllPlaceWt,1);
stdpsth = nanstd(AllPlaceWt,1)./sqrt(size(AllPlaceWt,1));
lArr  = meanpsth-stdpsth;
uArr = meanpsth+stdpsth;

fill([tBins; flipud(tBins)],[lArr'; flipud(uArr')],[0.2 0.2 0.2],'linestyle','none','FaceAlpha',0.5);                    
hold on
hi = line(tBins,meanpsth,'LineWidth',1,'Color',[0.2 0.2 0.2]);


meanpsth = nanmean(AllToneWt,1);
stdpsth = nanstd(AllToneWt,1)./sqrt(size(AllToneWt,1));
lArr  = meanpsth-stdpsth;
uArr = meanpsth+stdpsth;

fill([tBins; flipud(tBins)],[lArr'; flipud(uArr')],'m','linestyle','none','FaceAlpha',0.5);                    
hold on
hi = line(tBins,meanpsth,'LineWidth',1,'Color','m');
ylim([-0.4 1.3])
xlim([-2 0.2])
line([0 0],[-0.4 1.3],'Color','r','LineWidth',1.2)
box off

% Quantify each bin
data1 = AllPlaceWt(:,12);
data2 = AllToneWt(:,12);

subplot(2,5,8)
Summary.bin7assembly = groupStats([{data1} {data2}],[],'repeatedMeasures',true,'inAxis',true,'Color',[0 0 0;1 0 1]);
ylim([-1 3])
ylabel('Sum of weights')
title('T = -0.8')
box off

data1 = AllPlaceWt(:,19);
data2 = AllToneWt(:,19);
subplot(2,5,9)
Summary.bin14assembly = groupStats([{data1} {data2}],[],'repeatedMeasures',true,'inAxis',true,'Color',[0 0 0;1 0 1]);
ylim([-1 3])
ylabel('Sum of weights')
title('T = -0.1')
box off

% Save figure
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\MainFigures\';
saveas(gcf,strcat(expPath,'Figure2B_Assembly.png'));
saveas(gcf,strcat(expPath,'Figure2B_Assembly.eps'),'epsc');
saveas(gcf,strcat(expPath,'Figure2B_Assembly.fig'));
save(strcat(expPath,'Figure2B_Assembly.mat'),'Summary'); 
end
