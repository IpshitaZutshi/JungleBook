function Fig3B_lickResponses

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[680 42 975 962]);

%rows = 22, columns = 6
numrows = 13;
numcol = 12;
useMedian = 0;
usePGAM = 0;

if ~usePGAM
    Summary = examinePokePSTHs('plotfig',false);
    %load('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\LickPSTHSummary.mat')
else
    %Summary = examinePokePSTHs_usingPGAM('plotfig',false);    
    load('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\LickPSTHSummaryPGAM.mat')
end

col = [83/255 0/255 0/255;...
    184/255 15/255 10/255;...
    241/255 114/255 42/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];


%% Panel A: Example lick tuned cell 1

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220828_sess4';
cd(sessloc)
cellNum = 28;
file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);
file = dir(['*.spikes.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file.name);

subplot(numrows, numcol, [1:1:10 (1:1:10)+numcol (1:1:10)+(2*numcol) (1:1:10)+(3*numcol)])
plot(tracking.timestamps,tracking.position.y,'Color',[160/243 160/243 160/243], 'LineWidth',1)
hold on
box off
%plot(tracking.timestamps,tracking.position.v)
scatter(tracking.timestamps(spikeData.posIdx{cellNum}),tracking.position.y(spikeData.posIdx{cellNum}),5,'r','filled')
scatter(behavTrials.timestamps((behavTrials.linTrial==1),2),ones(1,sum(behavTrials.linTrial==1))*120,10,[0.5 0.5 0.5],'filled')
scatter(behavTrials.timestamps((behavTrials.lickLoc==0),2),ones(1,sum(behavTrials.lickLoc==0))*130,10,[83/255 0/255 0/255],'filled')
scatter(behavTrials.timestamps((behavTrials.lickLoc==1),2),ones(1,sum(behavTrials.lickLoc==1))*130,10,[184/255 15/255 10/255],'filled')
scatter(behavTrials.timestamps((behavTrials.lickLoc==2),2),ones(1,sum(behavTrials.lickLoc==2))*130,10,[241/255 114/255 42/255],'filled')
scatter(behavTrials.timestamps((behavTrials.lickLoc==3),2),ones(1,sum(behavTrials.lickLoc==3))*130,10,[249/255 197/255 81/255],'filled')
scatter(behavTrials.timestamps((behavTrials.lickLoc==4),2),ones(1,sum(behavTrials.lickLoc==4))*130,10,[143/255 189/255 107/255],'filled')
scatter(behavTrials.timestamps((behavTrials.lickLoc==5),2),ones(1,sum(behavTrials.lickLoc==5))*130,10,[87/255 116/255 144/255],'filled')   
scatter(behavTrials.timestamps((behavTrials.correct==1),2),ones(1,sum(behavTrials.correct==1))*140,10,[70/243 148/243 73/243],'filled')
scatter(behavTrials.timestamps((behavTrials.correct==0 & behavTrials.linTrial==0),2),ones(1,sum(behavTrials.correct==0 & behavTrials.linTrial==0))*140,10,[187/243 86/243 149/243],'filled')
ylim([0 150])
xlim([5850 6700])
xlabel('Time(s)')
ylabel('Position on track (cm)')

dataMat = [];
dataMatTone = [];

b = linspace(0,125,50);
a = linspace(2000,22000,50);

for kk = 1:6
    subplot(numrows, numcol, [11 12]);
    hold on        
    plot(b,firingMaps.forward.rateMaps{cellNum}{1+kk},'Color',col(kk,:),'LineWidth',1);
    xlim([0 122])
    set(gca,'xtick',[])
    box off    
    
%     subplot(numrows, numcol, [11+numcol 12+numcol]);
%     hold on
%     plot(a,firingMaps.tone.rateMaps{cellNum}{1+kk},'Color',col(kk,:),'LineWidth',1);      
%     xlim([1000 25000])
%     set(gca,'xtick',[])
%     box off
    
    dataMat = [dataMat;firingMaps.forward.rateMaps{cellNum}{1+kk}];
    atm = fillmissing(firingMaps.tone.rateMaps{cellNum}{1+kk},'linear');
    dataMatTone = [dataMatTone;atm];
    
end
    
YlGnBu=cbrewer('seq', 'YlGnBu', 11);
ax1 = subplot(numrows, numcol, [11+numcol 12+numcol]);
imagesc(1,b,nanmean(dataMat,1))
colormap(ax1,YlGnBu)
box off
set(gca,'xtick',[],'ytick',[])
xlabel(num2str(max(nanmean(dataMat,1))))
% 
% BuPu=cbrewer('seq', 'BuPu', 11);
% ax2 = subplot(numrows, numcol,[11+3*numcol 12+3*numcol]);
% imagesc(1,a,nanmean(dataMatTone,1))
% colormap(ax2, BuPu)
% box off
% set(gca,'xtick',[],'ytick',[])
% xlabel(num2str(max(nanmean(dataMatTone,1))))

% PSTH 
for ii = 1:6
    idx = behavTrials.lickLoc==(ii-1) & behavTrials.linTrial ==0; 
    st = behavTrials.timestamps(idx,2);
    [stccg, tPSTH] = CCG({spikes.times{cellNum} st},[],'binSize',0.1,'duration',2,'norm','rate');
    subplot(numrows, numcol, [11+2*numcol 12+2*numcol]);
    hold on
    plot(tPSTH,stccg(:,2,1)','Color',col(ii,:),'LineWidth',1);   
end

% Home port, end trial
idx = behavTrials.lickLoc==5 & behavTrials.linTrial ==1; 
st = behavTrials.timestamps(idx,2);
[stccg, tPSTH] = CCG({spikes.times{cellNum} st},[],'binSize',0.1,'duration',2,'norm','rate');
subplot(numrows, numcol, [11+3*numcol 12+3*numcol]);
hold on
plot(tPSTH,stccg(:,2,1)','Color',[52/243 52/243 52/243],'LineWidth',1);   

st = behavTrials.timestamps(:,1);
[stccg, tPSTH] = CCG({spikes.times{cellNum} st},[],'binSize',0.1,'duration',2,'norm','rate');
subplot(numrows, numcol, [11+3*numcol 12+3*numcol]);
plot(tPSTH,stccg(:,2,1)','Color',[187/243 86/243 149/243],'LineWidth',1);   
ylim([0 50])

%% Panel B: Average lick tuning
col = [56/243 61/243 150/243;...
    52/243 52/243 52/243;...
    8/243 133/243 161/243;...
    193/243 90/243 99/243;...
    187/243 86/243 149/243];

timeaxis = linspace(-1,1,21);
idxT = timeaxis<=0.3 & timeaxis>=-0.3;
spec=cbrewer('seq', 'Blues', 20);
spec(spec>1) = 1;
spec(spec<0) = 0;
% spec = spec(end:-1:1,:);
plotloc = [11+5*numcol 12+5*numcol 11+6*numcol 12+6*numcol];

ax1 = subplot(numrows,numcol,[1+5*numcol 2+5*numcol 1+6*numcol 2+6*numcol]);
selectedPSTH = Summary.psthReward.lickTypes{2};
[maxRate,maxRateIdx] = max(selectedPSTH,[],2);
[~,idxmax] = sort(maxRateIdx,'ascend'); 
norm = zscore(selectedPSTH,[],2);
imagesc(timeaxis,1:length(idxmax),norm(idxmax,:))
%imagesc(timeaxis,1:length(idxmax),selectedPSTH(idxmax,:))
colormap(ax1,spec)
caxis([-1 2])
colorbar
xlim([-1 1])
ylabel(strcat('cellID',num2str(length(idxmax))))
title('choice', 'Color',col(1,:))
xticks([])
a = selectedPSTH(:,idxT);
avgRate(:,1) = mean(a,2);%
if ~usePGAM
    plotAvgStd(selectedPSTH,numrows,numcol,plotloc,fig2,timeaxis',col(1,:),useMedian)
else
    plotAvgStd(norm,numrows,numcol,plotloc,fig2,timeaxis',col(1,:),useMedian)
end

ax1 = subplot(numrows,numcol,[3+5*numcol 4+5*numcol 3+6*numcol 4+6*numcol]);
selectedPSTH = Summary.psthReward.lickTypes{1};
norm = zscore(selectedPSTH,[],2);
imagesc(timeaxis,1:length(idxmax),norm(idxmax,:))
%imagesc(timeaxis,1:length(idxmax),selectedPSTH(idxmax,:))
colormap(ax1,spec)
caxis([-1 2])
xlim([-1 1])
xticks([])
yticks([])
title('no-tone end', 'Color',col(2,:))
a = selectedPSTH(:,idxT);
avgRate(:,2) = mean(a,2);%
if ~usePGAM
    plotAvgStd(selectedPSTH,numrows,numcol,plotloc,fig2,timeaxis',col(2,:),useMedian)
else
    plotAvgStd(norm,numrows,numcol,plotloc,fig2,timeaxis',col(2,:),useMedian)
end

ax1 = subplot(numrows,numcol,[5+5*numcol 6+5*numcol 5+6*numcol 6+6*numcol]);
selectedPSTH = Summary.psthReward.lickTypes{3};
norm = zscore(selectedPSTH,[],2);
imagesc(timeaxis,1:length(idxmax),norm(idxmax,:))
%imagesc(timeaxis,1:length(idxmax),selectedPSTH(idxmax,:))
colormap(ax1,spec)
caxis([-1 2])
xlim([-1 1])
xticks([])
yticks([])
title('Spont fwd', 'Color',col(3,:))
a = selectedPSTH(:,idxT);
avgRate(:,3) = mean(a,2);%
if ~usePGAM
    plotAvgStd(selectedPSTH,numrows,numcol,plotloc,fig2,timeaxis',col(3,:),useMedian)
else
    plotAvgStd(norm,numrows,numcol,plotloc,fig2,timeaxis',col(3,:),useMedian)
end

ax1 = subplot(numrows,numcol,[7+5*numcol 8+5*numcol 7+6*numcol 8+6*numcol]);
selectedPSTH = Summary.psthReward.lickTypes{4};
norm = zscore(selectedPSTH,[],2);
imagesc(timeaxis,1:length(idxmax),norm(idxmax,:))
%imagesc(timeaxis,1:length(idxmax),selectedPSTH(idxmax,:))
colormap(ax1,spec)
caxis([-1 2])
xlim([-1 1])
xticks([])
yticks([])
title('Spont rev', 'Color',col(4,:))
a = selectedPSTH(:,idxT);
avgRate(:,4) = mean(a,2);%
if ~usePGAM
    plotAvgStd(selectedPSTH,numrows,numcol,plotloc,fig2,timeaxis',col(4,:),useMedian)
else
    plotAvgStd(norm,numrows,numcol,plotloc,fig2,timeaxis',col(4,:),useMedian)
end

ax1 = subplot(numrows,numcol,[9+5*numcol 10+5*numcol 9+6*numcol 10+6*numcol]);
selectedPSTH = Summary.psthReward.lickTypes{5};
norm = zscore(selectedPSTH,[],2);
imagesc(timeaxis,1:length(idxmax),norm(idxmax,:))
%imagesc(timeaxis,1:length(idxmax),selectedPSTH(idxmax,:))
colormap(ax1,spec)
caxis([-1 2])
xlim([-1 1])
xticks([])
yticks([])
title('Home', 'Color',col(5,:))
a = selectedPSTH(:,idxT);
avgRate(:,5) = mean(a,2);%
if ~usePGAM
    plotAvgStd(selectedPSTH,numrows,numcol,plotloc,fig2,timeaxis',col(5,:),useMedian)
    ylim([2 8.5])
else
    plotAvgStd(norm,numrows,numcol,plotloc,fig2,timeaxis',col(5,:),useMedian)
    ylim([-0.5 1])
end

plotloc = [7+8*numcol 8+8*numcol 7+9*numcol 8+9*numcol];
subplot(numrows,numcol,plotloc)
Stats.lickType = groupStats([{avgRate(:,1)},{avgRate(:,2)},{avgRate(:,3)},{avgRate(:,4)},{avgRate(:,5)}],[],'Color',col,'inAxis',true,'labelSummary',false);

%% Forward correct vs incorrect
avgRate = [];

plotloc = [1+8*numcol 2+8*numcol 1+9*numcol 2+9*numcol];
ax1 = subplot(numrows,numcol,plotloc);
selectedPSTH = Summary.psthReward.lickTypes{6};
norm = zscore(selectedPSTH,[],2);
a = selectedPSTH(:,idxT);
avgRate(:,1) = mean(a,2);%
[maxRate,maxRateIdx] = max(selectedPSTH,[],2);
[~,idxmax] = sort(maxRateIdx,'ascend'); 
norm = zscore(selectedPSTH,[],2);
imagesc(timeaxis,1:length(idxmax),norm(idxmax,:))
colormap(ax1,spec)
caxis([-1 2])
xlim([-1 1])
ylabel(strcat('cellID',num2str(length(idxmax))))
title('Correct choice')
xticks([])

plotloc = [3+8*numcol 4+8*numcol 3+9*numcol 4+9*numcol];
ax1 = subplot(numrows,numcol,plotloc);
selectedPSTH = Summary.psthReward.lickTypes{7};
norm = zscore(selectedPSTH,[],2);
a = selectedPSTH(:,idxT);
avgRate(:,2) = mean(a,2);%
[maxRate,maxRateIdx] = max(selectedPSTH,[],2);
[~,idxmax] = sort(maxRateIdx,'ascend'); 
norm = zscore(selectedPSTH,[],2);
imagesc(timeaxis,1:length(idxmax),norm(idxmax,:))
colormap(ax1,spec)
caxis([-1 2])
xlim([-1 1])
ylabel(strcat('cellID',num2str(length(idxmax))))
title('Incorrect choice')
xticks([])

plotloc = [5+8*numcol 6+8*numcol 5+9*numcol 6+9*numcol];
selectedPSTH = Summary.psthReward.lickTypes{6};
norm = zscore(selectedPSTH,[],2);
if ~usePGAM
    plotAvgStd(selectedPSTH,numrows,numcol,plotloc,fig2,timeaxis',[70/255 148/255 73/255],useMedian)
else
    plotAvgStd(norm,numrows,numcol,plotloc,fig2,timeaxis',[70/255 148/255 73/255],useMedian)
end
hold on
selectedPSTH = Summary.psthReward.lickTypes{7};
norm = zscore(selectedPSTH,[],2);
if ~usePGAM
    plotAvgStd(selectedPSTH,numrows,numcol,plotloc,fig2,timeaxis',[175/243 54/243 60/243],useMedian)
    line([0 0],[2 7],'Color','r')
    ylim([2 8.5])
else
    plotAvgStd(norm,numrows,numcol,plotloc,fig2,timeaxis',[175/243 54/243 60/243],useMedian)
    line([0 0],[-1 1.5],'Color','r')
    ylim([-0.5 1]) 
end
box off
title('Corr/incorr choice')

Stats.correctIncorrect = groupStats([{avgRate(:,1)},{avgRate(:,2)}],[],'doPlot',false);

% Return correct vs incorrect
%plotloc = [7+8*numcol 8+8*numcol 7+9*numcol 8+9*numcol];
% selectedPSTH = Summary.psthReward.lickTypes{8};
% plotAvgStd(selectedPSTH,numrows,numcol,plotloc,fig2,timeaxis',[70/255 148/255 73/255])
% selectedPSTH = Summary.psthReward.lickTypes{9};
% plotAvgStd(selectedPSTH,numrows,numcol,plotloc,fig2,timeaxis',[175/243 54/243 60/243])
% line([0 0],[0 15],'Color','r')
% ylim([0 10])
% ylim([0 12])  
% box off
% title('Corr/incorr home')


%% Distribution across choice ports for choice cells

if ~usePGAM
    col = [83/255 0/255 0/255;...
        184/255 15/255 10/255;...
        241/255 114/255 42/255;...
        249/255 197/255 81/255;...
        143/255 189/255 107/255;...
        87/255 116/255 144/255];
else
    col = [56/243 61/243 150/243;...
        52/243 52/243 52/243;...
        187/243 86/243 149/243];
end

idxT = timeaxis<=0.3 & timeaxis>=-0.3;

selectedPSTH = Summary.psthReward.lickTypes{2};
[maxRate,maxRateIdx] = max(selectedPSTH,[],2);
[~,idxmax] = sort(maxRateIdx,'ascend');  

if ~usePGAM
    lickT = [10 11 12 13 14 15];
else   
    lickT = [2 1 5];
end

for tt = 1:length(lickT)
    normPSTH = zscore(Summary.psthReward.lickTypes{lickT(tt)},[],2);
    plotloc = [2*(tt-1)+1+11*numcol 2*(tt-1)+2+11*numcol 2*(tt-1)+1+12*numcol 2*(tt-1)+2+12*numcol];
    ax1 = subplot(numrows,numcol,plotloc);
    imagesc(timeaxis,1:length(idxmax),normPSTH(idxmax,:));
    ylabel(num2str(length(idxmax)))
    caxis([-1 2])
    xlim([-1 1])  
    colormap(ax1,spec)
    xticks([])
    if tt>1
        axis off
    end
    
    plotloc = [9+8*numcol 10+8*numcol 9+9*numcol 10+9*numcol];
    subplot(numrows,numcol,plotloc);
    temp = Summary.psthReward.lickTypes{lickT(tt)};
    norm = zscore(temp,[],2);
    
    if ~usePGAM
        plotAvgStd(temp,numrows,numcol,plotloc,fig2,timeaxis',col(tt,:),useMedian)    
        line([0 0],[2 7],'Color','r')
        ylim([2 8.5])        
    else
        plotAvgStd(norm,numrows,numcol,plotloc,fig2,timeaxis',col(tt,:),useMedian)    
        line([0 0],[-1 1.5],'Color','r')
        ylim([-0.5 1]) 
    end
    xlim([-1 1])  
    box off
    line([0 0],[2 7],'Color','r')
    
    norm = normPSTH(:,idxT);
    tempRate  = temp(:,idxT);
    avgRateTemp{tt} = log10(median(tempRate,2));
    
    for ss = 1:max(Summary.sessID)
        idxss = Summary.sessID==ss;       
        avgRate = nanmedian(norm(idxss,:),2);
        numProp(ss,tt) = sum(avgRate>0.1)./length(avgRate);
    end   
end

plotloc = [11+8*numcol 12+8*numcol 11+9*numcol 12+9*numcol];
subplot(numrows,numcol,plotloc)
if ~usePGAM
    Stats.portFract = groupStats([{numProp(:,1)},{numProp(:,2)},{numProp(:,3)},{numProp(:,4)},{numProp(:,5)},{numProp(:,6)}],[],'Color',col,'inAxis',true,'labelSummary',false,'repeatedMeasures',true);
    Stats.portRate= groupStats(avgRateTemp,[],'Color',col,'doPlot',false,'labelSummary',false,'repeatedMeasures',true);
else
  %  Stats.portFract = groupStats([{numProp(:,1)},{numProp(:,2)},{numProp(:,3)},{numProp(:,4)},{numProp(:,5)},{numProp(:,6)},{numProp(:,7)}],[],'Color',col,'inAxis',true,'labelSummary',false);
   Stats.portFract = groupStats([{numProp(:,1)},{numProp(:,2)},{numProp(:,3)}],[],'Color',col,'inAxis',true,'labelSummary',false);
end

box off

% Save figure

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\MainFigures\';
if ~usePGAM    
    saveas(gcf,strcat(expPath,'Figure3B_lickResponses.png'));
    saveas(gcf,strcat(expPath,'Figure3B_lickResponses.eps'),'epsc');
    saveas(gcf,strcat(expPath,'Figure3B_lickResponses.fig'));
    save(strcat(expPath,'Figure3B_lickResponses.mat'),'Stats'); 
else
    saveas(gcf,strcat(expPath,'Figure3B_lickResponsesPGAM.png'));
    saveas(gcf,strcat(expPath,'Figure3B_lickResponsesPGAM.eps'),'epsc');
    saveas(gcf,strcat(expPath,'Figure3B_lickResponsesPGAM.fig'));
    save(strcat(expPath,'Figure3B_lickResponsesPGAM.mat'),'Stats'); 
end

end

function plotAvgStd(array,numrows,numcol,subplotlocation,figureHandle,xAxis,col, useMedian)

    subplot(numrows, numcol, subplotlocation, 'Parent', figureHandle);
    if ~useMedian
        meanpsth = nanmean(array,1);
        stdpsth = nanstd(array,1)./sqrt(size(array,1));
        lArr  = meanpsth-stdpsth;
        uArr = meanpsth+stdpsth;
    else
        meanpsth = nanmedian(array,1);
        a = quantile(array,4,1);
        lArr  = meanpsth-a(2,:);
        uArr = meanpsth+a(3,:);
    end
    fill([xAxis; flipud(xAxis)],[lArr'; flipud(uArr')],col,'linestyle','none','FaceAlpha',0.5);                    
    hold on
    hi = line(xAxis,meanpsth,'LineWidth',1,'Color',col);

end