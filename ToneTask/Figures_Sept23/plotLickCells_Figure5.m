function plotLickCells_Figure5

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[680 42 975 962]);

%rows = 22, columns = 6
numrows = 14;
numcol = 12;

%Summary = examinePokePSTHs('plotfig',false);
%Summary = examinePokePSTHs_usingPGAM('plotfig',false);
load('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\LickPSTHSummary.mat')

col = [119/255 221/255 229/255;...
    122/255 201/255 229/255;...
    38/255 169/255 224/255;...
    73/255 136/255 189/255;...
    17/255 55/255 174/255;...
    0/255 0/255 134/255];


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
scatter(behavTrials.timestamps((behavTrials.lickLoc==0),2),ones(1,sum(behavTrials.lickLoc==0))*130,10,[176/243 223/243 229/243],'filled')
scatter(behavTrials.timestamps((behavTrials.lickLoc==1),2),ones(1,sum(behavTrials.lickLoc==1))*130,10,[149/243 200/243 216/243],'filled')
scatter(behavTrials.timestamps((behavTrials.lickLoc==2),2),ones(1,sum(behavTrials.lickLoc==2))*130,10,[137/243 207/243 240/243],'filled')
scatter(behavTrials.timestamps((behavTrials.lickLoc==3),2),ones(1,sum(behavTrials.lickLoc==3))*130,10,[70/243 130/243 180/243],'filled')
scatter(behavTrials.timestamps((behavTrials.lickLoc==4),2),ones(1,sum(behavTrials.lickLoc==4))*130,10,[16/243 52/243 166/243],'filled')
scatter(behavTrials.timestamps((behavTrials.lickLoc==5),2),ones(1,sum(behavTrials.lickLoc==5))*130,10,[0/243 0/243 128/243],'filled')   
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
colormap(ax1,spec)
caxis([-1 2])
colorbar
xlim([-1 1])
ylabel(strcat('cellID',num2str(length(idxmax))))
title('choice', 'Color',col(1,:))
xticks([])
a = norm(:,idxT);
avgRate(:,1) = max(a,[],2);%
plotAvgStd(selectedPSTH,numrows,numcol,plotloc,fig2,timeaxis',col(1,:))


ax1 = subplot(numrows,numcol,[3+5*numcol 4+5*numcol 3+6*numcol 4+6*numcol]);
selectedPSTH = Summary.psthReward.lickTypes{1};
norm = zscore(selectedPSTH,[],2);
imagesc(timeaxis,1:length(idxmax),norm(idxmax,:))
colormap(ax1,spec)
caxis([-1 2])
xlim([-1 1])
xticks([])
yticks([])
title('no-tone end', 'Color',col(2,:))
a = norm(:,idxT);
avgRate(:,2) = max(a,[],2);%
plotAvgStd(selectedPSTH,numrows,numcol,plotloc,fig2,timeaxis',col(2,:))

ax1 = subplot(numrows,numcol,[5+5*numcol 6+5*numcol 5+6*numcol 6+6*numcol]);
selectedPSTH = Summary.psthReward.lickTypes{3};
norm = zscore(selectedPSTH,[],2);
imagesc(timeaxis,1:length(idxmax),norm(idxmax,:))
colormap(ax1,spec)
caxis([-1 2])
xlim([-1 1])
xticks([])
yticks([])
title('Spont fwd', 'Color',col(3,:))
a = norm(:,idxT);
avgRate(:,3) = max(a,[],2);%
plotAvgStd(selectedPSTH,numrows,numcol,plotloc,fig2,timeaxis',col(3,:))

ax1 = subplot(numrows,numcol,[7+5*numcol 8+5*numcol 7+6*numcol 8+6*numcol]);
selectedPSTH = Summary.psthReward.lickTypes{4};
norm = zscore(selectedPSTH,[],2);
imagesc(timeaxis,1:length(idxmax),norm(idxmax,:))
colormap(ax1,spec)
caxis([-1 2])
xlim([-1 1])
xticks([])
yticks([])
title('Spont rev', 'Color',col(4,:))
a = norm(:,idxT);
avgRate(:,4) = max(a,[],2);%
plotAvgStd(selectedPSTH,numrows,numcol,plotloc,fig2,timeaxis',col(4,:))

ax1 = subplot(numrows,numcol,[9+5*numcol 10+5*numcol 9+6*numcol 10+6*numcol]);
selectedPSTH = Summary.psthReward.lickTypes{5};
norm = zscore(selectedPSTH,[],2);
imagesc(timeaxis,1:length(idxmax),norm(idxmax,:))
colormap(ax1,spec)
caxis([-1 2])
xlim([-1 1])
xticks([])
yticks([])
title('Home', 'Color',col(5,:))
a = norm(:,idxT);
avgRate(:,5) = max(a,[],2);%
plotAvgStd(selectedPSTH,numrows,numcol,plotloc,fig2,timeaxis',col(5,:))
ylim([0 12])

%% Forward correct vs incorrect
plotloc = [1+8*numcol 2+8*numcol 1+9*numcol 2+9*numcol];
ax1 = subplot(numrows,numcol,plotloc);
selectedPSTH = Summary.psthReward.lickTypes{6};
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
plotAvgStd(selectedPSTH,numrows,numcol,plotloc,fig2,timeaxis',[70/255 148/255 73/255])
selectedPSTH = Summary.psthReward.lickTypes{7};
plotAvgStd(selectedPSTH,numrows,numcol,plotloc,fig2,timeaxis',[175/243 54/243 60/243])
line([0 0],[0 15],'Color','r')
ylim([0 10])
ylim([0 12])  
box off
title('Corr/incorr choice')

% Return correct vs incorrect
plotloc = [7+8*numcol 8+8*numcol 7+9*numcol 8+9*numcol];
selectedPSTH = Summary.psthReward.lickTypes{8};
plotAvgStd(selectedPSTH,numrows,numcol,plotloc,fig2,timeaxis',[70/255 148/255 73/255])
selectedPSTH = Summary.psthReward.lickTypes{9};
plotAvgStd(selectedPSTH,numrows,numcol,plotloc,fig2,timeaxis',[175/243 54/243 60/243])
line([0 0],[0 15],'Color','r')
ylim([0 10])
ylim([0 12])  
box off
title('Corr/incorr home')


%% Distribution across choice ports for choice cells

col = [119/255 221/255 229/255;...
    122/255 201/255 229/255;...
    38/255 169/255 224/255;...
    73/255 136/255 189/255;...
    17/255 55/255 174/255;...
    0/255 0/255 134/255];

idxT = timeaxis<=0.3 & timeaxis>=-0.3;

selectedPSTH = Summary.psthReward.lickTypes{2};
[maxRate,maxRateIdx] = max(selectedPSTH,[],2);
[~,idxmax] = sort(maxRateIdx,'ascend');  

for tt = 1:6
%     if tt == 1
%         temp = Summary.psthReward.lickTypes{tt+9};
%         [maxRate,maxRateIdx] = max(temp,[],2);
%         [~,idxmax] = sort(maxRateIdx,'ascend');
%     end
    normPSTH = zscore(Summary.psthReward.lickTypes{tt+9},[],2);
    plotloc = [2*(tt-1)+1+11*numcol 2*(tt-1)+2+12*numcol 2*(tt-1)+1+11*numcol 2*(tt-1)+2+12*numcol];
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
    temp = Summary.psthReward.lickTypes{tt+9};
    plotAvgStd(temp,numrows,numcol,plotloc,fig2,timeaxis',col(tt,:))    
    ylim([0 12])
    xlim([-1 1])  
    box off
    line([0 0],[0 15],'Color','r')
    
    norm = normPSTH(:,idxT);
    avgRate = nanmedian(norm,2);
    numProp(tt) = sum(avgRate>0.5)./length(avgRate);
end

plotloc = [11+8*numcol 12+8*numcol 11+9*numcol 12+9*numcol];
subplot(numrows,numcol,plotloc)
b = bar(numProp, 'FaceColor','flat','EdgeColor','none');

for k = 1:6
    b.CData(k,:) = col(k,:);
end

box off

%% Save figure

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\Figure5_lickCells.png'));
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\Figure5_lickCells.eps'),'epsc');
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\Figure5_lickCells.fig'));
% save(strcat(expPath,'\Compiled\Figures\Figure4_taskTuning2.mat'),'Fig4Stats'); 

end

function plotAvgStd(array,numrows,numcol,subplotlocation,figureHandle,xAxis,col)

    subplot(numrows, numcol, subplotlocation, 'Parent', figureHandle);
    meanpsth = nanmean(array,1);
    stdpsth = nanstd(array,1)./sqrt(size(array,1));         
    fill([xAxis; flipud(xAxis)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],col,'linestyle','none','FaceAlpha',0.5);                    
    hold on
    hi = line(xAxis,meanpsth,'LineWidth',1,'Color',col);

end