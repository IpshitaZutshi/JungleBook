function plotLickCells_Figure5

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[680 42 975 962]);

%rows = 22, columns = 6
numrows = 22;
numcol = 16;

%Summary = examinePokePSTHs('plotfig',false);
%Summary = examinePokePSTHs_usingPGAM('plotfig',false);
load('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\LickPSTHSummary.mat')

%% Panel A: Example lick tuned cell 1

sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\IZ43_220828_sess4';
cd(sessloc)
cellNum = 28;
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
plotExampleCellErrorTrials(1,3, spikeData, tracking, behavTrials, firingMaps,errorMaps, fig2,numrows,numcol,cellNum,0)


%% Panel B: Average lick tuning
timeaxis = linspace(-2,2,41);
idxT = timeaxis<=0.3 & timeaxis>=-0.3;
spec=cbrewer('seq', 'Greys', 20);
spec(spec>1) = 1;
spec(spec<0) = 0;
% spec2 = spec(end:-1:1,:);
for mm = 1:2
    if mm == 1
        selectedPSTH = Summary.psthReward.lickTypes{1};       
        norm = zscore(selectedPSTH,[],2);
        a = norm(:,idxT);        
        avgRate1 = max(a,[],2);%mean(a,2);%
        
        selectedPSTH = Summary.psthReward.lickTypes{2};
        [maxRate,maxRateIdx] = max(selectedPSTH,[],2);
        [~,idxmax] = sort(maxRateIdx,'ascend'); 
        norm = zscore(selectedPSTH,[],2);
        a = norm(:,idxT);       
        avgRate2 = max(a,[],2);%mean(a,2);%
    else
        selectedPSTH = Summary.psthReward.lickTypes{5};
        [maxRate,maxRateIdx] = max(selectedPSTH,[],2);
        [~,idxmax] = sort(maxRateIdx,'ascend');     
        norm = zscore(selectedPSTH,[],2);        
        a = norm(:,idxT);
        avgRate5 = max(a,[],2);%mean(a,2);%
    end
    
    ax1 = subplot(numrows,numcol,[(numcol*2*(mm-1))+5 (numcol*2*(mm-1))+6 (numcol*2*(mm-1))+numcol+5 (numcol*2*(mm-1))+numcol+6]);
    normPSTH2 = zscore(Summary.psthReward.lickTypes{2},[],2);
    imagesc(timeaxis,1:length(idxmax),normPSTH2(idxmax,:))
    colormap(ax1,spec)
    caxis([-1 2])
    xlim([-1 1])
    ylabel(strcat('cellID',num2str(length(idxmax))))
    xticks([])
    
    ax1 = subplot(numrows,numcol,[(numcol*2*(mm-1))+7 (numcol*2*(mm-1))+8 (numcol*2*(mm-1))+numcol+7 (numcol*2*(mm-1))+numcol+8]);
    normPSTH5 = zscore(Summary.psthReward.lickTypes{5},[],2);
    imagesc(timeaxis,1:length(idxmax),normPSTH5(idxmax,:))
    colormap(ax1,spec)
    axis off
    caxis([-1 2])
    xlim([-1 1])
    
%     ax1 = subplot(numrows,numcol,[(numcol*2*(mm-1))+9 (numcol*2*(mm-1))+10 (numcol*2*(mm-1))+numcol+9 (numcol*2*(mm-1))+numcol+10]);
%     normPSTH5 = zscore(Summary.psthReward.lickTypes{5},[],2);
%     imagesc(timeaxis,1:length(idxmax),normPSTH5(idxmax,:))
%     colormap(ax1,spec)
%     axis off
%     caxis([-1 2])
%     xlim([-1 1])  
end

idxK = kmeans([avgRate2 avgRate5], 3);
mean1 = mean(avgRate2(idxK==1));
mean2 = mean(avgRate2(idxK==2));
mean3 = mean(avgRate2(idxK==3));

[~,idx] = min([mean1,mean2,mean3]);

if idx == 1
    homeLickidx = idxK == 1;
    meanA = mean(avgRate5(idxK==2));
    meanB = mean(avgRate5(idxK==3));
    [~,idx1] = min([meanA,meanB]);
    if idx1==1
       choiceLickidx = idxK == 2;
       bothLickidx = idxK == 3;
    else
       choiceLickidx = idxK == 3;
       bothLickidx = idxK == 2;
    end
elseif idx==2
    homeLickidx = idxK == 2;
    meanA = mean(avgRate5(idxK==1));
    meanB = mean(avgRate5(idxK==3));
    [~,idx1] = min([meanA,meanB]);
    if idx1==1
       choiceLickidx = idxK == 1;
       bothLickidx = idxK == 3;
    else
       choiceLickidx = idxK == 3;
       bothLickidx = idxK == 1;
    end
else
    homeLickidx = idxK == 3;
    meanA = mean(avgRate5(idxK==1));
    meanB = mean(avgRate5(idxK==2));
    [~,idx1] = min([meanA,meanB]);
    if idx1==1
       choiceLickidx = idxK == 1;
       bothLickidx = idxK== 2;
    else
       choiceLickidx = idxK == 2;
       bothLickidx = idxK == 1;
    end    
end

% choiceLickidx = avgRate2>=1 & avgRate5<1 & avgRate1<1;
% homeLickidx = avgRate2<1& avgRate5>=1 & avgRate1<1;
% bothLickidx = avgRate2>=1 & avgRate1>=1 & avgRate5<1;


%% Panel C: Choice vs no tone end scatter
bins = -1:0.25:4;%-1.5:0.1:3; %
subplot(numrows,numcol,[11 11+numcol 11+(2*numcol)])
histogram(avgRate1,bins,'Normalization','probability');
view([270 270])
set(gca,'xDir','reverse')
%ylim([0 0.15])
box off

subplot(numrows,numcol,[12+(3*numcol) 13+(3*numcol) 14+(3*numcol)])
histogram(avgRate2,bins,'Normalization','probability')
% ylim([0 0.15])
box off

subplot(numrows,numcol, [12 13 14 12+(1*numcol) 13+(1*numcol) 14+(1*numcol) 12+(2*numcol) 13+(2*numcol) 14+(2*numcol)])
scatter(avgRate2,avgRate1,5,[220/243 220/243 220/243],'filled')
hold on
scatter(avgRate2(choiceLickidx),avgRate1(choiceLickidx),5,[56/255 61/255 150/255],'filled')
scatter(avgRate2(homeLickidx),avgRate1(homeLickidx),5,[152/243 152/243 152/243],'filled')
scatter(avgRate2(bothLickidx),avgRate1(bothLickidx),5,[52/243 52/243 52/243],'filled')
refline(1)
ylim([bins(1) bins(end)])
xlim([bins(1) bins(end)])


%% Panel D: Scatter of licks - choice versus home
subplot(numrows,numcol,[11+(5*numcol) 11+(6*numcol) 11+(7*numcol)])
histogram(avgRate5,bins,'Normalization','probability');
view([270 270])
set(gca,'xDir','reverse')
%ylim([0 0.15])
box off

subplot(numrows,numcol,[12+(8*numcol) 13+(8*numcol) 14+(8*numcol)])
histogram(avgRate2,bins,'Normalization','probability')
% ylim([0 0.15])
box off

subplot(numrows,numcol, [12+(5*numcol) 13+(5*numcol) 14+(5*numcol) 12+(6*numcol) 13+(6*numcol) 14+(6*numcol) 12+(7*numcol) 13+(7*numcol) 14+(7*numcol)])
scatter(avgRate2,avgRate5,5,[220/243 220/243 220/243],'filled')
hold on
scatter(avgRate2(choiceLickidx),avgRate5(choiceLickidx),5,[56/255 61/255 150/255],'filled')
scatter(avgRate2(homeLickidx),avgRate5(homeLickidx),5,[152/243 152/243 152/243],'filled')
scatter(avgRate2(bothLickidx),avgRate5(bothLickidx),5,[52/243 52/243 52/243],'filled')
hold on
refline(1)
ylim([bins(1) bins(end)])
xlim([bins(1) bins(end)])


%% Panel D: Average rate of each type of cell
plotAvgStd(normPSTH2(choiceLickidx,:),numrows,numcol,[15 16 15+(1*numcol) 16+(1*numcol)],fig2,timeaxis',[56/255 61/255 150/255])
hold on
plotAvgStd(normPSTH2(homeLickidx,:),numrows,numcol,[15 16 15+(1*numcol) 16+(1*numcol)],fig2,timeaxis',[152/243 152/243 152/243])
plotAvgStd(normPSTH2(bothLickidx,:),numrows,numcol,[15 16 15+(1*numcol) 16+(1*numcol)],fig2,timeaxis',[52/243 52/243 52/243])
box off
line([0 0],[-1 3],'Color','r')
ylim([-1 3])
xlim([-1 1])  
yticklabels([])
xticklabels([])
title(strcat(num2str(sum(choiceLickidx)./length(choiceLickidx)),'...',...
   num2str(sum(bothLickidx)./length(bothLickidx))))

plotAvgStd(normPSTH5(choiceLickidx,:),numrows,numcol,[15+(2*numcol) 16+(2*numcol) 15+(3*numcol) 16+(3*numcol)],fig2,timeaxis',[56/255 61/255 150/255])
hold on
plotAvgStd(normPSTH5(homeLickidx,:),numrows,numcol,[15+(2*numcol) 16+(2*numcol) 15+(3*numcol) 16+(3*numcol)],fig2,timeaxis',[152/243 152/243 152/243])
plotAvgStd(normPSTH5(bothLickidx,:),numrows,numcol,[15+(2*numcol) 16+(2*numcol) 15+(3*numcol) 16+(3*numcol)],fig2,timeaxis',[52/243 52/243 52/243])
box off
line([0 0],[-1 3],'Color','r')
ylim([-1 3])
xlim([-1 1])  
yticklabels([])
xticklabels([])
title(num2str(sum(homeLickidx)./length(homeLickidx)))

% plotAvgStd(normPSTH5(choiceLickidx,:),numrows,numcol,[15+(4*numcol) 16+(4*numcol) 15+(5*numcol) 16+(5*numcol)],fig2,timeaxis',[56/255 61/255 150/255])
% hold on
% plotAvgStd(normPSTH5(homeLickidx,:),numrows,numcol,[15+(4*numcol) 16+(4*numcol) 15+(5*numcol) 16+(5*numcol)],fig2,timeaxis',[152/243 152/243 152/243])
% plotAvgStd(normPSTH5(bothLickidx,:),numrows,numcol,[15+(4*numcol) 16+(4*numcol) 15+(5*numcol) 16+(5*numcol)],fig2,timeaxis',[52/243 52/243 52/243])
% box off
% line([0 0],[-1 3],'Color','r')
% ylim([-1 3])
% xlim([-1 1])  
% yticklabels([])
% title(num2str(sum(bothLickidx)./length(bothLickidx)))


%% Panel E: Correct vs incorrect
    
for ii = 1:15
    normPSTH{ii} = [];
    %[maxRate] = max(Summary.psthReward.lickTypes{ii},[],2);
    normPSTH{ii} = Summary.psthReward.lickTypes{ii};%./maxRate;
end

%Forward correct vs incorrect
idx = choiceLickidx;% | bothLickidx;
plotloc = [5+(numcol*8) 6+(numcol*8) 5+(numcol*9) 6+(numcol*9)];
plotAvgStd(normPSTH{6}(idx,:),numrows,numcol,plotloc,fig2,timeaxis',[56/255 61/255 150/255])
plotAvgStd(normPSTH{7}(idx,:),numrows,numcol,plotloc,fig2,timeaxis',[8/243 133/243 161/243])
line([0 0],[0 15],'Color','r')
ylim([0 10])
xlim([-1 1])  
box off
title('Correct/incorrect')

idx = homeLickidx;
plotloc = [7+(numcol*8) 8+(numcol*8) 7+(numcol*9) 8+(numcol*9)];
plotAvgStd(normPSTH{8}(idx,:),numrows,numcol,plotloc,fig2,timeaxis',[187/243 86/243 149/243])
plotAvgStd(normPSTH{9}(idx,:),numrows,numcol,plotloc,fig2,timeaxis',[193/243 90/243 99/243])
line([0 0],[0 15],'Color','r')
ylim([0 15])
xlim([-1 1])  
box off

%% Panel F: Do the choice responsive cells also fire for end ports, and spontaneous licks?
col = [56/255 61/255 150/255;...
    152/243 152/243 152/243;...
    52/255 52/255 52/255];

%End ports
plotloc = [5+(numcol*11) 6+(numcol*11) 5+(numcol*12) 6+(numcol*12)];
plotAvgStd(normPSTH{1}(choiceLickidx,:),numrows,numcol,plotloc,fig2,timeaxis',col(1,:))
plotAvgStd(normPSTH{1}(homeLickidx,:),numrows,numcol,plotloc,fig2,timeaxis',col(2,:))
plotAvgStd(normPSTH{1}(bothLickidx,:),numrows,numcol,plotloc,fig2,timeaxis',col(3,:))
line([0 0],[0 15],'Color','r')
ylim([0 10])
xlim([-1 1])  
box off
title('No tone end')

%Spontaneous forward
plotloc = [7+(numcol*11) 8+(numcol*11) 7+(numcol*12) 8+(numcol*12)];
plotAvgStd(normPSTH{3}(choiceLickidx,:),numrows,numcol,plotloc,fig2,timeaxis',col(1,:))
plotAvgStd(normPSTH{3}(homeLickidx,:),numrows,numcol,plotloc,fig2,timeaxis',col(2,:))
plotAvgStd(normPSTH{3}(bothLickidx,:),numrows,numcol,plotloc,fig2,timeaxis',col(3,:))
line([0 0],[0 15],'Color','r')
ylim([0 10])
xlim([-1 1])  
box off
title('Spont fwd')

%Spontaneous reverse
plotloc = [9+(numcol*11) 10+(numcol*11) 9+(numcol*12) 10+(numcol*12)];
plotAvgStd(normPSTH{4}(choiceLickidx,:),numrows,numcol,plotloc,fig2,timeaxis',col(1,:))
plotAvgStd(normPSTH{4}(homeLickidx,:),numrows,numcol,plotloc,fig2,timeaxis',col(2,:))
plotAvgStd(normPSTH{4}(bothLickidx,:),numrows,numcol,plotloc,fig2,timeaxis',col(3,:))
line([0 0],[0 15],'Color','r')
ylim([0 10])
xlim([-1 1])  
box off
title('Spont rev')

%% G: Distribution across choice ports for choice cells

col = [119/255 221/255 229/255;...
    122/255 201/255 229/255;...
    38/255 169/255 224/255;...
    73/255 136/255 189/255;...
    17/255 55/255 174/255;...
    0/255 0/255 134/255];

idx = choiceLickidx | bothLickidx;
idxT = timeaxis<=0.3 & timeaxis>=-0.3;

selectedPSTH = Summary.psthReward.lickTypes{2}(idx,:);
[maxRate,maxRateIdx] = max(selectedPSTH,[],2);
[~,idxmax] = sort(maxRateIdx,'ascend');  

for tt = 1:6
    if tt == 1
        temp = Summary.psthReward.lickTypes{tt+9}(idx,:);
        [maxRate,maxRateIdx] = max(temp,[],2);
        [~,idxmax] = sort(maxRateIdx,'ascend');
    end
    normPSTH = zscore(Summary.psthReward.lickTypes{tt+9}(idx,:),[],2);
    plotloc = [4+tt+(numcol*14) 4+tt+(numcol*15)];
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
    
    plotloc = [4+tt+(numcol*16) 4+tt+(numcol*17)];
    subplot(numrows,numcol,plotloc);
    temp = Summary.psthReward.lickTypes{tt+9}(idx,:);
    plotAvgStd(temp,numrows,numcol,plotloc,fig2,timeaxis',col(tt,:))    
    ylim([0 10])
    xlim([-1 1])  
    box off
    line([0 0],[0 15],'Color','r')
    
    norm = normPSTH(:,idxT);
    avgRate{tt} = nanmedian(norm,2);
    numProp(tt) = sum(avgRate{tt}>=1)./length(avgRate{tt});
    newpsth{tt} = Summary.psthReward.lickTypes{tt+9};
end

plotloc = [11+(numcol*14) 12+(numcol*15) 11+(numcol*14) 12+(numcol*15)];
subplot(numrows,numcol,plotloc)
b = bar(numProp, 'FaceColor','flat','EdgeColor','none');

for k = 1:6
    b.CData(k,:) = col(k,:);
end

box off

%% G: Distribution across choice ports for both cells

% idx = bothLickidx;
% idxT = timeaxis<=0.3 & timeaxis>=-0.3;
% 
% selectedPSTH = Summary.psthReward.lickTypes{2}(idx,:);
% [maxRate,maxRateIdx] = max(selectedPSTH,[],2);
% [~,idxmax] = sort(maxRateIdx,'ascend');  
% 
% for tt = 1:6
%     normPSTH = zscore(Summary.psthReward.lickTypes{tt+9}(idx,:),[],2);
%     plotloc = [4+tt+(numcol*18) 4+tt+(numcol*19)];
%     ax1 = subplot(numrows,numcol,plotloc);
%     imagesc(timeaxis,1:length(idx),normPSTH(idxmax,:));
%     caxis([-1 2])
%     xlim([-1 1])  
%     colormap(ax1,spec)
%     axis off
%     
%     plotloc = [4+tt+(numcol*20) 4+tt+(numcol*21)];
%     subplot(numrows,numcol,plotloc);
%     temp = Summary.psthReward.lickTypes{tt+9}(idx,:);
%     plotAvgStd(temp,numrows,numcol,plotloc,fig2,timeaxis',col(tt,:))    
%     ylim([0 16])
%     box off
%     xlim([-1 1])  
%     line([0 0],[0 16],'Color','r')
%     
%     norm = normPSTH(:,idxT);
%     avgRate{tt} = nanmedian(norm,2);
%     numProp(tt) = sum(avgRate{tt}>=0.4)./length(avgRate{tt});
%     newpsth{tt} = Summary.psthReward.lickTypes{tt+9};
% end
% 
% plotloc = [11+(numcol*18) 12+(numcol*19) 11+(numcol*18) 12+(numcol*19)];
% subplot(numrows,numcol,plotloc)
% b = bar(numProp, 'FaceColor','flat','EdgeColor','none');
% 
% for k = 1:6
%     b.CData(k,:) = col(k,:);
% end
% 
% box off

%% Save figure

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';
saveas(gcf,strcat(expPath,'\Compiled\Figures\Figure5_lickCells.png'));
saveas(gcf,strcat(expPath,'\Compiled\Figures\Figure5_lickCells.eps'),'epsc');
saveas(gcf,strcat(expPath,'\Compiled\Figures\Figure5_lickCells.fig'));
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