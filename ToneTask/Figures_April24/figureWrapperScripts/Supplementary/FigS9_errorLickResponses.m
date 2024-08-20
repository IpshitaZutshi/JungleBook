function FigS9_errorLickResponses

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[204 27 1500 800]);

BuPu=cbrewer('seq', 'BuPu', 11);
YlGnBu=cbrewer('seq', 'YlGnBu', 11);

numrows = 2;
numcol = 4;

col = [83/255 0/255 0/255;...
    184/255 15/255 10/255;...
    241/255 114/255 42/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];

Summary = latencytoReturn;

%% Latency of returning to home
data.Correct = Summary.latencyCorrect;
data.Incorrect = Summary.latencyIncorrect;

subplot(numrows, numcol, 1)
nhist(data,'proportion','samebins')
subplot(numrows, numcol, 2)
Stats.Latency = groupStats([{data.Correct} {data.Incorrect}],[],'inAxis',true);

%% Number of spontaneous explorations
data.Correct = Summary.numberSamplesCorrect ;
data.Incorrect = Summary.numberSamplesIncorrect;

subplot(numrows, numcol, 3)
nhist(data,'proportion','samebins')
subplot(numrows, numcol, 4)
Stats.numSample = groupStats([{data.Correct} {data.Incorrect}],[],'inAxis',true);

%% Spike responses of cells
load('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\LickPSTHSummary.mat')
timeaxis = linspace(-2,2,41);
idxT = timeaxis<=0.3 & timeaxis>=-0.3;
spec=cbrewer('seq', 'Blues', 20);
spec(spec>1) = 1;
spec(spec<0) = 0;

selectedPSTH = Summary.psthReward.lickTypes{6};
[maxRate,maxRateIdx] = max(selectedPSTH,[],2);
[~,idxmax] = sort(maxRateIdx,'ascend'); 


ax1 = subplot(numrows, numcol, 5);
selectedPSTH = Summary.psthReward.lickTypes{7};
norm = zscore(selectedPSTH,[],2);
imagesc(timeaxis,1:length(idxmax),norm(idxmax,:))
colormap(ax1,spec)
caxis([-1 3])
colorbar
xlim([-2 2])
ylabel(strcat('cellID',num2str(length(idxmax))))
xticks([])
a = selectedPSTH(:,idxT);
avgRate(:,1) = mean(a,2);%
plotAvgStd(norm,numrows,numcol,7,fig2,timeaxis','m',1)
ylim([-1 1.5])

ax1 = subplot(numrows, numcol, 6);
selectedPSTH = Summary.psthReward.lickTypes{16};
% [maxRate,maxRateIdx] = max(selectedPSTH,[],2);
% [~,idxmax] = sort(maxRateIdx,'ascend'); 
norm = zscore(selectedPSTH,[],2);
imagesc(timeaxis,1:length(idxmax),norm(idxmax,:))
colormap(ax1,spec)
caxis([-1 3])
colorbar
xlim([-2 2])
ylabel(strcat('cellID',num2str(length(idxmax))))
xticks([])
a = selectedPSTH(:,idxT);
avgRate(:,2) = mean(a,2);%
plotAvgStd(norm,numrows,numcol,7,fig2,timeaxis','k',1)
ylim([-1 1.5])

subplot(numrows,numcol,8)
Stats.firstVssampleLicks = groupStats([{avgRate(:,1)},{avgRate(:,2)}],[],'inAxis',true,'repeatedMeasures',true);

%% Save figure
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
saveas(gcf,strcat(expPath,'SupFigure9_explorationLicks.png'));
saveas(gcf,strcat(expPath,'SupFigure9_explorationLicks.eps'),'epsc');
saveas(gcf,strcat(expPath,'SupFigure9_explorationLicks.fig'));
save(strcat(expPath,'SupFigure9_explorationLicks.mat'),'Stats'); 

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
        % Bootstrapping to estimate variability of the median
        n_bootstraps = 1000;
        bootstrap_medians = zeros(n_bootstraps, size(array, 2));
        
        for i = 1:n_bootstraps
            resample_indices = randi([1, size(array, 1)], size(array, 1), 1);  % Generate random indices with replacement
            bootstrap_sample = array(resample_indices, :);  % Create bootstrap sample
            bootstrap_medians(i, :) = nanmedian(bootstrap_sample);  % Compute median of the bootstrap sample
        end
        
        % Compute the 2.5th and 97.5th percentiles for the bounds
        lArr = prctile(bootstrap_medians, 2.5);
        uArr = prctile(bootstrap_medians, 97.5);

    end

    fill([xAxis; flipud(xAxis)],[lArr'; flipud(uArr')],col,'linestyle','none','FaceAlpha',0.5);                    
    hold on
    hi = line(xAxis,meanpsth,'LineWidth',1,'Color',col);
   % yscale log

end

