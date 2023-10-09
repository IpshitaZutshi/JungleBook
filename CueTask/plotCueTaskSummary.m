function plotCueTaskSummary

filename = 'C:\Users\ipshi\Dropbox (NYU Langone Health)\Cue Task\OptoSummary';

% Assume sheet 1 has hippocampal data
dataHip = xlsread(filename,'Hip');

% Assume sheet 2 has mPFC data
datamPFC = xlsread(filename,'mPFC');

%% Order of the data
labels = {'Cue','Delay','Cue+Delay','Choice (2s)', 'Choice (after 1.6s)','Delay+Choice', 'Choice (after 800ms)','Longer delay'};

% Skipping choice (after 800ms), we have 7 conditions
figure
set(gcf,'Renderer','painters')
set(gcf,'Position',[1921 41 1920 970])
set(gcf,'Color','w')

col = [0.5 0.5 0.5;
    8/243 133/243 161/243];

for ii = 1:8
    subplot(2,10,ii)
    data{1} = dataHip(:,(2*(ii-1)+1));
    data{2} = dataHip(:,(2*(ii-1)+2));
    stats.Hip{ii} = groupStats(data,[],'sigStarTest','anova','plotType','BoxLinesSEM','inAxis',true,'labelSummary',false,'repeatedMeasures',true,'Color',col);
    ylim([0 100])
    line([0.5 2.5],[50 50],'Color',[175/243 54/243 60/243],'LineStyle','--')
    title(labels{ii})
    xticks([])
    
%     subplot(2,8,ii+8)
%     data{1} = datamPFC(:,(2*(ii-1)+1));
%     data{2} = datamPFC(:,(2*(ii-1)+2));
%     stats.mPFC{ii} = groupStats(data,[],'sigStarTest','anova','plotType','BoxLinesSEM','inAxis',true,'labelSummary',false,'repeatedMeasures',true,'Color',col);
%     ylim([0 100])
%     line([0.5 2.5],[50 50],'Color',[175/243 54/243 60/243],'LineStyle','--')
%     title(labels{ii})
end

%% Plot the averages
avgHip = nanmean(dataHip,1);
avgHip1 = avgHip(1:2:end);
avgHip2 = avgHip(2:2:end);
subplot(2,9,9)
plot([avgHip1; avgHip2],'LineWidth',1.5)
box off
ylim([0 100])
legend(labels,'Location','southeast','NumColumns',2)
xticks([])

avgmPFC= nanmean(datamPFC,1);
avgmPFC1 = avgmPFC(1:2:end);
avgmPFC2 = avgmPFC(2:2:end);
subplot(2,9,18)
plot([avgmPFC1; avgmPFC2],'LineWidth',1.5)
box off
ylim([0 100])
legend(labels,'Location','southeast','NumColumns',2)
end