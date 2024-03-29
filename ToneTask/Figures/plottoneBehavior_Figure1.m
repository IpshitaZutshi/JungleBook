function plottoneBehavior_Figure1

figure
set(gcf,'Renderer','painters')
set(gcf,'Color','w')
set(gcf,'Position',[680 42 975 962]);

%% Panel A: Task schematic - ACgN [1,2]
subplot(4,4,1)
title('A. ACgN task schematic','FontName', 'Arial','FontSize', 9)
axis off

%% Panel B: ACgN Task trajectory example [3], C, ACgN licks heatmap[4]
subplot(4,4,2)
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ44\Final\IZ44_220913_sess11';
cd(sessloc)
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*Tracking.Behavior.mat']);
load(file.name);
plot(tracking.position.y,tracking.timestamps, 'Color',[0.5 0.5 0.5])
hold on
for tt  = 1:size(behavTrials.timestamps,1)
    [~,idx] = min(abs(tracking.timestamps-behavTrials.timestamps(tt,1)));
    scatter(tracking.position.y(idx),tracking.timestamps(idx),5,'k','filled')
    
    [~,idx] = min(abs(tracking.timestamps-behavTrials.timestamps(tt,2)));
    scatter(tracking.position.y(idx),tracking.timestamps(idx),5,'r','filled')
end
ylim([behavTrials.timestamps(1,1) behavTrials.timestamps(end,2)]);
set(gca,'YDir','reverse')
xlim([0 125])
box off

subplot(4,4,3)
colMap = cbrewer('seq','Greys',3);
imagesc(behavTrials.numLicks(:,2:7))
colormap(colMap)
caxis([0 1])
xlim([0.5 7.5])
ylabel(strcat('Trials (',num2str(1),'-',num2str(size(behavTrials.timestamps,1)),')'))
%xticklabels({'Lick1','Lick2','Lick3','Lick4','Lick5','Lick6'})
h=gca; h.XAxis.TickLength = [0 0];
hold on
%colMap = cbrewer('seq','Blues',18);
col = [119/255 221/255 229/255;...
    122/255 201/255 229/255;...
    38/255 169/255 224/255;...
    73/255 136/255 189/255;...
    17/255 55/255 174/255;...
    0/255 0/255 134/255];
scatter(ones(1,sum(behavTrials.correct))*6.5,find(behavTrials.correct),1,[70/243 148/243 73/243],'filled');
scatter(ones(1,sum(behavTrials.correct==0&behavTrials.linTrial==0))*6.75,find(behavTrials.correct==0&behavTrials.linTrial==0),1,[175/243 54/243 60/243],'filled');
scatter(ones(1,sum(behavTrials.linTrial==1))*6.9,find(behavTrials.linTrial==1),1,[50/243 50/243 50/243],'filled');
scatter(ones(1,sum(behavTrials.toneGain==0&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==0&behavTrials.linTrial==0),1,col(1,:),'filled');
scatter(ones(1,sum(behavTrials.toneGain==1&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==1&behavTrials.linTrial==0),1,col(2,:),'filled');
scatter(ones(1,sum(behavTrials.toneGain==2&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==2&behavTrials.linTrial==0),1,col(3,:),'filled');
scatter(ones(1,sum(behavTrials.toneGain==3&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==3&behavTrials.linTrial==0),1,col(4,:),'filled');
scatter(ones(1,sum(behavTrials.toneGain==4&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==4&behavTrials.linTrial==0),1,col(5,:),'filled');
scatter(ones(1,sum(behavTrials.toneGain==5&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==5&behavTrials.linTrial==0),1,col(6,:),'filled');
box off

%% Panel D: Performance across mice [5]
Summary = plotToneBehaviorStatistics;

subplot(4,4,4)
for ii = 1:6
    data = Summary.performance(Summary.mouseID==ii);
    scatter(ones*ii, data,10,'MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',.5)
    hold on
    line([ii-0.3 ii+0.3],[mean(data) mean(data)], 'Color','b','LineWidth',4)
end
xlim([0 7])
ylim([0.4 1])
ylabel('Performance')
xticks([1:6])
xticklabels(["39", "40", "43","44", "47", "48"])
xtickangle(0)

%% Panel E: Performance by port [7]% Trial probability [6]
subplot(4,4,5)
% colMap = cbrewer('seq','Blues',18);
% col = [colMap(5,:);colMap(8,:);colMap(10,:);colMap(13,:);colMap(16,:);colMap(18,:)];
col = [119/255 221/255 229/255;...
    122/255 201/255 229/255;...
    38/255 169/255 224/255;...
    73/255 136/255 189/255;...
    17/255 55/255 174/255;...
    0/255 0/255 134/255];

% Fig1Stats.trialType = groupStats([{Summary.trialType(:,1)},{Summary.trialType(:,2)},{Summary.trialType(:,3)},{Summary.trialType(:,4)},{Summary.trialType(:,5)},{Summary.trialType(:,6)}],[],'repeatedMeasures',true,'inAxis',true,'plotType','boxplot','color',col,'labelSummary',false);
% title('Trial probability','FontName', 'Arial','FontSize', 9)
% for ii = 1:6
%     hold on
%     scatter((ones*ii)-0.4, Summary.trialType(:,ii),10,'MarkerFaceColor',col(ii,:),'MarkerEdgeColor',col(ii,:))    
% end
% xticks([1:6])
% xticklabels(["1", "2", "3","4", "5", "6"])
Fig1Stats.performance = groupStats([{Summary.portCorrect(:,1)},{Summary.portCorrect(:,2)},{Summary.portCorrect(:,3)},{Summary.portCorrect(:,4)},{Summary.portCorrect(:,5)},{Summary.portCorrect(:,6)}],[],'repeatedMeasures',true,'inAxis',true,'plotType','boxplot','color',col,'labelSummary',false);
for ii = 1:6
    hold on
    scatter((ones*ii)-0.4, Summary.portCorrect(:,ii),10,'MarkerFaceColor',col(ii,:),'MarkerEdgeColor','none')    
end
ylim([0.4 1.05])
xticks([1:6])
xticklabels(["1", "2", "3","4", "5", "6"])
title('Performance by trial type','FontName', 'Arial','FontSize', 9)

%% Panel F: False alarm rate [7]
subplot(4,4,6)
Fig1Stats.falseAlarm = groupStats([{Summary.falseAlarm(:,1)},{Summary.falseAlarm(:,2)},{Summary.falseAlarm(:,3)},{Summary.falseAlarm(:,4)},{Summary.falseAlarm(:,5)},{Summary.falseAlarm(:,6)}],[],'repeatedMeasures',true,'inAxis',true,'plotType','boxplot','color',col,'labelSummary',false);
for ii = 1:6
    hold on
    scatter((ones*ii)-0.4, Summary.falseAlarm(:,ii),10,'MarkerFaceColor',col(ii,:),'MarkerEdgeColor',col(ii,:))    
end
%ylim([0.4 1.2])
xticks([1:6])
xticklabels(["1", "2", "3","4", "5", "6"])
title('False alarm rate','FontName', 'Arial','FontSize', 9)

%% Panel F: D prime 
subplot(4,4,7)
for kk = 1:6
    dPrime(kk)  = norminv(nanmean(Summary.portCorrect(:,kk)))-norminv(nanmean(Summary.falseAlarm(:,kk)));
end
b = bar(dPrime, 'FaceColor','flat','EdgeColor','none');

for k = 1:6
    b.CData(k,:) = col(k,:);
end
ylabel('D prime')
xticks([1:6])
xticklabels(["1", "2", "3","4", "5", "6"])
title('D prime','FontName', 'Arial','FontSize', 9)
box off

%% Panel G: Distance from target [8]
subplot(4,4,8)
histogram(-Summary.PortDiffDist,'Normalization','probability','FaceColor',[0.5 0.5 0.5],'EdgeColor','none')
title('Distance from target (incorrect)','FontName', 'Arial','FontSize', 9)
ylabel('Probability')
xlim([-5 5])
%ylim([0 0.7])
box off

%% Panel H: Schematic probe trials [9]
subplot(4,4,9)
Summary = extractProbeTrialSessions;
boxplot(Summary.fractProbe,'BoxStyle','filled','Colors',[0.5 0.5 0.5],'Symbol','')
title('H. Probe trials schematic','FontName', 'Arial','FontSize', 9)
box off

%% Panel I: Performance probe trials [10]
subplot(4,4,10)
freq = [4 8 16];
data = [];
for ii = 1:length(freq)
    % Stack the results from all animals together for each frequency    
    data{ii} = Summary.perf(Summary.freq==freq(ii),:);
    scatter(ones*(2*(ii-1)+1),data{ii}(:,1),10,[0.5 0.5 0.5],'filled')
    hold on
    scatter(ones*(2*(ii-1)+2),data{ii}(:,2),10,[96/243 60/243 108/243],'filled')
    plot([ones*(2*(ii-1)+1) ones*(2*(ii-1)+2)],[data{ii}(:,1) data{ii}(:,2)],'Color',[0.5 0.5 0.5])
    line([2*(ii-1)+1-0.3 2*(ii-1)+1+0.3],[median(data{ii}(:,1)) median(data{ii}(:,1))],'Color',[0.5 0.5 0.5],'LineWidth',4)
    line([2*(ii-1)+2-0.3 2*(ii-1)+2+0.3],[median(data{ii}(:,2)) median(data{ii}(:,2))],'Color',[96/243 60/243 108/243],'LineWidth',4)
end
ylabel('Performance')
box off
xticklabels([" ","4 kHz", "8 kHz", "16 kHz"])
xtickangle(0)

Fig1Stats.Probe4 = groupStats({data{1}(:,1) data{1}(:,2)},[],'repeatedMeasures',true,'doPlot',false);
Fig1Stats.Probe8 = groupStats({data{2}(:,1) data{2}(:,2)},[],'repeatedMeasures',true,'doPlot',false);
Fig1Stats.Probe16 = groupStats({data{3}(:,1) data{3}(:,2)},[],'repeatedMeasures',true,'doPlot',false);

%% Panel J: Distance to stop, 4 kHz [11]
subplot(4,4,11)
histogram(-Summary.PortDiffDist4{1},'Normalization','probability','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeColor','none')
hold on
histogram(-Summary.PortDiffDist4{2},'Normalization','probability','FaceColor',[96/243 60/243 108/243],'FaceAlpha',0.5,'EdgeColor','none')
xlim([-5 5])
ylabel('Probability')
title('4 Hz Probe trials','FontName', 'Arial','FontSize', 9)
Fig1Stats.PortDiffDist4 =  groupStats({Summary.PortDiffDist4{1} Summary.PortDiffDist4{2}},[],'doPlot',false);
box off
Fig1Stats.PortDiffDist4 =  groupStats({Summary.PortDiffDist4{1} Summary.PortDiffDist4{2}},[],'doPlot',false);

%% Panel K: Distance to stop, 8 kHz [12]
subplot(4,4,12)
histogram(-Summary.PortDiffDist8{1},'Normalization','probability','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeColor','none')
hold on
histogram(-Summary.PortDiffDist8{2},'Normalization','probability','FaceColor',[96/243 60/243 108/243],'FaceAlpha',0.5,'EdgeColor','none')
box off
xlim([-5 5])
ylabel('Probability')
title('8 Hz Probe trials','FontName', 'Arial','FontSize', 9)
box off
Fig1Stats.PortDiffDist8 =  groupStats({Summary.PortDiffDist8{1} Summary.PortDiffDist8{2}},[],'doPlot',false);

%% Panel L: Task schematic  - Control [13 14]
subplot(4,4,[13 14])
title('L. Control task schematic','FontName', 'Arial','FontSize', 9)
axis off

%% Panel M: Control Task trajectory example [15], N, Control licks heatmap[16]

subplot(4,4,15)
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ46\IZ46_230410_sess11';
cd(sessloc)
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*Tracking.Behavior.mat']);
load(file.name);
plot(tracking.position.y,tracking.timestamps, 'Color',[0.5 0.5 0.5])
hold on
for tt  = 1:size(behavTrials.timestamps,1)
    [~,idx] = min(abs(tracking.timestamps-behavTrials.timestamps(tt,1)));
    scatter(tracking.position.y(idx),tracking.timestamps(idx),5,'k','filled')
    
    [~,idx] = min(abs(tracking.timestamps-behavTrials.timestamps(tt,2)));
    scatter(tracking.position.y(idx),tracking.timestamps(idx),5,'r','filled')
end
ylim([behavTrials.timestamps(1,1) behavTrials.timestamps(end,2)]);
set(gca,'YDir','reverse')
xlim([0 125])
box off

subplot(4,4,16)
colMap = cbrewer('seq','Greys',3);
imagesc(behavTrials.numLicks(:,2:7))
colormap(colMap)
caxis([0 1])
xlim([0.5 7.5])
ylabel(strcat('Trials (',num2str(1),'-',num2str(size(behavTrials.timestamps,1)),')'))
%xticklabels({'Lick1','Lick2','Lick3','Lick4','Lick5','Lick6'})
h=gca; h.XAxis.TickLength = [0 0];
hold on
%colMap = cbrewer('seq','Blues',18);
scatter(ones(1,sum(behavTrials.linTrial==1))*6.9,find(behavTrials.linTrial==1),1,[50/243 50/243 50/243],'filled');
scatter(ones(1,sum(behavTrials.toneGain==0&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==0&behavTrials.linTrial==0),1,col(1,:),'filled');
scatter(ones(1,sum(behavTrials.toneGain==1&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==1&behavTrials.linTrial==0),1,col(2,:),'filled');
scatter(ones(1,sum(behavTrials.toneGain==2&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==2&behavTrials.linTrial==0),1,col(3,:),'filled');
scatter(ones(1,sum(behavTrials.toneGain==3&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==3&behavTrials.linTrial==0),1,col(4,:),'filled');
scatter(ones(1,sum(behavTrials.toneGain==4&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==4&behavTrials.linTrial==0),1,col(5,:),'filled');
scatter(ones(1,sum(behavTrials.toneGain==5&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==5&behavTrials.linTrial==0),1,col(6,:),'filled');
box off

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';
saveas(gcf,strcat(expPath,'\Compiled\Figures\Figure1_behaviorSummary.png'));
saveas(gcf,strcat(expPath,'\Compiled\Figures\Figure1_behaviorSummary.eps'),'epsc');
saveas(gcf,strcat(expPath,'\Compiled\Figures\Figure1_behaviorSummary.fig'));
save(strcat(expPath,'\Compiled\Figures\Figure1_behaviorStats.mat'),'Fig1Stats'); 
end