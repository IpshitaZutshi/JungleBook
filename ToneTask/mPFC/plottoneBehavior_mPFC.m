function plottoneBehavior_mPFC

fig2  = figure;
set(gcf,'Renderer','painters')
set(gcf,'Color','w')
set(gcf,'Position',[200 138 1675 745]);
numrows = 2;
numcol = 4;

%% Panel C: ACgN Task trajectory example [3], D, ACgN licks heatmap[4]
subplot(numrows,numcol,[1 2] ,'Parent', fig2)
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ44\Final\IZ44_220920_sess15';
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

subplot(numrows,numcol,3,'Parent', fig2)
colMap = [235/243 235/243 235/243;...
    0/243 0/243 0/243;...
    150/243 150/243 150/243;...
    0.2 0.6 0.2;...
    1 0 0.8];
lickHM = behavTrials.numLicks(:,2:7);
lickHM(lickHM>=1) = 2; %First make all licks the same

for ii = 1:length(behavTrials.linTrial)
    if behavTrials.linTrial(ii) == 0 
        lickLoc = behavTrials.lickLoc(ii)+1;
        if behavTrials.correct(ii) == 1
            lickHM(ii,lickLoc) = 3;
        else
            lickHM(ii,lickLoc) = 4;
        end
    else
        lickHM(ii,6) = 1;
    end
end
% 
imagesc(lickHM)
colormap(colMap)
caxis([0 4])
xlim([0.5 8.5])
ylabel(strcat('Trials (',num2str(1),'-',num2str(size(behavTrials.timestamps,1)),')'))
%xticklabels({'Lick1','Lick2','Lick3','Lick4','Lick5','Lick6'})
h=gca; h.XAxis.TickLength = [0 0];
hold on

col = [238/255 67/255 69/255;...
    241/255 114/255 42/255;...
    247/255 149/255 33/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];

% col = [119/255 221/255 229/255;...
%     122/255 201/255 229/255;...
%     38/255 169/255 224/255;...
%     73/255 136/255 189/255;...
%     17/255 55/255 174/255;...
%     0/255 0/255 134/255];
scatter(ones(1,sum(behavTrials.toneGain==0&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==0&behavTrials.linTrial==0),1,col(1,:),'filled');
scatter(ones(1,sum(behavTrials.toneGain==1&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==1&behavTrials.linTrial==0),1,col(2,:),'filled');
scatter(ones(1,sum(behavTrials.toneGain==2&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==2&behavTrials.linTrial==0),1,col(3,:),'filled');
scatter(ones(1,sum(behavTrials.toneGain==3&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==3&behavTrials.linTrial==0),1,col(4,:),'filled');
scatter(ones(1,sum(behavTrials.toneGain==4&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==4&behavTrials.linTrial==0),1,col(5,:),'filled');
scatter(ones(1,sum(behavTrials.toneGain==5&behavTrials.linTrial==0))*7.1,find(behavTrials.toneGain==5&behavTrials.linTrial==0),1,col(6,:),'filled');
scatter(ones(1,sum(behavTrials.stim==1))*8.1,find(behavTrials.stim==1),1,'b','filled');
box off

subplot(numrows,numcol, 4)

plot([25 45 65 85])

%% Panel E: Performance across mice [5]
Summary = plotToneBehaviorStatistics;

subplot(numrows,numcol,5,'Parent', fig2)
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

%% Panel F: Performance by port [7]% Trial probability [6]
subplot(numrows,numcol,6,'Parent', fig2)
% col = [119/255 221/255 229/255;...
%     122/255 201/255 229/255;...
%     38/255 169/255 224/255;...
%     73/255 136/255 189/255;...
%     17/255 55/255 174/255;...
%     0/255 0/255 134/255];

col = [238/255 67/255 69/255;...
    241/255 114/255 42/255;...
    247/255 149/255 33/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];

Fig1Stats.performance = groupStats([{Summary.portCorrect(:,1)},{Summary.portCorrect(:,2)},{Summary.portCorrect(:,3)},{Summary.portCorrect(:,4)},{Summary.portCorrect(:,5)},{Summary.portCorrect(:,6)}],[],'repeatedMeasures',true,'inAxis',true,'plotType','boxplot','color',col,'labelSummary',false);
for ii = 1:6
    hold on
    scatter((ones*ii)-0.4, Summary.portCorrect(:,ii),10,'MarkerFaceColor',col(ii,:),'MarkerEdgeColor','none')    
end
ylim([0.4 1.05])
xticks([1:6])
xticklabels(["1", "2", "3","4", "5", "6"])
title('Performance by trial type','FontName', 'Arial','FontSize', 9)

%% Panel G: False alarm rate [7]
subplot(numrows,numcol,7,'Parent', fig2)
Fig1Stats.falseAlarm = groupStats([{Summary.falseAlarm(:,1)},{Summary.falseAlarm(:,2)},{Summary.falseAlarm(:,3)},{Summary.falseAlarm(:,4)},{Summary.falseAlarm(:,5)},{Summary.falseAlarm(:,6)}],[],'repeatedMeasures',true,'inAxis',true,'plotType','boxplot','color',col,'labelSummary',false);
for ii = 1:6
    hold on
    scatter((ones*ii)-0.4, Summary.falseAlarm(:,ii),10,'MarkerFaceColor',col(ii,:),'MarkerEdgeColor',col(ii,:))    
end
%ylim([0.4 1.2])
xticks([1:6])
xticklabels(["1", "2", "3","4", "5", "6"])
title('False alarm rate','FontName', 'Arial','FontSize', 9)

%% Panel H: Distance from target [8]
subplot(numrows,numcol,8,'Parent', fig2)
histogram(-Summary.PortDiffDist,'Normalization','probability','FaceColor',[0.5 0.5 0.5],'EdgeColor','none')
title('Distance from target (incorrect)','FontName', 'Arial','FontSize', 9)
ylabel('Probability')
xlim([-5 5])
%ylim([0 0.7])
box off

%% Panel I - Time frequency for correct trials
plotFrequencyoverTime('fighandle',fig2,'numrows',numrows,'numcol',numcol,'rowloc',2,'colloc',4','correct',true)
ylim([0 25000])
box off

%% Panel J - Time frequency for incorrect trials
plotFrequencyoverTime('fighandle',fig2,'numrows',numrows,'numcol',numcol,'rowloc',2,'colloc',5','correct',false)
ylim([0 40000])
box off

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\MainFigures\';
saveas(gcf,strcat(expPath,'Figure1_behaviorSummary.png'));
saveas(gcf,strcat(expPath,'Figure1_behaviorSummary.eps'),'epsc');
saveas(gcf,strcat(expPath,'Figure1_behaviorSummary.fig'));
save(strcat(expPath,'Figure1_behaviorStats.mat'),'Fig1Stats'); 

end