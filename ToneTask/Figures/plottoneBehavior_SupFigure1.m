function plottoneBehavior_SupFigure1

figure
set(gcf,'Renderer','painters')
set(gcf,'Color','w')
set(gcf,'Position',[680 42 975 962]);

Summary = plotToneBehaviorStatistics;

%% Panel A: Trial probability [6]
subplot(4,4,1)
colMap = cbrewer('seq','Blues',18);
col = [colMap(5,:);colMap(8,:);colMap(10,:);colMap(13,:);colMap(16,:);colMap(18,:)];

SupFig1Stats.trialType = groupStats([{Summary.trialType(:,1)},{Summary.trialType(:,2)},{Summary.trialType(:,3)},{Summary.trialType(:,4)},{Summary.trialType(:,5)},{Summary.trialType(:,6)}],[],'repeatedMeasures',true,'inAxis',true,'plotType','boxplot','color',col,'labelSummary',false);
title('Trial probability','FontName', 'Arial','FontSize', 9)
for ii = 1:6
    hold on
    scatter((ones*ii)-0.4, Summary.trialType(:,ii),10,'MarkerFaceColor',col(ii,:),'MarkerEdgeColor',col(ii,:))    
end
xticks([1:6])
xticklabels(["1", "2", "3","4", "5", "6"])


%% Panel B: False alarm rate [2]
subplot(4,4,2)
SupFig1Stats.falseAlarm = groupStats([{Summary.falseAlarm(:,1)},{Summary.falseAlarm(:,2)},{Summary.falseAlarm(:,3)},{Summary.falseAlarm(:,4)},{Summary.falseAlarm(:,5)},{Summary.falseAlarm(:,6)}],[],'repeatedMeasures',true,'inAxis',true,'plotType','boxplot','color',col,'labelSummary',false);
for ii = 1:6
    hold on
    scatter((ones*ii)-0.4, Summary.falseAlarm(:,ii),10,'MarkerFaceColor',col(ii,:),'MarkerEdgeColor',col(ii,:))    
end
xticks([1:6])
xticklabels(["1", "2", "3","4", "5", "6"])
title('False alarm rate','FontName', 'Arial','FontSize', 9)

%% False alarm per mouse
for mm = 1:4
    subplot(4,4,2+mm)
    mouseidx = Summary.mouseID==mm;
    SupFig1Stats.falseAlarmMouse{mm} = groupStats([{Summary.falseAlarm(mouseidx,1)},{Summary.falseAlarm(mouseidx,2)},{Summary.falseAlarm(mouseidx,3)},{Summary.falseAlarm(mouseidx,4)},{Summary.falseAlarm(mouseidx,5)},{Summary.falseAlarm(mouseidx,6)}],[],'repeatedMeasures',true,'inAxis',true,'plotType','boxplot','color',col,'labelSummary',false);
    for ii = 1:6
        hold on
        scatter((ones*ii)-0.4, Summary.falseAlarm(mouseidx,ii),10,'MarkerFaceColor',col(ii,:),'MarkerEdgeColor',col(ii,:))   
    end
    ylim([0 0.2])
end


expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';
saveas(gcf,strcat(expPath,'\Compiled\Figures\SupFigure1_behaviorSummary.png'));
saveas(gcf,strcat(expPath,'\Compiled\Figures\SupFigure1_behaviorSummary.eps'),'epsc');
saveas(gcf,strcat(expPath,'\Compiled\Figures\SupFigure1_behaviorSummary.fig'));
save(strcat(expPath,'\Compiled\Figures\SupFigure1_behaviorStats.mat'),'SupFig1Stats'); 
end