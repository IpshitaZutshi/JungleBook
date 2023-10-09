function plottoneBehavior_SupFigure1

figure
set(gcf,'Renderer','painters')
set(gcf,'Color','w')
set(gcf,'Position',[2145 244 1085 600]);

numrows = 2;
numcol = 3;

Summary = plotToneBehaviorStatistics;

%% Panel A: Trial probability 
subplot(numrows,numcol,1)
col = [119/255 221/255 229/255;...
    122/255 201/255 229/255;...
    38/255 169/255 224/255;...
    73/255 136/255 189/255;...
    17/255 55/255 174/255;...
    0/255 0/255 134/255];

SupFig1Stats.trialType = groupStats([{Summary.trialType(:,1)},{Summary.trialType(:,2)},{Summary.trialType(:,3)},{Summary.trialType(:,4)},{Summary.trialType(:,5)},{Summary.trialType(:,6)}],[],'repeatedMeasures',true,'inAxis',true,'plotType','boxplot','color',col,'labelSummary',false);
title('Trial probability','FontName', 'Arial','FontSize', 9)
for ii = 1:6
    hold on
    scatter((ones*ii)-0.4, Summary.trialType(:,ii),10,'MarkerFaceColor',col(ii,:),'MarkerEdgeColor',col(ii,:))    
end
xticks([1:6])
xticklabels(["1", "2", "3","4", "5", "6"])


%% Panel B: D prime 
subplot(numrows,numcol,2)
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

%% Panel C: Schematic probe trials 
subplot(numrows,numcol,3)
Summary = extractProbeTrialSessions;
boxplot(Summary.fractProbe,'BoxStyle','filled','Colors',[0.5 0.5 0.5],'Symbol','')
title('Probe trials schematic','FontName', 'Arial','FontSize', 9)
box off

%% Panel D: Performance probe trials 
subplot(numrows,numcol,4)
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

SupFig1Stats.Probe4 = groupStats({data{1}(:,1) data{1}(:,2)},[],'repeatedMeasures',true,'doPlot',false);
SupFig1Stats.Probe8 = groupStats({data{2}(:,1) data{2}(:,2)},[],'repeatedMeasures',true,'doPlot',false);
SupFig1Stats.Probe16 = groupStats({data{3}(:,1) data{3}(:,2)},[],'repeatedMeasures',true,'doPlot',false);

%% Panel E: Distance to stop, 4 kHz 
subplot(numrows,numcol,5)
histogram(-Summary.PortDiffDist4{1},'Normalization','probability','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeColor','none')
hold on
histogram(-Summary.PortDiffDist4{2},'Normalization','probability','FaceColor',[96/243 60/243 108/243],'FaceAlpha',0.5,'EdgeColor','none')
xlim([-5 5])
ylabel('Probability')
title('4 Hz Probe trials','FontName', 'Arial','FontSize', 9)
SupFig1Stats.PortDiffDist4 =  groupStats({Summary.PortDiffDist4{1} Summary.PortDiffDist4{2}},[],'doPlot',false);
box off

%% Panel F: Distance to stop, 8 kHz
subplot(numrows,numcol,6)
histogram(-Summary.PortDiffDist8{1},'Normalization','probability','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeColor','none')
hold on
histogram(-Summary.PortDiffDist8{2},'Normalization','probability','FaceColor',[96/243 60/243 108/243],'FaceAlpha',0.5,'EdgeColor','none')
box off
xlim([-5 5])
ylabel('Probability')
title('8 Hz Probe trials','FontName', 'Arial','FontSize', 9)
box off
SupFig1Stats.PortDiffDist8 =  groupStats({Summary.PortDiffDist8{1} Summary.PortDiffDist8{2}},[],'doPlot',false);


expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure1_behaviorSummary.png'));
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure1_behaviorSummary.eps'),'epsc');
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure1_behaviorSummary.fig'));
save(strcat(expPath,'\Compiled\Figures_Sep23\SupFigure1_behaviorStats.mat'),'SupFig1Stats'); 
end