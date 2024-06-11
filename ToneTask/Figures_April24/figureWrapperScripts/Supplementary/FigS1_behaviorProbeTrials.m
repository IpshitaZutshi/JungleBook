function FigS1_behaviorProbeTrials

figure
set(gcf,'Renderer','painters')
set(gcf,'Color','w')
set(gcf,'Position',[2145 244 1085 600]);

numrows = 3;
numcol = 3;

Summary = plotToneBehaviorStatistics;

%% Panel A: Trial probability 
subplot(numrows,numcol,1)
col = [238/255 67/255 69/255;...
    241/255 114/255 42/255;...
    247/255 149/255 33/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];

SupFig1Stats.trialType = groupStats([{Summary.trialType(:,1)},{Summary.trialType(:,2)},{Summary.trialType(:,3)},{Summary.trialType(:,4)},{Summary.trialType(:,5)},{Summary.trialType(:,6)}],[],'repeatedMeasures',true,'inAxis',true,'plotType','boxplot','color',col,'labelSummary',false);
title('Trial probability','FontName', 'Arial','FontSize', 9)
for ii = 1:6
    hold on
    scatter((ones*ii)-0.4, Summary.trialType(:,ii),10,'MarkerFaceColor',col(ii,:),'MarkerEdgeColor',col(ii,:))    
end
xticks([1:6])
xticklabels(["1", "2", "3","4", "5", "6"])


%% Panel B: Schematic probe trials 
subplot(numrows,numcol,2)
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
%histogram(-Summary.PortDiffDist4{1},'Normalization','probability','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeColor','none')
%hold on
%histogram(-Summary.PortDiffDist4{2},'Normalization','probability','FaceColor',[96/243 60/243 108/243],'FaceAlpha',0.5,'EdgeColor','none')
edges = [-4.5 -3.5 -2.5 -1.5 -0.5 0.5 1.5 2.5 3.5];
[h(:,1)] = histcounts(-Summary.PortDiffDist4{1},edges,'Normalization','probability');
[h(:,2)] = histcounts(-Summary.PortDiffDist4{2},edges,'Normalization','probability');
M = movmean(edges,2);
b = bar(M(2:end),h,'FaceColor','flat');
b(1).BarWidth = 1.25;
b(1).CData = [0.5 0.5 0.5];
b(2).CData = [96/243 60/243 108/243];
xlim([-5 3])
ylabel('Probability')
title('4 Hz Probe trials','FontName', 'Arial','FontSize', 9)
SupFig1Stats.PortDiffDist4 =  groupStats({Summary.PortDiffDist4{1} Summary.PortDiffDist4{2}},[],'doPlot',false);
box off

%% Panel F: Distance to stop, 8 kHz
subplot(numrows,numcol,6)
% histogram(-Summary.PortDiffDist8{1},'Normalization','probability','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeColor','none')
% hold on
% histogram(-Summary.PortDiffDist8{2},'Normalization','probability','FaceColor',[96/243 60/243 108/243],'FaceAlpha',0.5,'EdgeColor','none')
% box off
% xlim([-5 5])
edges = [-4.5 -3.5 -2.5 -1.5 -0.5 0.5 1.5 2.5 3.5];
[h(:,1)] = histcounts(-Summary.PortDiffDist8{1},edges,'Normalization','probability');
[h(:,2)] = histcounts(-Summary.PortDiffDist8{2},edges,'Normalization','probability');
M = movmean(edges,2);
b = bar(M(2:end),h,'FaceColor','flat');
b(1).BarWidth = 1.25;
b(1).CData = [0.5 0.5 0.5];
b(2).CData = [96/243 60/243 108/243];
xlim([-5 3])
ylabel('Probability')
title('8 Hz Probe trials','FontName', 'Arial','FontSize', 9)
box off
SupFig1Stats.PortDiffDist8 =  groupStats({Summary.PortDiffDist8{1} Summary.PortDiffDist8{2}},[],'doPlot',false);

%% Performance in probe trials split up by target port

freqs = [4 8];
for ff = 1:length(freqs)
    subplot(numrows,numcol,6+ff)
    for ii = 1:6
        data = Summary.perfTrial(Summary.freq==freqs(ff),[ii ii+6]);
        %scatter(ones*(2*(ii-1)+1),data(:,1),10,[0.5 0.5 0.5],'filled')
        hold on
        %scatter(ones*(2*(ii-1)+2),data(:,2),10,[96/243 60/243 108/243],'filled')
        plot([ones*(2*(ii-1)+1) ones*(2*(ii-1)+2)],[data(:,1) data(:,2)],'Color',[0.5 0.5 0.5])
        line([2*(ii-1)+1-0.3 2*(ii-1)+1+0.3],[nanmedian(data(:,1)) nanmedian(data(:,1))],'Color',[0.5 0.5 0.5],'LineWidth',4)
        line([2*(ii-1)+2-0.3 2*(ii-1)+2+0.3],[nanmedian(data(:,2)) nanmedian(data(:,2))],'Color',[96/243 60/243 108/243],'LineWidth',4)
        ylim([0 1.2])

        x = data(:,1);
        y = data(:,2);
        str = strcat('ProbeTrials',num2str(freqs(ff)));
        [SupFig1Stats.(str){ii}.h, SupFig1Stats.(str){ii}.p,SupFig1Stats.(str){ii}.stats] = ttest(x(~isnan(x)&~isnan(y)), y(~isnan(x)&~isnan(y)));
        [SupFig1Stats.(str){ii}.n] = [sum(~isnan(x)&~isnan(y)), sum(~isnan(x)&~isnan(y))];
    end

    title(strcat(num2str(SupFig1Stats.(str){1}.p),'|',num2str(SupFig1Stats.(str){2}.p),'|',num2str(SupFig1Stats.(str){3}.p),'|',...
        num2str(SupFig1Stats.(str){4}.p),'|',num2str(SupFig1Stats.(str){5}.p),'|',num2str(SupFig1Stats.(str){6}.p)));
end


expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
saveas(gcf,strcat(expPath,'SupFigure1_probebehaviorSummary.png'));
saveas(gcf,strcat(expPath,'SupFigure1_probebehaviorSummary.eps'),'epsc');
saveas(gcf,strcat(expPath,'SupFigure1_probebehaviorSummary.fig'));
save(strcat(expPath,'SupFigure1_probebehaviorSummary.mat'),'SupFig1Stats'); 
end