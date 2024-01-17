function plotDecoding_SupFigure3

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[2025 411 1360 390]);

numrows = 2;
numcol = 8;

%% Calculate place cells
Summary = compileAvgTonePlaceMaps('plotfig',false,'savefig',false);


%% Fraction of tone tuned cells
ss=2;
sessIndices = unique(Summary.AllsessID(Summary.AllsessType==(ss-1)));
for ll = 1:length(sessIndices)
    sumActive(ll) = sum(Summary.AllsessType==(ss-1) & Summary.AllcellType == 1 & nanmean(Summary.AllspaceMapAvg,2)>0.2 & Summary.AllsessID==sessIndices(ll)); 
    sumMaps(ll) = sum(Summary.AllsessType==(ss-1) & Summary.AllcellType == 1 & Summary.AlltoneField & Summary.AlltoneCorr>0.1 & Summary.AllsessID==sessIndices(ll)); 
end
sessFractTuned = [sumMaps./sumActive];
sessActive = sumActive;
subplot(numrows,numcol,[1])
bar(sessFractTuned)
line([0 37],[median(sessFractTuned) median(sessFractTuned)])
title('Non spatial cells')

subplot(numrows,numcol,[2])
bar(sessActive)
line([0 37],[median(sessActive) median(sessActive)])
title('Active cells')

%% Panel A: Decoding Accuracy example
subplot(numrows,numcol,[9 10])
% sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ48\Final\IZ48_230705_sess22';
sessloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230710_sess25';
cd(sessloc)
result_names{1} = strcat(sessloc,'\Decoding\decoding_results_0');
result_names{2} = strcat(sessloc,'\Decoding\decoding_results_2');
plot_obj = plot_standard_results_object(result_names);

pval_dir_name{1} = strcat(sessloc,'\Decoding\shuff_results\');
pval_dir_name{2} = strcat(sessloc,'\Decoding\shuff_results\');
plot_obj.p_values = pval_dir_name;
    
% % use data from all time bins when creating the null distribution
plot_obj.collapse_all_times_when_estimating_pvals = 1;
plot_obj.significant_event_times = 5001;
plot_obj.plot_results;  

subplot(numrows,numcol,[12 11])
load(strcat(result_names{1},'.mat'))  
col = [119/255 221/255 229/255;...
    122/255 201/255 229/255;...
    38/255 169/255 224/255;...
    73/255 136/255 189/255;...
    17/255 55/255 174/255;...
    0/255 0/255 134/255];

for ii = 1:6
    a = (DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.confusion_matrix_results.confusion_matrix(ii,ii,:));
    b = reshape(a,[size(a,3),1,1])./60;
    plot(b,'Color',col(ii,:),'LineWidth',1.5)
    hold on
end

plot(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results,'Color','k')
ylim([0 1])

%% Panel B: Average decoding accuracy
[meanAccuracy, mouseID, Accuracy,trialDur] = compileDecodedTrialChoice;
subplot(numrows,numcol,[13 14])
plot(Accuracy','Color',[0.5 0.5 0.5])
hold on
plot(nanmean(Accuracy,1),'Color','k','LineWidth',1.5)

subplot(numrows,numcol,[15 16])
histogram(trialDur)

%% Panel C: Correlation between decoding accuracy and number of active cells
[meanAccuracy, mouseID] = compileDecodedTrialChoice;
mouseID2 = cell2mat(mouseID);
[a,~] = unique(mouseID2);
for ii = 1:length(a)
    mouseID_re(mouseID2==a(ii)) = ii;
end
a = ~isnan(meanAccuracy);
[SupFig3Stats.sessActiveCorr.r,SupFig3Stats.sessActiveCorr.p] = corrcoef(sessActive(a),meanAccuracy(a));
subplot(numrows,numcol,[3 4])
scatter(sessActive(a)',meanAccuracy(a)','.')
hold on
lsline
gscatter(sessActive(a)',meanAccuracy(a)',mouseID_re')
box off
ylabel('Decoding Accuracy')
xlabel('Number of active cells')
title(num2str(SupFig3Stats.sessActiveCorr.r(1,2)))


%% Panel D: Correlation between decoding accuracy and fraction tuned
subplot(numrows,numcol,[5 6])
[SupFig3Stats.sessTunedCorr.r,SupFig3Stats.sessTunedCorr.p] = corrcoef(sessFractTuned(a),meanAccuracy(a));
scatter(sessFractTuned,meanAccuracy,5,'filled')
hold on
lsline
box off
ylabel('Decoding Accuracy')
xlabel('Fraction non-spatial')
title(num2str(SupFig3Stats.sessTunedCorr.r(1,2)))

%% Panel E: Correlation between decoding accuracy and performance
subplot(numrows,numcol,[7 8])
Summary = plotToneBehaviorStatistics;
[SupFig3Stats.performanceCorr.r,SupFig3Stats.performanceCorr.p] = corrcoef(sessFractTuned,Summary.performance);
scatter(sessFractTuned,Summary.performance,5,'filled')
hold on
lsline
box off
xlabel('Fraction non-spatial')
ylabel('Behavior performance')
title(num2str(SupFig3Stats.performanceCorr.r(1,2)))

%% Save figure
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure3_Decoding.png'));
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure3_Decoding.eps'),'epsc');
saveas(gcf,strcat(expPath,'\Compiled\Figures_Sep23\SupFigure3_Decoding.fig'));
save(strcat(expPath,'\Compiled\Figures_Sep23\SupFigure3_Decoding.mat'),'SupFig3Stats'); 

end