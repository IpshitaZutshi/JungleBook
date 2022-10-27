function compilemPFCSilencing

sesstoAnalyze = {'IZ44\Final\IZ44_220920_sess15\IZ44_220920_130753','IZ44\Final\IZ44_220915_sess13\IZ44_220915_163652',...
    'IZ43\Final\IZ43_220920_sess15\IZ43_220920_122226','IZ43\Final\IZ43_220915_sess13\IZ43_220915_172502',...
    'IZ40\IZ40_220707_sess16\IZ40_220707_193451','IZ40\IZ40_220705_sess15\IZ40_220705_124255',...
    'IZ39\IZ39_220707_sess17\IZ39_220707_172123','IZ39\IZ39_220705_sess16\IZ39_220705_133945'};
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';
behavMat = [];
baselinePortDist = [];
baselinePortDiffDist = [];
stimPortDist = [];
stimPortDiffDist = [];

for ii = 1:length(sesstoAnalyze)
    cd(strcat(expPath,'\',sesstoAnalyze{ii}))
    file = dir('*.TrialBehavior.Events.mat');
    load(file.name);
    behavMat = catpad(1,behavMat,behavTrials.performance);
    
    %Also look at port behavior
    idxNoStim = behavTrials.linTrial==0&behavTrials.stim==0;
    idxStim = behavTrials.linTrial==0&behavTrials.stim==1;
    baselinePortDist = [baselinePortDist; behavTrials.lickLoc(idxNoStim)];
    stimPortDist = [stimPortDist;behavTrials.lickLoc(idxStim)];

    baselinePortDiffDist = [baselinePortDiffDist; (behavTrials.toneGain(idxNoStim)-behavTrials.lickLoc(idxNoStim))];
    stimPortDiffDist = [stimPortDiffDist;(behavTrials.toneGain(idxStim)-behavTrials.lickLoc(idxStim))];    
end

behavMat(behavMat>1)  = 1;
figure
set(gcf, 'renderer','painters')
set(gcf,'Color','w')
subplot(2,3,[1 2])
plot(behavMat(:,2:(end-2))','Color',[0.5 0.5 0.5])
hold on
plot(nanmean(behavMat(:,2:(end-2)),1),'LineWidth',2,'Color','k')
ylim([0 1.1])
ylabel('Performance');
avgBase = nanmean(behavMat(:,2:2:(end-2)),2);
avgStim = nanmean(behavMat(:,3:2:(end-2)),2);
subplot(2,3,3)
plot([avgBase avgStim]','Color',[0.5 0.5 0.5])
xlim([0.5 2.5])
hold on
plot([nanmean(avgBase) nanmean(avgStim)],'LineWidth',2,'Color','k')
pVal = signrank(avgBase,avgStim);
title(strcat('n:',num2str(length(avgBase)),',p:',num2str(pVal)))

subplot(2,3,4)
data.base = baselinePortDist;
data.stim = stimPortDist;
nhist(data,'samebins','proportion')
pVal = ranksum(baselinePortDist,stimPortDist);
title(strcat('n:',num2str(length(baselinePortDist)),',p:',num2str(pVal)))
xlabel('Distribution of licks')

subplot(2,3,5)
data.base = baselinePortDiffDist;
data.stim = stimPortDiffDist;
nhist(data,'samebins','proportion')
pVal = ranksum(baselinePortDiffDist,stimPortDiffDist);
title(strcat('n:',num2str(length(baselinePortDiffDist)),',p:',num2str(pVal)))
xlabel('Distribution of lick difference')

end