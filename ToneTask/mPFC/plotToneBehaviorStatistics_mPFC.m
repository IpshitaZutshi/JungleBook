function Summary = plotToneBehaviorStatistics_mPFC

sesstoAnalyze = {'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   23
        'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...        
        'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    36        
        'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15'};  

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';

Summary.PortDiffDistBase = [];
Summary.PortDiffDistStim = [];

for ii = 1:length(sesstoAnalyze)
    cd(strcat(expPath,'\',sesstoAnalyze{ii}))
    file = dir('*.TrialBehavior.Behavior.mat');
    load(file.name);
    
    if strcmp(sesstoAnalyze{ii}(1:4),'IZ39')==1
        Summary.mouseID(ii) = 1;
    elseif strcmp(sesstoAnalyze{ii}(1:4),'IZ40')==1
        Summary.mouseID(ii) = 2;
    elseif strcmp(sesstoAnalyze{ii}(1:4),'IZ43')==1
        Summary.mouseID(ii) = 3;
    elseif strcmp(sesstoAnalyze{ii}(1:4),'IZ44')==1
        Summary.mouseID(ii) = 4;
    end
        
    numBaseToneTrials = sum(behavTrials.linTrial == 0 & behavTrials.stim == 0);
    numStimToneTrials = sum(behavTrials.linTrial == 0 & behavTrials.stim == 1);
    
    for kk = 1:6
        % Collect distribution of tone gains
        Summary.trialTypeBase(ii,kk)  = sum(behavTrials.toneGain == (kk-1) & behavTrials.linTrial == 0 & behavTrials.stim == 0)./numBaseToneTrials; 
        Summary.trialTypeStim(ii,kk)  = sum(behavTrials.toneGain == (kk-1) & behavTrials.linTrial == 0 & behavTrials.stim == 1)./numStimToneTrials; 
        
        % Collect distribution of lick locations
        Summary.lickChoiceBase(ii,kk)  = sum(behavTrials.lickLoc == (kk-1) & behavTrials.linTrial == 0 & behavTrials.stim == 0)./numBaseToneTrials; 
        Summary.lickChoiceStim(ii,kk)  = sum(behavTrials.lickLoc == (kk-1) & behavTrials.linTrial == 0 & behavTrials.stim == 1)./numStimToneTrials; 
        
    end
    
    Summary.PortDiffDistBase = [Summary.PortDiffDistBase; (behavTrials.toneGain(behavTrials.linTrial == 0 & behavTrials.stim == 0)-behavTrials.lickLoc(behavTrials.linTrial == 0 & behavTrials.stim == 0))];
    Summary.PortDiffDistStim = [Summary.PortDiffDistStim; (behavTrials.toneGain(behavTrials.linTrial == 0 & behavTrials.stim == 1)-behavTrials.lickLoc(behavTrials.linTrial == 0 & behavTrials.stim == 1))];
    
    % Collect performance 
    Summary.performance(ii,:) = [sum(behavTrials.correct==1 & behavTrials.linTrial == 0 & behavTrials.stim == 0)./numBaseToneTrials ...
        sum(behavTrials.correct==1 & behavTrials.linTrial == 0 & behavTrials.stim == 1)./numStimToneTrials];
    
    
end

% saveas(gcf,strcat(expPath,'\Compiled\behaviorSummary.png'));
% saveas(gcf,strcat(expPath,'\Compiled\behaviorSummary.eps'),'epsc');
% saveas(gcf,strcat(expPath,'\Compiled\behaviorSummary.fig'));
% save(strcat(expPath,'\Compiled\behaviorSummary.mat'),'stats'); 

end

