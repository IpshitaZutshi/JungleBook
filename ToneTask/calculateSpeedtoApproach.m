function Summary = calculateSpeedtoApproach

sess= {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
        'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
        'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   
        'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...
        'IZ40\Final\IZ40_220708_sess17','IZ40\Final\IZ40_220714_sess18',...
        'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
        'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
        'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
        'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    
        'IZ44\Final\IZ44_220827_sess4', 'IZ44\Final\IZ44_220828_sess5',...
        'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...
        'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
        'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',...
        };
    
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';

timetoApproach = 1;
velCorr = [];
velInCorr = [];

for ii = 1:length(sess)
    cd(strcat(expPath,sess{ii}))    
    file = dir(['*.Tracking.Behavior.mat']);
    load(file(1).name);
    file = dir(['*TrialBehavior.Behavior.mat']);
    load(file(1).name);
    [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts', true);    
    
    if ~isfield(behavTrials,'probe')
        behavTrials.probe(1:length(behavTrials.linTrial),1) = 0;
    end
    
    idxCorrect = find(behavTrials.linTrial==0 & behavTrials.correct==1 & behavTrials.probe==0 & behavTrials.stim==0 & behavTrials.toneGain~=5);
    
    idxIncorrect = find(behavTrials.linTrial==0 & behavTrials.correct==0 & behavTrials.probe==0 & behavTrials.stim==0 & behavTrials.toneGain~=5);
    
    for kk = 1:(length(idxCorrect)-1)
        [~,idxTS] = min(abs(tracking.timestamps-behavTrials.timestamps(idxCorrect(kk),2)));
        velCorr = [velCorr;tracking.position.v(idxTS-(timetoApproach*30):idxTS+(0.5*30))'];
    end
    
    for kk = 1:(length(idxIncorrect)-1)
        [~,idxTS] = min(abs(tracking.timestamps-behavTrials.timestamps(idxIncorrect(kk),2)));
        velInCorr = [velInCorr;tracking.position.v(idxTS-(timetoApproach*30):idxTS+(0.5*30))'];
    end
end
end