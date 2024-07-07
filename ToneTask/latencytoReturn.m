function Summary = latencytoReturn

sesstoAnalyze = {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
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
    'IZ47\Final\IZ47_230626_sess15','IZ47\Final\IZ47_230707_sess24',...
    'IZ47\Final\IZ47_230710_sess25','IZ47\Final\IZ47_230712_sess27',...
    'IZ48\Final\IZ48_230628_sess17','IZ48\Final\IZ48_230703_sess21',...
    'IZ48\Final\IZ48_230705_sess22','IZ48\Final\IZ48_230714_sess28',... 
    }; 

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';

Summary.latencyCorrect = [];
Summary.latencyIncorrect = [];

Summary.numberSamplesCorrect = [];
Summary.lickedPortCorrect = [];
Summary.numberSamplesIncorrect = [];
Summary.lickedPortIncorrect = [];

for ii = 1:length(sesstoAnalyze)
    cd(strcat(expPath,'\',sesstoAnalyze{ii}))
    file = dir('*.TrialBehavior.Behavior.mat');
    load(file.name);
    file = dir('*.DigitalIn.Events.mat');
    load(file(1).name);

    licktimes = [];
    lickport1 = [];
    lickport = [1 2 3 4 5 6 7];
    for ll = 1:7
        licktimes = [licktimes; digitalIn.timestampsOn{lickport(ll)+2}]; 
        lickport1 = [lickport1; ones(length(digitalIn.timestampsOn{lickport(ll)+2}),1)*lickport(ll)];
    end  

    %Sort licktimes
    [licktimesSorted,lickIdx] = sort(licktimes,'ascend');
    lickportSorted = lickport1(lickIdx);
    idxNew = [1;find(diff(lickportSorted)~=0)+1];
    
    tsLick2 = licktimesSorted(idxNew);
    lickPortLoc = lickportSorted(idxNew);

    % Exclude the home port licks
    tsLick = tsLick2(lickPortLoc~=1);
    lickLoc = lickPortLoc(lickPortLoc~=1);

    %% Latency after a correct trial
    idxLin = find(behavTrials.linTrial(1:end-1) ==0 & behavTrials.correct(1:end-1) ==1);
    tsFwd = behavTrials.timestamps(idxLin+1,1);
    tsRet = behavTrials.timestamps(idxLin,2);
    ints = [tsRet tsFwd];
    Summary.latencyCorrect = [Summary.latencyCorrect; ints(:,2)-ints(:,1)];

    %% Latency after an incorrect trial
    idxLin = find(behavTrials.linTrial(1:end-1) ==0 & behavTrials.correct(1:end-1) ==0);
    tsFwd = behavTrials.timestamps(idxLin+1,1);
    tsRet = behavTrials.timestamps(idxLin,2);
    ints = [tsRet tsFwd];
    Summary.latencyIncorrect = [Summary.latencyIncorrect; ints(:,2)-ints(:,1)];

    %% Number of samples after a correct trial
    idxLin = find(behavTrials.linTrial(1:end-1) ==0 & behavTrials.correct(1:end-1) ==1);
    tsFwd = behavTrials.timestamps(idxLin+1,1);
    tsRet = behavTrials.timestamps(idxLin,2);
    ints = [tsRet tsFwd];
    for yy = 1:size(ints,1)
        intPer = InIntervals(tsLick(:,1),ints(yy,:));  
        Summary.numberSamplesCorrect = [Summary.numberSamplesCorrect;sum(intPer)];   
        Summary.lickedPortCorrect = [Summary.lickedPortCorrect;lickLoc(intPer)];
    end

    %% Number of samples after an incorrect trial
    idxLin = find(behavTrials.linTrial(1:end-1) ==0 & behavTrials.correct(1:end-1) ==0);
    tsFwd = behavTrials.timestamps(idxLin+1,1);
    tsRet = behavTrials.timestamps(idxLin,2);
    ints = [tsRet tsFwd];
    for yy = 1:size(ints,1)
        intPer = InIntervals(tsLick(:,1),ints(yy,:));
        Summary.numberSamplesIncorrect = [Summary.numberSamplesIncorrect;sum(intPer)];   
        Summary.lickedPortIncorrect = [Summary.lickedPortIncorrect;lickLoc(intPer)];
    end
end
end