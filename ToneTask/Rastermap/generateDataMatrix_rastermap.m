%% Variables of interest: 
% y position on the track 
% Licking events 
% TrialEnd 
%
% Spike counts 30 ms
% 
% Spike matrix

function generateDataMatrix_rastermap(varargin)

%% Defaults and Parms - Change for each session
p = inputParser;
addParameter(p,'root','Z:\Homes\zutshi01\Recordings\Auditory_Task\',@isstr);

parse(p,varargin{:});
root = p.Results.root;

sess = {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
    'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
    'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   23
    'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...
    'IZ40\Final\IZ40_220708_sess17','IZ40\Final\IZ40_220714_sess18',...27
    'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
    'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
    'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
    'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    36
    'IZ44\Final\IZ44_220827_sess4', 'IZ44\Final\IZ44_220828_sess5',...
    'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...
    'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
    'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',... 45
    'IZ47\Final\IZ47_230626_sess15','IZ47\Final\IZ47_230707_sess24',...
    'IZ47\Final\IZ47_230710_sess25','IZ47\Final\IZ47_230712_sess27',...49
    'IZ48\Final\IZ48_230628_sess17','IZ48\Final\IZ48_230703_sess21',...
    'IZ48\Final\IZ48_230705_sess22','IZ48\Final\IZ48_230714_sess28'};    

for s = 1:length(sess)
    %close all

    disp(sess{s})
    basepath = fullfile(root, sess{s});
    destination = 'C:\Data\Rastermap\preProcessedData';

    cd(fullfile(root,sess{s}))
    
    %% (1) Copy all files in a different directory and load all the .mat files
    if ~isempty(dir([basepath filesep '*.Tracking.Behavior.mat'])) 
        file = dir([basepath filesep '*.Tracking.Behavior.mat']);
        load(file(1).name)
    end
    
    if ~isempty(dir([basepath filesep '*TrialBehavior.Behavior.mat'])) 
        file = dir([basepath filesep '*TrialBehavior.Behavior.mat']);
        load(file(1).name)
    end
    
    if ~isempty(dir([basepath filesep '*.spikes.cellinfo.mat']))
        file = dir([basepath filesep '*.spikes.cellinfo.mat']);
        load(file(1).name)
    end
    
    if ~isempty(dir([basepath filesep '*.cell_metrics.cellinfo.mat']))
        file = dir([basepath filesep '*.cell_metrics.cellinfo.mat']);
        load(file(1).name)
    end
    

    [sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);
    
    % Load digitalIn
    if ~isempty(dir([basepath filesep '*.DigitalIn.events.mat']))
        file = dir([basepath filesep '*.DigitalIn.events.mat']);
        load(file(1).name)
    end

    %% (2) Change what defines the beginning and the end of the trials. 
    % Until now,it was defined based on DigitalIn {10} and {2} respectively. 
    % Now, we want to define them based on the licking, as there is a ~33 ms 
    % delay between the two signals. 

    indStart = [];
    for i = 1:length(digitalIn.timestampsOn{3})
        [~,indStart(i)] = min(abs(digitalIn.timestampsOn{3}(i).' - digitalIn.timestampsOn{10}));
    end
    
    indStart_new = [];
    indStart_new(1) = 1;
    indStart_new = [indStart_new, find(diff(indStart) == 1) + 1];

    behavTrials.timestampsNew(:,1) = digitalIn.timestampsOn{3}(indStart_new);
    
    % Find the median delay between the lick and solenoid
    delay = median(digitalIn.timestampsOn{10}(1:end) - ...
        behavTrials.timestampsNew(:));
    
    % Adjust and shift the timestamps by the median delay between a lick and
    % the solenoid
    behavTrials.timestamps = behavTrials.timestamps - abs(delay);
    
    %% (3) Generate spike matrix 
    % To do that, we use the dt as defined by the tracking.timestamps. 
    
    dtime = mean(diff(tracking.timestamps));
    
    % Excluding the last trial
    win = [behavTrials.timestamps(1,1) behavTrials.timestamps(end,2)];
    spkData = bz_SpktToSpkmat(spikes,'dt',dtime,'win',win);
    
    spkMat = spkData.data';
    timestamps = spkData.timestamps';

    %Only select pyramidal cells
    logicalVector = cellfun(@(x) strcmp(x, 'Pyramidal Cell'), cell_metrics.putativeCellType);
    %Only select cells with a rate> 0.1 Hz
    rate = sum(spkMat,2)./(length(timestamps)*(1/30));
    logVector2 = rate>0.1;
    keepCells = logicalVector'& logVector2;

    spkMat = spkMat(keepCells,:);

    %% (4) Generate matrix for continuous variables x, vel, vy, cyclicY, cyclicYLin.
    % Note, that we later change the length of these variables to match the
    % length of the behavior variables in case frames have been dropped. 
    
    % Find indices of start of first trial and end of last trial 
    [~,idxStart] = min(abs(tracking.timestamps-timestamps(1)));
    [~,idxEnd] = min(abs(tracking.timestamps-timestamps(end)));
    
    % Interpolate
    timeS = tracking.timestamps(idxStart:idxEnd);
    constVar.x = interp1(timeS,tracking.position.x(idxStart:idxEnd),spkData.timestamps)';
    constVar.y = interp1(timeS,tracking.position.y(idxStart:idxEnd),spkData.timestamps)';
    constVar.vy = interp1(timeS,tracking.position.vy(idxStart:idxEnd),spkData.timestamps)';

    %% Find index in matrix to start and end points of each trial
    numTrials = size(behavTrials.timestamps,1);
    
    startIdx = [];
    endIdx = [];
    for tt = 1:numTrials
        [~,startIdx(tt)] = min(abs(spkData.timestamps - behavTrials.timestamps(tt,1))); 
        [~,endIdx(tt)] = min(abs(spkData.timestamps - behavTrials.timestamps(tt,2))); 
    end
    
    %% Generate matrix for events - licks, trialStart, trialEnd
    
    lickport = [1 2 3 4 5 6 7];
    for ll = 1:7
        event.times{ll} = digitalIn.timestampsOn{lickport(ll)+2}; 
        event.UID(ll) = ll;
    end
    
    event.times{8} = behavTrials.timestamps(:,1);
    event.times{9} = behavTrials.timestamps(:,2);
    
    choiceIncorrect = behavTrials.correct == 0 & behavTrials.linTrial == 0;
    event.times{10} = behavTrials.timestamps(choiceIncorrect,2);
    
    choiceCorrect = behavTrials.correct == 1;
    event.times{11} = behavTrials.timestamps(choiceIncorrect,2);
    event.UID(8:11)  = [8 9 10 11];
    
    eventMat = bz_SpktToSpkmat(event,'dt',dtime,'win',win);
    
    % Add the first and last trials
    eventMat.data(1,8) = 1;
    eventMat.data(end,9) = 1;
    
    % Store other trial relevant information.
    eventVar.trialStart = eventMat.data(:,8)';
    eventVar.trialEnd = eventMat.data(:,9)';

    eventVar.licks = eventMat.data(:,1:7)';

    % Create a single vector with all licking events - logical.
    eventVar.licksAll = eventMat.data(:,1)' | eventMat.data(:,2)' | ...
        eventMat.data(:,3)' | eventMat.data(:,4)' | eventMat.data(:,5)' | ...
        eventMat.data(:,6)' | eventMat.data(:,7)';
    
    % Create a vector with licking events according to port.
    eventVar.licksPorts = eventMat.data(:,1:7)';
    for i = 1:7
        indices = find(eventVar.licks(i,:) == 1);
        eventVar.licksPorts(i,indices) = i;
    end
    
    eventVar.licksPortsAll = zeros(size(eventVar.licksAll));
    for i = 1:7
        indices = find(eventVar.licks(i,:) == 1);
        eventVar.licksPortsAll(indices) = ...
            eventVar.licksPorts(i,indices);
    end
    
    % Create a vector with the id of the lick port only at the index where it
    % first appears anew. 
    eventVar.licksPortsAllSingle = zeros(size(eventVar.licksAll));
    
    licksIdx = find(eventVar.licksPortsAll ~= 0);
    difference = diff(eventVar.licksPortsAll(licksIdx));
    singlesIdx = licksIdx(find(difference ~=0) + 1);
    
    eventVar.licksPortsAllSingle(singlesIdx) = eventVar.licksPortsAll(singlesIdx);
    eventVar.licksPortsAllSingle(1) = 1;
    
    % Create a vector with single licks when they first appear - logical. 
    eventVar.lickEvents = zeros(size(eventVar.licksAll));
    eventVar.lickEvents(singlesIdx) = 1;
   
    %% Generate logical variables - 
    % trial number, forward/reverse, toneOn/Off, trial gain, current choice, 
    % correct/incorrect, past choice, past choice correct/incorrect,
    
    logVar.correct = zeros(1,length(spkData.timestamps));
    logVar.incorrect = zeros(1,length(spkData.timestamps));
    logVar.probe = zeros(1,length(spkData.timestamps));
    logVar.lickloc = zeros(1,length(spkData.timestamps));
    logVar.trialNum = zeros(1, length(length(spkData.timestamps)));


    %last lin track trial 
    choiceCorrect = behavTrials.correct;
    choiceInCorrect = ~behavTrials.correct;
    if isfield(behavTrials,'probe')
        probe = behavTrials.probe;
    end
    
    for tt = 1:(length(startIdx)-1)
        logVar.correct(startIdx(tt):startIdx(tt+1)) = choiceCorrect(tt);
        logVar.incorrect(startIdx(tt):startIdx(tt+1)) = choiceInCorrect(tt);
        logVar.lickloc(startIdx(tt):startIdx(tt+1)) = (behavTrials.lickLoc(tt)+1)/6;
        logVar.trialNum(startIdx(tt):startIdx(tt+1)) = tt;
        if isfield(behavTrials,'probe')
            logVar.probe(startIdx(tt):startIdx(tt+1)) = probe(tt);
        end        
    end

    %% Output for the rastermap
    eventVariables.trialStart = eventVar.trialStart';
    eventVariables.trialEnd = eventVar.trialEnd';    
    eventVariables.licks = eventVar.lickEvents';
    eventVariables.correct = logVar.correct';
    eventVariables.incorrect = logVar.incorrect';  
    eventVariables.lickloc = logVar.lickloc'; 
    eventVariables.trialNum = logVar.trialNum;
    if isfield(behavTrials,'probe')
        eventVariables.probe = logVar.probe';    
    end
    constVariables.x = constVar.x';
    constVariables.y = constVar.y';
    constVariables.vy = constVar.vy';
    
    save(fullfile(destination, strcat(extractAfter(sess{s}, 'Final\'), ...
        '.rastermapData.mat')), 'eventVariables', 'constVariables','spkMat','timestamps'); 

end

end