%% Variables of interest: 
% y position on the track (circular) - tone and no-tone trials
% Distance to stop: normalized and rescaled, equivalent to path length
% Frequency: normalized and rescaled 
% Licking type: tone/choice, tone/spontaneous, no-tone/choice,
%   no-tone/spontaneous, home 
% Licking events 
% Theta LFP 
%
% Spike counts 30 ms
% 
% Spike matrix

function [spkMat,constVar,logVar,eventVar,timestamps] = generateDataMatrix_AAIZ_new(varargin)

%% Defaults and Parms - Change for each session
p = inputParser;
addParameter(p,'root','Z:\Homes\zutshi01\Recordings\Auditory_Task\',@isstr);

parse(p,varargin{:});
root = p.Results.root;

sess = {'IZ48\Final\IZ48_230628_sess17','IZ48\Final\IZ48_230714_sess28'};
% sess= {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
%     'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
%     'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   
%     'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...
%     'IZ40\Final\IZ40_220708_sess17','IZ40\Final\IZ40_220714_sess18',...
%     'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
%     'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
%     'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
%     'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    
%     'IZ44\Final\IZ44_220827_sess4', 'IZ44\Final\IZ44_220828_sess5',...
%     'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...
%     'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
%     'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',...
%     'IZ47\Final\IZ47_230626_sess15','IZ47\Final\IZ47_230707_sess24','IZ47\Final\IZ47_230710_sess25',...
%     'IZ47\Final\IZ47_230712_sess27'}; 

for s = 1:length(sess)
    %close all

    disp(sess{s})
    basepath = fullfile(root, sess{s});
    destination = 'C:\Data\PGAMAnalysis\data';
    %destination = strcat('C:\Data\PGAMAnalysis\sessData\', sess{s}(11:end));
%     if ~exist(destination)
%         mkdir(destination)
%     end
%     cd(sessions{s})

    cd(fullfile(root,sess{s}))
    
    %% (1) Copy all files in a different directory and load all the .mat files
    if ~isempty(dir([basepath filesep '*.Tracking.Behavior.mat'])) 
        disp('Tracking already detected! Loading file.');
        file = dir([basepath filesep '*.Tracking.Behavior.mat']);
%         copyfile(file(1).name, destination);
        load(file(1).name)
    end
    
    if ~isempty(dir([basepath filesep '*TrialBehavior.Behavior.mat'])) 
        disp('Behavior already detected! Loading file.');
        file = dir([basepath filesep '*TrialBehavior.Behavior.mat']);
%         copyfile(file(1).name, destination);
        load(file(1).name)
    end
    
    if ~isempty(dir([basepath filesep '*.spikes.cellinfo.mat']))
        disp('Spikes already detected! Loading file.');
        file = dir([basepath filesep '*.spikes.cellinfo.mat']);
%         copyfile(file(1).name, destination);
        load(file(1).name)
    end
    
    [sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);
    
    % Load digitalIn
    if ~isempty(dir([basepath filesep '*.DigitalIn.events.mat']))
        disp('DigitalIn already detected! Loading file.');
        file = dir([basepath filesep '*.DigitalIn.events.mat']);
%         copyfile(file(1).name, destination);
        load(file(1).name)
    end
%     filename = [basepath filesep sessionInfo.session.name '.DigitalIn.events.mat'];
%     load();
    
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
    
%     indStart = [];
%     indStart(1) = 1; 
%     
%     % Note that the threshold might need to be different in some cases. Check
%     % that the sizes of all variables match in the end. 
%     for i = 2:length(digitalIn.timestampsOn{3})
%         if abs(digitalIn.timestampsOn{3}(i)-digitalIn.timestampsOn{3}(i-1))>=4
%             indStart = [indStart, i];
%         end
%     end
%     
%     % Ensure that one licking timestamp is kept for each trial. Every timestamp 
%     % that is far away in time from either the start or the end of that trial 
%     % is discarded.  
%     j = 1;
%     indStart_new = indStart;
%     
%     for i = 1:(length(indStart))
%         if abs(digitalIn.timestampsOn{3}(indStart(i))-digitalIn.timestampsOn{10}(j))>0.1 & ...
%             abs(digitalIn.timestampsOn{3}(indStart(i))-digitalIn.timestampsOn{2}(j+1))>0.1
%     
%             indStart_new(indStart_new == indStart(i)) = [];
%         else
%             j = j + 1;
%             continue
%         end  
%     end
%     
%     % Remove the last element, as this would mark the end of the last trial and
%     % not the beginning of a new. 
%     % indStart(end) = [];
%     % indStart_new(end) = [];
    
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

%     constVar.x = tracking.position.x(idxStart:idxEnd)'; % [cm]
%     constVar.y = tracking.position.y(idxStart:idxEnd)'; % [cm]
%     constVar.vy = tracking.position.vy(idxStart:idxEnd)'; % [cm/s]
%     
    constVar.yFwd = constVar.y;
    constVar.yRev = constVar.y;
    constVar.yFwd(constVar.vy<-1) = nan;
    constVar.yLin = constVar.yFwd;
    constVar.yRev(constVar.vy>0) = nan;

%     constVar.hd = tracking.position.angle(idxStart:idxEnd)'; % [degrees]
    
    % Convert the y position to cover both directions in the range [0,240].
%     constVar.ylong = constVar.y;
% %     indices = find(constVar.vy < 0 & (constVar.hd > 85 | constVar.hd < -85));
%     indices = find(constVar.vy < 0);
%     constVar.ylong(indices) = ...
%         240 - constVar.ylong(indices);
%     constVar.ylong = smooth(constVar.ylong, 20, 'rlowess')';
%     
%     constVar.cyclicY = mod(constVar.ylong, 240); % Tone trials
%     constVar.cyclicYLin = constVar.cyclicY'; % No-tone trials
    
    %% Find index in matrix to start and end points of each trial
    numTrials = size(behavTrials.timestamps,1);
    
    startIdx = [];
    endIdx = [];
    for tt = 1:numTrials
        [~,startIdx(tt)] = min(abs(spkData.timestamps - behavTrials.timestamps(tt,1))); 
        [~,endIdx(tt)] = min(abs(spkData.timestamps - behavTrials.timestamps(tt,2))); 
    end
    
    %% (5) Generate logical variables - 
    % trial number, forward/reverse, toneOn/Off, trial gain, current choice, 
    % correct/incorrect, past choice, past choice correct/incorrect
    
    logVar.forward = zeros(1, length(length(spkData.timestamps)));
    logVar.toneOn = zeros(1, length(length(spkData.timestamps)));
    logVar.trialType = zeros(1, length(length(spkData.timestamps)));
    logVar.currChoice = zeros(1, length(length(spkData.timestamps)));
    logVar.currCorrect = zeros(1, length(length(spkData.timestamps)));
    logVar.trialNum = zeros(1, length(length(spkData.timestamps)));
    
    %last lin track trial 
    lastLin = behavTrials.linTrial;
    toneGain = (behavTrials.toneGain);
    toneGain(~lastLin) = toneGain(~lastLin)+1;
    toneGain(lastLin==1) = 0;
    choiceCorrect = behavTrials.correct;
    choiceCorrect(lastLin==1) = NaN;
    
    % for tt = 1:(length(startIdx)-1)        
    for tt = 1:length(startIdx)        
        if tt <= length(startIdx)-1
            logVar.toneOn(startIdx(tt):startIdx(tt+1)) = ~behavTrials.linTrial(tt);
            logVar.trialNum(startIdx(tt):startIdx(tt+1)) = tt;
            logVar.trialType(startIdx(tt):startIdx(tt+1)) = toneGain(tt);
            logVar.currChoice(startIdx(tt):startIdx(tt+1)) = behavTrials.lickLoc(tt)+1;
            logVar.currCorrect(startIdx(tt):startIdx(tt+1)) = choiceCorrect(tt);
        elseif tt == length(startIdx)
            logVar.toneOn(startIdx(tt):endIdx(tt)) = ~behavTrials.linTrial(tt);
            logVar.trialNum(startIdx(tt):endIdx(tt)) = tt;
            logVar.trialType(startIdx(tt):endIdx(tt)) = toneGain(tt);
            logVar.currChoice(startIdx(tt):endIdx(tt)) = behavTrials.lickLoc(tt)+1;
            logVar.currCorrect(startIdx(tt):endIdx(tt)) = choiceCorrect(tt);
        end
    end
    
    
    %% (7) Now we can generate the remaining continuous variables:
    %Calculate port locations flexibly according to the session
    for pp=1:6
        ts = behavTrials.timestamps(behavTrials.lickLoc==(pp-1) & behavTrials.linTrial==0,2);
        %Find median position at that timestamp
        pos = [];
        for tsidx = 1:length(ts)
            [~,posIdx] = min(abs(tracking.timestamps-ts(tsidx)));
            pos(tsidx) = tracking.position.y(posIdx);
        end
        port_locations(pp) = nanmedian(pos);
    end
    %port_locations = [11.6, 32.27, 55.53, 79.62, 102.79, 120];
    gain = 120./port_locations;%[120/11.6, 120/32.27 120/55.53 120/79.62 120/102.79 120/120];
    freqExp = log10(22000/1000);
    
    constVar.relDistStop = NaN(size(logVar.currChoice));
    constVar.freq = NaN(size(logVar.currChoice));
    
    for i = 1:numTrials
        if logVar.trialType(startIdx(i):endIdx(i)) > 0
            
            % a. Distance to stop location (relative distance): normalized to 
            % 0-1 and then rescaled to 1-120. 
            constVar.relDistStop(startIdx(i):endIdx(i)) =  ...
                port_locations(logVar.currChoice(startIdx(i):endIdx(i))) - ...
                constVar.y(startIdx(i):endIdx(i));
    
            constVar.relDistStop(startIdx(i):endIdx(i)) = ...
                constVar.relDistStop(startIdx(i):endIdx(i))./ ...
                port_locations(logVar.currChoice(startIdx(i):endIdx(i))) * 120;
    
            % b. Frequency: normalized to 0-1 and then rescaled to 1-120. 
            freq = (constVar.y(startIdx(i):endIdx(i)).* ...
                gain(logVar.trialType(startIdx(i))))/120;
            constVar.freq(startIdx(i):endIdx(i)) = ...
                ((1000*(10.^(freqExp*freq)))/22000)*120;
        end
    end
    
    correct = find(constVar.relDistStop < 0 & ~isnan(constVar.relDistStop));
    constVar.relDistStop(correct) = 0;
    
    % Update relDistStop and trialType to include the return runs with goal = home port
    constVar.relDistStop_2d = NaN(size(logVar.currChoice));
    logVar.trialType_2d = logVar.currChoice;
    
    % forward runs
    for i = 1:numTrials
        constVar.relDistStop_2d(startIdx(i):endIdx(i)) =  ...
                port_locations(logVar.currChoice(startIdx(i):endIdx(i))) - ...
                constVar.y(startIdx(i):endIdx(i));

        constVar.relDistStop_2d(startIdx(i):endIdx(i)) = ...
                constVar.relDistStop_2d(startIdx(i):endIdx(i))./ ...
                port_locations(logVar.currChoice(startIdx(i):endIdx(i))) * 120;
    end
    correct = find(constVar.relDistStop_2d < 0 & ~isnan(constVar.relDistStop_2d));
    constVar.relDistStop_2d(correct) = 0;

    % return runs
    for i = 1:numTrials-1
        constVar.relDistStop_2d(endIdx(i):startIdx(i+1)) = ...
            abs(0 - constVar.y(endIdx(i):startIdx(i+1)));

        relDistStop_return = constVar.relDistStop_2d(endIdx(i):startIdx(i+1));
        relDistStop_return(constVar.vy(endIdx(i):startIdx(i+1)) > 0) = NaN;

        constVar.relDistStop_2d(endIdx(i):startIdx(i+1)) = relDistStop_return/ ...
            max(constVar.y(endIdx(i):startIdx(i+1))) * 120;

        trialType_return = logVar.trialType_2d(endIdx(i):startIdx(i+1));
        trialType_return(constVar.vy(endIdx(i):startIdx(i+1)) < 0) = 7;
        logVar.trialType_2d(endIdx(i):startIdx(i+1)) = trialType_return;
    end
    
    constVar.freq(constVar.freq > 120) = nan;  % no tone after correct lick
    
    % c. Update the y and ylin variables
    constVar.yLin(logVar.trialType>0) = nan; % No-tone trials
    constVar.yFwd(logVar.trialType==0) = nan; % Tone trials
    
    %% (8) Generate matrix for events - licks, trialStart, trialEnd
    
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
    eventVar.incorrect = eventMat.data(:,10)';
    eventVar.correct = eventMat.data(:,11)';
    
    % Note that these licks now are only the ones within the boundaries of the
    % first and last trial (win) - logical.
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
    
    %% (9) Separate the types of licks: 
    % 1. tone choice lick, 2. tone spontaneous lick, 3. no tone end lick, 
    % 4. no-tone spontaneous lick, 5. home lick 
    
    % 1. All tone trial choice licks
    % ind1 = find(eventVar.licksPortsAllSingle == (logVar.currChoice+1) & ...
    %     eventVar.licksPortsAllSingle ~= 0 & logVar.toneOn==1); 
    ind1 = find(eventVar.licksPortsAllSingle == (logVar.currChoice+1) & ...
        eventVar.licksPortsAllSingle ~= 0 & ...
        eventVar.licksPortsAllSingle ~= 1 & logVar.toneOn==1); 
    
    % currChoice goes from 1-6, while licks are from port 1 (home) to 7. add 1
    
    % 2. All tone trial spontaneous licks
    ind2 = find(eventVar.licksPortsAllSingle ~= (logVar.currChoice+1) & ...
        eventVar.licksPortsAllSingle ~= 0 & ...
        eventVar.licksPortsAllSingle ~= 1 & logVar.toneOn==1);
    
    % 3. All no tone trial end licks
    ind3 = find(eventVar.licksPortsAllSingle == 7 & logVar.toneOn==0);
    
    % 4. All no tone trial middle port licks
    ind4 = find(eventVar.licksPortsAllSingle ~= 7 & ...
        eventVar.licksPortsAllSingle ~= 1 & ...
        eventVar.licksPortsAllSingle ~= 0 & logVar.toneOn==0);
    
    % 5. All home ports
    ind5 = find(eventVar.licksPortsAllSingle == 1);
    
    % Create a single vector with a different type of lick indicated.
    eventVar.licksType = zeros(size(eventVar.licksAll));
    eventVar.licksType(ind1) = 1;
    eventVar.licksType(ind2) = 2;
    eventVar.licksType(ind3) = 3;
    eventVar.licksType(ind4) = 4;
    eventVar.licksType(ind5) = 5;
    
    % Create 5 different binary vectors for each type of lick. 
    eventVar.('licksType1') = zeros(size(eventVar.licksAll));
    eventVar.('licksType2') = zeros(size(eventVar.licksAll));
    eventVar.('licksType3') = zeros(size(eventVar.licksAll));
    eventVar.('licksType4') = zeros(size(eventVar.licksAll));
    eventVar.('licksType5') = zeros(size(eventVar.licksAll));
    
    eventVar.('licksType1')(ind1) = 1;
    eventVar.('licksType2')(ind2) = 1;
    eventVar.('licksType3')(ind3) = 1;
    eventVar.('licksType4')(ind4) = 1;
    eventVar.('licksType5')(ind5) = 1;
    
    %% (10) Generate LFP data 
    % % cd 'IZ44_220830_sess7'
    % % cd(string(basepath))
    % 
    % lfp = bz_GetLFP(69,'noprompts',true,'interval',win);
    % signal.data = bz_Filter(double(lfp.data),'filter','butter','passband',[6 12],'order',3);
    % signal.timestamps = lfp.timestamps;
    % signal.samplingRate = lfp.samplingRate;
    % 
    % downFactor = 20;
    % lfpdown1 = bz_DownsampleLFP(signal, downFactor);
    % 
    % [~,matchIdx] = min(abs(lfpdown1.timestamps-timestamps));
    % 
    % lfpdown = lfpdown1.data(matchIdx);
    % constVar.lfp = lfpdown';
    % 
    % % cd ..
    
    %% Plot distributions and time series of the variables and save data. 
    
    % all_variables.x = constVar.x;
    all_variables.ylin = constVar.yLin;
    all_variables.y = constVar.yFwd;
    all_variables.yRev = constVar.yRev;
    all_variables.freq = constVar.freq;
    all_variables.relDistStop = constVar.relDistStop;
    all_variables.licks = eventVar.lickEvents;
    all_variables.('licksType1') = eventVar.('licksType1');
    all_variables.('licksType2') = eventVar.('licksType2');
    all_variables.('licksType3') = eventVar.('licksType3');
    all_variables.('licksType4') = eventVar.('licksType4');
    all_variables.('licksType5') = eventVar.('licksType5');
    % all_variables.lfp = constVar.lfp;
%     
%     save(fullfile(destination, 'sessionData.mat'),'eventVar','constVar', ...
%         'timestamps','spkMat','tracking','all_variables'); 
    
    % Create plots
    fields = fieldnames(all_variables);
    
    figure;
    
    nc = 3;
    nr = ceil(length(fields)/3);
    
    for i = 1:length(fields)
        subplot(nr, nc, i);
    
        fieldname = fields{i};
        field = all_variables.(fieldname);
    
        if strcmp(fieldname, 'licksType1') || strcmp(fieldname, 'licksType2')...
                || strcmp(fieldname, 'licksType3') || strcmp(fieldname, 'licksType4')...
                || strcmp(fieldname, 'licksType5') || strcmp(fieldname, 'licks')
            field = nonzeros(field);
        end
        
        histogram(field);    
        xlabel(fieldname);
        ylabel('density');
    end
    
%     saveas(gcf, fullfile(destination, 'Variables_distribution'), 'fig');
%     saveas(gcf, fullfile(destination, 'Variables_distribution'), 'png');
%     
    figure;
    set(gcf,'Color','w')
    nr = length(fields);
    
    for i = 1:length(fields)
        subplot(nr, 1, i);
    
        fieldname = fields{i};
        field = all_variables.(fieldname);
    
        if strcmp(fieldname, 'lfp')
            plot(timestamps(1:length(all_variables.lfp)),field)
        else
            plot(timestamps,field);
        end
    
        if strcmp(fieldname, 'licksType1')
            title(strcat(fieldname,': Tone choice')) 
        elseif strcmp(fieldname, 'licksType2')
            title(strcat(fieldname,': Tone spont.')) 
        elseif strcmp(fieldname, 'licksType3')
            title(strcat(fieldname,': No tone choice')) 
        elseif strcmp(fieldname, 'licksType4')
            title(strcat(fieldname,': No-tone spont.'))
        elseif strcmp(fieldname, 'licksType5')
            title(strcat(fieldname,': home')) 
        else
            title(fieldname)
        end
        ylabel(field);
        xlim([timestamps(1) timestamps(end)])
    
        if i < length(fields)
            set(gca,'xtick',[])
        end
    
    %     axis off
    end
    xlabel('time (s)');
    
%     saveas(gcf, fullfile(destination, 'Variables_time_series'), 'fig');
%     saveas(gcf, fullfile(destination, 'Variables_time_series'), 'png');
    
    %% Output for the PGAM
    eventVariables.trialStart = eventVar.trialStart;
    eventVariables.trialEnd = eventVar.trialEnd;    
    eventVariables.('licksType1') = eventVar.('licksType1');
    eventVariables.('licksType2') = eventVar.('licksType2');
    eventVariables.('licksType3') = eventVar.('licksType3');
    eventVariables.('licksType4') = eventVar.('licksType4');
    eventVariables.('licksType5') = eventVar.('licksType5');
    eventVariables.licks = eventVar.lickEvents;
    constVariables.x = constVar.x;
    constVariables.yFwd = constVar.yFwd';
    constVariables.yRev = constVar.yRev';
    constVariables.yLin = constVar.yLin';
    constVariables.freq = constVar.freq;
    constVariables.relDistStop = constVar.relDistStop;
    % constVariables.lfp = constVar.lfp;
    
    save(fullfile(destination, strcat(extractAfter(sess{s}, 'Final\'), ...
        '.sessionDataPGAM.mat')), 'eventVariables', 'constVariables','spkMat'); 
    
    %% Output for CEBRA
%     x = constVar.x;
%     
%     % Forward direction - 1, backward direction - 0.
%     direction = ones(size(constVar.y));
%     direction(constVar.ylong > 120) = 0;
%     
%     y = zeros(2, length(constVar.y));
%     y(1,:) = constVar.y;
%     y(2,:) = direction;
%     frequency = constVar.freq;
%     relDistStop = constVar.relDistStop;
%     relDistStop_2d = constVar.relDistStop_2d;
%     % lfp = constVar.lfp;
%     licksType1 = eventVar.('licksType1');
%     licksType2 = eventVar.('licksType2');
%     licksType3 = eventVar.('licksType3');
%     licksType4 = eventVar.('licksType4');
%     licksType5 = eventVar.('licksType5');
%     licks = eventVar.lickEvents;
%     trialType = logVar.currChoice;
%     trialType_2d = logVar.trialType_2d;
%     toneTrial = logVar.toneOn;
%     
%     save(fullfile(destination, strcat(extractAfter(destination, 'Final\'), '.labelsCEBRA.mat')), ...
%         'x', 'y', 'direction', 'frequency', 'relDistStop', 'relDistStop_2d', 'toneTrial', ...
%         'licksType1', 'licksType2', 'licksType3', 'licksType4', 'licksType5', ...
%         'licks', 'trialType', 'trialType_2d', 'spkMat');
%     cd(root)
     
end

end