%% Variables of interest: 
% x position on the track
% y position on the track (circular) - tone and no-tone trials
% Absolute distance to the stop location
% Licking type: tone/choice, tone/spontaneous, no-tone/spontaneous, home 
% Theta LFP 
% Tone frequency?
%
% Spike counts 33.3 ms
% 
% Spike matrix

function [spkMat,constVar,logVar,eventVar,timestamps] = generateDataMatrix_AAIZ(varargin)

%addpath(genpath('F:\github\buzcode'));

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
%addParameter(p,'dt',1/30,@isnumeric); % Change this to 1/30 s

parse(p,varargin{:});
basepath = p.Results.basepath;
%dt = p.Results.dt;

%% (1) Load all the .mat files
if ~isempty(dir([basepath filesep '*.Tracking.Behavior.mat'])) 
    disp('Loading tracking');
    file = dir([basepath filesep '*.Tracking.Behavior.mat']);
    load(file(1).name);
end

if ~isempty(dir([basepath filesep '*TrialBehavior.Behavior.mat'])) 
    disp('Behavior already detected! Loading file.');
    file = dir([basepath filesep '*TrialBehavior.Behavior.mat']);
    load(file(1).name);
end

if ~isempty(dir([basepath filesep '*.spikes.cellinfo.mat']))
    disp('Spikes already detected! Loading file.');
    file = dir([basepath filesep '*.spikes.cellinfo.mat']);
    load(file.name);
end

[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);

% Load digitalIn
filename = [basepath filesep sessionInfo.session.name '.DigitalIn.events.mat'];
load(filename);

%% (2) Change what defines the beginning and the end of the trials. 
% Until now,it was defined based on DigitalIn {10} and {2} respectively. 
% Now, we want to define them based on the licking, as there is a ~33 ms 
% delay between the two signals. 

indStart = [];
indStart(1) = 1; 

for i = 2:length(digitalIn.timestampsOn{3})
    if abs(digitalIn.timestampsOn{3}(i)-digitalIn.timestampsOn{3}(i-1))>=10
        indStart = [indStart, i];
    end
end

% Remove the last element, as this would mark the end of the last trial and
% not the beginning of a new. 
indStart(end) = [];

behavTrials.timestampsNew(:,1) = digitalIn.timestampsOn{3}(indStart);

% Find the median delay between the lick and solenoid
delay = median(digitalIn.timestampsOn{10}(1:30) - ...
    behavTrials.timestampsNew(1:30));

% Adjust and shift the timestamps by the median delay between a lick and
% the solenoid
behavTrials.timestamps = behavTrials.timestamps - delay;

%% (3) Generate spike matrix 
% To do that, we use the dt as defined by the tracking.timestamps. 

dtime = mean(diff(tracking.timestamps));

win = [behavTrials.timestamps(1,1) behavTrials.timestamps(end,2)];%behavTrials.timestamps((end-1),2)];
spkData = bz_SpktToSpkmat(spikes,'dt',dtime,'win',win);

spkMat = spkData.data';
timestamps = spkData.timestamps';

%% (4) Generate matrix for continuous variables x, y, vel, velY, cyclicY.
% Note, that we later change the length of these variables to match the
% length of the behavior variables in case frames have been dropped. 

% Find indices of start of first trial and end of last trial 
[~,idxStart] = min(abs(tracking.timestamps-timestamps(1)));
[~,idxEnd] = min(abs(tracking.timestamps-timestamps(end)));


constVar.x = tracking.position.x(idxStart:idxEnd)'; % [cm]
constVar.y = tracking.position.y(idxStart:idxEnd)'; % [cm]
constVar.velY = tracking.position.vy(idxStart:idxEnd)';

% Convert the y position to cover both directions in the range [0,240].
constVar.ylong = constVar.y;
constVar.ylong(constVar.velY < -1) = 240 - constVar.ylong(constVar.velY < -1);
constVar.cyclicY = mod(constVar.ylong, 240); 
constVar.cyclicYLin = constVar.cyclicY;

%% Find index in matrix to start and end points of each trial
numTrials = size(behavTrials.timestamps,1);

for tt = 1:numTrials
    [~,startIdx(tt)] = min(abs(spkData.timestamps - behavTrials.timestamps(tt,1))); 
    [~,endIdx(tt)] = min(abs(spkData.timestamps - behavTrials.timestamps(tt,2))); 
end


%% (5) Generate logical variables - 
% trial number, forward/reverse, toneOn/Off, trial gain, current choice, 
% correct/incorrect, past choice, past choice correct/incorrect

logVar.forward(1:length(spkData.timestamps)) = 0;
logVar.toneOn(1:length(spkData.timestamps)) = 0;
logVar.trialType(1:length(spkData.timestamps)) = 0;
logVar.currChoice(1:length(spkData.timestamps)) = 0;
logVar.currCorrect(1:length(spkData.timestamps)) = 0;
logVar.trialNum(1:length(spkData.timestamps)) = 0;

%last lin track trial 
lastLin = behavTrials.linTrial;
toneGain = (behavTrials.toneGain);
toneGain(~lastLin) = toneGain(~lastLin)+1;
toneGain(lastLin==1) = 0;
choiceCorrect = behavTrials.correct;
choiceCorrect(lastLin==1) = NaN;

for tt = 1:(length(startIdx)-1)        
    if tt<length(startIdx)-1
        logVar.toneOn(startIdx(tt):startIdx(tt+1)) = ~behavTrials.linTrial(tt);
        logVar.trialNum(startIdx(tt):startIdx(tt+1)) = tt;
        logVar.trialType(startIdx(tt):startIdx(tt+1)) = toneGain(tt);
        logVar.currChoice(startIdx(tt):startIdx(tt+1)) = behavTrials.lickLoc(tt)+1;
        logVar.currCorrect(startIdx(tt):startIdx(tt+1)) = choiceCorrect(tt);
    elseif tt == length(startIdx)-1
        logVar.toneOn(startIdx(tt):endIdx(tt)) = ~behavTrials.linTrial(tt);
        logVar.trialNum(startIdx(tt):endIdx(tt)) = tt;
        logVar.trialType(startIdx(tt):endIdx(tt)) = toneGain(tt);
        logVar.currChoice(startIdx(tt):endIdx(tt)) = behavTrials.lickLoc(tt)+1;
        logVar.currCorrect(startIdx(tt):endIdx(tt)) = choiceCorrect(tt);
    end
end

%% (7) Now can generate the remaining continuous variables:

% a. Distance to stop location (absolute distance). It is only defined in 
% the boundaries of the trials and only for tone trials.  
port_locations = [11.6, 32.27, 55.53, 79.62, 102.79, 120];
constVar.distStop = NaN(size(logVar.currChoice));

for i = 1:numTrials
    if logVar.trialType(startIdx(i):endIdx(i)) > 0
        constVar.distStop(startIdx(i):endIdx(i)) =  ...
        port_locations(logVar.currChoice(startIdx(i):endIdx(i))) - ...
        constVar.y(startIdx(i):endIdx(i));
    end
end

correct = find(constVar.distStop < 0 & ~isnan(constVar.distStop));
constVar.distStop(correct) = 0;

constVar.cyclicYLin(logVar.trialType>0) = nan; % No-tone trials
constVar.cyclicY(logVar.trialType==0) = nan; % Tone trials

%% (8) Generate matrix for events - licks, trialStart, trialEnd

lickport = [1 2 3 4 5 6 7];
for ll = 1:7
    event.times{ll} = digitalIn.timestampsOn{lickport(ll)+2}; 
    event.UID(ll) = ll;
end

event.times{8} = behavTrials.timestamps(:,1);%((1:(end-1)),1);
event.times{9} = behavTrials.timestamps(:,2);%((1:(end-1)),2);

choiceIncorrect = behavTrials.correct == 0 & behavTrials.linTrial == 0;
event.times{10} = behavTrials.timestamps(choiceIncorrect,2);%(choiceIncorrect(1:(end-1)),2);

choiceCorrect = behavTrials.correct == 1;
event.times{11} = behavTrials.timestamps(choiceIncorrect,2);%(choiceCorrect(1:(end-1)),2);
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

%% (9) Separate the types of licks: 
% 1. tone choice lick, 2. tone spontaneous lick, 3. tone end lick, 4. no-tone spontaneous lick, 5. home lick 

% idxNoToneTrials = find(toneGain == 0);
% idxToneTrials = find(toneGain ~= 0);
% 
% ind1 = [];
% ind2 = [];
% ind3 = [];

% tone_end = end_indices(idxToneTrials(end));
% noTone_start = start_indices(idxNoToneTrials(find(diff(idxNoToneTrials)>1)+1));
% 
% toneTrials = any(idxToneTrials(:) == find(start_indices)); 
% 
% for j = 1:(length(idxToneTrials) + length(idxNoToneTrials))
%     % Tone trials
%     if toneTrials(j) == 1
%         if j == idxToneTrials(end)
%             break
%         end
%         ind1 = [ind1, end_indices(j)];
% 
%         indices = end_indices(j) + find(eventVar.licksPortsAllSingle(...
%             (end_indices(j)+1):(start_indices(j+1)-1)) > 1);
%         ind2 = [ind2, indices];
%     
%     % No-tone trials
%     else
%         if j == idxNoToneTrials(end) | j == idxToneTrials(1)
%             break
%         end
%         indices = end_indices(j) + find(eventVar.licksPortsAllSingle(...
%             (end_indices(j)+1):(start_indices(j+1)-1)) > 1);
%         ind3 = [ind3, indices];
%     end
% end
% 
% ind2 = [ind2, tone_end + find(eventVar.licksPortsAllSingle(...
%     tone_end+1:noTone_start-1) > 1)];
% 
% ind3 = [ind3, end_indices(end) + find(eventVar.licksPortsAllSingle(...
%     end_indices(end)+1:end) > 1)];

% 1. All tone trial choice licks
ind1 = find(eventVar.licksPortsAllSingle == (logVar.currChoice+1) & ...
    eventVar.licksPortsAllSingle ~= 0 & logVar.toneOn==1); % currChoice goes from 1-6, while licks are from port 1 (home) to 7. add 1

% 2. All tone trial spontaneous licks
ind2 = find(eventVar.licksPortsAllSingle ~= (logVar.currChoice+1) & ...
    eventVar.licksPortsAllSingle ~= 0 & logVar.toneOn==1);

% 3. All no tone trial end licks
ind3 = find(eventVar.licksPortsAllSingle == 7 & logVar.toneOn==0);

% 4. All no tone trial middle port licks
ind4 = find(eventVar.licksPortsAllSingle ~= 7 & ...
    eventVar.licksPortsAllSingle ~= 1 & ...
    eventVar.licksPortsAllSingle ~= 0 & logVar.toneOn==0);

% 5. All home ports
ind5 = find(eventVar.licksPortsAllSingle == 1);

eventVar.licksType = zeros(size(eventVar.licksAll));
eventVar.licksType(ind1) = 1;
eventVar.licksType(ind2) = 2;
eventVar.licksType(ind3) = 3;
eventVar.licksType(ind4) = 4;
eventVar.licksType(ind5) = 5;


%% (10) Generate LFP data 
lfp = bz_GetLFP(69,'noprompts',true);%'basepath','I:\Videos\IZ43_220828_sess4'
signal = bz_Filter(double(lfp.data),'filter','butter','passband',[6 12],'order',3);
downFactor = floor(length(signal)/length(timestamps));
lfpdown = downsample(signal,downFactor);
lfpdown = lfpdown(1:length(timestamps));


%% Plot distributions and time series of the variables. 
% Variables of interest are: eventVar.lickType, eventVar.licks, constVar.cyclicYLin, 
% constVar.cyclicY, constVar.distStop, constVar.x, lfp (theta)
all_variables = struct('x', [], 'cyclicYLin', [], 'cyclicY', [], 'distStop', ...
    [], 'licksType', []);
all_variables.x = constVar.x;
all_variables.cyclicYLin = constVar.cyclicYLin;
all_variables.cyclicY = constVar.cyclicY;
all_variables.distStop = constVar.distStop;
all_variables.licks = eventVar.licksPortsAllSingle;
all_variables.licksType = eventVar.licksType;
all_variables.lfp = lfpdown';

save('sessionData.mat','eventVar','constVar','timestamps','spkMat',...
    'tracking','lfpdown','all_variables'); 

fields = fieldnames(all_variables);

figure;

nc = 3;
nr = ceil(length(fields)/3);
np = 1;

for i = 1:length(fields)
    subplot(nr, nc, np);

    fieldname = fields{i};
    field = all_variables.(fieldname);

    if strcmp(fieldname, 'licksType') || strcmp(fieldname, 'licks') 
        field = nonzeros(field);
    end
    
    histogram(field);    
    xlabel(fieldname);
    ylabel('density');

    np = np + 1;
end
% 
% saveas(gcf, fullfile(basepath, 'Variables distribution'), '.fig');
% saveas(gcf, fullfile(basepath, 'Variables distribution'), '.png');

figure;
set(gcf,'Color','w')
nr = length(fields);

for i = 1:length(fields)
    subplot(nr, 1, i);

    fieldname = fields{i};
    field = all_variables.(fieldname);
    plot(timestamps,field);

    xlabel('time (ms)');
    ylabel(field);
    if strcmp(fieldname, 'licksType')
        title(strcat(fieldname,': 1-Tone choice,2-Tone spont., 3-no tone choice, 4- no-tone spont., 5- home'))
    else
        title(fieldname)
    end
    xlim([timestamps(1) timestamps(end)])
    axis off
end

% saveas(gcf, fullfile(basepath, 'Variables time series'), '.fig');
% saveas(gcf, fullfile(basepath, 'Variables time series'), '.png');
end