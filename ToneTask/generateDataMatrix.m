% Variables of interest- 
% Velocity
% X ? y location
% Tone/no tone
% Reward location (timestamp of delivery)
% Velocity threshold
% Tone frequency
% Reward history?
% Firing rate of other cells?
% 
% Spike counts 6 ms
% 
% Spike matrix

function [spkMat,constVar,logVar,eventVar,timestamps] = generateDataMatrix(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'dt',0.006,@isnumeric);

parse(p,varargin{:});
basepath = p.Results.basepath;
dt = p.Results.dt;

%% Deal with inputs
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
file = dir([basepath filesep '*.MergePoints.events.mat']);
load(file.name); 

% Load digitalIn
%filename = [MergePoints.foldernames{2} filesep 'amplifier.DigitalIn.events.mat'];
filename = [sessionInfo.session.name '.DigitalIn.events.mat'];
load(filename)
%% First generate spike matrix 
win = [behavTrials.timestamps(1,1) behavTrials.timestamps((end-1),2)];
spkData = bz_SpktToSpkmat(spikes,'dt',dt, 'win',win);

spkMat = spkData.data';
timestamps = spkData.timestamps';

%% Next, generate matrix for constant variables x, y, v

[~,idxStart] = min(abs(tracking.timestamps-win(1)));
[~,idxEnd] = min(abs(tracking.timestamps-win(2)));

% Interpolate
timeS = tracking.timestamps(idxStart:idxEnd);
constVar.x = interp1(timeS,tracking.position.x(idxStart:idxEnd),spkData.timestamps)';
constVar.y = interp1(timeS,tracking.position.y(idxStart:idxEnd),spkData.timestamps)';
constVar.vel = interp1(timeS,tracking.position.v(idxStart:idxEnd),spkData.timestamps)';
constVar.velY = interp1(timeS,tracking.position.vy(idxStart:idxEnd),spkData.timestamps)';

% Also generate separate y variables for forward versus reverse directions.
constVar.yFwd = constVar.y;
constVar.yRev = constVar.y;
constVar.yFwd(constVar.velY<-1) = nan;
constVar.yLin = constVar.yFwd;
constVar.yRev(constVar.velY>0) = nan;

%% Find index in matrix to start and end points of each trial
for tt = 1:size(behavTrials.timestamps,1)
    [~,startIdx(tt)] = min(abs(spkData.timestamps - behavTrials.timestamps(tt,1))); 
    [~,endIdx(tt)] = min(abs(spkData.timestamps - behavTrials.timestamps(tt,2))); 
end

%% Generate logical variables - 
% trial number, forward/reverse, toneOn/Off, trial gain, current choice, 
% correct/incorrect, past choice, past choice correct/incorrect,

logVar.forward(1:length(spkData.timestamps)) = 0;
logVar.toneOn(1:length(spkData.timestamps)) = 0;

%last lin track trial 
lastLin = behavTrials.linTrial;
toneGain = (behavTrials.toneGain);
toneGain(~lastLin) = toneGain(~lastLin)+1;
toneGain(lastLin==1) = 0;
choiceCorrect = behavTrials.correct;
choiceCorrect(lastLin==1) = NaN;

for tt = 1:(length(startIdx)-1)
    logVar.trialNum(startIdx(tt):startIdx(tt+1)) = tt;
    logVar.forward(startIdx(tt):endIdx(tt)) = 1;
    logVar.toneOn(startIdx(tt):endIdx(tt)) = ~behavTrials.linTrial(tt);
    logVar.trialType(startIdx(tt):startIdx(tt+1)) = toneGain(tt);
    logVar.currChoice(startIdx(tt):startIdx(tt+1)) = behavTrials.lickLoc(tt)+1;
    logVar.currCorrect(startIdx(tt):startIdx(tt+1)) = choiceCorrect(tt);
    if tt==1
        logVar.prevChoice(startIdx(tt):startIdx(tt+1)) = NaN;
        logVar.prevCorrect(startIdx(tt):startIdx(tt+1)) = NaN; 
    else
        logVar.prevChoice(startIdx(tt):startIdx(tt+1)) = behavTrials.lickLoc(tt-1)+1;
        logVar.prevCorrect(startIdx(tt):startIdx(tt+1)) = choiceCorrect(tt-1);    
    end
end

% Now can generate cont freq variable - transformed to Hz
%gain =[435/72, 435/144, 435/216, 435/390, 435/370, 435/435];
gain = [120/11.6, 120/32.27 120/55.53 120/79.62 120/102.79 120/120];
%gain = [120/10.3, 120/31.27 120/53.53 120/78.62 120/102.79 120/120];

constVar.freq(1:length(spkData.timestamps)) = nan;

freqExp = log10(22000/1000);
for ii = 1:length(logVar.trialType)
    if logVar.trialType(ii)>0 && logVar.toneOn(ii)>0
        freq = (constVar.y(ii)*gain(logVar.trialType(ii)))/120;
%         tonepos(ii) = 1000*(10.^(freqExp*freq));        
%          freq = (constVar.y(ii)*gain(logVar.trialType(ii)))/120;
         constVar.freq(ii) = ((1000*(10.^(freqExp*freq)))/22000)*120;
    end
end

constVar.yLin(logVar.trialType>0) = nan;
constVar.yFwd(logVar.trialType==0) = nan;
%% Generate matrix for events - licks, trialStart, trialEnd

lickport = [1 2 3 4 5 6 7];
for ll = 1:7
    event.times{ll} = digitalIn.timestampsOn{lickport(ll)+2};%+MergePoints.timestamps(2,1);
    event.UID(ll) = ll;
end

event.times{8} = behavTrials.timestamps((1:(end-1)),1);
event.times{9} = behavTrials.timestamps((1:(end-1)),2);

choiceIncorrect = behavTrials.correct == 0 & behavTrials.linTrial == 0;
event.times{10} = behavTrials.timestamps(choiceIncorrect(1:(end-1)),2);

choiceCorrect = behavTrials.correct == 1;
event.times{11} = behavTrials.timestamps(choiceCorrect(1:(end-1)),2);
event.UID(8:11)  = [8 9 10 11];

eventMat = bz_SpktToSpkmat(event,'dt',dt, 'win',win);
% Add the first and last trials
eventMat.data(1,8) = 1;
eventMat.data(end,9) = 1;

% add a variable that is shifted 0.5 sec previous to trialEnd
trialRamp = zeros(size(eventMat.data(:,9)));
n = round(500/6);
trialRamp(1:(end-n)) = eventMat.data(n+1:end,9)';

eventVar.licks = eventMat.data(:,1:7)';
eventVar.trialStart = eventMat.data(:,8)';
eventVar.trialEnd = eventMat.data(:,9)';
eventVar.incorrect = eventMat.data(:,10)';
eventVar.correct = eventMat.data(:,11)';
eventVar.trialRamp = trialRamp;

save('sessionData.mat','eventVar', 'constVar','timestamps','spkMat'); 
end

% 
% a = [];
% %b = [];
% figure
% for t = 1:size(behavTrials.timestamps,1)
%     [~,idx] = min(abs(tracking.timestamps-behavTrials.timestamps(t,1)));
%     %[~,idx2] = min(abs(tracking.timestamps-behavTrials.timestamps(t,2)));
%    
%     a(t,:) = tracking.position.vy(idx-120:idx+500);
%     hold on
%     plot(a(t,:),'Color',[0.5 0.5 0.5])
%     
% end