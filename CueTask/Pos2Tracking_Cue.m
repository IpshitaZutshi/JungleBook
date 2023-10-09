
function [tracking] = Pos2Tracking_Cue(varargin)
% Get position tracking
%
% USAGE
%
%   [behavior] = Pos2Tracking_Cue(varagin)
%
% INPUTS
%   excelfile      Look for it in the basepath
%   It requires an excel file (deeplabcut Output) and a digitalin.dat file in the
%   basepath folder.
% 
%   (OPTIONAL)
%   basePath       -(default: pwd) basePath for the recording file, in
%                   buzcode format.
%   convFact       - Spatial conversion factor (cm/px). If not provide,
%                   normalize maze size.
%   saveFrames     - Creates mat file containin all frames (default false).
%   forceReload    - default false.
%   bazlerTTL      - Rx1 with timestamps from bazler ttl pulses. If not
%                   provided try to extract ttl pulses from digitalIn
%                   channel 1. If it fails, gives video time.
%   saveMat        - default true
%   artifactThreshold - max allow movements per frame (in cm, default 3).
%                   Disabled if not convFact is provided.
% 
% OUTPUT
%       - tracking.behaviour output structure, with the fields:
%   x               - x position in cm/ normalize
%   y               - y position in cm/ normalize
%   timestamps      - in seconds, if Bazler ttl detected, sync by them
%   folder          - 
%   sync.sync       - Rx1 LED luminance.
%   sync.timestamps - 2xC with start stops of sync LED.
%   samplingRate    - in Hz
%   averageFrame    - 
%   
%
%   Manu Valero 2019
% TO DO: Generalize for non-T maze behaviour
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'convFact',[],@isnumeric); % 0.1149
addParameter(p,'roiLED',[],@ismatrix);
addParameter(p,'ledSync',true,@islogical);
addParameter(p,'forceReload',false,@islogical)
addParameter(p,'saveFrames',true,@islogical)
addParameter(p,'thresh',220.^2,@isnumeric)
addParameter(p,'bazlerTTL',[],@isnumeric)
addParameter(p,'saveMat',true,@islogical);
% addParameter(p,'RGBChannel',[],@isstr);

parse(p,varargin{:});
basepath = p.Results.basepath;
fs = p.Results.fs;
artifactThreshold = p.Results.artifactThreshold;
roiLED = p.Results.roiLED;
convFact = p.Results.convFact;
forceReload = p.Results.forceReload;
ledSync = p.Results.ledSync;
saveFrames = p.Results.saveFrames;
verbose = p.Results.verbose;
thresh = p.Results.thresh;
bazlerTtl = p.Results.bazlerTTL;
saveMat = p.Results.saveMat;
% RGBChannel = p.Results.RGBChannel;

%% Deal with inputs
if ~isempty(dir([basepath filesep '*Tracking.Behavior.mat'])) && ~forceReload
    disp('Trajectory already detected! Loading file.');
    file = dir([basepath filesep '*Tracking.Behavior.mat']);
    load(file.name);
    return
end

    if ~isempty(dir([basepath filesep '*avi']))
        aviFile = dir([basepath filesep '*avi']); 
        aviFile = erase(aviFile.name,'.avi');
    else
        warning('No video file!!');
        tracking = [];
        return
    end
    
    if ~isempty(dir([basepath filesep 'CueTask*csv']))
        csvFile = dir([basepath filesep 'CueTask*csv']); 
        csvTS = load(csvFile.name);
    else
        warning('No csv file!!');
        tracking = [];
        return
    end
end

if ~isempty(dir([basepath filesep '*TrialBehavior.Events.mat']))
    file = dir([basepath filesep '*TrialBehavior.Events.mat']);
    load(file.name)
end

% xMaze = [0 size(frames,2) * convFact];
% yMaze = [0 size(frames,1) * convFact];

%% Figure out LED for sync
if ledSync    
    disp('Get average frame...');
    videoObj = VideoReader([aviFile '.avi']);   % get video
    numFrames = get(videoObj, 'NumFrames');
    frames = [];
    tic
    f = waitbar(0,'Getting frames...');
    for ii = 1:numFrames
        waitbar(ii/numFrames,f)
        temp_frames = read(videoObj,ii);        % get all frames
        frames_red(:,:,ii) = temp_frames(:,:,1);      % convert to grayscale
    end
    close(f)
    toc      
end

disp('Draw ROI for LED...');
h1 = figure;
idx = ONidx;
imshow(frames_red(:,:,idx))

roi = drawpolygon;
roiLED = [roi.Position; roi.Position(1,:)];

%% detect LED pulses for sync
if ~isempty(roiLED)
    disp('Detect LED for sync...');
    bwLED = uint8(poly2mask(roiLED(:,1),roiLED(:,2),size(frames_red,1),size(frames_red,2)));
    parfor ii = 1:size(frames_red,3)
        fr = double(frames_red(:,:,ii).*bwLED);
        fr(fr==0) = NaN;
        sync(ii) = nanmedian(fr(:)); 
    end
    
    sync2 = sync.^2;
    syncBin = (sync2>thresh); % binarize signal
    locsA = find(diff(syncBin)==1); % start of pulses
    locsB = find(diff(syncBin)==-1); % end of pulses
    pul = locsA(1:min([length(locsA) length(locsB)]));
    for ii = 1 : size(pul,2) % pair begining and end of the pulse
        if sum(locsB > pul(1,ii)) > 0
            pul(2,ii) =  locsB(find(locsB - pul(1,ii) ==...
                min(locsB(locsB > pul(1,ii)) - pul(1,ii))));
        else
            pul(2,ii) = nan;
        end
    end
    % if a jump happened only for 1-2 frames, ignore
    a = pul(2,:)-pul(1,:);
    
else
    sync = []; pul = [];
end

if ~isempty(dir('*digitalIn.events.mat'))
    file = dir('*digitalIn.events.mat');
    load(file.name)
end

%% Syncs should align with base solenoid
timestampsEnd = digitalIn.timestampsOn{10}(1:end);
sync_ts  = (sync_ts-sync_ts(1))+timestampsEnd(1);

%% Align csvTS
csvTS_corr = (csvTS-csvTS(sync_fin(1)))./1000 + timestampsEnd(1);

%% Add tracking time data
tracking.position.x = [];
tracking.position.y = [];
tracking.logical.csvTS_corr = csvTS_corr;
tracking.logical.end = timestampsEnd;
tracking.timestamps = csvTS_corr(csvTS_corr>0 & csvTS_corr<timestampsEnd(end));

intervals = [behavTrials.timestamps(1:(end-1),1) behavTrials.timestamps(2:(end),1)]
[status,interval,index] = InIntervals(tracking.timestamps,intervals);
tracking.masks.trial = interval;
tracking.masks.trial(tracking.masks.trial==0) = nan;

tracking.masks.forward(1:length(tracking.timestamps)) = nan;
intervals = [behavTrials.timestamps(:,1) behavTrials.choiceTS];
[status,interval,index] = InIntervals(tracking.timestamps,intervals);
tracking.masks.forward(status==1) = 1;

tracking.masks.choice(1:length(tracking.timestamps)) = nan
lTrials = find(behavTrials.choice==1);
Lia = ismember(tracking.masks.trial,lTrials);
tracking.masks.choice(Lia==1) = 1;
lTrials = find(behavTrials.choice==0);
Lia = ismember(tracking.masks.trial,lTrials);
tracking.masks.choice(Lia==1) = 0;

tracking.masks.cue(1:length(tracking.timestamps)) = nan
lTrials = find(behavTrials.cue==1);
Lia = ismember(tracking.masks.trial,lTrials);
tracking.masks.cue(Lia==1) = 1;
lTrials = find(behavTrials.cue==0);
Lia = ismember(tracking.masks.trial,lTrials);
tracking.masks.cue(Lia==1) = 0;

tracking.masks.stim(1:length(tracking.timestamps)) = nan
lTrials = find(behavTrials.stim==1);
Lia = ismember(tracking.masks.trial,lTrials);
tracking.masks.stim(Lia==1) = 1;
lTrials = find(behavTrials.stim==0);
Lia = ismember(tracking.masks.trial,lTrials);
tracking.masks.stim(Lia==1) = 0;

if saveMat
    save([basepath filesep fbasename '.Tracking.Behavior.mat'],'tracking');
end

end

figure
set(gcf,'Renderer','painters')
subplot(1,3,1)
hold on
for ii = 2:(length(behavTrials.cue)-1)
    if (behavTrials.stim(ii)==0 && behavTrials.choice(ii) == 0 && behavTrials.correct(ii) == 1)        
        idx = tracking.position.conf'>0.9 & tracking.masks.forward == 1 & tracking.masks.trial' == ii;
        plot(tracking.position.y(idx),tracking.position.x(idx),'k')
    elseif (behavTrials.stim(ii)==0 && behavTrials.choice(ii) == 1 && behavTrials.correct(ii) == 1)
        idx = tracking.position.conf'>0.9 & tracking.masks.forward == 1 & tracking.masks.trial' == ii;
        plot(tracking.position.y(idx),tracking.position.x(idx),'b')    
    end
end
ylim([80 650])

subplot(1,3,2)
hold on
for ii = 2:(length(behavTrials.cue)-1)
    if (behavTrials.stim(ii)==0 && behavTrials.correct(ii) == 1)       
        idx = tracking.position.conf'>0.9 & tracking.masks.forward == 1 & tracking.masks.trial' == ii;
        plot(tracking.position.y(idx),tracking.position.x(idx),'k')
    elseif (behavTrials.stim(ii)==1 && behavTrials.correct(ii) == 1)
        idx = tracking.position.conf'>0.9 & tracking.masks.forward == 1 & tracking.masks.trial' == ii;
        plot(tracking.position.y(idx),tracking.position.x(idx),'b')    
    end
end
ylim([80 650])

subplot(1,3,3)
hold on
for ii = 2:(length(behavTrials.cue)-1)
    if (behavTrials.stim(ii)==0 && behavTrials.correct(ii) == 1)       
        idx = tracking.position.conf'>0.9 & tracking.masks.forward == 1 & tracking.masks.trial' == ii;
        plot(tracking.position.y(idx),tracking.position.x(idx),'k')
    elseif (behavTrials.stim(ii)==1 && behavTrials.correct(ii) == 1)
        idx = tracking.position.conf'>0.9 & tracking.masks.forward == 1 & tracking.masks.trial' == ii;
        plot(tracking.position.y(idx),tracking.position.x(idx),'b')    
    end
end
ylim([88 160])

% maze is 48 x 67 (exterior wall to exterior wall). Virtual maze is 580 x
% 420 pixels. Conv factor is ~ 0.1143 - 0.1155 cm/pix 
