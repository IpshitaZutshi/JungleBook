function [tracking] = getCueTracking(varargin)

% Gets position tracking for each sub-session and concatenate all of them so they are 
% aligned with LFP and spikes. Default is using the analog input from Bonsai detected position. 
%
%% USAGE
%
%   [tracking] = getToneTracking(varargin)
%
% INPUTS
%   basePath        -(default: pwd) basePath for the recording file, in buzcode format:
%   analogInputPos  - Analog channel that has tracking information.
%   analogInputTone - Analog input that has tone information.
%   fs              - sampling rate for behavior. default 1250
%   trackImgLength  - Distance of track in the video file, default, 410.
%   trackLength     - Actual length of the track (in cm)
%   freqRange       - Frequency range of the Tone, default 1000 to 22000
%   saveMat         - default true
%   forceReload     - default false
%
% OUTPUT
%       - tracking.behaviour output structure, with the fields:
%   position.x               - x position in cm/ normalize
%   position.tone
%   timestamps      - in seconds, if Basler ttl detected, sync by them

%   HISTORY:
%     - Ipshita Zutshi 2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'forceReload',false,@islogical)

parse(p,varargin{:});
basepath = p.Results.basepath;
saveMat = p.Results.saveMat;
forceReload = p.Results.forceReload;

%% In case tracking already exists 
if ~isempty(dir([basepath filesep '*Tracking.Behavior.mat'])) && ~forceReload
    disp('Trajectory already detected! Loading file.');
    file = dir([basepath filesep '*Tracking.Behavior.mat']);
    load(file.name);
    return
end

 %% Find subfolder recordings
cd(basepath);
[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);
%Load concatenated digitalIn
dig_np = getDigitalIn_Neuropixels;

if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
    load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
    count = 1;
    for ii = 1:size(MergePoints.foldernames,2)
         if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep 'Top*']))   
            cd([basepath filesep MergePoints.foldernames{ii}]); 
            fprintf('Computing tracking in %s folder \n',MergePoints.foldernames{ii});            
            tempTracking{count}= Pos2Tracking_Cue('forceReload',forceReload); % computing trajectory

            items = dir;

            % Step 2: Filter for directories whose names start with 'T'
            folders = items([items.isdir] & startsWith({items.name}, 'T'));
            folderNames = {folders.name};
            cd(folderNames{1})
            dig_Intan = bz_getDigitalIn;

            fprintf('Computing behavior in %s folder \n',folderNames{1});  
            tempBehav{count} = getCueBehavior_newSetup('forceRun',forceReload);
            trackFolder(count) = ii;             

            %% Adjust the tracking and behavior timestamps to match neuropixels (i.e., the merge points timestamps)
            ts_npx = dig_np.timestampsOn(dig_np.timestampsOn>MergePoints.timestamps(ii,1) & dig_np.timestampsOn<MergePoints.timestamps(ii,2));
            
            master_timestamps = ts_npx - MergePoints.timestamps(ii,1);
            device_timestamps = dig_Intan.timestampsOn{13};

            diff_npx = diff(master_timestamps);
            diff_intan = diff(device_timestamps);

            [cross_corr, lags] = xcorr(diff_npx, diff_intan);
            [~,idx] = max(cross_corr);
            shift = abs(lags(idx))+1;

            % So ts_npx(1) = ts_intan(shift);
            tsdiff = device_timestamps(shift)-master_timestamps(1);

            %% Remove that much time from the digital and tracking timestamps;

            tempTracking{count}.originalTimestamps = tempTracking{count}.timestamps;
            tempTracking{count}.timestamps = tempTracking{count}.timestamps-tsdiff;
            % Save positive idx 
            posIdx = tempTracking{count}.timestamps>0;
            tempTracking{count}.timestamps = tempTracking{count}.timestamps(posIdx);
            tempTracking{count}.position.x = tempTracking{count}.position.x(posIdx);
            tempTracking{count}.position.y = tempTracking{count}.position.y(posIdx);
            tempTracking{count}.position.vx = tempTracking{count}.position.vx(posIdx);
            tempTracking{count}.position.vy = tempTracking{count}.position.vy(posIdx);
            tempTracking{count}.position.v = tempTracking{count}.position.v(posIdx);
            tempTracking{count}.tsshift = tsdiff;

            % Also behavior
            tempBehav{count}.timestamps = tempBehav{count}.timestamps-tsdiff;
            tempBehav{count}.choiceTS = tempBehav{count}.choiceTS-tsdiff;
            tempBehav{count}.tsshift = tsdiff;

            % For a sanity check, visually confirm that the arduino barcode
            % pulses are aligned
            time_base = 0:1:round(max([master_timestamps; device_timestamps]*1000)); % 1ms resolution
            device_timestamps2 = device_timestamps-tsdiff;
            device_timestamps2 = device_timestamps2(device_timestamps2>0);
            
            npx_event_series = ismember(time_base, round(master_timestamps*1000));
            intan_event_series = ismember(time_base, round(device_timestamps2*1000));
            h = figure
            subplot(2,2,[1 2])            
            plot(npx_event_series)
            hold on
            plot(intan_event_series)

            subplot(2,2,3)
            plot(npx_event_series(1:100000))
            hold on
            plot(intan_event_series(1:100000))

            subplot(2,2,4)
            plot(npx_event_series(end-100000:end))
            hold on
            plot(intan_event_series(end-100000:end))
            meanmismatch = mean(master_timestamps-device_timestamps2(1:length(master_timestamps)));
            title(strcat('avg mismatch = ',num2str(meanmismatch)))
            cd(basepath)
            saveas(h,'Pulses\barcodeAlign.png');
            count = count + 1;
        end
    end
    cd(basepath);
else
    error('missing MergePoints, quitting...');
end

%% Concatenate and sync timestamps
ts = []; subSessions = []; maskSessions = [];
tsBehav= []; maskSessionsBehav = []; tsBehavCue = [];
if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
    load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
    for ii = 1:length(trackFolder)
        if strcmpi(MergePoints.foldernames{trackFolder(ii)},tempTracking{ii}.folder)
            sumTs = tempTracking{ii}.timestamps + MergePoints.timestamps(trackFolder(ii),1);
            subSessions = [subSessions; MergePoints.timestamps(trackFolder(ii),1:2)];
            maskSessions = [maskSessions; ones(size(sumTs))*ii];
            ts = [ts; sumTs];

            sumTs = tempBehav{ii}.timestamps + MergePoints.timestamps(trackFolder(ii),1);
            sumTsCue = tempBehav{ii}.choiceTS + MergePoints.timestamps(trackFolder(ii),1);
            maskSessionsBehav = [maskSessionsBehav; ones(size(sumTs))*ii];
            tsBehav = [tsBehav; sumTs];
            tsBehavCue = [tsBehavCue;sumTsCue];
            
        else
            error('Folders name does not match!!');
        end
    end
else
    warning('No MergePoints file found. Concatenating timestamps...');
    for ii = 1:length(trackFolder)
        sumTs = max(ts)+ tempTracking{ii}.timestamps;
        subSessions = [subSessions; [sumTs(1) sumTs(end)]];
        ts = [ts; sumTs];
    end
end

% Concatenating tracking fields...
x = []; y = []; vx = []; vy = []; v = []; folder = []; samplingRate = []; description = [];framesDropped = [];
for ii = 1:length(tempTracking)
    x = [x; tempTracking{ii}.position.x];     
    y = [y; tempTracking{ii}.position.y]; 

    vx = [vx; tempTracking{ii}.position.vx];     
    vy = [vy; tempTracking{ii}.position.vy]; 
    v = [v; tempTracking{ii}.position.v];     

    folder{ii} = tempTracking{ii}.folder; 
    tshift{ii} = tempTracking{ii}.tsshift;
    samplingRate = [samplingRate; tempTracking{ii}.samplingRate];  
    description{ii} = tempTracking{ii}.description;  
    framesDropped{ii} = tempTracking{ii}.framesDropped;  
end

tracking.position.x = x;
tracking.position.y = y;
tracking.position.vx = vx;
tracking.position.vy = vy;
tracking.position.v = v;
tracking.tsShift = tshift;

tracking.folders = folder;
tracking.samplingRate = samplingRate;
tracking.timestamps = ts;
tracking.framesDropped = framesDropped;
tracking.events.subSessions =  subSessions;
tracking.events.subSessionsMask = maskSessions;

behavTrials = tempBehav{1};
behavTrials.timestamps = tsBehav;
behavTrials.choiceTS = tsBehavCue;

%% save tracking
cd(basepath);
[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
if saveMat
    save([basepath filesep sessionInfo.FileName '.Tracking.Behavior.mat'],'tracking');
    save([basepath filesep sessionInfo.FileName '.TrialBehavior.Behavior.mat'],'behavTrials');
end

% figure
% scatter(tracking.position.x,tracking.position.y,2,[0.9 0.9 0.9])
% hold on
% for ii = 1:length(behavTrials.choiceTS)
%     [~,idx1] = min(abs(tracking.timestamps-behavTrials.choiceTS(ii)));
%     [~,idx2] = min(abs(tracking.timestamps-behavTrials.timestamps(ii,1)));
%     [~,idx3] = min(abs(tracking.timestamps-behavTrials.timestamps(ii,2)));
% 
%     scatter(tracking.position.x(idx1),tracking.position.y(idx1),3,[1 0 0])
%     scatter(tracking.position.x(idx2),tracking.position.y(idx2),3,[0 1 0])
%     scatter(tracking.position.x(idx3),tracking.position.y(idx3),3,[0 0 1])
% end
end

