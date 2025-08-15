
% Preprocessing and folder summary 
% 
% 1. extract behavior/ tracking?
% 2. extract LFP
% 3. extract spikes
% 4. extract SWR
% 5. sleep score
% 6. cell Explorer

disp('Get components...');

%% 1. Patch Behavior
getPatchTracking('basePath',pwd)

%% 2. Extract LFP
[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true); sessionInfo.rates.lfp = 1250;  save(strcat(sessionInfo.session.name,'.sessionInfo.mat'),'sessionInfo');
if isempty(dir('*.lfp'))
    try 
        bz_LFPfromDat(pwd,'outFs',1250); % generating lfp
    catch
        disp('Problems with bz_LFPfromDat, resampling...');
        ResampleBinary(strcat(sessionInfo.session.name,'.dat'),strcat(sessionInfo.session.name,'.lfp'),...
            sessionInfo.nChannels,1,sessionInfo.rates.wideband/sessionInfo.rates.lfp);
    end
end

%% 3. Extract spikes
spikes = loadSpikes('getWaveformsFromDat', false);

%% 4. Extract sharp wave ripples
pyrCh = 121; % 75 for n11 115 67 for n11
%67 90 66
noiseCh = 111;
[ripples] = bz_FindRipples(pwd,pyrCh,'noise',noiseCh,'savemat',true,'durations',[30 100],'passband',[130 200]);

%% 5. Sleep score
%badChannels = [24:38 48:63]; %N7
badChannels = [0:3 15:18 21:30 43 50 95 97]; %N9
%badChannels = [0:3 15:18 21:30 41 43 46 47 50 52 95 97]; %N15
badChannels = [42 48 56:59 61 70:72]; % N17
SleepScoreMaster(pwd,'stickytrigger',true,'rejectChannels',badChannels); % try to sleep score

%% 6. Cell Explorer
cell_metrics = ProcessCellMetrics('manualAdjustMonoSyn',false,'forceReload',true,'submitToDatabase',false,'showGUI',false);


