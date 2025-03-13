
% Preprocessing and folder summary 
% 
% 1. extract LFP
% 2. extract spikes
% 3. extract SWR
% 4. sleep score
% 5. cell Explorer
% 6. extract behavior/ tracking?
%clear 

disp('Get components...');

%% 1. Extract LFP
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

%% 2. Extract spikes
spikes = loadSpikes('getWaveformsFromDat', false);

%% 3. Extract sharp wave ripples
pyrCh = 7;
noiseCh = 111;
[ripples] = bz_FindRipples(pwd,pyrCh,'noise',noiseCh,'savemat',true);

%% 4. Sleep score
badChannels = [24:38 48:63]; %N7
% badChannels = [0:3 15:18 21:30 43 50 95 97]; %N9
SleepScoreMaster(pwd,'stickytrigger',true,'rejectChannels',badChannels); % try to sleep score


%% 5. Cell Explorer
cell_metrics = ProcessCellMetrics('manualAdjustMonoSyn',false,'forceReload',true,'submitToDatabase',false,'showGUI',false);

