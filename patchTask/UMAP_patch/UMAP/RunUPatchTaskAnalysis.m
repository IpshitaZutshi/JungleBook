
basepath = pwd;

%% 1. Plot behavior
getPatchBehavior()

%% 2. Run Kilosort
AnalysisPreprocessing()
% Manual curation with Phy before running Step 2

%% 3. LFP, sleep score, ripples, spikes
AnalysisBatchPatch()

%% Photometry processing
% first, use python script iso_processing.py on each ppd file from the session
% (or processing.py if there are two fluorophores, i.e. not isosbestic)

%% 4. Get Tracking and Behavior
getPatchTracking() 

% This function returns the tracking 
% syncs intan, tracking, and photometry together

%% 5. Sleep photometry
sleepPhotometry()
% aligns photometry from sleep session to intan data. also analyzes
% photometry from sleep sessions at times of ripples

%% 6. Photometry during behavior
photometry_behav.m

% aligns photometry from behavior session to intan data. also analyzes
% dopamine levels around licks and different types of events

%% 7. Analyzing photometry across sessions
averageSleepPhotometryAcrossSessions.m
averageBehaviorPhotometryAcrossSessions.m

% these scripts require saving the photometry matrices around specified
% events generated in photometry_behav.m or sleepPhotometry.m as variables
% which you can then load into these scripts.

%% UMAP
%% 8. Load spikes
kilosort_folder = dir(fullfile(basepath, 'Kilosort*'));
cd(fullfile(basepath, kilosort_folder.name));
[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true); sessionInfo.rates.lfp = 1250;  save(strcat(sessionInfo.session.name,'.sessionInfo.mat'),'sessionInfo');
spikes = loadSpikes('getWaveformsFromDat', false);

%% 9. Run preprocessing for the manifold
cd(basepath)
preprocess_manifold_patch()

%% 10. Run Umap with python
%  Run umap_patch code

%% 11. Plot UMAP
plotPatchUMAP()

















