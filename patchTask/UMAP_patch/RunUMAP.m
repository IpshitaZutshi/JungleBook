
basepath = pwd;

%% 1. Plot behavior
getPatchBehavior()

%patchStatistics()

%% 2. Run Kilosort
AnalysisPreprocessing()
% Manual curation with Phy before running Step 2

%% 3. LFP, sleep score, ripples, spikes
AnalysisBatchPatch()

%% 4. Get Tracking and Behavior
getPatchTracking() %This function returns the tracking and runs getPatchBehavior
% syncs intan, tracking, and photometry together

%% 5. Load spikes
kilosort_folder = dir(fullfile(basepath, 'Kilosort*'));
cd(fullfile(basepath, kilosort_folder.name));
[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true); sessionInfo.rates.lfp = 1250;  save(strcat(sessionInfo.session.name,'.sessionInfo.mat'),'sessionInfo');
spikes = loadSpikes('getWaveformsFromDat', false);

%% 6. Run preprocessing for the manifold
cd(basepath)
preprocess_manifold_patch()

%% 7. Run Umap with python
%  Run umap_patch code

%% 8. Plot UMAP
plotPatchUMAP()
