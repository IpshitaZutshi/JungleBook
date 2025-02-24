function preprocess_manifold

SPIKEbin_beh = 0.1;

file = dir('*.spikes.cellinfo.mat');
load(file.name)

file = dir('*.cell_metrics.cellinfo.mat');
load(file.name)

file = dir('*.session.mat');
load(file.name)

file = dir('*.Tracking.Behavior.mat');
load(file.name)

file = dir('*.TrialBehavior.Behavior.mat');
load(file.name)

beh_interval = [behavTrials.timestamps(1,1) behavTrials.timestamps(end,2)];

SPIKEMAT = bz_SpktToSpkmat_manifold(spikes, 'dt',SPIKEbin_beh,'win',beh_interval,'units','counts');
timestamp = [];
timestamp = [timestamp, SPIKEMAT.timestamps'];

%% tracking data
pos_x = smoothdata(tracking.position.x(InIntervals(tracking.timestamps,beh_interval)),'movmean',5);
pos_y = smoothdata(tracking.position.y(InIntervals(tracking.timestamps,beh_interval)),'movmean',5);

%% Make masks for the behavior
trialMask = zeros(size(tracking.timestamps));
cueMask = zeros(size(tracking.timestamps));
correctMask = zeros(size(tracking.timestamps));
choiceMask = zeros(size(tracking.timestamps));
stimMask = zeros(size(tracking.timestamps));
directionMask = zeros(size(tracking.timestamps));

for ii = 1:length(behavTrials.cue)
    posTrials = tracking.timestamps >= behavTrials.timestamps(ii,1) & tracking.timestamps <= behavTrials.timestamps(ii,2);
    trialMask(posTrials) = ii;
   
    cueMask(posTrials) = behavTrials.cue(ii)+1;
    correctMask(posTrials) = behavTrials.correct(ii)+1;
    choiceMask(posTrials) = behavTrials.choice(ii)+1;
    stimMask(posTrials) = behavTrials.stim(ii)+1;
    posTrials = tracking.timestamps >= behavTrials.timestamps(ii,1) & tracking.timestamps <= behavTrials.choiceTS(ii);
    directionMask(posTrials) = 1;

    posTrials = tracking.timestamps > behavTrials.choiceTS(ii) & tracking.timestamps <= behavTrials.timestamps(ii,2);
    directionMask(posTrials) = 2;
end


trialMask = trialMask(InIntervals(tracking.timestamps,beh_interval));
cueMask = cueMask(InIntervals(tracking.timestamps,beh_interval));
correctMask = correctMask(InIntervals(tracking.timestamps,beh_interval));
choiceMask = choiceMask(InIntervals(tracking.timestamps,beh_interval));
stimMask = stimMask(InIntervals(tracking.timestamps,beh_interval));
directionMask = directionMask(InIntervals(tracking.timestamps,beh_interval));


% interpolate
trial_num_ds = interp1(tracking.timestamps(InIntervals(tracking.timestamps,beh_interval)),trialMask,timestamp,'nearest');
correct_ds = interp1(tracking.timestamps(InIntervals(tracking.timestamps,beh_interval)),correctMask,timestamp,'nearest');
choice_ds = interp1(tracking.timestamps(InIntervals(tracking.timestamps,beh_interval)),choiceMask,timestamp,'nearest');
cue_ds = interp1(tracking.timestamps(InIntervals(tracking.timestamps,beh_interval)),cueMask,timestamp,'nearest');
stim_ds = interp1(tracking.timestamps(InIntervals(tracking.timestamps,beh_interval)),stimMask,timestamp,'nearest');
direction_ds = interp1(tracking.timestamps(InIntervals(tracking.timestamps,beh_interval)),directionMask,timestamp,'nearest');

pos_x_ds = interp1(tracking.timestamps(InIntervals(tracking.timestamps,beh_interval)),pos_x,timestamp,'linear'); 
pos_y_ds = interp1(tracking.timestamps(InIntervals(tracking.timestamps,beh_interval)),pos_y,timestamp,'linear'); 
timestamp_beh = timestamp;

%% Separate CA1 and V1 cells
spike_count_beh = [];
spike_count_beh = [spike_count_beh,SPIKEMAT.data(:,cell_metrics.brainRegion==0)'];

spike_counts = spike_count_beh;
data = normalize(spike_counts,1,'zscore');
data(isnan(data))=0;
data = data';
save('T17_240421_sess6_V1.data.mat','data','timestamp')

spike_count_beh = [];
spike_count_beh = [spike_count_beh,SPIKEMAT.data(:,cell_metrics.brainRegion==1)'];
spike_counts = spike_count_beh;
data = normalize(spike_counts,1,'zscore');
data(isnan(data))=0;
data = data';
save('T17_240421_sess6_CA1.data.mat','data','timestamp')

save('Umap_behavior.mat','timestamp_beh','trial_num_ds','correct_ds',...
    'choice_ds','cue_ds','choice_ds','pos_x_ds','pos_y_ds','direction_ds','stim_ds');


end