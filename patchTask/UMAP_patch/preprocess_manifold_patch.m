function preprocess_manifold_patch

SPIKEbin_beh = 0.1; % size of time bin
basepath = pwd;
[~, currentFolderName, ~] = fileparts(pwd);


file = dir(fullfile(pwd, '*', '*.spikes.cellinfo.mat'));
cd(file.folder)
load(file.name)
cd(basepath)

% file = dir('*.cell_metrics.cellinfo.mat');
% load(file.name)

file = dir('*.session.mat');
load(file.name)

file = dir('*.Tracking.Behavior.mat');
load(file.name)

file = dir('*.TrialBehavior.mat'); 
load(file.name)

file = dir('*.PhotometryBehav.mat'); 
load(file.name)

beh_interval = [tracking.timestamps(1), tracking.timestamps(end)]; 
speed_lim = 1;


SPIKEMAT = bz_SpktToSpkmat_manifold(spikes, 'dt',SPIKEbin_beh,'win',beh_interval,'units','counts');
timestamp = [];
timestamp = [timestamp, SPIKEMAT.timestamps'];

%Get lick trajectories
[trajectory_numbers_direction, trajectory_numbers_no_direction] = getTrajectories(behavTrials.port);


%% smooth tracking data
pos_x = smoothdata(tracking.position.x(InIntervals(tracking.timestamps,beh_interval)),'movmean',5);
pos_y = smoothdata(tracking.position.y(InIntervals(tracking.timestamps,beh_interval)),'movmean',5);
pos_v = smoothdata(tracking.position.v(InIntervals(tracking.timestamps,beh_interval)),'movmean',5);
pos_vy = smoothdata(tracking.position.vy(InIntervals(tracking.timestamps,beh_interval)),'movmean',5);

%% Make masks for the behavior
trialnumberMask = zeros(size(tracking.timestamps)); % total trials within session
lickedPortMask = zeros(size(tracking.timestamps)); % chosen port
outcomeMask = zeros(size(tracking.timestamps)); % rewarded or not rewarded
chosenPortProbMask = zeros(size(tracking.timestamps));
patchNumMask = zeros(size(tracking.timestamps)); % 0 or 1
patchTrialMask = zeros(size(tracking.timestamps)); %patch number within trial
patchTypeMask = zeros(size(tracking.timestamps)); %High or low prob patch
staySwitchMask = zeros(size(tracking.timestamps));
trajectoryMask = zeros(size(tracking.timestamps));
trajectoryDirMask = zeros(size(tracking.timestamps));
licksMask = zeros(size(tracking.timestamps));
firstLickMask = zeros(size(tracking.timestamps));



trial_begin = 0;
for ii = 1:(length(behavTrials.timestamps))
    posTrials = tracking.timestamps >= trial_begin & ...
                tracking.timestamps < behavTrials.timestamps(ii); % maybe i should keep all timestamps in behavTrials
    posLicks = tracking.timestamps >= trial_begin & ...
                tracking.timestamps < behavTrials.timestamps(ii); % irrelevant 
    

    trialnumberMask(posTrials) = ii;
    lickedPortMask(posTrials) = behavTrials.port(ii);
    outcomeMask(posTrials) = behavTrials.reward_outcome(ii)+1; % not rewarded is 1, rewarded is 2
    chosenPortProbMask(posTrials) = behavTrials.ports_probability(ii, behavTrials.port(ii));
    patchNumMask(posTrials) = behavTrials.patch_number(ii);
    patchTrialMask(posTrials) = behavTrials.patch_trials(ii);
    %patchTypeMask(posTrials) = behavTrials.patch_type(ii)+1;
    staySwitchMask(posTrials) = behavTrials.stay_switch(ii)+1; 
    trajectoryDirMask(posTrials) = trajectory_numbers_direction(ii);
    trajectoryMask(posTrials) = trajectory_numbers_no_direction(ii);
    licksMask(posLicks) = behavTrials.port(ii); % irrelevant

    trial_begin = behavTrials.timestamps(ii);
end


% interpolate - this downsamples the data to match intan sampling
trial_num_ds = interp1(tracking.timestamps,trialnumberMask,timestamp,'nearest');
licked_port_ds = interp1(tracking.timestamps,lickedPortMask,timestamp,'nearest');
outcome_ds = interp1(tracking.timestamps,outcomeMask,timestamp,'nearest');
chosen_port_prob_ds = interp1(tracking.timestamps,chosenPortProbMask,timestamp,'nearest');
patch_num_ds = interp1(tracking.timestamps,patchNumMask,timestamp,'nearest');
patch_trial_num_ds = interp1(tracking.timestamps,patchTrialMask,timestamp,'nearest');
patch_type_ds = interp1(tracking.timestamps,patchTypeMask,timestamp,'nearest');
stay_switch_ds = interp1(tracking.timestamps,staySwitchMask,timestamp,'nearest');
trajectory_dir_ds =interp1(tracking.timestamps,trajectoryDirMask,timestamp,'nearest');
trajectory_ds = interp1(tracking.timestamps,trajectoryMask,timestamp,'nearest');
licks_ds = interp1(tracking.timestamps,licksMask,timestamp,'nearest');

photometry_ds.sampling_rate = photometry.sampling_rate;
photometry_ds.timestamps = interp1(photometry.timestamps, photometry.timestamps, timestamp, 'nearest');
photometry_ds.grabDA_z = interp1(photometry.timestamps, photometry.grabDA_z, timestamp, 'nearest');
photometry_ds.grabDA_df = interp1(photometry.timestamps, photometry.grabDA_df, timestamp, 'nearest');
photometry_ds.grabDA_raw = interp1(photometry.timestamps, photometry.grabDA_raw, timestamp, 'nearest');


% maybe i should put the tracking back in a structure
position_x_all = interp1(tracking.timestamps,pos_x,timestamp,'linear'); 
position_y_all = interp1(tracking.timestamps,pos_y,timestamp,'linear'); 
speed_ds = interp1(tracking.timestamps,pos_v,timestamp,'linear'); 
speed_all = speed_ds';
speed_dir_ds = interp1(tracking.timestamps,pos_vy,timestamp,'linear');
speed_dir = speed_dir_ds';
timestamp_beh = timestamp;


%% Speed limit - why do i not want the data points when the mouse is stopped/moving slowly?
trial_num_ds = trial_num_ds(speed_ds >= speed_lim);
licked_port_ds = licked_port_ds(speed_ds >= speed_lim);
outcome_ds = outcome_ds(speed_ds >= speed_lim);
chosen_port_prob_ds = chosen_port_prob_ds(speed_ds >= speed_lim);
patch_num_ds = patch_num_ds(speed_ds >= speed_lim);
patch_trial_num_ds = patch_trial_num_ds(speed_ds >= speed_lim);
patch_type_ds = patch_type_ds(speed_ds >= speed_lim);
stay_switch_ds = stay_switch_ds(speed_ds >= speed_lim);
trajectory_dir_ds = trajectory_dir_ds(speed_ds >= speed_lim);
trajectory_ds = trajectory_ds(speed_ds >= speed_lim);
licks_ds = licks_ds(speed_ds >= speed_lim);

position_x_all = position_x_all(speed_ds >= speed_lim); 
position_y_all = position_y_all(speed_ds >= speed_lim); 
speed_dsa = speed_ds(speed_ds>= speed_lim);
speed_dir_lim = speed_dir_ds(speed_ds>= speed_lim);
speed_dir = speed_dir_lim';
direction = speed_dir > 0;
speed_all = speed_dsa';
timestamp_beh = timestamp_beh(speed_ds >= speed_lim);

photometry_ds.timestamps = photometry_ds.timestamps(speed_ds >= speed_lim);
photometry_ds.grabDA_z = photometry_ds.grabDA_z(speed_ds >= speed_lim);
photometry_ds.grabDA_df = photometry_ds.grabDA_df(speed_ds >= speed_lim);
photometry_ds.grabDA_raw = photometry_ds.grabDA_raw(speed_ds >= speed_lim);


spike_counts = [];
spike_counts = [spike_counts,SPIKEMAT.data(:,:)'];

data = normalize(spike_counts,1,'zscore');
data(isnan(data))=0;
data = data';

data = data(speed_ds >= speed_lim, :);
timestamp = timestamp(speed_ds >= speed_lim);

tracking_ds.x = position_x_all; 
tracking_ds.y = position_y_all; 
tracking_ds.speed = speed_all;
tracking_ds.speed_dir = speed_dir;
tracking_ds.direction = direction;
tracking_ds.timestamp = timestamp_beh;

%x = 1;

%% Save
save([currentFolderName, '.data.mat'], 'data', 'timestamp')
save([currentFolderName, '.PhotometryBehavDS.mat'], 'photometry_ds')
save([currentFolderName, 'TrackingDS'], 'tracking_ds')

save([currentFolderName,'.position_behavior_speed.mat'],'timestamp_beh','trial_num_ds','licked_port_ds','outcome_ds',...
    'chosen_port_prob_ds','position_x_all','position_y_all','licks_ds','speed_all','speed_dir','direction', 'trajectory_ds', 'trajectory_dir_ds','patch_num_ds','patch_trial_num_ds', 'patch_type_ds', 'stay_switch_ds');

end
