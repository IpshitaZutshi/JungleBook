 

function manifold_saveData_zscore_ipshita(varargin)
%%
% Jan 2023, Winnie Yang
% Feb 2023, adding acceleration
%%
p = inputParser;
addParameter(p,'basepath',pwd);
addParameter(p,'save_name','linear');
addParameter(p,'save_ripple',false);
addParameter(p,'save_rippleHSE',false);
addParameter(p,'save_HSE',false);
addParameter(p,'save_theta',false);
addParameter(p,'save_REM',false);

addParameter(p,'rippleID',[]);
addParameter(p,'thetaID',[]);

addParameter(p,'behavior_epoch_index',[]);
addParameter(p,'beh_interval',[]);
addParameter(p,'REM_interval',[]);

addParameter(p,'SPIKEbin_beh',0.1);
addParameter(p,'SPIKEbin_ripple',0.02);
addParameter(p,'SPIKEbin_theta',0.02);

addParameter(p,'maze_type','linear');
addParameter(p,'sub_cell','all');
addParameter(p,'smooth_win',5);
addParameter(p,'brain_region',[]);
addParameter(p,'save_spikemat',true);

addParameter(p,'shuffle_theta',false);
addParameter(p,'shuffle_REM',false); % 
addParameter(p,'shuffle_rippleHSE',false);
addParameter(p,'shuffle_HSE',false);
addParameter(p,'shuffle_type','ir'); % 

addParameter(p,'num_shuffle',1); % 

addParameter(p,'rippleID_shuffle',[]);
addParameter(p,'thetaID_shuffle',[]);


addParameter(p,'shuffle_beh',false);

addParameter(p,'multiply',[]);

addParameter(p,'downsample_cell',false);
addParameter(p,'down_cell_num',60);



addParameter(p,'speed_lim',0);
addParameter(p,'theta_source','HSE'); % 

addParameter(p,'example_events',false); % 

addParameter(p,'leastNumTimeBin',3); % 
addParameter(p,'leastNumCell',3); % 
addParameter(p,'exclude_tuned_cell',0); % 
addParameter(p,'tuned_cell',0); % 
addParameter(p,'seed',0); % 
addParameter(p,'correct_only',0); % 
addParameter(p,'error_only',0); % 




% example usage:
% manifold_saveData_zscore('ripple',false,'theta',true, 'behavior_epoch_index',[2,4],'save_name', 'Tmaze1Linear2_zscore_goodTheta', 'maze_type', {'Tmaze','Linear'}, 'thetaID',[])
% Winnie Yang, JUNE 2022
%%
%ripple_idx = [1803, 1997,1794,1844,2021];
%beh_interval = [behavior.timestamps(1) behavior.timestamps(end)];
%SPIKEbin = 0.1;

parse(p,varargin{:})
error_only = p.Results.error_only;
correct_only = p.Results.correct_only;
basepath = p.Results.basepath;
save_name = p.Results.save_name;
save_ripple = p.Results.save_ripple;
save_HSE = p.Results.save_HSE;
save_rippleHSE = p.Results.save_rippleHSE;

save_theta = p.Results.save_theta;
save_REM = p.Results.save_REM;
rippleID = p.Results.rippleID;
thetaID = p.Results.thetaID;

behavior_epoch_index = p.Results.behavior_epoch_index;

beh_interval = p.Results.beh_interval;
REM_interval = p.Results.REM_interval;

SPIKEbin_beh = p.Results.SPIKEbin_beh;
SPIKEbin_ripple = p.Results.SPIKEbin_ripple;
SPIKEbin_theta = p.Results.SPIKEbin_theta;

maze_type = p.Results.maze_type;
sub_cell = p.Results.sub_cell;
smooth_win = p.Results.smooth_win;
brain_region = p.Results.brain_region;
save_spikemat= p.Results.save_spikemat;
shuffle_theta = p.Results.shuffle_theta;
shuffle_REM = p.Results.shuffle_REM;
shuffle_rippleHSE = p.Results.shuffle_rippleHSE;
shuffle_HSE = p.Results.shuffle_HSE;
rippleID_shuffle = p.Results.rippleID_shuffle;
thetaID_shuffle = p.Results.thetaID_shuffle;
shuffle_beh = p.Results.shuffle_beh;
downsample_cell =  p.Results.downsample_cell;
down_cell_num = p.Results.down_cell_num;
speed_lim = p.Results.speed_lim;
theta_source = p.Results.theta_source;
%ripple_source =p.Results.ripple_source;
example_events = p.Results.example_events;
shuffle_type = p.Results.shuffle_type;
num_shuffle = p.Results.num_shuffle;
leastNumTimeBin = p.Results.leastNumTimeBin;
leastNumCell = p.Results.leastNumCell;
multiply = p.Results.multiply;
exclude_tuned_cell = p.Results.exclude_tuned_cell;
tuned_cell = p.Results.tuned_cell;
seed = p.Results.seed;

%% load
basename = basenameFromBasepath(basepath);
save_path= [basepath, '\manifold'];
try 
    load([basepath, filesep, basename,'.Behavior.mat']);
catch
    load([basepath, filesep,basename,'.Tracking.Behavior.mat']);
    load([basepath, filesep,basename,'.TrialBehavior.Behavior.mat'])

end
load([basepath, filesep, basename,'.spikes.cellinfo.mat']);
try
    load([basepath, filesep,basename,'.session.mat']);
end

%% save path
if ~exist([basepath, '\manifold'])
    mkdir([basepath, '\manifold']);
end
%% cell type and brain region
switch sub_cell

    case 'all'
        sub_cells = 1:length(spikes.UID);
    case 'interneuron'
        load([basename,'.cell_metrics.cellinfo.mat'])
        Intcells = find(strcmp(cell_metrics.putativeCellType,'Narrow Interneuron') | strcmp(cell_metrics.putativeCellType,'Wide Interneuron'));
        sub_cells = Intcells;
    case 'pyramidal'       
        load([basename,'.cell_metrics.cellinfo.mat'])
        Extcells = find(strcmp(cell_metrics.putativeCellType,'Pyramidal Cell'));
        sub_cells = Extcells;
    case 'rHPC'
        load([basename,'.cell_metrics.cellinfo.mat'])
        rightCells = cell_metrics.rHPC;
        sub_cells = rightCells;
    case 'lHPC'
        load([basename,'.cell_metrics.cellinfo.mat'])
        leftCells = cell_metrics.lHPC;
        sub_cells = leftCells;
end

if exclude_tuned_cell
    load([basepath,filesep,basename,'.resultsPGAM_UMAP.mat'])
    switch tuned_cell
        case 'distStop'
            sub_cells = find(resultsPGAM_UMAP.tuned_distStop==0);
        case 'tone'
            sub_cells = find(resultsPGAM_UMAP.tuned_tone==0);  
    end
end
%%

if ~isempty(brain_region)

    switch brain_region
        case 'CA1'
            %CA1 cells only
            region_ind = [];
            shankID = session.brainRegions.CA1.electrodeGroups;
            for shank= 1:length(shankID)
                region_ind = [region_ind,find(spikes.shankID ==shankID(shank))];
            end
        case 'SUB'
            region_ind = [];
            shankID = session.brainRegions.SUB.electrodeGroups;
            for shank= 1:length(shankID)
                region_ind = [region_ind,find(spikes.shankID ==shankID(shank))];
            end
        case 'MEC'
            region_ind = [];
            region_ind = session.brainRegions.ENTm.channels;

    end
    sub_cells = intersect(sub_cells,region_ind);
end

%%
processInfo = {};
if exist('session','var')
    processInfo.epoch_index = 1:length(session.epochs);
else
    processInfo.epoch_index = 1;
end
processInfo.smooth_win = smooth_win;
%% (1) interval 
if exist('session','var')
    if isempty(beh_interval) 
        if ~isempty(behavior_epoch_index)
            beh_interval = zeros(length(behavior_epoch_index),2);
            for ep = 1:length(behavior_epoch_index)
                epoch = behavior_epoch_index(ep);
                if exist('behavior','var')
                    [mm1,I1]= min(abs(session.epochs{1,epoch}.startTime-behavior.timestamps));
                    [mm2,I2]= min(abs(session.epochs{1,epoch}.stopTime-behavior.timestamps));
                    beh_interval(ep,1) = behavior.timestamps(I1);
                    beh_interval(ep,2) = behavior.timestamps(I2);
                else
                    [mm1,I1]= min(abs(session.epochs{1,epoch}.startTime-tracking.timestamps));
                    [mm2,I2]= min(abs(session.epochs{1,epoch}.stopTime-tracking.timestamps));
                    beh_interval(ep,1) = tracking.timestamps(I1);
                    beh_interval(ep,2) = tracking.timestamps(I2);
                end
        
            end
        else
            beh_interval = [behavior.trials.interval(1,1), behavior.trials.interval(end,2)];
            behavior_epoch_index = unique(behavior.masks.EPOCH);
        end
            
    end
else
    beh_interval = [tracking.timestamps(1), tracking.timestamps(end)];
end


%%

spike_count_beh = [];
timestamp = [];
epoch_ind = [];
for int = 1:size(beh_interval,1)
    SPIKEMAT = bz_SpktToSpkmat(spikes, 'dt',SPIKEbin_beh,'win',beh_interval(int,:),'units','counts', 'smooth_win',smooth_win);
    
    if downsample_cell
%         SPIKEMAT = downsample_matrix(SPIKEMAT,'down_cell_num',down_cell_num,'down_cell',true,'sub_cell',sub_cells);
        rng(seed);
        downsample_ID = randsample(sub_cells,down_cell_num);
        spike_count_beh = [spike_count_beh,SPIKEMAT.data(:,downsample_ID)'];
        sub_cells =  1:size(spike_count_beh,1);
    else
        spike_count_beh = [spike_count_beh,SPIKEMAT.data(:,:)'];

    end
    timestamp = [timestamp, SPIKEMAT.timestamps'];
    epoch_ind = [epoch_ind,behavior_epoch_index(int)*ones(1,length(SPIKEMAT.timestamps))];
end


%% direction/ arm
direction = [];
for int = 1:length(maze_type)
    
    switch maze_type{int}
        case 'Tmaze'
                DIRECTION = behavior.masks.ARM( InIntervals(behavior.timestamps,beh_interval(int,:)));
                direction = [direction;DIRECTION];
        case 'Linear'
                DIRECTION = behavior.position.direction(InIntervals(behavior.timestamps,beh_interval(int,:)));
                direction = [direction;DIRECTION];
        case 'Radial'
                DIRECTION = behavior.masks.ARM(InIntervals(behavior.timestamps,beh_interval(int,:)));
                direction = [direction;DIRECTION];

    end
end
%direction = smooth(direction, smooth_bin);
%% behavior pos 
% smooth
if exist('behavior','var')
    trial_num = behavior.masks.TRIALS(InIntervals(behavior.timestamps,beh_interval));

    
    pos = behavior.position.lin(InIntervals(behavior.timestamps,beh_interval));
    %lin_norm = behavior.position.lin_norm(InIntervals(behavior.timestamps,beh_interval));
    pos_lin_norm =  behavior.position.lin_norm(InIntervals(behavior.timestamps,beh_interval));
    pos_x = smoothdata(behavior.position.x(InIntervals(behavior.timestamps,beh_interval)),'movmean',smooth_win);
    pos_y = smoothdata(behavior.position.y(InIntervals(behavior.timestamps,beh_interval)),'movmean',smooth_win);
    pos_x_bin = behavior.position.bin_x(InIntervals(behavior.timestamps,beh_interval));
    pos_y_bin = behavior.position.bin_y(InIntervals(behavior.timestamps,beh_interval));
    pos_rad = behavior.position.pos_rad(InIntervals(behavior.timestamps,beh_interval));
    pos_deg = behavior.position.pos_deg(InIntervals(behavior.timestamps,beh_interval));
    pos_rad_bin = behavior.position.pos_rad_bin(InIntervals(behavior.timestamps,beh_interval));
    pos_deg_bin = behavior.position.pos_deg_bin(InIntervals(behavior.timestamps,beh_interval));
    bin_size = behavior.position.bin_size;
    bin_size_rad = behavior.position.bin_size_rad;
else
    % for Ipshita's data
    trial_num = tracking.masks.TRIALS(InIntervals(tracking.timestamps,beh_interval));    
    trial_type = tracking.masks.TRIAL_TYPE(InIntervals(tracking.timestamps,beh_interval));  
    lick_loc = tracking.masks.LICK_LOC(InIntervals(tracking.timestamps,beh_interval));  
    correct = tracking.masks.CORRECT(InIntervals(tracking.timestamps,beh_interval));
    probe = tracking.masks.PROBE(InIntervals(tracking.timestamps,beh_interval));
    % Add stim as a variable
    

    
    pos_x = smoothdata(tracking.position.x(InIntervals(tracking.timestamps,beh_interval)),'movmean',smooth_win);
    pos_y = smoothdata(tracking.position.y(InIntervals(tracking.timestamps,beh_interval)),'movmean',smooth_win);
    speed = tracking.position.v;
    %acl = abs(diff(speed)./diff(tracking.timestamps));
    acl = diff(speed)./diff(tracking.timestamps);

    acceleration = [];
    acceleration(1) = acl(1);
    acceleration(2:length(speed)) = acl;
    bin_size = [];
    bin_size_rad =[];
end

%% interpolate
if exist('behavior','var')

    trial_num_ds = interp1(behavior.timestamps(InIntervals(behavior.timestamps,beh_interval)),trial_num,timestamp,'nearest');

    pos_ds = interp1(behavior.timestamps(InIntervals(behavior.timestamps,beh_interval)),pos,timestamp,'linear'); 
    pos_lin_norm_ds = interp1(behavior.timestamps(InIntervals(behavior.timestamps,beh_interval)),pos_lin_norm,timestamp,'linear'); 
    pos_x_ds = interp1(behavior.timestamps(InIntervals(behavior.timestamps,beh_interval)),pos_x,timestamp,'linear'); 
    pos_y_ds = interp1(behavior.timestamps(InIntervals(behavior.timestamps,beh_interval)),pos_y,timestamp,'linear'); 
    %pos_arm_ds = interp1(behavior.timestamps(InIntervals(behavior.timestamps,beh_interval)),pos_arm,timestamp,'linear'); 
    direction = double(direction);
    direction = direction(1:length(behavior.timestamps(InIntervals(behavior.timestamps,beh_interval))));
    direction_ds =  interp1(behavior.timestamps(InIntervals(behavior.timestamps,beh_interval)),double(direction),timestamp,'nearest'); 
    
    pos_x_bin_ds = interp1(behavior.timestamps(InIntervals(behavior.timestamps,beh_interval)),pos_x_bin,timestamp,'linear'); 
    pos_y_bin_ds = interp1(behavior.timestamps(InIntervals(behavior.timestamps,beh_interval)),pos_y_bin,timestamp,'linear');
    pos_rad_ds = interp1(behavior.timestamps(InIntervals(behavior.timestamps,beh_interval)),pos_rad,timestamp,'linear');
    pos_deg_ds = interp1(behavior.timestamps(InIntervals(behavior.timestamps,beh_interval)),pos_deg,timestamp,'linear');
    pos_rad_bin_ds = interp1(behavior.timestamps(InIntervals(behavior.timestamps,beh_interval)),pos_rad_bin,timestamp,'linear');
    pos_deg_bin_ds = interp1(behavior.timestamps(InIntervals(behavior.timestamps,beh_interval)),pos_deg_bin,timestamp,'linear');
    
    %speed
    t = 1:length(pos_ds) ;
    speed = zeros(length(t)-1,1) ;
    for i = 1:length(t)-1
        speed(i) = abs((pos_ds(i+1)-pos_ds(i))/(t(i+1)-t(i))) ;
    end
    speed = speed* mean(diff(timestamp));
    lick_loc_ds = [];
else
    trial_num_ds = interp1(tracking.timestamps(InIntervals(tracking.timestamps,beh_interval)),trial_num,timestamp,'nearest'); 
    correct_ds = interp1(tracking.timestamps(InIntervals(tracking.timestamps,beh_interval)),correct,timestamp,'nearest');
    probe_ds = interp1(tracking.timestamps(InIntervals(tracking.timestamps,beh_interval)),probe,timestamp,'nearest');

    trial_type_ds = interp1(tracking.timestamps(InIntervals(tracking.timestamps,beh_interval)),trial_type,timestamp,'nearest');
    lick_loc_ds = interp1(tracking.timestamps(InIntervals(tracking.timestamps,beh_interval)),lick_loc,timestamp,'nearest');
    pos_x_ds = interp1(tracking.timestamps(InIntervals(tracking.timestamps,beh_interval)),pos_x,timestamp,'linear'); 
    pos_y_ds = interp1(tracking.timestamps(InIntervals(tracking.timestamps,beh_interval)),pos_y,timestamp,'linear'); 
    speed_ds = interp1(tracking.timestamps(InIntervals(tracking.timestamps,beh_interval)),speed,timestamp,'linear'); 
    acceleration_ds =  interp1(tracking.timestamps(InIntervals(tracking.timestamps,beh_interval)),acceleration,timestamp,'linear'); 
end
%% 
figure()
plot(acceleration)
%% add speed limit

figure; plot(speed_ds);ylim([0,10])
if ~isempty(speed_lim) && (~correct_only && ~error_only)
    high_speed_ind = find(speed_ds>speed_lim);
    speed_ds = speed_ds(high_speed_ind);
    acceleration_ds = acceleration_ds(high_speed_ind);
    if exist('behavior','var')

        trial_num_ds = trial_num_ds(high_speed_ind);
        direction_ds =direction_ds(high_speed_ind);
        trial_type_ds = direction_ds;
        pos_ds = pos_ds(high_speed_ind) ;
        pos_lin_norm_ds = pos_lin_norm_ds(high_speed_ind);
        %lin_norm_ds = lin_norm_ds(high_speed_ind) ;
        pos_x_ds = pos_x_ds(high_speed_ind);
        pos_y_ds = pos_y_ds(high_speed_ind);
    
        pos_x_bin_ds = pos_x_bin_ds(high_speed_ind);
        pos_y_bin_ds = pos_y_bin_ds(high_speed_ind);
    
        pos_rad_ds = pos_rad_ds(high_speed_ind);
        pos_deg_ds = pos_deg_ds(high_speed_ind);
        pos_rad_bin_ds = pos_rad_bin_ds(high_speed_ind);
        pos_deg_bin_ds = pos_deg_bin_ds(high_speed_ind);
        %pos_arm_ds = pos_arm_ds(high_speed_ind);
        


        figure;
        plot(pos_ds);
        title('after speed limit');

    elseif exist('behavTrials','var') % for Ipshita's data

        pos_x_ds = pos_x_ds(high_speed_ind);
        pos_y_ds = pos_y_ds(high_speed_ind);
        trial_num_ds = trial_num_ds(high_speed_ind);
        correct_ds = correct_ds(high_speed_ind);
        probe_ds = probe_ds(high_speed_ind);
  
        trial_type_ds = trial_type_ds(high_speed_ind);
        lick_loc_ds = lick_loc_ds(high_speed_ind);


    end

    timestamp =timestamp(high_speed_ind);
    epoch_ind = epoch_ind(high_speed_ind);
    spike_count_beh = spike_count_beh(:,high_speed_ind);

end


epoch_ind_beh=epoch_ind;
timestamp_beh = timestamp;
%% filter correct or error trials
if correct_only
    high_speed_ind = find(correct_ds==1 & speed_ds>speed_lim);
    speed_ds = speed_ds(high_speed_ind);
    acceleration_ds = acceleration_ds(high_speed_ind);
    if exist('behavTrials','var') % for Ipshita's data

        pos_x_ds = pos_x_ds(high_speed_ind);
        pos_y_ds = pos_y_ds(high_speed_ind);
        trial_num_ds = trial_num_ds(high_speed_ind);
        correct_ds = correct_ds(high_speed_ind);
        probe_ds = probe_ds(high_speed_ind);
  
        trial_type_ds = trial_type_ds(high_speed_ind);
        lick_loc_ds = lick_loc_ds(high_speed_ind);

        timestamp =timestamp(high_speed_ind);
        epoch_ind = epoch_ind(high_speed_ind);
        spike_count_beh = spike_count_beh(:,high_speed_ind);

    end
elseif error_only
    high_speed_ind = find(correct_ds==0 & speed_ds>speed_lim);
    speed_ds = speed_ds(high_speed_ind);
    acceleration_ds = acceleration_ds(high_speed_ind);
    if exist('behavTrials','var') % for Ipshita's data

        pos_x_ds = pos_x_ds(high_speed_ind);
        pos_y_ds = pos_y_ds(high_speed_ind);
        trial_num_ds = trial_num_ds(high_speed_ind);
        correct_ds = correct_ds(high_speed_ind);
        probe_ds = probe_ds(high_speed_ind);
  
        trial_type_ds = trial_type_ds(high_speed_ind);
        lick_loc_ds = lick_loc_ds(high_speed_ind);

        timestamp =timestamp(high_speed_ind);
        epoch_ind = epoch_ind(high_speed_ind);
        spike_count_beh = spike_count_beh(:,high_speed_ind);

    end

end

epoch_ind_beh=epoch_ind;
timestamp_beh = timestamp;
%% shuffle behavior
if shuffle_beh
    [spike_count_shuffle,~] = shuffle_spikemat(spike_count_beh,'shuffle_type',shuffle_type,'num_shuffle',num_shuffle);
    spike_count_beh_shuffle = [spike_count_beh,spike_count_shuffle];
end
%% (2) ripple spiking data

%%
spike_eventIDs_ripple = [];
%spike_epochIDs_ripple = [];
spike_count_ripple=[];
spike_ts_ripple = [];
%ripple_idx = [];
if save_HSE
        disp('preprocessing HSE events...')
    
        load([basepath, '\',basename,'.HSE.events.mat']);
        evtIDs_all = 1:length(HSE.timestamps);

        if ~exist([basepath,'\',basename, '.spikemat_HSE.mat'])
            [SPIKEMAT,eventIDs,epochIDs] = getRipThetaSpikeMat(spikes,HSE, evtIDs_all, SPIKEbin_ripple,'event_type','rippleHSE','basepath',basepath,...
                'leastNumTimeBin',leastNumTimeBin,'leastNumCell',leastNumCell,'processInfo',processInfo);
            if save_spikemat 

                save([basepath, '\',basename,'.spikemat_HSE.mat'],'SPIKEMAT');
            end
        else
            load([basepath,'\',basename, '.spikemat_HSE.mat']);
        end
        spike_count_ripple = SPIKEMAT.spikemat_all_event;
        spike_eventIDs_ripple = SPIKEMAT.spikemat_all_eventID;
        spike_epochIDs_ripple = SPIKEMAT.spikemat_all_epochID;
        spike_ts_ripple = SPIKEMAT.spikemat_all_ts;
        epoch_ind = [epoch_ind,spike_epochIDs_ripple];
    
elseif save_rippleHSE
        disp('preprocessing rippleHSE events...')
    
        load([basepath, '\',basename,'.rippleHSE.events.mat']);
        evtIDs_all = 1:length(ripple_HSE.timestamps);

        if ~exist([basepath,'\',basename, '.spikemat_RSE.mat'])
            [SPIKEMAT,eventIDs,epochIDs] = getRipThetaSpikeMat(spikes,ripple_HSE, evtIDs_all, SPIKEbin_ripple,'event_type','rippleHSE','basepath',basepath,...
                'leastNumTimeBin',leastNumTimeBin,'leastNumCell',leastNumCell,'processInfo',processInfo);
            if save_spikemat 

                save([basepath, '\',basename,'.spikemat_RSE.mat'],'SPIKEMAT');
            end
        else
            load([basepath,'\',basename, '.spikemat_RSE.mat']);
        end
        spike_count_ripple = SPIKEMAT.spikemat_all_event;
        spike_eventIDs_ripple = SPIKEMAT.spikemat_all_eventID;
        spike_epochIDs_ripple = SPIKEMAT.spikemat_all_epochID;
        spike_ts_ripple = SPIKEMAT.spikemat_all_ts;
        epoch_ind = [epoch_ind,spike_epochIDs_ripple];
    


    

    
elseif save_ripple
    %if ~exist([basepath,'\',basename,'.spikemat_ripple_good_bad.mat']) || example_events
        disp('preprocessing ripple events...')
       
        load([basepath, '\',basename,'.ripples.events.mat']);
        evtIDs_all = 1:length(ripples.timestamps);

        [SPIKEMAT_ripple,eventIDs_ripple,epochIDs_ripple] = getRipThetaSpikeMat(spikes,ripples, evtIDs_all, SPIKEbin_ripple,'event_type','rippleHSE','basepath',basepath,...
            'leastNumTimeBin',leastNumTimeBin,'leastNumCell',leastNumCell);
        spike_count_ripple = SPIKEMAT_ripple.spikemat_all_event;
        spike_eventIDs_ripple = SPIKEMAT_ripple.spikemat_all_eventID;
        spike_epochIDs_ripple = SPIKEMAT_ripple.spikemat_all_epochID;
        spike_ts_ripple = SPIKEMAT_ripple.spikemat_all_ts;
        epoch_ind = [epoch_ind,spike_epochIDs_ripple];
    
        %% save
        if save_spikemat 
            SPIKEMAT.processInfo.SPIKEbin = SPIKEbin_ripple;
            SPIKEMAT.processInfo.epoch_index = 1:length(session.epochs);
            SPIKEMAT.processInfo.units = 'counts';
            SPIKEMAT.processInfo.smooth_win = smooth_win;
            SPIKEMAT.eventIDs =eventIDs_ripple;
            SPIKEMAT.epochIDs = epochIDs_ripple;
            save([basepath, '\',basename,'.spikemat_ripple.mat'],'SPIKEMAT');
        end

end

%% (3) theta spike data


spike_eventIDs_theta = [];
%spike_epochIDs_theta =[];
%theta_idx = [];
spike_count_theta = [];
spike_ts_theta = [];
if save_theta
        load([basepath, '\',basename,'.sleepState.states.mat']);

        %if ~exist([basepath,'\',basename,'.spikemat_theta.mat'],'file')|| example_events
        disp('preprocessing theta events...')
           
        switch theta_source
                case 'thetaEpochs'

                    load([basepath, '\',basename,'.thetaEpochs.states.mat']);
                    [SPIKEMAT,eventIDs,epochIDs] = getRipThetaSpikeMat(spikes,thetaEpochs, thetaID, SPIKEbin_theta,'event_type','theta','basepath',basepath,...
                        'leastNumTimeBin',leastNumTimeBin,'leastNumCell',leastNumCell,'theta_source',theta_source);
                 case 'REM_cycles'   
                    if ~exist([basepath,'\',basename, '.spikemat_REM_cycles.mat'])

                        load([basepath, '\',basename,'.REM_cycles.events.mat']);
                        [SPIKEMAT,eventIDs,epochIDs] = getRipThetaSpikeMat(spikes,REM_cycles, thetaID, SPIKEbin_theta,'event_type','theta','basepath',basepath,...
                            'leastNumTimeBin',leastNumTimeBin,'leastNumCell',leastNumCell,'theta_source',theta_source);
                    else
                        load([basepath,'\',basename, '.spikemat_REM_cycles.mat'])
                    end
                case 'RUN_cycles'
                    if ~exist([basepath,'\',basename, '.spikemat_RUN_cycles.mat'])
                        load([basepath, '\',basename,'.RUN_cycles.events.mat']);
                        [SPIKEMAT,eventIDs,epochIDs] = getRipThetaSpikeMat(spikes,RUN_cycles, thetaID, SPIKEbin_theta,'event_type','theta','basepath',basepath,...
                            'leastNumTimeBin',leastNumTimeBin,'leastNumCell',leastNumCell,'theta_source',theta_source);
                    else
                        load([basepath,'\',basename, '.spikemat_RUN_cycles.mat'])
                    end
                case 'HSE'
                    load([basepath, '\',basename,'.HSE.events.mat']);
                    evtIDs_all = 1:length(HSE.timestamps);
                    theta_idx = zeros(length(thetaID),1);
                    for evt = 1:length(thetaID)
                        theta_idx(evt)=find(evtIDs_all==thetaID(evt));
                    end
                    [SPIKEMAT,eventIDs,epochIDs] = getRipThetaSpikeMat(spikes,HSE, theta_idx, SPIKEbin_ripple,'event_type','ripple','basepath',basepath,...
                        'leastNumTimeBin',leastNumTimeBin,'leastNumCell',leastNumCell);
              
        end           
            spike_count_theta = SPIKEMAT.spikemat_all_event;
            spike_eventIDs_theta = SPIKEMAT.spikemat_all_eventID;
            spike_epochIDs_theta = SPIKEMAT.spikemat_all_epochID;
            spike_ts_theta = SPIKEMAT.spikemat_all_ts;
            epoch_ind = [epoch_ind,spike_epochIDs_theta];
            %% save
            if save_spikemat 
                SPIKEMAT.processInfo.SPIKEbin = SPIKEbin_theta;
                SPIKEMAT.processInfo.epoch_index = 1:length(session.epochs);
                SPIKEMAT.processsInfo.units = 'counts';
                SPIKEMAT.processInfo.smooth_win = smooth_win;
                SPIKEMAT.eventIDs =SPIKEMAT.spikemat_all_eventID;
                SPIKEMAT.epochIDs = SPIKEMAT.spikemat_all_epochID;
           
                save([basepath, '\',basename,'.spikemat_',theta_source,'.mat'],'SPIKEMAT');
            end

end

if save_REM
        disp('preprocessing REM events...')
            load([basepath,'\',basename,'.sleepStateEpisodes.states.mat'])
            REM_ints = SleepStateEpisodes.ints.REMepisode;
     
            [SPIKEMAT,eventIDs,epochIDs] = getRipThetaSpikeMat(spikes,SleepState, 1:length(REM_ints), SPIKEbin_theta,'event_type','REM','basepath',basepath,...
                'leastNumTimeBin',leastNumTimeBin,'leastNumCell',leastNumCell,'theta_source',theta_source);
    
            spike_count_theta = SPIKEMAT.spikemat_all_event;
            spike_eventIDs_theta = SPIKEMAT.spikemat_all_eventID;
            spike_epochIDs_theta = SPIKEMAT.spikemat_all_epochID;
            spike_ts_theta = SPIKEMAT.spikemat_all_ts;
            epoch_ind = [epoch_ind,spike_epochIDs_theta];
            %% save
            if save_spikemat && ~example_events
                SPIKEMAT.processInfo.SPIKEbin = SPIKEbin_theta;
                SPIKEMAT.processInfo.epoch_index = 1:length(session.epochs);
                SPIKEMAT.processsInfo.units = 'counts';
                SPIKEMAT.processInfo.smooth_win = smooth_win;
                SPIKEMAT.eventIDs =eventIDs;
                SPIKEMAT.epochIDs = epochIDs;
           
                save([basepath, '\',basename,'.spikemat_REM.mat'],'SPIKEMAT');
            end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (5) shuffle ripple events
if shuffle_rippleHSE

%     if ~exist([basepath,'\',basename,'.spikemat_rippleHSE_shuffle.mat'],'file') || example_events
        disp('preprocessing rippleHSE shuffle events...')
       
        load([basepath, '\',basename,'.rippleHSE.events.mat']);
        disp('loading spikemat...')

        load([basepath,'\',basename,'.spikemat_RSE.mat']);
        SPIKEMAT_ripple = SPIKEMAT;

        [spike_count_shuffle_ripple,shuffle_num_ripple] = shuffle_spikemat( SPIKEMAT_ripple.spikemat_all_event,'shuffle_type',shuffle_type,'num_shuffle',num_shuffle);
        SPIKEMAT.spikemat_all_event = spike_count_shuffle_ripple;
        spike_eventIDs_shuffle_ripple = repmat(SPIKEMAT.spikemat_all_eventID,1,num_shuffle);
        spike_epochIDs_shuffle_ripple = repmat(SPIKEMAT.spikemat_all_epochID,1,num_shuffle);
        spike_ts_shuffle_ripple = repmat(SPIKEMAT.spikemat_all_ts,1,num_shuffle);
        epoch_ind = [epoch_ind,spike_epochIDs_shuffle_ripple];
        eventIDs = repmat(SPIKEMAT.eventIDs,1,num_shuffle);
        epochIDs = repmat(SPIKEMAT.epochIDs,1,num_shuffle);
        %% save
        if save_spikemat && ~example_events
            SPIKEMAT.processInfo.SPIKEbin = SPIKEbin_ripple;
            %SPIKEMAT.processInfo.epoch_index = 1:length(session.epochs);
            SPIKEMAT.processsInfo.units = 'counts';
            SPIKEMAT.processInfo.smooth_win = smooth_win;
            SPIKEMAT.eventIDs =eventIDs;
            SPIKEMAT.epochIDs = epochIDs;
            save([basepath, '\',basename,'.spikemat_',save_name,'.mat'],'SPIKEMAT');
        end

else
    spike_count_shuffle_ripple = [];
    spike_ts_shuffle_ripple = [];
    spike_eventIDs_shuffle_ripple = [];
    shuffle_num_ripple = [];
end
%%
if shuffle_HSE

    %if ~exist([basepath,'\',basename,'.spikemat_HSE_shuffle.mat'],'file') || example_events
        disp('preprocessing HSE shuffle events...')
        
        load([basepath, '\',basename,'.HSE.events.mat']);


        load([basepath,'\',basename,'.spikemat_HSE.mat']);

        [spike_count_shuffle_ripple,shuffle_num_ripple] = shuffle_spikemat( SPIKEMAT.spikemat_all_event,'shuffle_type',shuffle_type,'num_shuffle',num_shuffle);
        SPIKEMAT.spikemat_all_event = spike_count_shuffle_ripple;
        spike_eventIDs_shuffle_ripple = repmat(SPIKEMAT.spikemat_all_eventID,1,num_shuffle);
        spike_epochIDs_shuffle_ripple = repmat(SPIKEMAT.spikemat_all_epochID,1,num_shuffle);
        spike_ts_shuffle_ripple = repmat(SPIKEMAT.spikemat_all_ts,1,num_shuffle);
        epoch_ind = [epoch_ind,spike_epochIDs_shuffle_ripple];
        eventIDs = repmat(SPIKEMAT.eventIDs,1,num_shuffle);
        epochIDs = repmat(SPIKEMAT.epochIDs,1,num_shuffle);
        %% save
        if save_spikemat && ~example_events
            SPIKEMAT.spikemat_all_event = spike_count_shuffle_ripple;
            SPIKEMAT.processInfo.SPIKEbin = SPIKEbin_ripple;
            SPIKEMAT.processInfo.epoch_index = 1:length(session.epochs);
            SPIKEMAT.processsInfo.units = 'counts';
            SPIKEMAT.processInfo.smooth_win = smooth_win;
            SPIKEMAT.eventIDs =eventIDs;
            SPIKEMAT.epochIDs = epochIDs;
            save([basepath, '\',basename,'.',save_name,'.mat'],'SPIKEMAT');
        end

end
%% (6) shuffle theta events
if shuffle_theta

        disp('preprocessing theta shuffle events...')

        load([basepath,'\',basename,'.spikemat_theta.mat']);
        [spike_count_shuffle_theta,shuffle_num_theta] = shuffle_spikemat( SPIKEMAT.spikemat_all_event,'shuffle_type',shuffle_type,'num_shuffle',num_shuffle);
        SPIKEMAT.spikemat_all_event = spike_count_shuffle_theta;
        spike_eventIDs_shuffle_theta = repmat(SPIKEMAT.spikemat_all_eventID,1,num_shuffle);
        spike_epochIDs_shuffle_theta = repmat(SPIKEMAT.spikemat_all_epochID,1,num_shuffle);
        spike_ts_shuffle_theta = repmat(SPIKEMAT.spikemat_all_ts,1,num_shuffle);
        epoch_ind = [epoch_ind,spike_epochIDs_shuffle_theta];
        eventIDs = repmat(eventIDs,1,num_shuffle);
        epochIDs = repmat(epochIDs,1,num_shuffle);
        %% save
        if save_spikemat && ~example_events
            SPIKEMAT.processInfo.SPIKEbin = SPIKEbin_theta;
            SPIKEMAT.processInfo.epoch_index = 1:length(session.epochs);
            SPIKEMAT.processsInfo.units = 'counts';
            SPIKEMAT.processInfo.smooth_win = smooth_win;
            SPIKEMAT.eventIDs =eventIDs;
            SPIKEMAT.epochIDs = epochIDs;
            save([basepath, '\',basename,'.spikemat_theta_shuffle.mat'],'SPIKEMAT');
        end

else
    spike_count_shuffle_theta = [];
    spike_ts_shuffle_theta = [];
    spike_eventIDs_shuffle_theta = [];
    shuffle_num_theta = [];
end

if shuffle_REM

        disp('preprocessing REM shuffle events...')

        load([basepath,'\',basename,'.spikemat_REM.mat']);
        [spike_count_shuffle_theta,shuffle_num_theta] = shuffle_spikemat( SPIKEMAT.spikemat_all_event,'shuffle_type',shuffle_type,'num_shuffle',num_shuffle);
        SPIKEMAT.spikemat_all_event = spike_count_shuffle_theta;
        spike_eventIDs_shuffle_theta = repmat(SPIKEMAT.spikemat_all_eventID,1,num_shuffle);
        spike_epochIDs_shuffle_theta = repmat(SPIKEMAT.spikemat_all_epochID,1,num_shuffle);
        spike_ts_shuffle_theta = repmat(SPIKEMAT.spikemat_all_ts,1,num_shuffle);
        epoch_ind = [epoch_ind,spike_epochIDs_shuffle_theta];
        eventIDs = repmat(eventIDs,1,num_shuffle);
        epochIDs = repmat(epochIDs,1,num_shuffle);
        %% save
        if save_spikemat && ~example_events
            SPIKEMAT.processInfo.SPIKEbin = SPIKEbin_theta;
            SPIKEMAT.processInfo.epoch_index = 1:length(session.epochs);
            SPIKEMAT.processsInfo.units = 'counts';
            SPIKEMAT.processInfo.smooth_win = smooth_win;
            SPIKEMAT.eventIDs =eventIDs;
            SPIKEMAT.epochIDs = epochIDs;
            save([basepath, '\',basename,save_name,'.mat'],'SPIKEMAT');
        end

else
    spike_count_shuffle_theta = [];
    spike_ts_shuffle_theta = [];
    spike_eventIDs_shuffle_theta = [];
    shuffle_num_theta = [];
end
%% (6) concatenate all


if shuffle_beh
    spike_counts = [spike_count_beh_shuffle,spike_count_ripple,spike_count_theta,spike_count_shuffle_ripple,spike_count_shuffle_theta];
else
    spike_counts = [spike_count_beh,spike_count_ripple,spike_count_theta,spike_count_shuffle_ripple,spike_count_shuffle_theta];

end
timestamps = [timestamp,spike_ts_ripple,spike_ts_theta,spike_ts_shuffle_ripple,spike_ts_shuffle_theta];

%position_arm_ds = nan(1,length(spike_counts));
%position_arm_ds(1:length(spike_count_beh)) = pos_arm_ds;

beh_ind = zeros(1,length(spike_counts));
beh_ind(1:size(spike_count_beh,2)) = ones(1,size(spike_count_beh,2));

ripple_ind = zeros(1,length(spike_counts));
ripple_ind(size(spike_count_beh,2)+1:size(spike_count_beh,2)+size(spike_count_ripple,2)) = spike_eventIDs_ripple;


theta_ind = zeros(1,length(spike_counts));
theta_ind(size(spike_count_beh,2)+size(spike_count_ripple,2)+1: size(spike_count_beh,2)+size(spike_count_ripple,2)+size(spike_count_theta,2)) = spike_eventIDs_theta;



ripple_shuffle_ind = zeros(1,length(spike_counts));
ripple_shuffle_ind(size(spike_count_beh,2)+size(spike_count_ripple,2)+size(spike_count_theta,2)+1: size(spike_count_beh,2)+size(spike_count_ripple,2)+size(spike_count_theta,2)+size(spike_count_shuffle_ripple,2)) = spike_eventIDs_shuffle_ripple;


theta_shuffle_ind = zeros(1,length(spike_counts));
theta_shuffle_ind(size(spike_count_beh,2)+size(spike_count_ripple,2)+size(spike_count_theta,2)+size(spike_count_shuffle_ripple,2)+1: size(spike_count_beh,2)+size(spike_count_ripple,2)+size(spike_count_theta,2)+size(spike_count_shuffle_ripple,2)+size(spike_count_shuffle_theta,2)) = spike_eventIDs_shuffle_theta;

shuffle_num_ripple_all = zeros(1,length(spike_counts));
shuffle_num_ripple_all(size(spike_count_beh,2)+size(spike_count_ripple,2)+size(spike_count_theta,2)+1: size(spike_count_beh,2)+size(spike_count_ripple,2)+size(spike_count_theta,2)+size(spike_count_shuffle_ripple,2)) = shuffle_num_ripple;

shuffle_num_theta_all = zeros(1,length(spike_counts));
shuffle_num_theta_all(size(spike_count_beh,2)+size(spike_count_ripple,2)+size(spike_count_theta,2)+size(spike_count_shuffle_ripple,2)+1: size(spike_count_beh,2)+size(spike_count_ripple,2)+size(spike_count_theta,2)+size(spike_count_shuffle_ripple,2)+size(spike_count_shuffle_theta,2)) = shuffle_num_theta;


%%
%data = softmax(spike_counts(sub_cells,:)');
data = normalize(spike_counts(sub_cells,:),1,'zscore');

data(isnan(data))=0;
data = data';


%% save spiking info

if downsample_cell
    save([save_path, '\', basename, '.data_',save_name,'_cell_',num2str(down_cell_num),'_seed_',num2str(seed),'.mat'],'data','timestamps','beh_ind','ripple_ind','theta_ind','ripple_shuffle_ind','theta_shuffle_ind','shuffle_num_theta_all','shuffle_num_ripple_all','epoch_ind','-v7.3'); 

else
    save([save_path, '\', basename, '.data_',save_name,'.mat'],'data','timestamps','beh_ind','ripple_ind','theta_ind','ripple_shuffle_ind','theta_shuffle_ind','shuffle_num_theta_all','shuffle_num_ripple_all','epoch_ind','-v7.3'); 
end
%% save behavior info
position_x_all = pos_x_ds;
position_y_all = pos_y_ds;
speed_all = speed_ds';
acceleration_all = acceleration_ds';
if exist('behavior','var')
    position_all = pos_ds;
    pos_lin_norm_all = pos_lin_norm_ds;
    
    direction_all = direction_ds;
    pos_x_bin = pos_x_bin_ds;
    pos_y_bin = pos_y_bin_ds;
    pos_rad = pos_rad_ds;
    pos_deg = pos_deg_ds;
    
    pos_rad_bin = pos_rad_bin_ds;
    pos_deg_bin = pos_deg_bin_ds;
else
    position_all = [];
    pos_lin_norm_all = [];
    
    direction_all = [];
    pos_x_bin = [];
    pos_y_bin = [];
    pos_rad = [];
    pos_deg = [];
    
    pos_rad_bin = [];
    pos_deg_bin = [];
end
save([save_path, '\', basename, '.position_',save_name,'.mat'],'timestamp_beh','epoch_ind_beh','position_all','pos_lin_norm_all','position_x_all','position_y_all',...
    'pos_x_bin','pos_y_bin','pos_rad','pos_deg','pos_rad_bin','pos_deg_bin',...
    'direction_all','speed_all','bin_size','bin_size_rad','trial_num_ds','trial_type_ds',...
    'lick_loc_ds','acceleration_all',...
    'correct_ds','probe_ds');

