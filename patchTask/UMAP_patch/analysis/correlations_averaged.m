% correlations_2 but averaged across sessions


direc = '\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\patchTask';

sessions = {'N17\N17_250519_sess21', ...
    'N17\N17_250520_sess22'};

% Pre-define states and color palette
state_names = {'WAKEstate', 'NREMstate', 'REMstate'};
%colors = {[0.7412 0.8784 0.9176], [0.2431 0.7451 0.8706], [0.1922 0.4000 0.5804]};
colors = {[0.5451    0.9412    0.9686], [0.2706    0.1020    0.9412], [0.2431    0.7020    0.9294]}; % colors for plotting
sampling_rate = 130;
max_lag_sec = 2;
max_lag_samp = max_lag_sec * sampling_rate;

% Initialize xcorr storage per state
xcorr_all_sessions = struct();
for s = 1:length(state_names)
    xcorr_all_sessions.(state_names{s}) = {};  % each cell: one session’s xcorrs for this state
end

for sess = 1:length(sessions)
    session_path = fullfile(direc, sessions{sess});
    cd(session_path);
    basepath = pwd;
    disp(['Processing session: ', session_path]);

    % Load necessary files
    merge_file = dir(fullfile(basepath,'*.MergePoints.events.mat'));
    load(merge_file.name);

    sleep_state_file = dir(fullfile(basepath,'*.SleepState.states.mat'));
    load(sleep_state_file.name);

    % Determine which sleep block this is (e.g. pre or post)
    for ii = 1:size(MergePoints.foldernames,2)
        if isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep 'top*']))
            cd([basepath filesep MergePoints.foldernames{ii}]);
            curr_folder = pwd;
            % set pre_post to be 1 for pre behavior sleep, 2 for post behavior sleep  
            if ii == 1
                pre_post = 1;
            else
                pre_post = 2;
            end
            
            if pre_post == 1 
                sleep_start = MergePoints.timestamps(1,1); % start of sleep session
                sleep_end = MergePoints.timestamps(1,2); % time stamp of end of sleep session 
            elseif pre_post == 2
                sleep_start = MergePoints.timestamps(3,1); 
                sleep_end = MergePoints.timestamps(3,2);  
            end
            
            if ~isempty (dir(fullfile(curr_folder, '*HPCSleepPhotomSynced.mat'))) 
                load(dir(fullfile(curr_folder,'*HPCSleepPhotomSynced.mat')).name);
            end
    
            if ~isempty (dir(fullfile(curr_folder, '*StriatumSleepPhotomSynced.mat')))
                load(dir(fullfile(curr_folder,'*StriatumSleepPhotomSynced.mat')).name);
            end
            curr_folder = pwd;
       

            adjusted_ts_hpc = hpc_sleep_sync.timestamps + sleep_start;

            for s = 1:length(state_names)
                curr_state_name = state_names{s};
                if ~isfield(SleepState.ints, curr_state_name)
                    continue;
                end
        
                curr_state_epochs = SleepState.ints.(curr_state_name);
                % Only include epochs within this MergePoints window
                valid_epochs = curr_state_epochs(:,1) > sleep_start & ...
                               curr_state_epochs(:,2) < sleep_end;
                curr_state_epochs = curr_state_epochs(valid_epochs, :);
        
                if isempty(curr_state_epochs)
                    continue;
                end
        
                session_xcorrs = [];
        
                for i = 1:size(curr_state_epochs,1)
                    start_time = curr_state_epochs(i,1);
                    end_time   = curr_state_epochs(i,2);
        
                    if end_time > adjusted_ts_hpc(end)
                        continue;
                    end
        
                    [~, start_idx] = min(abs(adjusted_ts_hpc - start_time));
                    [~, end_idx] = min(abs(adjusted_ts_hpc - end_time));
        
                    if end_idx - start_idx < 2 * max_lag_samp || end_idx > length(adjusted_ts_hpc)
                        continue;
                    end
        
                    hpc_win = hpc_sleep_sync.grabDA_z(start_idx:end_idx);
                    str_win = striatum_sleep_sync.grabDA_z(start_idx:end_idx);
        
                    hpc_win = (hpc_win - mean(hpc_win)) / std(hpc_win);
                    str_win = (str_win - mean(str_win)) / std(str_win);
        
                    [xc, ~] = xcorr(hpc_win, str_win, max_lag_samp, 'coeff');
                    session_xcorrs = [session_xcorrs; xc'];
                end
        
                if ~isempty(session_xcorrs)
                    % Store for later averaging
                    xcorr_all_sessions.(curr_state_name){end+1} = mean(session_xcorrs, 1);
                end
            end
        end
    end
end


lags = -max_lag_samp:max_lag_samp;
lags_sec = lags / sampling_rate;

figure('color', 'white');
hold on;
legend_entries = {};

for s = 1:length(state_names)
    curr_state = state_names{s};
    state_xcorrs = xcorr_all_sessions.(curr_state);

    if isempty(state_xcorrs)
        continue;
    end

    state_matrix = cat(1, state_xcorrs{:});  % sessions × timepoints
    mean_xc = mean(state_matrix, 1);
    sem_xc = std(state_matrix, 0, 1) / sqrt(size(state_matrix,1));

    fill([lags_sec, fliplr(lags_sec)], ...
         [mean_xc + sem_xc, fliplr(mean_xc - sem_xc)], ...
         colors{s}, 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'HandleVisibility', 'off');
    plot(lags_sec, mean_xc, 'Color', colors{s} * 0.6, 'LineWidth', 2);
    legend_entries{end+1} = curr_state(1:end-5);
end

xline(0, '--k');
xlabel('Lag (s)');
ylabel('Cross-correlation coefficient');
title('Avg HPC–STR DA cross-correlation across sessions');
legend(legend_entries);

