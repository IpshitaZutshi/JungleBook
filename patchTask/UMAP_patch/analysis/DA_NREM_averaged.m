%% AVERAGED ACROSS SESSIONS DA throughout NREM
% NEED TO FIX SAVING SO IT HAS PROPER TITLE


%% Inputs
basepath = pwd;

channels = [0:41 43:47 49:55 60 62:69 73:127]; %N17
%lfp = bz_GetLFP(channels);
%lfp = bz_GetLFP(channels, 'restrict', [0 3600]);


%% Step 1: Extract NREM packets


save_path = '\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\patchTask\figures';
saveLoc = strcat(save_path,'\sleep_phot');

% hpc sessions
% sessions_N11 = {'N11\Final\N11_250314_sess18', ...
%     'N11\Final\N11_250318_sess19' ...
%     'N11\Final\N11_250321_sess22'};
% 
% % Set up session list
% sessions_N17 = {'N17\N17_250501_sess10', ...
%     'N17\N17_250509_sess15', ...
%     'N17\N17_250513_sess19', ...
%     'N17\N17_250519_sess21', ...
%     'N17\N17_250520_sess22'};


% striatum sessions
sessions_N11 = {'N11\Final\N11_250312_sess16', ...
    'N11\Final\N11_250313_sess17' ...
    'N11\Final\N11_250401_sess25'};

% Set up session list
sessions_N17 = {'N17\N17_250509_sess15', ...
    'N17\N17_250512_sess18', ...
    'N17\N17_250519_sess21', ...
    'N17\N17_250520_sess22'};


sessions = sessions_N17;  
%sessions = cat(2, sessions_N11, sessions_N17);

direc = '\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\patchTask';

% Settings
fs_da = 130;
nPoints = 100;

% Initialize storage
all_sessions_pre = {};  % interpolated DA traces for pre-sleep
all_sessions_post = {}; % interpolated DA traces for post-sleep

for s = 1:length(sessions)
    sessionpath = fullfile(direc, sessions{s});
    cd(sessionpath);
    disp(['Processing session: ' sessions{s}]);

    % Load session files
    load(dir(fullfile(sessionpath, '*.MergePoints.events.mat')).name);
    load(dir(fullfile(sessionpath, '*.SleepState.states.mat')).name);

    for ii = 1:size(MergePoints.foldernames, 2)
        sleep_dir = fullfile(sessionpath, MergePoints.foldernames{ii});

        if isempty(dir(fullfile(sleep_dir, 'top*')))
            cd(sleep_dir);

            % Determine pre or post sleep
            if ii == 1
                pre_post = 1; % pre sleep
            else
                pre_post = 2; % post sleep
            end

            % Load photometry data
            % if ~isempty(dir(fullfile(sleep_dir, '*HPCSleepPhotomSynced.mat')))
            %     load(dir(fullfile(sleep_dir, '*HPCSleepPhotomSynced.mat')).name);
            %     daSignal = hpc_sleep_sync.grabDA_z;
            %     timestamps_da = hpc_sleep_sync.timestamps;
            %     color_pre = [0.588235294117647   0.800000000000000   0.345098039215686]; % green
            %     color_post = [0.301960784313725   0.458823529411765   0.192156862745098]; % darker green
            if ~isempty(dir(fullfile(sleep_dir, '*StriatumSleepPhotomSynced.mat')))
                load(dir(fullfile(sleep_dir, '*StriatumSleepPhotomSynced.mat')).name);
                daSignal = striatum_sleep_sync.grabDA_z;
                timestamps_da = striatum_sleep_sync.timestamps;
                color_pre = [0.909803921568627   0.290196078431373   0.454901960784314]; % bright pink
                color_post = [0.647058823529412   0.196078431372549   0.313725490196078]; % darker pink
            % else
            %     warning('No photometry data found in %s', sleep_dir);
            %     continue;
            end

            % Get sleep start time for offset
            if pre_post == 1
                sleep_start = MergePoints.timestamps(1, 1);
            elseif pre_post == 2
                sleep_start = MergePoints.timestamps(3, 1);
            else
                continue;
            end
            timestamps_da = timestamps_da + sleep_start;

            % Get NREM intervals
            nrem_intervals = SleepState.ints.NREMstate;
            durations = nrem_intervals(:, 2) - nrem_intervals(:, 1);
            valid_idx = durations > 10;
            nrem_intervals = nrem_intervals(valid_idx, :);

            % Interpolate DA for each NREM
            all_da_interp = nan(length(nrem_intervals), nPoints);

            for i = 1:size(nrem_intervals, 1)
                t_start = nrem_intervals(i, 1);
                t_end = nrem_intervals(i, 2);

                idx = find(timestamps_da >= t_start & timestamps_da <= t_end);
                if length(idx) < 2
                    continue;
                end

                da_seg = daSignal(idx);
                t_seg = timestamps_da(idx);
                
                % Normalize time to 0-1 scale
                t_norm = (t_seg - t_start) / (t_end - t_start);
                
                % Normalize da_seg within this NREM interval (mean 0, std 1)
                da_seg = (da_seg - mean(da_seg)) / std(da_seg);
                
                % Interpolate to fixed number of points
                tq = linspace(0, 1, nPoints);
                da_interp = interp1(t_norm, da_seg, tq, 'linear');
                
                all_da_interp(i, :) = da_interp;

            end

            % Store per session
            if pre_post == 1
                all_sessions_pre{end+1} = all_da_interp;
            elseif pre_post == 2
                all_sessions_post{end+1} = all_da_interp;
            end
        end
    end
end

%% Grand averaging across all sessions
% Combine all rows from all sessions into a single matrix
concat_pre = vertcat(all_sessions_pre{:});
concat_post = vertcat(all_sessions_post{:});

% Mean and SEM
mean_pre = nanmean(concat_pre, 1);
sem_pre = nanstd(concat_pre, 0, 1) ./ sqrt(sum(~isnan(concat_pre), 1));

mean_post = nanmean(concat_post, 1);
sem_post = nanstd(concat_post, 0, 1) ./ sqrt(sum(~isnan(concat_post), 1));

% Plot
t_plot = linspace(0, 1, nPoints);
figure('Color','w'); 
hold on;

% Pre-sleep
% fill([t_plot fliplr(t_plot)], [mean_pre+sem_pre fliplr(mean_pre-sem_pre)], ...
%     [0.5 0.7 0.9], 'EdgeColor','none', 'FaceAlpha', 0.3);
plot(t_plot, mean_pre, 'Color', color_pre, 'LineWidth', 2);

% Post-sleep
% fill([t_plot fliplr(t_plot)], [mean_post+sem_post fliplr(mean_post-sem_post)], ...
%     [0.9 0.5 0.5], 'EdgeColor','none', 'FaceAlpha', 0.3);
plot(t_plot, mean_post, 'Color', color_post, 'LineWidth', 2);

xlabel('Normalized NREM interval duration');
ylabel('DA signal (z-score)');
legend('Pre sleep', 'Post sleep');
title('Average Dopamine Across NREM Intervals (All Sessions)');
grid on;
set(gcf, 'Position', [300, 200, 1000, 600]);

% Save
saveas(gcf,[saveLoc,filesep ,'str_avgDA_NREM_allSessions.png','.png'],'png');

