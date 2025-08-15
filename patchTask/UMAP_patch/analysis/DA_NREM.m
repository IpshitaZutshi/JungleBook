%% DA throughout NREM
% NEED TO FIX SAVING SO IT HAS PROPER TITLE


%% Inputs
basepath = pwd;

fs_lfp = 1250;
fs_da = 130;

channels = [0:41 43:47 49:55 60 62:69 73:127]; %N17
%lfp = bz_GetLFP(channels);
%lfp = bz_GetLFP(channels, 'restrict', [0 3600]);


%% Step 1: Extract NREM packets
% nremIdx = strcmp(sleepStates.state, 'NREM');
% nremStarts = sleepStates.start(nremIdx);
% nremStops = sleepStates.stop(nremIdx);

save_path = '\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\patchTask\figures';
saveLoc = strcat(save_path,'\sleep_phot');

merge_file = dir(fullfile(basepath,'*.MergePoints.events.mat'));
load(merge_file.name);
sleep_state_file = dir(fullfile(basepath,'*.SleepState.states.mat'));
load(sleep_state_file.name);

for ii = 1:size(MergePoints.foldernames,2)
    cd(basepath);
    if isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep 'top*']))
        cd([basepath filesep MergePoints.foldernames{ii}]); 
        % set pre_post to be 1 for pre behavior sleep, 2 for post behavior sleep  
        if ii == 1
            pre_post = 1;
        else
            pre_post = 2;
        end
        
        curr_folder = pwd;
        if ~isempty (dir(fullfile(curr_folder, '*HPCSleepPhotomSynced.mat'))) 
            load(dir(fullfile(curr_folder,'*HPCSleepPhotomSynced.mat')).name);
        end

        if ~isempty (dir(fullfile(curr_folder, '*StriatumSleepPhotomSynced.mat')))
            load(dir(fullfile(curr_folder,'*StriatumSleepPhotomSynced.mat')).name);
        end

        if pre_post == 1 
            sleep_start = MergePoints.timestamps(1,1); % start of sleep session
            sleep_end = MergePoints.timestamps(1,2); % time stamp of end of sleep session 
        elseif pre_post == 2
            sleep_start = MergePoints.timestamps(3,1); 
            sleep_end = MergePoints.timestamps(3,2);  
        end


        % Extract NREM intervals
        nrem_intervals = SleepState.ints.NREMstate; % Nx2 matrix [start, stop]
        
        % Filter for minimum duration, e.g. 10 seconds
        min_duration = 10;
        durations = nrem_intervals(:,2) - nrem_intervals(:,1);
        valid_idx = durations > min_duration;
        nrem_intervals = nrem_intervals(valid_idx, :);
        
        % Initialize to store dopamine traces
        % We'll normalize each NREM interval to a fixed number of points for averaging
        nPoints = 100; % Number of points to interpolate each NREM interval to
        
        % Assuming daSignal is your dopamine trace sampled at da_fs (e.g. 130 Hz)
        % And timestamps_da is the time vector for daSignal, in seconds
        
        all_da_interp = nan(length(nrem_intervals), nPoints);
        
        % daSignal = hpc_sleep_sync.grabDA_z;
        % plot_color = [0.588235294117647   0.800000000000000   0.345098039215686]; % green
        % timestamps_da = hpc_sleep_sync.timestamps + sleep_start;

        daSignal = striatum_sleep_sync.grabDA_z;
        plot_color = [0.909803921568627   0.290196078431373   0.454901960784314]; % pink
        timestamps_da = striatum_sleep_sync.timestamps + sleep_start;


        for i = 1:length(nrem_intervals)
            t_start = nrem_intervals(i,1);
            t_end = nrem_intervals(i,2);
            
            % Find indices of DA timestamps within this interval
            idx = find(timestamps_da >= t_start & timestamps_da <= t_end);
            if length(idx) < 2
                continue; % skip if too few samples
            end
            
            % Extract DA segment and timestamps
            da_seg = daSignal(idx);
            t_seg = timestamps_da(idx);
            
            % Normalize time to 0-1 scale for interpolation
            t_norm = (t_seg - t_start) / (t_end - t_start);
            
            % Interpolate DA segment to nPoints evenly spaced between 0 and 1
            tq = linspace(0,1,nPoints);
            da_interp = interp1(t_norm, da_seg, tq, 'linear');
            
            all_da_interp(i,:) = da_interp;
        end
        
        % Compute mean and SEM across all NREM intervals
        mean_da = nanmean(all_da_interp,1);
        sem_da = nanstd(all_da_interp,0,1) ./ sqrt(sum(~isnan(all_da_interp),1));
        
        if pre_post == 1 
            pp = 'pre sleep';
        elseif pre_post == 2
            pp = 'post sleep';  
        end

        % Plot
        figure('color', 'white'); 
        hold on;
        t_plot = linspace(0,1,nPoints);
        fill([t_plot, fliplr(t_plot)], [mean_da+sem_da, fliplr(mean_da-sem_da)], ...
            plot_color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        plot(t_plot, mean_da, 'Color', plot_color, 'LineWidth', 2);
        xlabel('Normalized NREM duration');
        ylabel('Dopamine trace z score');
        title('Average dopamine trace across NREM intervals');
        grid on;
        set(gcf,'position',[500,200,1120,840])
        hold off;
        saveas(gcf,[saveLoc,filesep ,'striatum_NREM_packets_22',pp,'.png'],'png');


        % %% Step 2: Define sliding windows within NREM packets
        % winSize = 5;    % seconds
        % stepSize = 1;   % seconds
        % 
        % % Resample DA signal
        % da_resampled = resample(daSignal, fs_lfp, fs_da);
        % 
        % % Frequency bands
        % freqBands = {
        %     'Delta', [1 4];
        %     'Theta', [6 10];
        %     'Beta', [15 30];
        %     'Gamma', [30 80];
        %     'Ripple', [100 250];
        % };
        % nBands = size(freqBands,1);
        % 
        % % Prepare band filters
        % filters = cell(nBands,1);
        % for b = 1:nBands
        %     filters{b} = designfilt('bandpassiir', 'FilterOrder', 4, ...
        %         'HalfPowerFrequency1', freqBands{b,2}(1), ...
        %         'HalfPowerFrequency2', freqBands{b,2}(2), ...
        %         'SampleRate', fs_lfp);
        % end
        % 
        % %% Step 3: Analyze each window
        % allBandPower = [];
        % allDA = [];
        % 
        % for i = 1:length(nremStarts)
        %     t_start = nremStarts(i);
        %     t_stop = nremStops(i);
        % 
        %     t = t_start:stepSize:(t_stop - winSize);
        % 
        %     for j = 1:length(t)
        %         win_t = [t(j), t(j) + winSize];
        %         idx_lfp = round(win_t(1) * fs_lfp):round(win_t(2) * fs_lfp) - 1;
        % 
        %         % Check bounds
        %         if idx_lfp(end) > length(lfp)
        %             continue;
        %         end
        % 
        %         lfp_win = lfp(idx_lfp);
        %         da_win = da_resampled(idx_lfp); % resampled to match
        % 
        %         % Compute DA mean
        %         allDA(end+1) = mean(da_win);
        % 
        %         % Compute band power
        %         for b = 1:nBands
        %             lfp_filt = filtfilt(filters{b}, lfp_win);
        %             power = mean(abs(hilbert(lfp_filt)).^2);
        %             allBandPower(b, end+1) = power;
        %         end
        %     end
        % end
        % 
        % 
        % %% Additional Step: Average DA trace across normalized NREM packets
        % 
        % nPackets = length(nremStarts);
        % nPoints = 100;  % number of points to normalize all packets to
        % 
        % da_packets = nan(nPackets, nPoints);
        % 
        % for i = 1:nPackets
        %     t_start = nremStarts(i);
        %     t_stop = nremStops(i);
        % 
        %     idx_start = round(t_start * fs_da);
        %     idx_stop = round(t_stop * fs_da);
        % 
        %     if idx_stop > length(daSignal)
        %         idx_stop = length(daSignal);
        %     end
        % 
        %     da_segment = daSignal(idx_start:idx_stop);
        % 
        %     % Normalize length by interpolation to nPoints
        %     da_interp = interp1(linspace(0,1,length(da_segment)), da_segment, linspace(0,1,nPoints));
        % 
        %     % Z-score normalization per packet (optional, but recommended)
        %     da_interp = (da_interp - mean(da_interp)) / std(da_interp);
        % 
        %     da_packets(i,:) = da_interp;
        % end
        % 
        % % Average across all packets
        % mean_da = nanmean(da_packets, 1);
        % sem_da = nanstd(da_packets, 0, 1) / sqrt(nPackets);
        % 
        % % Plot average dopamine trace normalized over NREM packets
        % figure;
        % hold on;
        % xvals = linspace(0, 100, nPoints); % % normalized NREM duration (0-100%)
        % fill([xvals fliplr(xvals)], [mean_da+sem_da fliplr(mean_da-sem_da)], ...
        %      [0.6 0.6 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        % plot(xvals, mean_da, 'b-', 'LineWidth', 2);
        % xlabel('Normalized NREM Packet Duration (%)');
        % ylabel('Z-scored Dopamine Signal');
        % title('Average Dopamine Trace Over NREM Packets');
        % grid on;
        % 
        % %% Step 4: Correlate band power and DA
        % rvals = nan(nBands,1);
        % pvals = nan(nBands,1);
        % for b = 1:nBands
        %     [r, p] = corr(allBandPower(b,:)', allDA');
        %     rvals(b) = r;
        %     pvals(b) = p;
        % end
        % 
        % %% Step 5: Plot
        % figure;
        % bar(rvals);
        % hold on;
        % xticks(1:nBands);
        % xticklabels(freqBands(:,1));
        % ylabel('Correlation with DA');
        % title('Band power–DA correlation during NREM packets');
        % ylim([-1 1]);
        % grid on;
        % 
        % % Optional: mark significant
        % for b = 1:nBands
        %     if pvals(b) < 0.05
        %         text(b, rvals(b) + 0.05*sign(rvals(b)), '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
        %     end
        % end
    end
end
