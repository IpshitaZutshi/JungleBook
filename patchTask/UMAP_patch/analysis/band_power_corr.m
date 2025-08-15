
% correlate band powers with da during sleep
basepath = pwd;

ripple_file = dir(fullfile(basepath, '*ripples.events.mat'));
load(ripple_file.name)

session_file = dir(fullfile(basepath, '*sessionInfo.mat'));
load(session_file.name)

sleep_state_file = dir(fullfile(basepath, '*SleepState.states.mat'));
load(sleep_state_file.name)


if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
    load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
    
    tokens = regexp(sessionInfo.FileName, 'sess(\d+)$', 'tokens');
    sessionStr = tokens{1}{1};
    sleep_photometry.session = sessionStr; 
    
    for ii = 1:size(MergePoints.foldernames,2)
        cd(basepath);
        channels = [0:41 43:47 49:55 60 62:69 73:127]; %N17
        channels = [70];
        %lfp = bz_GetLFP(channels);
        lfpStruct = bz_GetLFP(channels, 'restrict', [0 3600]);

        if isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep 'top*']))
            cd([basepath filesep MergePoints.foldernames{ii}]); 
            % set pre_post to be 1 for pre behavior sleep, 2 for post behavior sleep  
            if ii == 1
                pre_post = 1;
            else
                pre_post = 2;
            end
            
            curr_folder = pwd;
            if ~isempty (dir(fullfile(curr_folder, '*HPC*photometry.mat'))) 
                photometry_file = dir(fullfile(curr_folder, '*HPCSleepPhotomSynced.mat'));
                load(photometry_file.name);
                photometry_hpc = hpc_sleep_sync;
            end
    
            if ~isempty (dir(fullfile(curr_folder, '*striatum*photometry.mat')))
                photometry_file = dir(fullfile(curr_folder, '*StriatumSleepPhotomSynced.mat'));
                load(photometry_file.name);
                photometry_striatum = striatum_sleep_sync;
            end

            if pre_post == 1 
                sleep_start = MergePoints.timestamps(1,1); % start of sleep session
                sleep_end = MergePoints.timestamps(1,2); % time stamp of end of sleep session 
            elseif pre_post == 2
                sleep_start = MergePoints.timestamps(3,1); 
                sleep_end = MergePoints.timestamps(3,2);  
            end

            lfp = double(lfpStruct.data);
            fs = lfpStruct.samplingRate;
            
            [deltaPower, deltaTime] = computeBandPower(lfp, fs, [0.5 4], 2);
            [sigmaPower, ~] = computeBandPower(lfp, fs, [10 20], 2);
            [gammaPower, ~] = computeBandPower(lfp, fs, [30 100], 2);
            
            daT = photometry_hpc.timestamps + sleep_start;
            da = photometry_hpc.grabDA_z;  % or grabDA_df
            
            % Interpolate DA to match power time
            da_interp = interp1(daT, da, deltaTime, 'linear', 'extrap');
            
            
            corrDelta = corr(da_interp(:), deltaPower(:), 'rows', 'complete');
            corrSigma = corr(da_interp(:), sigmaPower(:), 'rows', 'complete');
            corrGamma = corr(da_interp(:), gammaPower(:), 'rows', 'complete');
            
            fprintf('Correlations:\nDelta: %.3f\nSigma: %.3f\nGamma: %.3f\n', ...
                    corrDelta, corrSigma, corrGamma);
            
            
            packets = SleepState.ints.NREMstate; % Nx2 start/end
            nBins = 100; % normalize to 100 bins
            
            allDA = [];
            allDelta = [];
            allSigma = [];
            
            for i = 1:size(packets,1)
                pktStart = packets(i,1);
                pktEnd = packets(i,2);
            
                t_norm = linspace(pktStart, pktEnd, nBins);
            
                % Interpolate signals onto normalized timebase
                daPkt = interp1(daT, da, t_norm, 'linear', NaN);
                deltaPkt = interp1(deltaTime, deltaPower, t_norm, 'linear', NaN);
                sigmaPkt = interp1(deltaTime, sigmaPower, t_norm, 'linear', NaN);
            
                allDA = [allDA; daPkt];
                allDelta = [allDelta; deltaPkt];
                allSigma = [allSigma; sigmaPkt];
            end
            
            % Compute mean trace
            meanDA = nanmean(allDA);
            meanDelta = nanmean(allDelta);
            meanSigma = nanmean(allSigma);
            
            % Plot
            figure('color', 'white');
            subplot(3,1,1); plot(meanDA); title('DA across NREM Packets');
            subplot(3,1,2); plot(meanDelta); title('Delta Power across NREM Packets');
            subplot(3,1,3); plot(meanSigma); title('Sigma Power across NREM Packets');
            
            
            [maxLagSec] = 10;
            Fs = 1 / mean(diff(deltaTime));
            maxLagSamples = round(maxLagSec * Fs);
            [xc,lags] = xcorr(zscore(da_interp), zscore(sigmaPower), maxLagSamples, 'coeff');
            timeLags = lags / Fs;
            
            figure('color', 'white'); plot(timeLags, xc);
            xlabel('Lag (s)'); ylabel('Cross-corr (DA vs Sigma)');


        end
    end
end



