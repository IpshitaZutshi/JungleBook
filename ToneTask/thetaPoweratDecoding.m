function thetaPoweratDecoding

sess = {'IZ48_230714_sess28'};
decodingSess = {'IZ48\Final\IZ48_230714_sess28'};

twin = 2;
nfreqs = 100;
ncyc = 7;
fBand = [2 150];

for ss = 1:length(sess)
    
    % Load behavior file
    load(strcat('Z:\Homes\zutshi01\Recordings\Auditory_Task\',decodingSess{ss},'\',sess{ss},'.TrialBehavior.Behavior.mat'))
    if ~isfield(behavTrials,'probe')
        behavTrials.probe(1:length(behavTrials.linTrial))=0;
    end

    % Load theta file
    cd(strcat('E:\',decodingSess{ss},'\'))
    [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
    channelNum = 7;
    lfp = bz_GetLFP(channelNum,'noPrompts',true);

    % Load decoding
    decodingPath = strcat('Z:\Homes\zz737\ipshita_data\Auditory_Task\',decodingSess{ss},'\py_data\theta_decoding_lickLoc_y\up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].nc');
    changePointName = strcat('Z:\Homes\zz737\ipshita_data\Auditory_Task\',decodingSess{ss},'\py_data\theta_decoding_lickLoc_y\change_point_posterior_up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].mat');

    posterior_goal =  ncread(decodingPath,'x_position') ;
    posterior_pos =  ncread(decodingPath,'y_position') ;
    post_time =  ncread(decodingPath,'time') ;
    post_pos =  ncread(decodingPath,'y_position_value') ;
    post_goal =  ncread(decodingPath,'x_position_value') ;

    load(changePointName);
    
    ts_dec = [];
    trial_dec = [];

    %% Within each trial, find the change point time that the goal is decoded 
     for tt = 1:length(behavTrials.lickLoc)
        if behavTrials.linTrial(tt)==0 && behavTrials.probe(tt) == 0 && behavTrials.correct(tt)==1
    
            [~,idxstart] = min(abs(post_time-behavTrials.timestamps(tt,1)));
            if post_time(idxstart)<behavTrials.timestamps(tt,1) %Take the next index
                idxstart = idxstart+1;
            end        
            [~,idxend] = min(abs(post_time-behavTrials.timestamps(tt,2)));
            if post_time(idxend)>behavTrials.timestamps(tt,2) %Take the previous index
                idxend = idxend-1;
            end   
    
            [~,decGoal] = max(posterior_goal(:,idxstart:idxend));
    
            %% Get last change point for that trial
            if sum(trial==(tt-1))>0
                curChanges = change_point{trial==(tt-1)};
    
                idxGoal = curChanges(end);
                trialDecGoal = mode(decGoal(curChanges(end)+1:end));
    
                if trialDecGoal==(behavTrials.lickLoc(tt)+1)             
                    ts_dec = [ts_dec post_time(idxGoal+idxstart)];
                    trial_dec = [trial_dec behavTrials.timestamps(tt,2)];
                end
    
            end
        end
     end

     %% Find theta power around these timepoints
     lfpSpect_cp = [];
     lfpSpect_trial = [];
    for iRip = 1:length(ts_dec)
        [~,idx] = min(abs(lfp.timestamps-ts_dec(iRip)));
        intervals = [lfp.timestamps(idx-(twin*1250)) lfp.timestamps(idx+(1*1250))];
        wavelet = bz_WaveSpec(lfp,'frange',fBand,'nfreqs',nfreqs,'ncyc',ncyc,'intervals',intervals);
        lfpSpect_cp(:,:,iRip) = log10(abs(wavelet.data));   
        freqs = wavelet.freqs;                
    end
    for iRip = 1:length(trial_dec)
        [~,idx] = min(abs(lfp.timestamps-trial_dec(iRip)));
        intervals = [lfp.timestamps(idx-(twin*1250)) lfp.timestamps(idx+(1*1250))];
        wavelet = bz_WaveSpec(lfp,'frange',fBand,'nfreqs',nfreqs,'ncyc',ncyc,'intervals',intervals);
        lfpSpect_trial(:,:,iRip) = log10(abs(wavelet.data));   
        freqs = wavelet.freqs;                
    end
    ts = wavelet.timestamps-wavelet.timestamps(1);

end

figure

subplot(2,1,1)
imagesc(ts,freqs,nanmedian(lfpSpect_cp,3)')
set(gca,'YDir','normal')
set(gca,'YScale','log')
ylim([2 150])
colormap jet
hold on
line([2 2 ],[2 150],'Color','w','LineWidth',1.5);

subplot(2,1,2)
imagesc(ts,freqs,nanmedian(lfpSpect_trial,3)')
set(gca,'YDir','normal')
set(gca,'YScale','log')
ylim([2 150])
colormap jet
hold on
line([2 2],[2 150],'Color','w','LineWidth',1.5);

%% ---- Theta (6–12 Hz) band: mean ± SEM across trials ----
thetaBand = [6 12];
thetaIdx  = freqs >= thetaBand(1) & freqs <= thetaBand(2);

% Average across theta freqs for each trial: [time x trials]
theta_cp_trials    = squeeze(nanmean(lfpSpect_cp(:, thetaIdx, :), 2));     % change-point aligned
theta_trial_trials = squeeze(nanmean(lfpSpect_trial(:, thetaIdx, :), 2));  % trial-end aligned

% Ensure 2D [time x nTrials] if only one trial
if isvector(theta_cp_trials),    theta_cp_trials    = theta_cp_trials(:);    end
if isvector(theta_trial_trials), theta_trial_trials = theta_trial_trials(:); end

% Mean and SEM across trials (handle NaNs per timepoint)
mu_cp  = nanmean(theta_cp_trials, 2);
sd_cp  = nanstd(theta_cp_trials, 0, 2);
n_cp   = sum(~isnan(theta_cp_trials), 2);
sem_cp = sd_cp ./ sqrt(max(n_cp,1));

mu_tr  = nanmean(theta_trial_trials, 2);
sd_tr  = nanstd(theta_trial_trials, 0, 2);
n_tr   = sum(~isnan(theta_trial_trials), 2);
sem_tr = sd_tr ./ sqrt(max(n_tr,1));

% Plot mean ± SEM
figure; hold on;
fill([ts; flipud(ts)], [mu_cp - sem_cp; flipud(mu_cp + sem_cp)], [0 0.4470 0.7410], ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(ts, mu_cp, 'LineWidth', 2);

fill([ts; flipud(ts)], [mu_tr - sem_tr; flipud(mu_tr + sem_tr)], [0.8500 0.3250 0.0980], ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(ts, mu_tr, 'LineWidth', 2);

xline(twin, '--k', 'LineWidth', 1); % event at t = twin seconds

xlabel('Time from event (s)');
ylabel('Theta power (6–12 Hz)');
title('Theta-band power around decoded-goal change point vs. trial end');
legend({'CP \pm SEM','CP mean','Trial end \pm SEM','Trial end mean'}, 'Location', 'best');
box on; grid on;

end