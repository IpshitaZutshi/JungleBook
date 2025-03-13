%% Analyze the stimulation with the photometry data
% 

basepath = pwd;
%photometry_file = dir(fullfile(basepath, '*N11_striatum.mat'));
%load(photometry_file.name);

sampling_rate = 130; % sampling rate of photometry set up - 130



%% Plot photometry with behavior
%{
hold on
figure(1)
plot(photom_var.timestamps, photom_var.grabDA_z, 'g'); %plot(photometry(:,1), photometry(:,2), 'g');
ylabel('z score');
xlabel('time');
xlim([photom_var.timestamps(1), photom_var.timestamps(length(photom_var.timestamps(~isnan(photom_var.timestamps))))]);
for j = 1:length(rewarded_times)
    xline(rewarded_times(j, 1), '-b');
end
for k = 1:length(nonrewarded_times)
    xline(nonrewarded_times(k, 1), '-r');
end
hold off

%}


%% Average photometry around stim pulses

window = 10; % window of time around lick to average
samples = window*sampling_rate;

zscore_matrix = nan(length(photometryData.pulse_times), (samples*2)+1); 
photometryData.pulse_times = photometryData.pulse_times';
% convert to seconds!
photometryData.pulse_times = photometryData.pulse_times / 1000;

% average photometry data within a specified time window around stim
for j = 1:length(photometryData.pulse_times)
    curr_stim_time = photometryData.pulse_times(j);
    [~, stim_idx] = min(abs(photometryData.timestamps - photometryData.pulse_times(j)));
    start_idx = stim_idx - samples;
    end_idx = stim_idx + samples;
    if start_idx >= 1 && end_idx <= length(photometryData.timestamps)
        zscore_matrix(j, :) = photometryData.grabDA_z(start_idx:end_idx);
    end
end

%new_zscore_matrix = zscore_matrix;
%new_zscore_matrix(any(isnan(new_zscore_matrix), 2), :) = [];

zscore_matrix(any(isnan(zscore_matrix), 2), :) = [];
med_z_stim = median(zscore_matrix, 1); % median at each timepoint
time = linspace(-window, window, ((samples*2)+1));

% calculate confidence intervals
N = height(zscore_matrix);                          % Number of eExperimentsn In Data Set
avg_z_stim = mean(zscore_matrix, 1);              % Mean Of All Experiments At Each Value Of ,x 
stim_SEM = std(zscore_matrix, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
stim_CI95 = bsxfun(@times, stim_SEM, CI95(:));  % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of exw


%% Plot  average photometry level around stim pulses
figure;
hold on
plot(time, med_z_stim, 'g', 'LineWidth', 2);
fill([time,fliplr(time)], [(stim_CI95(1,:)+avg_z_stim),fliplr((stim_CI95(2,:)+avg_z_stim))], 'b', 'EdgeColor','none', 'FaceAlpha',0.25)
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title('Average Z-score Around stims');
grid on;
hold off



figure;
customGreen = [0.2039 0.7294 0.2039];
plot(photometryData.timestamps, photometryData.grabDA_z, 'Color', customGreen)
xlabel('time (s)');
ylabel('z-score');
title('GRAB-DA in hippocampus');








