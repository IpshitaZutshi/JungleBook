% AGGREGATING BUT KEEPING SESSIONS SEPARATE ( OR MICE?)
% 
% z score across sleep sessions
%
% USAGE
%    z score photometry data from across sleep sessions and plot the 
%    dopamine traces around ripples, in the HPC and striatum
% 
%
% INPUTS 
%    
%
%    =========================================================================


direc = '\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\patchTask';

sessions_N11 = {'N11\Final\N11_250312_sess16', ...
    'N11\Final\N11_250313_sess17', ...
    'N11\Final\N11_250314_sess18', ...
    'N11\Final\N11_250318_sess19' ...
    'N11\Final\N11_250319_sess20', ...
    'N11\Final\N11_250321_sess22', ...
    'N11\Final\N11_250331_sess24', ...
    'N11\Final\N11_250401_sess25', ...
    'N11\Final\N11_250403_sess26', ...
    'N11\Final\N11_250407_sess27', ...
    'N11\Final\N11_250408_sess28', ...
    'N11\Final\N11_250410_sess30', ...
    'N11\Final\N11_250411_sess31'};

sessions_N14 = {'N14\N14_250429_sess8', ...
    'N14\N14_250501_sess10', ...
    'N14\N14_250509_sess15', ...
    'N14\N14_250510_sess16', ...
    'N14\N14_250511_sess17', ...
    'N14\N14_250512_sess18', ...
    'N14\N14_250519_sess20', ...
    'N14\N14_250520_sess21', ...
    'N14\N14_250521_sess22'};
%    'N14\N14_250507_sess13', ...

sessions_N17 = {'N17\N17_250430_sess9', ...
    'N17\N17_250501_sess10', ...
    'N17\N17_250509_sess15', ...
    'N17\N17_250510_sess16', ...
    'N17\N17_250511_sess17', ...
    'N17\N17_250512_sess18', ...
    'N17\N17_250513_sess19', ...
    'N17\N17_250519_sess21'};

% Create a structure to hold session groups and labels
mice = {'N11', 'N14', 'N17'};
sessions_grouped = {sessions_N11, sessions_N14, sessions_N17};

% Colors for plotting each subject group
colors = [0.2 0.6 0.8;  % N11 - blueish
          0.9 0.3 0.3;  % N14 - reddish
          0.3 0.8 0.3]; % N17 - greenish

% Initialize storage for each group
all_data = struct();

for m = 1:length(mice)
    group_name = mice{m};
    sessions = sessions_grouped{m};
    
    all_pre_ripples_hpc = {};
    all_post_ripples_hpc = {};
    all_pre_ripples_striatum = {};
    all_post_ripples_striatum = {};
    
    for s = 1:length(sessions)
        sessionpath = [direc filesep sessions{s}];
        cd(sessionpath);
        disp(['Collecting data from ' sessions{s}]);

        session_data = dualSleepPhotometry('showfig', false); 

        % HPC data
        if isfield(session_data, 'hpc')
            if isfield(session_data.hpc, 'pre_sleep')
                all_pre_ripples_hpc{end+1} = session_data.hpc.pre_sleep;
            end
            if isfield(session_data.hpc, 'post_sleep')
                all_post_ripples_hpc{end+1} = session_data.hpc.post_sleep;
            end
        end

        % Striatum data
        if isfield(session_data, 'striatum')
            if isfield(session_data.striatum, 'pre_sleep')
                all_pre_ripples_striatum{end+1} = session_data.striatum.pre_sleep;
            end
            if isfield(session_data.striatum, 'post_sleep')
                all_post_ripples_striatum{end+1} = session_data.striatum.post_sleep;
            end
        end
    end

    % Concatenate all ripples for this group
    pre_ripples_hpc_all = cat(1, all_pre_ripples_hpc{:});
    post_ripples_hpc_all = cat(1, all_post_ripples_hpc{:});
    pre_ripples_striatum_all = cat(1, all_pre_ripples_striatum{:});
    post_ripples_striatum_all = cat(1, all_post_ripples_striatum{:});

    % Calculate means, medians, SEM and 95% CI for HPC
    N_pre_hpc = size(pre_ripples_hpc_all, 1);
    N_post_hpc = size(post_ripples_hpc_all, 1);
    
    all_data.(group_name).pre_hpc.mean = mean(pre_ripples_hpc_all, 1);
    all_data.(group_name).pre_hpc.median = median(pre_ripples_hpc_all, 1);
    all_data.(group_name).pre_hpc.sem = std(pre_ripples_hpc_all, 1) / sqrt(N_pre_hpc);
    all_data.(group_name).pre_hpc.ci95 = tinv([0.025 0.975], N_pre_hpc-1);
    all_data.(group_name).pre_hpc.n = N_pre_hpc;

    all_data.(group_name).post_hpc.mean = mean(post_ripples_hpc_all, 1);
    all_data.(group_name).post_hpc.median = median(post_ripples_hpc_all, 1);
    all_data.(group_name).post_hpc.sem = std(post_ripples_hpc_all, 1) / sqrt(N_post_hpc);
    all_data.(group_name).post_hpc.ci95 = tinv([0.025 0.975], N_post_hpc-1);
    all_data.(group_name).post_hpc.n = N_post_hpc;

    % Calculate means, medians, SEM and 95% CI for Striatum
    N_pre_str = size(pre_ripples_striatum_all, 1);
    N_post_str = size(post_ripples_striatum_all, 1);
    
    all_data.(group_name).pre_striatum.mean = mean(pre_ripples_striatum_all, 1);
    all_data.(group_name).pre_striatum.median = median(pre_ripples_striatum_all, 1);
    all_data.(group_name).pre_striatum.sem = std(pre_ripples_striatum_all, 1) / sqrt(N_pre_str);
    all_data.(group_name).pre_striatum.ci95 = tinv([0.025 0.975], N_pre_str-1);
    all_data.(group_name).pre_striatum.n = N_pre_str;

    all_data.(group_name).post_striatum.mean = mean(post_ripples_striatum_all, 1);
    all_data.(group_name).post_striatum.median = median(post_ripples_striatum_all, 1);
    all_data.(group_name).post_striatum.sem = std(post_ripples_striatum_all, 1) / sqrt(N_post_str);
    all_data.(group_name).post_striatum.ci95 = tinv([0.025 0.975], N_post_str-1);
    all_data.(group_name).post_striatum.n = N_post_str;
end

%% Plotting

sampling_rate = 130;
window = 5;
samples = window * sampling_rate;
time = linspace(-window, window, 2*samples+1);

figure('color','white');

% HPC pre-sleep subplot
subplot(2,2,1);
hold on;
cmap = summer(length(mice)+1);
for m = 1:length(mice)
    group_name = mice{m};
    med_data = all_data.(group_name).pre_hpc.median;
    sem_data = all_data.(group_name).pre_hpc.sem;
    N = all_data.(group_name).pre_hpc.n;
    ci95 = all_data.(group_name).pre_hpc.ci95;

    CI_lower = med_data + ci95(1)*sem_data;
    CI_upper = med_data + ci95(2)*sem_data;

    plot(time, med_data, 'LineWidth', 2, 'Color', cmap(m, :));
    fill([time fliplr(time)], [CI_upper fliplr(CI_lower)], cmap(m, :), 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
xline(0, '--r', 'LineWidth', 1);
title('HPC Pre-Sleep');
xlabel('Time (s)');
ylabel('Avg Z-score');
legend(mice);
grid on;
hold off;

% HPC post-sleep subplot
subplot(2,2,2);
hold on;
cmap = summer(length(mice)+1);
for m = 1:length(mice)
    group_name = mice{m};
    med_data = all_data.(group_name).post_hpc.median;
    sem_data = all_data.(group_name).post_hpc.sem;
    N = all_data.(group_name).post_hpc.n;
    ci95 = all_data.(group_name).post_hpc.ci95;

    CI_lower = med_data + ci95(1)*sem_data;
    CI_upper = med_data + ci95(2)*sem_data;

    plot(time, med_data, 'LineWidth', 2, 'Color', cmap(m,:));
    fill([time fliplr(time)], [CI_upper fliplr(CI_lower)], cmap(m,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
xline(0, '--r', 'LineWidth', 1);
title('HPC Post-Sleep');
xlabel('Time (s)');
ylabel('Avg Z-score');
legend(mice);
grid on;
hold off;

% Striatum pre-sleep subplot
subplot(2,2,3);
hold on;
cmap = spring(length(mice)+1);
for m = 1:length(mice)
    group_name = mice{m};
    med_data = all_data.(group_name).pre_striatum.median;
    sem_data = all_data.(group_name).pre_striatum.sem;
    N = all_data.(group_name).pre_striatum.n;
    ci95 = all_data.(group_name).pre_striatum.ci95;

    CI_lower = med_data + ci95(1)*sem_data;
    CI_upper = med_data + ci95(2)*sem_data;

    plot(time, med_data, 'LineWidth', 2, 'Color', cmap(m,:));
    fill([time fliplr(time)], [CI_upper fliplr(CI_lower)], cmap(m,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
xline(0, '--r', 'LineWidth', 1);
title('Striatum Pre-Sleep');
xlabel('Time (s)');
ylabel('Avg Z-score');
legend(mice);
grid on;
hold off;

% Striatum post-sleep subplot
subplot(2,2,4);
hold on;
cmap = spring(length(mice)+1);
for m = 1:length(mice)
    group_name = mice{m};
    med_data = all_data.(group_name).post_striatum.median;
    sem_data = all_data.(group_name).post_striatum.sem;
    N = all_data.(group_name).post_striatum.n;
    ci95 = all_data.(group_name).post_striatum.ci95;

    CI_lower = med_data + ci95(1)*sem_data;
    CI_upper = med_data + ci95(2)*sem_data;

    plot(time, med_data, 'LineWidth', 2, 'Color', cmap(m,:));
    fill([time fliplr(time)], [CI_upper fliplr(CI_lower)], cmap(m,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
xline(0, '--r', 'LineWidth', 1);
title('Striatum Post-Sleep');
xlabel('Time (s)');
ylabel('Avg Z-score');
legend(mice);
grid on;
hold off;

% =========================================================================

% smoothed
figure('color','white');

% HPC pre-sleep subplot
subplot(2,2,1);
hold on;
cmap = summer(length(mice)+1);
for m = 1:length(mice)
    group_name = mice{m};
    med_data = all_data.(group_name).pre_hpc.median;
    sem_data = all_data.(group_name).pre_hpc.sem;
    N = all_data.(group_name).pre_hpc.n;
    ci95 = all_data.(group_name).pre_hpc.ci95;

    CI_lower = med_data + ci95(1)*sem_data;
    CI_upper = med_data + ci95(2)*sem_data;

    plot(time, smoothdata(med_data, "sgolay"), 'LineWidth', 2, 'Color', cmap(m, :));
    fill([time fliplr(time)], [smoothdata(CI_upper, 2, "sgolay") fliplr(smoothdata(CI_lower, 2, "sgolay"))], cmap(m, :), 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
xline(0, '--r', 'LineWidth', 1);
title('HPC Pre-Sleep');
xlabel('Time (s)');
ylabel('Avg Z-score');
legend(mice);
grid on;
hold off;

% HPC post-sleep subplot
subplot(2,2,2);
hold on;
cmap = summer(length(mice)+1);
for m = 1:length(mice)
    group_name = mice{m};
    med_data = all_data.(group_name).post_hpc.median;
    sem_data = all_data.(group_name).post_hpc.sem;
    N = all_data.(group_name).post_hpc.n;
    ci95 = all_data.(group_name).post_hpc.ci95;

    CI_lower = med_data + ci95(1)*sem_data;
    CI_upper = med_data + ci95(2)*sem_data;

    plot(time, smoothdata(med_data, "sgolay"), 'LineWidth', 2, 'Color', cmap(m,:));
    fill([time fliplr(time)], [smoothdata(CI_upper, 2, "sgolay") fliplr(smoothdata(CI_lower, 2, "sgolay"))], cmap(m,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
xline(0, '--r', 'LineWidth', 1);
title('HPC Post-Sleep');
xlabel('Time (s)');
ylabel('Avg Z-score');
legend(mice);
grid on;
hold off;

% Striatum pre-sleep subplot
subplot(2,2,3);
hold on;
cmap = spring(length(mice)+1);
for m = 1:length(mice)
    group_name = mice{m};
    med_data = all_data.(group_name).pre_striatum.median;
    sem_data = all_data.(group_name).pre_striatum.sem;
    N = all_data.(group_name).pre_striatum.n;
    ci95 = all_data.(group_name).pre_striatum.ci95;

    CI_lower = med_data + ci95(1)*sem_data;
    CI_upper = med_data + ci95(2)*sem_data;

    plot(time, smoothdata(med_data, "sgolay"), 'LineWidth', 2, 'Color', cmap(m,:));
    fill([time fliplr(time)], [smoothdata(CI_upper, 2, "sgolay") fliplr(smoothdata(CI_lower, 2, "sgolay"))], cmap(m,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
xline(0, '--r', 'LineWidth', 1);
title('Striatum Pre-Sleep');
xlabel('Time (s)');
ylabel('Avg Z-score');
legend(mice);
grid on;
hold off;

% Striatum post-sleep subplot
subplot(2,2,4);
hold on;
cmap = spring(length(mice)+1);
for m = 1:length(mice)
    group_name = mice{m};
    med_data = all_data.(group_name).post_striatum.median;
    sem_data = all_data.(group_name).post_striatum.sem;
    N = all_data.(group_name).post_striatum.n;
    ci95 = all_data.(group_name).post_striatum.ci95;

    CI_lower = med_data + ci95(1)*sem_data;
    CI_upper = med_data + ci95(2)*sem_data;

    plot(time, smoothdata(med_data, "sgolay"), 'LineWidth', 2, 'Color', cmap(m,:));
    fill([time fliplr(time)], [smoothdata(CI_upper, 2, "sgolay") fliplr(smoothdata(CI_lower, 2, "sgolay"))], cmap(m,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
xline(0, '--r', 'LineWidth', 1);
title('Striatum Post-Sleep');
xlabel('Time (s)');
ylabel('Avg Z-score');
legend(mice);
grid on;
hold off;


