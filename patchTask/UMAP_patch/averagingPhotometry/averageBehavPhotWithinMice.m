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
%     'N14\N14_250507_sess13', ...

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

% Initialize storage for each group
all_data = struct();

for m = 1:length(mice)
    group_name = mice{m};
    sessions = sessions_grouped{m};
    
    % Initialize containers
    types = {'rewarded', 'nonrewarded', ...
             'expected_reward', 'unexpected_reward', ...
             'expected_nonreward', 'unexpected_nonreward', ...
             'awake_ripples'};

    regions = {'hpc', 'striatum'};
    
    % Initialize dynamic storage
    for t = 1:length(types)
        for r = 1:length(regions)
            all_data.(group_name).([types{t} '_' regions{r}]) = {};
        end
    end

    % Loop through sessions
    for s = 1:length(sessions)
        sessionpath = [direc filesep sessions{s}];
        cd(sessionpath);
        disp(['Collecting data from ' sessions{s}]);

        session_data = photometry_behav('showfig', false); 

        for r = 1:length(regions)
            region = regions{r};
            if isfield(session_data, region)
                region_data = session_data.(region);
                for t = 1:length(types)
                    type = types{t};
                    if isfield(region_data, type)
                        all_data.(group_name).([type '_' region]){end+1} = region_data.(type);
                    end
                end
            end
        end
    end

    % Aggregate and compute stats
    for t = 1:length(types)
        for r = 1:length(regions)
            fieldname = [types{t} '_' regions{r}];
            if isempty(all_data.(group_name).(fieldname))
                continue;
            end
            data_all = cat(1, all_data.(group_name).(fieldname){:});
            N = size(data_all, 1);
            sem = std(data_all, 1) / sqrt(N);
            ci95 = tinv([0.025 0.975], N-1);
            all_data.(group_name).(fieldname) = struct( ...
                'mean', mean(data_all, 1), ...
                'median', median(data_all, 1), ...
                'sem', sem, ...
                'ci95', ci95, ...
                'n', N ...
            );
        end
    end
end

%% Plotting

save_path = '\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\patchTask\figures';
saveLoc = strcat(save_path,'\behavior');

sampling_rate = 130;
window = 3.1;
samples = window * sampling_rate;
time = linspace(-window, window, 2*samples+1);

% expected HPC reward
figure('color', 'white');
subplot(1,1,1);
hold on;
cmap = summer(length(mice)+1);

for m = 1:length(mice)
    group_name = mice{m};
    d = all_data.(group_name).expected_reward_hpc;
    if isempty(d), continue; end
    CI_lower = d.median + d.ci95(1)*d.sem;
    CI_upper = d.median + d.ci95(2)*d.sem;

    plot(time, d.median, 'Color', cmap(m,:), 'LineWidth', 2);
    fill([time fliplr(time)], [CI_upper fliplr(CI_lower)], cmap(m,:), ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

xline(0, '--r', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Avg Z-score');
title('HPC - Expected Reward');
legend(mice);
grid on;
hold off;


% expected striatum reward
figure('color', 'white');
subplot(1,1,1);
hold on;
cmap = spring(length(mice)+1);

for m = 1:length(mice)
    group_name = mice{m};
    d = all_data.(group_name).expected_reward_striatum;
    if isempty(d), continue; end
    CI_lower = d.median + d.ci95(1)*d.sem;
    CI_upper = d.median + d.ci95(2)*d.sem;

    plot(time, d.median, 'Color', cmap(m,:), 'LineWidth', 2);
    fill([time fliplr(time)], [CI_upper fliplr(CI_lower)], cmap(m,:), ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

xline(0, '--r', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Avg Z-score');
title('Striatum - Expected Reward');
legend(mice);
grid on;
hold off;


%% AWAKE RIPPLES
window = 5;
samples = window * sampling_rate;
time = linspace(-window, window, 2*samples+1);
% awake ripples HPC 
figure('color', 'white');
subplot(1,1,1);
hold on;
cmap = summer(length(mice)+1);

for m = 1:length(mice)
    group_name = mice{m};
    d = all_data.(group_name).awake_ripples_hpc;
    if isempty(d), continue; end
    CI_lower = d.median + d.ci95(1)*d.sem;
    CI_upper = d.median + d.ci95(2)*d.sem;

    plot(time, d.median, 'Color', cmap(m,:), 'LineWidth', 2);
    fill([time fliplr(time)], [CI_upper fliplr(CI_lower)], cmap(m,:), ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

xline(0, '--r', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Avg Z-score');
title('HPC - Awake ripples');
legend(mice);
grid on;
hold off;
set(gcf,'position',[500,200,1120,840])
saveas(gcf,[saveLoc,filesep ,'hpc_awake_ripples', '.png'],'png');


% awake ripples striatum 
figure('color', 'white');
subplot(1,1,1);
hold on;
cmap = spring(length(mice)+1);

for m = 1:length(mice)
    group_name = mice{m};
    d = all_data.(group_name).awake_ripples_striatum;
    if isempty(d), continue; end
    CI_lower = d.median + d.ci95(1)*d.sem;
    CI_upper = d.median + d.ci95(2)*d.sem;

    plot(time, d.median, 'Color', cmap(m,:), 'LineWidth', 2);
    fill([time fliplr(time)], [CI_upper fliplr(CI_lower)], cmap(m,:), ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

xline(0, '--r', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Avg Z-score');
title('Striatum - Awake Ripples');
legend(mice);
grid on;
hold off;

set(gcf,'position',[500,200,1120,840])
saveas(gcf,[saveLoc,filesep ,'striatum_awake_ripples', '.png'],'png');


%saveas(gcf,[saveLoc,filesep ,'striatum_behav_phot_sep_mouse', '.png'],'png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other way to do it
%{
for m = 1:length(mice)
    group_name = mice{m};
    sessions = sessions_grouped{m};
    
    all_rewarded_hpc = {};
    all_nonrewarded_hpc = {};
    all_expected_reward_hpc = {};
    all_unexpected_reward_hpc = {};
    all_expected_nonreward_hpc = {};
    all_unexpected_nonreward_hpc = {};
    all_awake_ripples_hpc = {};
    
    all_rewarded_striatum = {};
    all_nonrewarded_striatum = {};
    all_expected_reward_striatum = {};
    all_unexpected_reward_striatum = {};
    all_expected_nonreward_striatum = {};
    all_unexpected_nonreward_striatum = {};
    all_awake_ripples_striatum = {};
    
    for s = 1:length(sessions)
        sessionpath = [direc filesep sessions{s}];
        cd(sessionpath);
        disp(['Collecting data from ' sessions{s}]);

        session_data = photometry_behav('showfig', false); 

        % HPC data
        if isfield(session_data, 'hpc')
            if isfield(session_data.hpc, 'rewarded')
                all_rewarded_hpc{end+1} = session_data.hpc.rewarded;
            end
            if isfield(session_data.hpc, 'nonrewarded')
                all_nonrewarded_hpc{end+1} = session_data.hpc.nonrewarded;
            end
            if isfield(session_data.hpc, 'expected_reward')
                all_expected_reward_hpc{end+1} = session_data.hpc.expected_reward;
            end
            if isfield(session_data.hpc, 'expected_no_reward')
                all_expected_nonreward_hpc{end+1} = session_data.hpc.expected_no_reward;
            end
            if isfield(session_data.hpc, 'unexpected_reward')
                all_unexpected_reward_hpc{end+1} = session_data.hpc.unexpected_reward;
            end
            if isfield(session_data.hpc, 'unexpected_no_reward')
                all_unexpected_nonreward_hpc{end+1} = session_data.hpc.unexpected_no_reward;
            end
            if isfield(session_data.hpc, 'awake_ripples')
                all_awake_ripples_hpc{end+1} = session_data.hpc.awake_ripples;
            end
        end


        % Striatum data
        if isfield(session_data, 'striatum')
            if isfield(session_data.striatum, 'rewarded')
                all_rewarded_striatum{end+1} = session_data.striatum.rewarded;
            end
            if isfield(session_data.striatum, 'nonrewarded')
                all_nonrewarded_striatum{end+1} = session_data.striatum.nonrewarded;
            end
            if isfield(session_data.striatum, 'expected_reward')
                all_expected_reward_striatum{end+1} = session_data.striatum.expected_reward;
            end
            if isfield(session_data.striatum, 'expected_no_reward')
                all_expected_nonreward_striatum{end+1} = session_data.striatum.expected_no_reward;
            end
            if isfield(session_data.striatum, 'unexpected_reward')
                all_unexpected_reward_striatum{end+1} = session_data.striatum.unexpected_reward;
            end
            if isfield(session_data.striatum, 'unexpected_no_reward')
                all_unexpected_nonreward_striatum{end+1} = session_data.striatum.unexpected_no_reward;
            end
            if isfield(session_data.striatum, 'awake_ripples')
                all_awake_ripples_striatum{end+1} = session_data.striatum.awake_ripples;
            end
        end
    end
end

%}

