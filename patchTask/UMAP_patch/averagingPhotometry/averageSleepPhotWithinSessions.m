
%% z score across sleep sessions
%
% USAGE
%    z score photometry data from within sleep sessions and plot the 
%    dopamine traces around ripples, in the HPC and striatum
%
%    plots each session separately
% 
%
% INPUTS 
%    
%
%    =========================================================================


%% Settings

direc = '\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\patchTask';

sessions = {'N17\N17_250430_sess9', ...
    'N17\N17_250501_sess10', ...
    'N17\N17_250509_sess15', ...
    'N17\N17_250510_sess16', ...
    'N17\N17_250511_sess17', ...
    'N17\N17_250512_sess18', ...
    'N17\N17_250513_sess19', ...
    'N17\N17_250519_sess21'};

sessions = {'N11\Final\N11_250312_sess16', ...
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


sampling_rate = 130;
window = 5;
samples = window * sampling_rate;
time = linspace(-window, window, 2*samples+1);

%% Collect data
session_data_struct = struct();
for i = 1:length(sessions)
    session_path = fullfile(direc, sessions{i});
    cd(session_path);
    disp(['Loading session: ', sessions{i}]);

    data = dualSleepPhotometry('showfig', false);
    session_name = matlab.lang.makeValidName(strrep(sessions{i}, '\', '_'));

    if isfield(data, 'hpc')
        if isfield(data.hpc, 'pre_sleep')
            session_data_struct.(session_name).hpc_pre = data.hpc.pre_sleep;
        end
        if isfield(data.hpc, 'post_sleep')
            session_data_struct.(session_name).hpc_post = data.hpc.post_sleep;
        end
    end
    if isfield(data, 'striatum')
        if isfield(data.striatum, 'pre_sleep')
            session_data_struct.(session_name).str_pre = data.striatum.pre_sleep;
        end
        if isfield(data.striatum, 'post_sleep')
            session_data_struct.(session_name).str_post = data.striatum.post_sleep;
        end
    end
end

%% Plot each session trace (HPC & Striatum, pre & post sleep)
fields = fieldnames(session_data_struct);
colors = summer(length(fields)+1); % fix so it's different for hpc/striatum

figure('Color', 'w');
titles = {'HPC Pre-Sleep', 'HPC Post-Sleep', 'Striatum Pre-Sleep', 'Striatum Post-Sleep'};
field_targets = {'hpc_pre', 'hpc_post', 'str_pre', 'str_post'};

for subplot_idx = 1:4
    subplot(2,2,subplot_idx); hold on;
    title(titles{subplot_idx});
    xlabel('Time (s)');
    ylabel('Z-score');
    legend_labels = {};

    for i = 1:length(fields)
        this_session = session_data_struct.(fields{i});
        if isfield(this_session, field_targets{subplot_idx})
            data_matrix = this_session.(field_targets{subplot_idx});
            if ~isempty(data_matrix)
                med_trace = median(data_matrix, 1);
                smoothed = smoothdata(med_trace, 'sgolay');
                plot(time, smoothed, 'LineWidth', 2, 'Color', colors(i,:));
                legend_labels{end+1} = fields{i}; %#ok<SAGROW>
            end
        end
    end

    xline(0, '--r');
    legend(legend_labels, 'Interpreter', 'none', 'FontSize', 6);
    grid on;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


direc = '\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\patchTask';
%{
sessions_N17 = {'N17\N17_250430_sess9', ...
    'N17\N17_250501_sess10', ...
    'N17\N17_250509_sess15', ...
    'N17\N17_250510_sess16', ...
    'N17\N17_250511_sess17', ...
    'N17\N17_250512_sess18', ...
    'N17\N17_250513_sess19', ...
    'N17\N17_250519_sess21'};


% Initialize storage
session_data_struct = struct();

for i = 1:length(sessions_N17)
    session_path = fullfile(direc, sessions_N17{i});
    cd(session_path);
    disp(['Loading session: ', sessions_N17{i}]);

    data = dualSleepPhotometry('showfig', false);

    session_name = strrep(sessions_N17{i}, '\', '_');

    if isfield(data, 'hpc')
        if isfield(data.hpc, 'pre_sleep')
            session_data_struct.(session_name).hpc_pre = data.hpc.pre_sleep;
        end
        if isfield(data.hpc, 'post_sleep')
            session_data_struct.(session_name).hpc_post = data.hpc.post_sleep;
        end
    end

    if isfield(data, 'striatum')
        if isfield(data.striatum, 'pre_sleep')
            session_data_struct.(session_name).str_pre = data.striatum.pre_sleep;
        end
        if isfield(data.striatum, 'post_sleep')
            session_data_struct.(session_name).str_post = data.striatum.post_sleep;
        end
    end
end

%% Plotting session-wise traces
sampling_rate = 130;
window = 5;
samples = window * sampling_rate;
time = linspace(-window, window, 2*samples+1);

fields = fieldnames(session_data_struct);
colors = parula(length(fields));

figure('Color', 'white');
subplot(2,2,1); hold on;
title('HPC Pre-Sleep'); xlabel('Time (s)'); ylabel('Z-score');
for i = 1:length(fields)
    trace = session_data_struct.(fields{i}).hpc_pre;
    if ~isempty(trace)
        med_trace = median(trace,1);
        plot(time, med_trace, 'LineWidth', 2, 'Color', colors(i,:));
    end
end
xline(0, '--r'); legend(fields, 'Interpreter', 'none'); grid on;

subplot(2,2,2); hold on;
title('HPC Post-Sleep'); xlabel('Time (s)'); ylabel('Z-score');
for i = 1:length(fields)
    trace = session_data_struct.(fields{i}).hpc_post;
    if ~isempty(trace)
        med_trace = median(trace,1);
        plot(time, med_trace, 'LineWidth', 2, 'Color', colors(i,:));
    end
end
xline(0, '--r'); legend(fields, 'Interpreter', 'none'); grid on;

subplot(2,2,3); hold on;
title('Striatum Pre-Sleep'); xlabel('Time (s)'); ylabel('Z-score');
for i = 1:length(fields)
    trace = session_data_struct.(fields{i}).str_pre;
    if ~isempty(trace)
        med_trace = median(trace,1);
        plot(time, med_trace, 'LineWidth', 2, 'Color', colors(i,:));
    end
end
xline(0, '--r'); legend(fields, 'Interpreter', 'none'); grid on;

subplot(2,2,4); hold on;
title('Striatum Post-Sleep'); xlabel('Time (s)'); ylabel('Z-score');
for i = 1:length(fields)
    trace = session_data_struct.(fields{i}).str_post;
    if ~isempty(trace)
        med_trace = median(trace,1);
        plot(time, med_trace, 'LineWidth', 2, 'Color', colors(i,:));
    end
end
xline(0, '--r'); legend(fields, 'Interpreter', 'none'); grid on;


%}
