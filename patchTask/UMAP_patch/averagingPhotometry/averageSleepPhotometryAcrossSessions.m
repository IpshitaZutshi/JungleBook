
%% z score across sleep sessions
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

sessions = cat(2, sessions_N11, sessions_N17);

% Initialize arrays to accumulate data
% note - this is all ripple data. average photometry around ripples
all_pre_ripples_hpc = {};
all_post_ripples_hpc = {};

all_pre_ripples_striatum = {};
all_post_ripples_striatum = {};


for s = 1:length(sessions)
    sessionpath = [direc filesep sessions{s}];
    cd(sessionpath);
    disp(['Collecting data from ' sessions{s}]);

    %fprintf('Running session: %s\n', sessions(s).name);
    session_data = dualSleepPhotometry('showfig', false); 

    % Accumulate HPC
    if isfield(session_data, 'hpc')
        if isfield(session_data.hpc, 'pre_sleep')
            all_pre_ripples_hpc{end+1} = session_data.hpc.pre_sleep;
        end
        if isfield(session_data.hpc, 'post_sleep')
            all_post_ripples_hpc{end+1} = session_data.hpc.post_sleep;
        end
    end
    
    % Accumulate STRIATUM
    if isfield(session_data, 'striatum')
        if isfield(session_data.striatum, 'pre_sleep')
            all_pre_ripples_striatum{end+1} = session_data.striatum.pre_sleep;
        end
        if isfield(session_data.striatum, 'post_sleep')
            all_post_ripples_striatum{end+1} = session_data.striatum.post_sleep;
        end
    end
end

%% average and plot across sessions

% Concatenate and average for HPC
pre_ripples_hpc_all = cat(1, all_pre_ripples_hpc{:});
post_ripples_hpc_all = cat(1, all_post_ripples_hpc{:});

        % rewarded
avg_pre_ripples_hpc = mean(pre_ripples_hpc_all, 1);
med_pre_ripples_hpc = median(pre_ripples_hpc_all, 1); % median at each timepoint

            % confidence intervals
N = height(pre_ripples_hpc_all);                          % # of experiments
a_SEM = std(pre_ripples_hpc_all, 1)/sqrt(N);         % Standard Error Of The Mean 
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
a_CI95 = bsxfun(@times, a_SEM, CI95(:));  % Calculate 95% Confidence Intervals 

        % nonrewarded
avg_post_ripples_hpc = mean(post_ripples_hpc_all, 1);
med_post_ripples_hpc = median(post_ripples_hpc_all, 1);

            % confidence intervals
N = height(post_ripples_hpc_all);                          % Number of eExperimentsn In Data Set
b_SEM = std(post_ripples_hpc_all, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
b_CI95 = bsxfun(@times, b_SEM, CI95(:));  % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of exw


% PLOT

% Create time vector
sampling_rate = 130;
window = 5;
samples = window * sampling_rate;
time = linspace(-window, window, 2*samples+1);

plot_color = [0.588235294117647   0.800000000000000   0.345098039215686]; % green
all_axes = gobjects(2,1);  % Preallocate axes handles
y_min = inf;
y_max = -inf;

figure('color','white');
subplot(2,1,1)
hold on
plot(time, med_pre_ripples_hpc, 'color', plot_color, 'LineWidth', 2);
fill([time,fliplr(time)], [(a_CI95(1,:)+med_pre_ripples_hpc),fliplr((a_CI95(2,:)+med_pre_ripples_hpc))], plot_color, 'EdgeColor','none', 'FaceAlpha',0.25)
ax = gca;
ax.FontSize = 15;
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title(['Average Z-score Around Ripples'], [num2str(size(pre_ripples_hpc_all, 1)), ' ripples']);
grid on;
% Store axes handle
all_axes(1) = gca;  % for subplot 1
% Update Y-axis limits
current_ylim = ylim(gca);
y_min = min(y_min, current_ylim(1));
y_max = max(y_max, current_ylim(2));
hold off

subplot(2,1,2)
hold on
plot(time, med_post_ripples_hpc, 'color', plot_color, 'LineWidth', 2);
fill([time,fliplr(time)], [(b_CI95(1,:)+med_post_ripples_hpc),fliplr((b_CI95(2,:)+med_post_ripples_hpc))], plot_color, 'EdgeColor','none', 'FaceAlpha',0.25)
ax = gca;
ax.FontSize = 15;
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title(['Average Z-score Around Ripples'], [num2str(size(post_ripples_hpc_all, 1)), ' ripples']);
grid on;
all_axes(2) = gca; 
% Update Y-axis limits
current_ylim = ylim(gca);
y_min = min(y_min, current_ylim(1));
y_max = max(y_max, current_ylim(2));
hold off

% Set all to same y-limits
for i = 1:2
    ylim(all_axes(i), [y_min y_max]);
end

% Link all axes
linkaxes(all_axes, 'xy');



% Concatenate and average for STRIATUM
pre_ripples_striatum_all = cat(1, all_pre_ripples_striatum{:});
post_ripples_striatum_all = cat(1, all_post_ripples_striatum{:});

        % rewarded
avg_pre_ripples_striatum = mean(pre_ripples_striatum_all, 1);
med_pre_ripples_striatum = median(pre_ripples_striatum_all, 1);

            % confidence intervals
N = height(pre_ripples_striatum_all);                          % # of experiments
a_SEM = std(pre_ripples_striatum_all, 1)/sqrt(N);         % Standard Error Of The Mean 
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
a_CI95 = bsxfun(@times, a_SEM, CI95(:));  % Calculate 95% Confidence Intervals 

        % nonrewarded
avg_post_ripples_striatum = mean(post_ripples_striatum_all, 1);
med_post_ripples_striatum = median(post_ripples_striatum_all, 1);

            % confidence intervals
N = height(post_ripples_striatum_all);                          % Number of eExperimentsn In Data Set
b_SEM = std(post_ripples_striatum_all, 1)/sqrt(N);         % Compute rStandard Error Of The Meane Of All Experiments At Each Value Of wxa
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
b_CI95 = bsxfun(@times, b_SEM, CI95(:));  % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of exw


% PLOT

% Create time vector
sampling_rate = 130;
window = 5;
samples = window * sampling_rate;
time = linspace(-window, window, 2*samples+1);

plot_color = [0.909803921568627   0.290196078431373   0.454901960784314]; % pink
all_axes = gobjects(2,1);  % Preallocate axes handles
y_min = inf;
y_max = -inf;

figure('color','white');
subplot(2,1,1)
hold on
plot(time, med_pre_ripples_striatum, 'color', plot_color, 'LineWidth', 2);
fill([time,fliplr(time)], [(a_CI95(1,:)+med_pre_ripples_striatum),fliplr((a_CI95(2,:)+med_pre_ripples_striatum))], plot_color, 'EdgeColor','none', 'FaceAlpha',0.25)
ax = gca;
ax.FontSize = 15;
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title(['Average Z-score Around Ripples'], [num2str(size(pre_ripples_striatum_all, 1)), ' ripples']);
grid on;
% Store axes handle
all_axes(1) = gca;  % for subplot 1
% Update Y-axis limits
current_ylim = ylim(gca);
y_min = min(y_min, current_ylim(1));
y_max = max(y_max, current_ylim(2));
hold off

subplot(2,1,2)
hold on
plot(time, med_post_ripples_striatum, 'color', plot_color, 'LineWidth', 2);
fill([time,fliplr(time)], [(b_CI95(1,:)+med_post_ripples_striatum),fliplr((b_CI95(2,:)+med_post_ripples_striatum))], plot_color, 'EdgeColor','none', 'FaceAlpha',0.25)
ax = gca;
ax.FontSize = 15;
xline(0, '--r', 'LineWidth', 1)
xlabel('time (s)');
ylabel('avg z-score');
title(['Average Z-score Around Ripples'], [num2str(size(post_ripples_striatum_all, 1)), ' ripples']);
grid on;
all_axes(2) = gca; 
% Update Y-axis limits
current_ylim = ylim(gca);
y_min = min(y_min, current_ylim(1));
y_max = max(y_max, current_ylim(2));
hold off

% Set all to same y-limits
for i = 1:2
    ylim(all_axes(i), [y_min y_max]);
end

% Link all axes
linkaxes(all_axes, 'xy');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% aggregated pre and post - HPC

save_path = '\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\patchTask\figures';
saveLoc = strcat(save_path,'\sleep_phot');
% if ~isfolder('sleep_phot')
%     mkdir('sleep_phot')
% end 

figure('color','white');
hold on;

% Colors for pre and post
color_pre = [0.588235294117647   0.800000000000000   0.345098039215686]; % green
color_post = [0.301960784313725   0.458823529411765   0.192156862745098]; % darker green

% Plot pre sleep median and CI
h_pre = plot(time, med_pre_ripples_hpc, 'Color', color_pre, 'LineWidth', 2);
fill([time, fliplr(time)], [(a_CI95(1,:) + med_pre_ripples_hpc), fliplr(a_CI95(2,:) + med_pre_ripples_hpc)], ...
    color_pre, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Plot post sleep median and CI
h_post = plot(time, med_post_ripples_hpc, 'Color', color_post, 'LineWidth', 2);
fill([time, fliplr(time)], [(b_CI95(1,:) + med_post_ripples_hpc), fliplr(b_CI95(2,:) + med_post_ripples_hpc)], ...
    color_post, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');

xline(0, '--r', 'LineWidth', 1);
xlabel('time (s)');
ylabel('avg z-score');
%title('Combined HPC: Pre and Post Sleep Around Ripples');
legend([h_pre, h_post], {'Pre Sleep', 'Post Sleep'}, 'Location', 'Best');
grid on;
ax = gca;
ax.FontSize = 15;
set(gcf,'position',[500,200,1120,840])
hold off;
saveas(gcf,[saveLoc,filesep ,'HPC_separated','.png'],'png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aggregated pre and post - striatum

figure('color','white');
hold on;

% Colors for pre and post (pink shades)
color_pre = [0.909803921568627   0.290196078431373   0.454901960784314]; % bright pink
color_post = [0.647058823529412   0.196078431372549   0.313725490196078]; % darker pink

% Plot pre sleep median and CI
h_pre = plot(time, med_pre_ripples_striatum, 'Color', color_pre, 'LineWidth', 2);
fill([time, fliplr(time)], [(a_CI95(1,:) + med_pre_ripples_striatum), fliplr(a_CI95(2,:) + med_pre_ripples_striatum)], ...
    color_pre, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Plot post sleep median and CI
h_post = plot(time, med_post_ripples_striatum, 'Color', color_post, 'LineWidth', 2);
fill([time, fliplr(time)], [(b_CI95(1,:) + med_post_ripples_striatum), fliplr(b_CI95(2,:) + med_post_ripples_striatum)], ...
    color_post, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');

xline(0, '--r', 'LineWidth', 1);
xlabel('time (s)');
ylabel('avg z-score');
%title('Combined Striatum: Pre and Post Sleep Around Ripples');
legend([h_pre, h_post], {'Pre Sleep', 'Post Sleep'}, 'Location', 'Best');
grid on;
ax = gca;
ax.FontSize = 15;
set(gcf,'position',[500,200,1120,840])
hold off;
saveas(gcf,[saveLoc,filesep ,'striatum_separated','.png'],'png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMBINED HPC TRACE (pre + post pooled together)

% Combine pre and post data
all_ripples_hpc = cat(1, pre_ripples_hpc_all, post_ripples_hpc_all);

% Median and CI
med_all_hpc = median(all_ripples_hpc, 1);
N = size(all_ripples_hpc, 1);
SEM_all = std(all_ripples_hpc, 1) / sqrt(N);
CI95 = tinv([0.025 0.975], N-1);
CI_all = bsxfun(@times, SEM_all, CI95(:));

% Plot
figure('color','white');
hold on;
plot_color = [0.588235294117647   0.800000000000000   0.345098039215686]; % green
h = plot(time, med_all_hpc, 'Color', plot_color, 'LineWidth', 2);
fill([time, fliplr(time)], ...
    [(CI_all(1,:) + med_all_hpc), fliplr(CI_all(2,:) + med_all_hpc)], ...
    plot_color, 'FaceAlpha', 0.25, 'EdgeColor', 'none');

xline(0, '--r', 'LineWidth', 1);
xlabel('time (s)');
ylabel('avg z-score');
title(['Combined HPC Trace (Pre + Post), ' num2str(N) ' ripples']);
grid on;
ax = gca;
ax.FontSize = 15;
set(gcf,'position',[500,200,1120,840])
hold off;
saveas(gcf,[saveLoc,filesep ,'HPC_aggreg','.png'],'png');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMBINED STRIATUM TRACE (pre + post pooled together)

% Combine pre and post data
all_ripples_striatum = cat(1, pre_ripples_striatum_all, post_ripples_striatum_all);

% Median and CI
med_all_striatum = median(all_ripples_striatum, 1);
N = size(all_ripples_striatum, 1);
SEM_all = std(all_ripples_striatum, 1) / sqrt(N);
CI95 = tinv([0.025 0.975], N-1);
CI_all = bsxfun(@times, SEM_all, CI95(:));

% Plot
figure('color','white');
hold on;
plot_color = [0.909803921568627   0.290196078431373   0.454901960784314]; % pink
h = plot(time, med_all_striatum, 'Color', plot_color, 'LineWidth', 2);
fill([time, fliplr(time)], ...
    [(CI_all(1,:) + med_all_striatum), fliplr(CI_all(2,:) + med_all_striatum)], ...
    plot_color, 'FaceAlpha', 0.25, 'EdgeColor', 'none');

xline(0, '--r', 'LineWidth', 1);
xlabel('time (s)');
ylabel('avg z-score');
title(['Combined Striatum Trace (Pre + Post), ' num2str(N) ' ripples']);
grid on;
ax = gca;
ax.FontSize = 15;
set(gcf,'position',[500,200,1120,840])
hold off;
saveas(gcf,[saveLoc,filesep ,'striatum_aggreg','.png'],'png');







