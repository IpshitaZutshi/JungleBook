% Get patch behavior
% Make sure to set mode to 1 or 0


function [behavTrials] = getPatchBehavior(varargin)

p = inputParser;
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'plotfig',true,@islogical);
addParameter(p,'forceRun',true,@islogical);
addParameter(p,'updatedIntan',true,@islogical);

parse(p,varargin{:});
saveMat = p.Results.saveMat;
plotfig = p.Results.plotfig;
forceRun = p.Results.forceRun;
updatedIntan = p.Results.updatedIntan;

basepath = pwd;
[~, currentFolderName] = fileparts(basepath);

%% Deal with inputs
if ~isempty(dir([basepath filesep '*.TrialBehavior.Events.mat'])) && ~forceRun
    disp('Trial behavior already detected! Loading file.');
    file = dir([basepath filesep '*.TrialBehavior.Events.mat']);
    load(file.name);
    return
end

%% Get digital inputs
if exist('settings.xml')
    delete 'settings.xml'
end

disp('Loading digital In...');
digitalIn = bz_getDigitalIn; 
if isempty(digitalIn)
    patchBehav = [];
    return
end


%% trial data from intan

% add columns to each data with their port number
lick_1 = [digitalIn.timestampsOn{1,3}, ones(length(digitalIn.timestampsOn{1,3}),1)];
lick_2 = [digitalIn.timestampsOn{1,4}, 2*ones(length(digitalIn.timestampsOn{1,4}),1)];
lick_3 = [digitalIn.timestampsOn{1,5}, 3*ones(length(digitalIn.timestampsOn{1,5}),1)];
lick_4 = [digitalIn.timestampsOn{1,6}, 4*ones(length(digitalIn.timestampsOn{1,6}),1)];
lick_5 = [digitalIn.timestampsOn{1,7}, 5*ones(length(digitalIn.timestampsOn{1,7}),1)];
lick_6 = [digitalIn.timestampsOn{1,8}, 6*ones(length(digitalIn.timestampsOn{1,8}),1)];
lick_7 = [digitalIn.timestampsOn{1,9}, 7*ones(length(digitalIn.timestampsOn{1,9}),1)];

% creat matrix from individual arrays
all_starts = [lick_1; lick_2; lick_3; lick_4; lick_5; lick_6; lick_7];
all_starts = sortrows(all_starts);    

trial_licks = []; % array which stores the time of the lick indicating a trial. so it excludes multiple licks at one port, and only counts the first lick (which is the one that would be rewarded)
choice_port = [];
most_recent_lick = 100;

% determines times of licks which determine the start and end of a trial
% (the end of one trial is the same as the start of the next)
for j = 1:length(all_starts)
    start_time = all_starts(j, 1); % gets nonzero number in the row
    current_port = all_starts(j, 2); 
    if most_recent_lick == current_port
        continue
    else
        most_recent_lick = current_port;
        trial_licks = [trial_licks; start_time]; 
        choice_port = [choice_port; current_port];
    end
end


%% Initialize
num_trials = length(trial_licks);
behavTrials.num_trials = num_trials;
behavTrials.timestamps = trial_licks;
behavTrials.port =  choice_port;
behavTrials.reward_outcome = zeros(size(trial_licks,1),1); % rewarded (1), not rewarded (0)
behavTrials.ports_probability = zeros(size(trial_licks,1),7); % probabilities of all ports
behavTrials.patch_number = zeros(size(trial_licks,1),1); % which patch has higher rewards
behavTrials.patch_trials = zeros(size(trial_licks,1),1); % what number trial it is within the patch
behavTrials.stay_switch = zeros(size(trial_licks,1),1);% stay patch(0), switch patch(1)

%% Get solenoid data
sol = 0;
for i = 1:length(trial_licks)
    time = trial_licks(i);
    prt = choice_port(i);
    if prt == 1
        sol = 10;
    elseif prt == 2
        sol = 11;
    elseif prt == 3
        sol = 12;
    elseif prt == 5 % skipping middle port
        sol = 14;
    elseif prt == 6
        sol = 15;
    elseif prt == 7
        sol = 16;
    end
    vals = (digitalIn.timestampsOn{1,sol} > (time - 0.01)) & (digitalIn.timestampsOn{1,sol} < (time + 0.01));

    if ismember(1, vals)
        behavTrials.reward_outcome(i) = 1;
    else
        behavTrials.reward_outcome(i) = 0;
    end
end

%% Detect switch patch trials
%{
for i = 2:length(trial_licks)
    current_port = behavTrials.port(i);
    previous_port = behavTrials.port(i-1);
    
    if (ismember(current_port, [1, 2, 3, 0]) && ismember(previous_port, [0, 5, 6, 7])) || ...
       (ismember(current_port, [5, 6, 7, 0]) && ismember(previous_port, [0, 1, 2, 3]))
        behavTrials.stay_switch(i) = 1;
    else
        behavTrials.stay_switch(i) = 0;
    end
end
%}

%% Analyze text file
% this only gets info on the probabilities and patch, because reward and lick data is
% coming from intan. to read all info from the text file, use
% analyze_patch_data.m
% fileInfo = dir(fullfile(basepath, '*', 'PatchBehav2024*'));
fileInfo = dir('PatchBehav2024*');
if isempty(fileInfo)
        error('No file starting with "PatchBehav" found in the current directory.');
elseif numel(fileInfo) > 1
        error('Multiple files starting with "PatchBehav" found. Please specify the file.');
end
    
filename = fileInfo.name;
    
fid = fopen(filename, 'r');

% Initialize variables
trial_count = 0;
trial_numbers = [];
timestamps_licks = [];
   
mode = 0; % because the serial monitor got screwed up partway through the experiments. 
% if patchBehav file (saved from arduino) only saves two lines for each
% trial, set mode equal to 0. if arduino saves multiple lines (laura's
% change), set mode equal to 1.

if mode == 0
    while ~feof(fid)
        % Read first line (trial data)
        trial_line = fgetl(fid);
        trial_count = trial_count + 1;  
        trial_info = split(trial_line, '_');
        trial_numbers(trial_count) = str2double(trial_info{2});  
        current_rewarded_count = str2double(trial_info{3});     
        
        % get reward probabilities
        prob_str = trial_info{5};
        reward_probs = (split(prob_str, ','));
        behavTrials.ports_probability(trial_count, 1) = str2double(reward_probs{1});
        behavTrials.ports_probability(trial_count, 2) = str2double(reward_probs{2});
        behavTrials.ports_probability(trial_count, 3) = str2double(reward_probs{3});
        behavTrials.ports_probability(trial_count, 4) = 0;
        behavTrials.ports_probability(trial_count, 5) = str2double(reward_probs{4});
        behavTrials.ports_probability(trial_count, 6) = str2double(reward_probs{5});
        behavTrials.ports_probability(trial_count, 7) = str2double(reward_probs{6});
    
        
        % read next line (patch info)
        patch_line = fgetl(fid);
        patch_info = split(patch_line, '_');
        behavTrials.patch_trials(trial_count) = str2double(patch_info{1});
        behavTrials.patch_number(trial_count) = str2double(patch_info{2});
    
    end

elseif mode == 1
    while ~feof(fid)

    % Read the first line (last trials)
    % irrelevant data
    last_trials = fgetl(fid);

    % Skip the empty line
    empty_line = fgetl(fid);
    if ~isempty(strtrim(empty_line))
        % If not empty, warn and skip
        disp('Expected an empty line, skipping unexpected content...');
    end

    % Read first real line (trial data)
    trial_line = fgetl(fid);
    trial_count = trial_count + 1;  
    trial_info = split(trial_line, '_');
    timestamps = trial_info{1};
    timestamps_licks(trial_count) = str2double(timestamps);
    trial_numbers(trial_count) = str2double(trial_info{2});  
    current_rewarded_count = str2double(trial_info{3});     
    
    % get reward probabilities
    prob_str = trial_info{5};
    reward_probs = (split(prob_str, ','));
    behavTrials.ports_probability(trial_count, 1) = str2double(reward_probs{1});
    behavTrials.ports_probability(trial_count, 2) = str2double(reward_probs{2});
    behavTrials.ports_probability(trial_count, 3) = str2double(reward_probs{3});
    behavTrials.ports_probability(trial_count, 4) = 0;
    behavTrials.ports_probability(trial_count, 5) = str2double(reward_probs{4});
    behavTrials.ports_probability(trial_count, 6) = str2double(reward_probs{5});
    behavTrials.ports_probability(trial_count, 7) = str2double(reward_probs{6});

    
    % read next line (patch info)
    patch_line = fgetl(fid);
    patch_info = split(patch_line, '_');
    behavTrials.patch_trials(trial_count) = str2double(patch_info{1});
    behavTrials.patch_number(trial_count) = str2double(patch_info{2});

    end
end

fclose(fid);


%% Plot
lick_matrix = zeros(num_trials,7);

% fill matrix
for k = 1:num_trials
    pt = behavTrials.port(k);
    if behavTrials.reward_outcome(k) == 1
        lick_matrix(k, pt) = 2; % rewarded
    elseif behavTrials.reward_outcome(k) == 0
        lick_matrix(k, pt) = 1; % not rewarded
    end
    % unlicked ports remain 0
end
behavTrials.lick_matrix = lick_matrix;


patch_licks = zeros(3,1);
patch_licks(1) = nnz(behavTrials.port == 1) + nnz(behavTrials.port == 2) + nnz(behavTrials.port == 3);
patch_licks(2) = nnz(behavTrials.port == 4);
patch_licks(3) = nnz(behavTrials.port == 5) + nnz(behavTrials.port == 6) + nnz(behavTrials.port == 7);


figure
subplot(1,3,[1,2])
ylabel('Trial #');
xlabel('Reward Port');
heatmap(lick_matrix);

subplot(1,3,3)
title('Distribution of Licks across patches');
x = ["Patch 0" "Middle port" "Patch 1"];
bar(x, patch_licks);


%% Plot how often they lick in the correct patch

lick_outcome = zeros(3, 1);
% first row of lick outcome = rewarded patch, second is middle port,
% third is unrewarded patch

for i = 1:length(behavTrials.patch_number)
    % determine which patch has higher rewards
    if behavTrials.patch_number(i) == 0
        rewarded_ports = [1, 2, 3];
        nonrewarded_ports = [5, 6, 7];
    else 
        rewarded_ports = [5, 6, 7];
        nonrewarded_ports = [1, 2, 3];
    end

    % determine if mouse was licking high or low patch
    if ismember(behavTrials.port(i), rewarded_ports)
        lick_outcome(1) = lick_outcome(1)+1;
    elseif ismember(behavTrials.port(i), nonrewarded_ports)
        lick_outcome(3) = lick_outcome(3)+1;
    else 
        lick_outcome(2) = lick_outcome(2)+1;
    end
end

figure(2)
title('Distribution of Licks across patches');
z = ["Rewarded Patch" "Middle port" "Nonrewarded Patch"];
bar(z, lick_outcome, 'FaceAlpha',0.5);

%% Plot behavior over time

licked_ports = behavTrials.port; 

figure (3);
colormap(flipud(gray)); 
hold on;
customGreen = [152, 194, 9] / 255;
customRed = [238, 75, 43] / 255;
y_limits = [min(licked_ports), max(licked_ports)]; 
timestamps_minutes = timestamps_licks / 60000;
x_limits = [min(timestamps_minutes), max(timestamps_minutes)];
hold on

% plot patch with higher probability
first_in_patch = 1;
prev = behavTrials.patch_number(1);
for j = 2:length(behavTrials.patch_number)
    if behavTrials.patch_number(j) == behavTrials.patch_number(first_in_patch) && (j ~= length(behavTrials.patch_number))
        prev = j;
        continue
    else
        patch_end = prev;
        if behavTrials.patch_number(patch_end) == 0 % patch 0 is high prob
            y = [1, 1, 3, 3];
        else % patch 1 is high prob
            y = [5, 5, 7, 7];
        end
        x = [timestamps_minutes(first_in_patch), timestamps_minutes(patch_end), timestamps_minutes(patch_end), timestamps_minutes(first_in_patch)];
        first_in_patch = j;
        patch(x, y, [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeAlpha', 0)
    end
end

% plot gray trajectory lines
% fix if the two arrays are different sizes
plot(timestamps_minutes, licked_ports, 'Color', [0.1, 0.1, 0.1]);  

% plot lick points
rewarded_indices = find(behavTrials.reward_outcome == 1);
not_rewarded_indices = find(behavTrials.reward_outcome == 0);
scatter(timestamps_minutes(rewarded_indices), licked_ports(rewarded_indices), 36, customGreen, 'filled');
scatter(timestamps_minutes(not_rewarded_indices), licked_ports(not_rewarded_indices), 36, customRed, 'filled');

ylabel('Port');
xlabel('Time');
title(strjoin(string(currentFolderName), ' '));
hold off
%}

%% Save
if saveMat
    C = strsplit(pwd,'\');
    save([basepath filesep C{end} '.TrialBehaviorUnsynced.mat'],'behavTrials');
end

disp('done!');

