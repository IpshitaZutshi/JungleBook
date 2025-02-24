function [behavTrials] = getPatchBehaviorL(varargin)

p = inputParser;
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'plotfig',true,@islogical)
addParameter(p,'forceRun',true,@islogical)
addParameter(p,'updatedIntan',true,@islogical)

parse(p,varargin{:});
saveMat = p.Results.saveMat;
plotfig = p.Results.plotfig;
forceRun = p.Results.forceRun;
updatedIntan = p.Results.updatedIntan;

basepath = pwd;

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

% THIS METHOD SKIPS FIRST TRIAL

%Organize timestamps
ports = 3:9; % Port indices in the digitalIn file
all_On = [];
all_Off = [];

for port = ports
    % Create the On and Off matrices for each port
    port_On = [digitalIn.timestampsOn{1, port}, port - 2 * ones(length(digitalIn.timestampsOn{1, port}), 1)];
    port_Off = [digitalIn.timestampsOff{1, port}, port - 2 * ones(length(digitalIn.timestampsOff{1, port}), 1)];
    all_On = [all_On; port_On];
    all_Off = [all_Off; port_Off];
end

all_On = sortrows(all_On);
all_Off = sortrows(all_Off);


start_trial = []; %when the aimal leaves previous port
end_trial = []; %when the animal leaves chosen port
chosen_port = []; %chosen licked port
outcome_time = []; %timestamsp where the solenoid should be activated (or not if non rewarded)
start_lick = [];
end_lick = [];
% Trial onset: the animal leaves the port. Trial offset: the animal leaves
% the chosen port.

% Calculating trial onset and chosen port

for i = 1:size(all_Off, 1) - 1
    current_time = all_Off(i, 1);     % Actual time
    current_port = all_Off(i, 2);     % Actual port
    next_time = all_Off(i + 1, 1);    % Next time
    next_port = all_Off(i + 1, 2);    % Next port    

    if current_port ~= next_port
        start_trial = [start_trial; current_time];
        chosen_port = [chosen_port; next_port];   
    end
end

%Calculating trial offset (=onset of next trial)
for i = 1:length(start_trial)
    end_trial = [start_trial(2:end); all_Off(end, 1)];
end

%Calculating lick outcome timestamps
for i = 2:size(all_On, 1)
    current_time_on = all_On(i, 1);     % Actual time
    current_port_on = all_On(i, 2);     % Actual port
    pre_port_on = all_On(i - 1, 2);     % Previous port    

    if current_port_on ~= pre_port_on
        outcome_time = [outcome_time; current_time_on];
    end  
end

for i = 2:size(all_On, 1) 
    current_time_on = all_On(i, 1);     % Actual time
    current_port_on = all_On(i, 2);     % Actual port
    pre_port_on = all_On(i - 1, 2);
    if current_port_on ~= pre_port_on
       start_lick = [start_lick; all_On(i, 1)];
       end_lick = [end_lick; all_Off(i, 1)];
    end
end

%% Initialize
num_trials = length(start_trial);
behavTrials.timestamps = [start_trial,end_trial];
behavTrials.port =  chosen_port;
behavTrials.reward_outcome = zeros(size(start_trial,1),1); % rewarded (1), not rewarded (0)
behavTrials.ports_probability = zeros(size(start_trial,1),7); % probabilities of all ports THIS IS SCREWED UP
behavTrials.patch_number = zeros(size(start_trial,1),1); % number of patch switches SO IS THIS
behavTrials.patch_type = zeros(size(start_trial,1),1); % which patch has higher rewards
behavTrials.patch_trials = zeros(size(start_trial,1),1); % what number trial it is within the patch
behavTrials.stay_switch = zeros(size(start_trial,1),1);% stay patch(0), switch patch(1)
behavTrials.lick_timestamps = [start_lick, end_lick]; % only first lick


%% Get solenoid data
sol = 0;
for i = 1:length(outcome_time)
    time = outcome_time(i);
    prt = chosen_port(i);
    if prt == 1
        sol = 10;
    elseif prt == 2
        sol = 11;
    elseif prt == 3
        sol = 12;
    elseif prt == 4
        sol = 13;
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
for i = 2:length(start_trial)
    current_port = behavTrials.port(i);
    previous_port = behavTrials.port(i-1);
    
    if (ismember(current_port, [1, 2, 3, 0]) && ismember(previous_port, [0, 5, 6, 7])) || ...
       (ismember(current_port, [5, 6, 7, 0]) && ismember(previous_port, [0, 1, 2, 3]))
        behavTrials.stay_switch(i) = 1;
    else
        behavTrials.stay_switch(i) = 0;
    end
end

%% Analyze text file
% this only gets info on the probabilities and patch, because reward and lick data is
% coming from intan. to read all info from the text file, use
% analyze_patch_data.m
fileInfo = dir('PatchBehav*');
    
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
count = 1;    
while ~feof(fid)
    % Read the first line (last trials)
    trial_count = trial_count + 1; % Increment trial count
    lastTrials = fgetl(fid); 
    if startsWith(lastTrials, 'Updated lastTrials: ')
        % Remove the prefix 'Updated lastTrials: ' and keep only the numbers
        lastTrials = strrep(lastTrials, 'Updated lastTrials: ', '');
    end
    % Skip the empty line
    empty_line = fgetl(fid);
    if ~isempty(strtrim(empty_line))
        % If not empty, warn and skip
        disp('Expected an empty line, skipping unexpected content...');
    end
    
    % Read second line (trial data)
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
    behavTrials.patch_trials(trial_count) = str2double(patch_info{1})+1;
    behavTrials.patch_type(trial_count) = str2double(patch_info{2});
    
    if trial_count > 1 && behavTrials.patch_type(trial_count) ~= behavTrials.patch_type(trial_count - 1)
        count = count + 1;
    end
    behavTrials.patch_number(trial_count) = count;
end
fclose(fid);


%% Plot
% lick_matrix = zeros(num_trials,7);
% 
% % fill matrix
% for k = 1:num_trials
%     pt = behavTrials.port(k);
%     if behavTrials.reward_outcome(k) == 1
%         lick_matrix(k, pt) = 2; % rewarded
%     elseif behavTrials.reward_outcome(k) == 0
%         lick_matrix(k, pt) = 1; % not rewarded
%     end
%     % unlicked ports remain 0
% end
% behavTrials.lick_matrix = lick_matrix;
% 
% 
% patch_licks = zeros(3,1);
% patch_licks(1) = nnz(behavTrials.port == 1) + nnz(behavTrials.port == 2) + nnz(behavTrials.port == 3);
% patch_licks(2) = nnz(behavTrials.port == 4);
% patch_licks(3) = nnz(behavTrials.port == 5) + nnz(behavTrials.port == 6) + nnz(behavTrials.port == 7);
% 
% 
% figure
% subplot(1,3,[1,2])
% ylabel('Trial #');
% xlabel('Reward Port');
% heatmap(lick_matrix);
% 
% subplot(1,3,3)
% title('Distribution of Licks across patches');
% x = ["Patch 0" "Middle port" "Patch 1"];
% bar(x, patch_licks);

if saveMat
    C = strsplit(pwd,'\');
    save([basepath filesep C{end} '.TrialBehavior.Behavior.mat'],'behavTrials');
end

disp('done!');



