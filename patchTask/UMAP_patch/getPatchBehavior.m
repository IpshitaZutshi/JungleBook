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
    
    if prt == 4
        behavTrials.reward_outcome(i) = 0;
    else
        vals = (digitalIn.timestampsOn{1,sol} > (time - 0.01)) & (digitalIn.timestampsOn{1,sol} < (time + 0.01));

        if ismember(1, vals)
            behavTrials.reward_outcome(i) = 1;
        else
            behavTrials.reward_outcome(i) = 0;
        end
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
fileInfo = dir('PatchBehav20*');
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
   
mode = 1; % because the serial monitor got screwed up partway through the experiments. 
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
if plotfig
    h2 = figure;
    set(gcf,'Renderer','painters')
    set(gcf,'Color','w')

    if ~exist('Behavior','dir')
        mkdir('Behavior\')
    end

    % Plot licks - rewarded/unrewarded based on changing patch probability
    subplot(2,4,1:4)
    box off
    hold on
    col = [0/255 151/255 150/255; 128/255 0/255 128/255];
    for ii= 1:length(behavTrials.timestamps)
       % current high patch
       patch1 = mean(behavTrials.ports_probability(ii,1:3));
       patch2 = mean(behavTrials.ports_probability(ii,5:7));
       if patch1>patch2
           highPatch = 1;
       else
           highPatch = 2;
       end
       if behavTrials.reward_outcome(ii) == 1
            scatter(behavTrials.timestamps(ii),behavTrials.port(ii),45,col(highPatch,:),"filled") 
       else
            scatter(behavTrials.timestamps(ii),behavTrials.port(ii),45,col(highPatch,:)) 
       end
    end
    xlabel('Time')
    ylabel('Port')
    title(strcat('Total trials=',num2str(behavTrials.num_trials),' Rewarded =',num2str(sum(behavTrials.reward_outcome))));
    
    %% Running percentage of licks in a patch
    % Define the number of trials for the running window
    windowSize = 30;
    
    % Identify high-probability patch at each trial
    for ii = 1:behavTrials.num_trials
        highProbPatch(ii) = mean(behavTrials.ports_probability(ii,1:3)) > mean(behavTrials.ports_probability(ii,5:7));
    end
    
    % Determine whether each lick was in the high-probability patch
    licksInHighProbPatch = ismember(behavTrials.port, 1:3) & highProbPatch' | ...
                            ismember(behavTrials.port, 5:7) & ~highProbPatch';
    
    % Compute running percentage over a sliding window
    runningPercentage = zeros(size(licksInHighProbPatch));
    
    for t = 1:length(licksInHighProbPatch)
        if t >= windowSize
            windowLicks = licksInHighProbPatch(t-windowSize+1:t);
        else
            windowLicks = licksInHighProbPatch(1:t);
        end
        runningPercentage(t) = mean(windowLicks) * 100;
    end
    
    % Plot the running percentage
    subplot(2,4,5:8)
    plot(behavTrials.timestamps, runningPercentage, 'b', 'LineWidth', 1);
    hold on
    plot(behavTrials.timestamps,highProbPatch*100,'Color','r','LineWidth',1)
    xlabel('Trial Number');
    ylabel('Running %');

    saveas(h2,'Behavior\Behavior.png');
    saveas(h2,'Behavior\Behavior.fig');

end

%% Save
if saveMat
    C = strsplit(pwd,'\');
    save([basepath filesep C{end} '.TrialBehaviorUnsynced.mat'],'behavTrials');
end

disp('done!');

