%animalList = [N1, N2, N3].
behavModel = runBehavModel(varargin)

allEntries = dir(pwd);
allDirs = allEntries([allEntries.isdir]);
folderNames = {allDirs(~ismember({allDirs.name}, {'.', '..'})).name};
filteredNames = folderNames(~cellfun('isempty', regexp(folderNames, '^N\d+$')));
animalList = strjoin(filteredNames, ', ');

p = inputParser;
addParameter(p, 'basepath', pwd, @(x) ischar(x) || isstring(x));  
addParameter(p, 'animalList', animalList, @(x) iscell(x) && all(cellfun(@isfolder, x))); 

parse(p,varargin{:});
basepath = p.Results.basepath;
animalList = p.Results.animalList;


animalArray = strsplit(animalList, ', ');
for i = 1:length(animalArray)
    animalModelData = struct();
    animal = animalArray{i};
    disp(['Processing: ', animal]);
    cd(fullfile(basepath, animal));  
    behavTrials = getAnimalBehav();

    %Get variables for the model
    chosenPort = behavTrials.port;
    patchTrial = behavTrials.patch_trials;
    sessionTrial = 1:length(behavTrials.patch_trials);
    chosenPortProb = behavTrials.ports_probability(behavTrials.port);
    previousReward = [2, behavTrials.reward_outcome(1:end-1)];
    twoPreviousReward = [2, 2, behavTrials.reward_outcome(1:end-2)];
    consecutiveNonRew = zeros(size(behavTrials.reward_outcome))
    for i = 2:length(behavTrials.reward_outcome)
        if consecutiveNonRew(i) == 0
            consecutiveNonRew(i) = consecutiveNonRew(i-1) + 1;
        else
            consecutiveNonRew(i) = 0; % Reset to 0 when a 1 is encountered
        end
    end

    tau = 2.1; % Decay constant.
    num_trials = length(behavTrials.reward_outcome);
    n = 5; % Number of previous trials to include

    rewardTrace = zeros(size(behavTrials.reward_outcome));
    
    for t = 2:num_trials
        start_idx = max(1, t-n); % Start of the window
        past_outcomes = behavTrials.reward_outcome(start_idx:t-1); % ourcomes in the window
        time_diff = (t-1:-1:start_idx); % Time difference for weights
        weights = exp(-time_diff / tau);
        rewardTrace(t) = sum(weights .* past_outcomes);
    end


    patchTotalTrials = 


end
    