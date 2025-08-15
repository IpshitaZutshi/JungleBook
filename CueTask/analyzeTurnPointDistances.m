function analyzeTurnPointDistances(tracking, behavTrials)
% Analyze turn point distances (entry to max turn) across trials, comparing
% baseline vs stim on the preferred side (where mouse tends to go during stim)

% Settings
yThresh = 90;  % Y threshold to define T-junction entry

% Get all trial info
nTrials = length(behavTrials.choice);
turnData = struct('trial', {}, 'stim', {}, 'choice', {}, ...
                  'euclideanDist', {}, 'pathDist', {}, ...
                  'entryIdx', {}, 'turnIdx', {});

for i = 1:nTrials
    tStart = behavTrials.timestamps(i,1);
    tEnd   = behavTrials.choiceTS(i);

    idx = tracking.timestamps >= tStart & tracking.timestamps <= tEnd;
    x = tracking.position.x(idx);
    y = tracking.position.y(idx);
    
    if length(x) < 5, continue; end

    % Entry to junction
    entryIdx = find(y > yThresh, 1, 'first');
    if isempty(entryIdx) || entryIdx > length(x)-3, continue; end

    % Angular velocity
    dx = gradient(x);
    dy = gradient(y);
    angles = atan2(dy, dx);
    dAngles = abs(gradient(unwrap(angles)));

    searchIdx = entryIdx:length(x);
    [~, maxTurnRel] = max(dAngles(searchIdx));
    turnIdx = searchIdx(1) + maxTurnRel - 1;

    % Distances
    eucDist = sqrt((x(turnIdx) - x(entryIdx))^2 + (y(turnIdx) - y(entryIdx))^2);
    pathDist = sum(sqrt(diff(x(entryIdx:turnIdx)).^2 + diff(y(entryIdx:turnIdx)).^2));

    % Store
    turnData(end+1) = struct( ...
        'trial', i, ...
        'stim', behavTrials.stim(i), ...
        'choice', behavTrials.choice(i), ...
        'euclideanDist', eucDist, ...
        'pathDist', pathDist, ...
        'entryIdx', entryIdx, ...
        'turnIdx', turnIdx);
end

% Convert to table for easy filtering
T = struct2table(turnData);

% Determine preferred side during stim
stimTrials = T.stim == 1;
leftCount = sum(T.choice(stimTrials) == 1);
rightCount = sum(T.choice(stimTrials) == 0);
if leftCount > rightCount
    preferred = 1;  % left
    disp('Preferred side during stim: LEFT');
else
    preferred = 0;  % right
    disp('Preferred side during stim: RIGHT');
end

% Filter only preferred side
preferredTrials = T.choice == preferred;

% Get path distances for stim and baseline
dStim = T.pathDist(T.stim == 1 & preferredTrials);
dBase = T.pathDist(T.stim == 0 & preferredTrials);

% Remove NaNs
dStim = dStim(~isnan(dStim));
dBase = dBase(~isnan(dBase));

% Plot
figure; hold on;
boxplot([dBase; dStim], [zeros(size(dBase)); ones(size(dStim))], ...
    'Labels', {'Baseline', 'Stim'});
ylabel('Path Distance to Turn Point');
title('Turn Distance (Preferred Side Only)');
[h, p] = ttest2(dBase, dStim);
text(1.5, max([dBase; dStim]) * 1.05, sprintf('p = %.3f', p), 'HorizontalAlignment', 'center');

end
