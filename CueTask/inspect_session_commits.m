function inspect_session_commits

basepath = pwd;
%% --- Params ---
NBins      = 120;      % Y bins for template
dprimeThr  = 1.0;      % Template separation threshold
dotSize    = 20;
minCov     = 5;        % minimum trials per side per bin for template

clrBC = [0.20 0.60 0.80];  % BC correct
clrBE = [0.30 0.75 0.40];  % BC error
clrST = [0.85 0.40 0.40];  % Stim trials

%% --- Load data ---
noseFile = dir(fullfile(basepath, '*Tracking.Nose.mat'));
behFile  = dir(fullfile(basepath, '*TrialBehavior.Events.mat'));
assert(~isempty(noseFile) && ~isempty(behFile), 'Missing files in %s', basepath);

load(fullfile(noseFile.folder, noseFile.name), 'tracking');
load(fullfile(behFile.folder, behFile.name), 'behavTrials');

stim   = behavTrials.stim(:);
corr   = behavTrials.correct(:);
choice = behavTrials.choice(:);

idxBC = find(stim==0 & corr==1);
idxBE = find(stim==0 & corr==0);
idxST = find(stim==1);

%% --- Build baseline template from BC trials ---
L = idxBC(choice(idxBC)==1);
R = idxBC(choice(idxBC)==0);

% Full Y range from all BC trials
allYmin = []; allYmax = [];
for tr = [L(:); R(:)]'
    idx = InIntervals(tracking.timestamps, [behavTrials.timestamps(tr,1), behavTrials.choiceTS(tr)]);
    y   = tracking.position.y(idx);
    allYmin(end+1) = min(y); %#ok<AGROW>
    allYmax(end+1) = max(y); %#ok<AGROW>
end
yG = linspace(min(allYmin), max(allYmax), NBins);

% Interp without extrapolation
interp_no_extrap = @(x,y,yq) interp1(y, x, yq, 'linear', NaN);

Xl = nan(NBins, numel(L));
Xr = nan(NBins, numel(R));
for i=1:numel(L), [x,y]=xy(tracking,behavTrials,L(i)); Xl(:,i)=interp_no_extrap(x,y,yG); end
for i=1:numel(R), [x,y]=xy(tracking,behavTrials,R(i)); Xr(:,i)=interp_no_extrap(x,y,yG); end

% Coverage per bin
nL = sum(isfinite(Xl), 2);
nR = sum(isfinite(Xr), 2);
goodBin = (nL >= minCov) & (nR >= minCov);

% Means & SDs, masked for low coverage
muL = nan(NBins,1); muR = nan(NBins,1);
sdL = nan(NBins,1); sdR = nan(NBins,1);

muL(goodBin) = movmean(mean(Xl(goodBin,:), 2, 'omitnan'), 5, 'omitnan');
muR(goodBin) = movmean(mean(Xr(goodBin,:), 2, 'omitnan'), 5, 'omitnan');
sdL(goodBin) = max(movmean(std(Xl(goodBin,:),0,2,'omitnan'), 5, 'omitnan'), 1e-6);
sdR(goodBin) = max(movmean(std(Xr(goodBin,:),0,2,'omitnan'), 5, 'omitnan'), 1e-6);

%% --- Commit detection function ---
getCommitY = @(xq) ...
    find(abs(muL - muR) ./ sqrt(0.5*(sdL.^2 + sdR.^2)) >= dprimeThr & ...
         (abs( (-0.5*((xq - muL)./sdL).^2 - log(sdL)) ...
              - (-0.5*((xq - muR)./sdR).^2 - log(sdR)) ) ) > 0, 1, 'first');

%% --- Get commit points ---
S = struct('trials',[],'commitY',[],'commitX',[]);
conds = {idxBC, idxBE, idxST};
for c = 1:3
    ids = conds{c};
    commY = nan(size(ids));
    commX = nan(size(ids));
    for k=1:numel(ids)
        tr = ids(k);
        [x,y] = xy(tracking, behavTrials, tr);
        xq = interp_no_extrap(x, y, yG);
        ci = getCommitY(xq);
        if ~isempty(ci) && ci >= 1 && ci <= numel(yG)
            commY(k) = yG(ci);
            commX(k) = xq(ci);
        end
    end
    S(c).trials = ids;
    S(c).commitY = commY;
    S(c).commitX = commX;
end

%% --- FIGURE 1: BC template + BC trajectories ---
figure; hold on;
for k = 1:numel(idxBC)
    tr = idxBC(k);
    idx = InIntervals(tracking.timestamps, ...
        [behavTrials.timestamps(tr,1), behavTrials.choiceTS(tr)]);
    plot(tracking.position.x(idx), tracking.position.y(idx), ...
         '-', 'Color', [0.7 0.7 0.7 0.3]);
end
plot(muL, yG, 'b', 'LineWidth', 2);
plot(muR, yG, 'r', 'LineWidth', 2);
xlabel('X'); ylabel('Y'); title('Baseline template + BC trajectories');
axis equal; grid on;

%% --- FIGURE 2: All trials w/ commit points ---
figure; hold on;
colors = {clrBC, clrBE, clrST};
labels = {'BC','BE','Stim'};
for c = 1:3
    ids = conds{c};
    for k=1:numel(ids)
        tr = ids(k);
        idx = InIntervals(tracking.timestamps, ...
            [behavTrials.timestamps(tr,1), behavTrials.choiceTS(tr)]);
        plot(tracking.position.x(idx), tracking.position.y(idx), ...
             '-', 'Color', [colors{c} 0.2]);
        if isfinite(S(c).commitY(k))
            plot(S(c).commitX(k), S(c).commitY(k), 'o', 'MarkerSize', 5, ...
                'MarkerFaceColor', colors{c}, 'MarkerEdgeColor','k');
        end
    end
end
xlabel('X'); ylabel('Y'); title('Trajectories with commit points');
axis equal; grid on;

% Scatter of commit Y
figure; hold on;
for c=1:3
    yvals = S(c).commitY;
    scatter(c + 0.1*(rand(size(yvals))-0.5), yvals, dotSize, 'filled', ...
        'MarkerFaceColor', colors{c}, 'MarkerEdgeAlpha',0.8);
    yline(mean(yvals,'omitnan'), '--', 'Color', colors{c});
end
xlim([0.5 3.5]); xticks(1:3); xticklabels(labels);
ylabel('Commit Y'); title('Commit Y per condition');
grid on;

end

%% --- Helper ---
function [x,y] = xy(tracking, behavTrials, tr)
    idx = InIntervals(tracking.timestamps,[behavTrials.timestamps(tr,1), behavTrials.choiceTS(tr)]);
    x = tracking.position.x(idx);
    y = tracking.position.y(idx);
    [y,ord] = sort(y(:)); 
    x = x(ord);
    [y, ui] = unique(y); 
    x = x(ui);
end