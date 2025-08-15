function [y_diverge, stats] = findChoiceDivergenceY(tracking, behavTrials,varargin)
% findChoiceDivergenceY
% Detects the earliest Y-position along the T-maze stem where left vs right
% trajectories occupy significantly different X values.
%
% INPUTS
%   tracking: struct with fields .timestamps, .position.x, .position.y
%   behavTrials: table with fields:
%       .timestamps(:,1)  (trial start)
%       .choiceTS(:)      (choice time)
%       .choice(:)        (1=Left, 0=Right)
%       .stim(:)          (0=baseline, 1=stim)  [row or col]
%       .correct(:)       (1/0)  [only used if stim==0]
%   stim: 0 -> compare Correct vs Error within baseline (your original logic uses this),
%         1 -> compare Stim vs Baseline; BUT for divergence we specifically
%              compare Left vs Right choices within the selected condition(s).
%              So we follow the same filtering you use to define groups, then
%              pool *Left* trials vs *Right* trials from those groups.
%
% Name-Value options:
%   'NBins'          (default 100)  : number of Y bins
%   'Alpha'          (default 0.05) : FDR level
%   'MinRun'         (default 3)    : consecutive significant bins required
%   'Test'           (default 'ranksum') : 'ranksum' or 'ttest2'
%   'ShowPlot'       (default true) : overlay y_diverge on current figure
%
% OUTPUTS
%   y_diverge: scalar Y at which divergence first becomes significant (NaN if none)
%   stats: struct with fields:
%       .yGrid, .p_raw, .p_fdr, .dCohen, .meanLeft, .meanRight, .nLeft, .nRight

% ------------------- Parameters -------------------
p = inputParser;
addParameter(p,'NBins',100,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'Alpha',0.05,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'MinRun',3,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'MinTrials',3,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'Test','ranksum',@(s)ischar(s)||isstring(s));
addParameter(p,'ShowPlot',true,@islogical);
addParameter(p,'TrialFilter', []); % <-- n
parse(p,varargin{:});

NBins   = p.Results.NBins;
Alpha   = p.Results.Alpha;
MinRun  = p.Results.MinRun;
testStr = lower(string(p.Results.Test));
showPlot= p.Results.ShowPlot;
TrialFilter = p.Results.TrialFilter;
MinTrials = p.Results.MinTrials;

% ------------------- Trial selection -------------------

L_idx = find(TrialFilter & behavTrials.choice==1);
R_idx = find(TrialFilter & behavTrials.choice==0);

% Early exits
if numel(L_idx) < 3 || numel(R_idx) < 3
    warning('Not enough trials: Left=%d, Right=%d. Need >=3 each.', numel(L_idx), numel(R_idx));
end

% ------------------- Extract per-trial X(Y) -------------------
% First pass: compute per-trial Y ranges to build a common overlap grid
trialYmin = [];
trialYmax = [];

% Helper to pull (x,y) for a trial segment
    function [x,y] = getXY(tr)
        idx = InIntervals(tracking.timestamps, [behavTrials.timestamps(tr,1), behavTrials.choiceTS(tr)]);
        x = tracking.position.x(idx);
        y = tracking.position.y(idx);
        % Require at least 2 samples
        if numel(x) < 2
            x = []; y = [];
            return;
        end
        % Sort by Y to use Y as independent variable; drop duplicate Y
        [y, ord] = sort(y(:));
        x = x(ord);
        [y, uniqIdx] = unique(y); %#ok<ASGLU>
        x = x(uniqIdx);
    end

% Collect trial lists
allTrials = [L_idx(:); R_idx(:)];
groupFlag = [ones(numel(L_idx),1); zeros(numel(R_idx),1)]; % 1=Left, 0=Right

% Pre-scan for overlap
valid = false(size(allTrials));
for k = 1:numel(allTrials)
    tr = allTrials(k);
    [x,y] = getXY(tr);
    if isempty(x), continue; end
    % Require Y to span some range
    if numel(y) < 5, continue; end
    trialYmin(end+1,1) = min(y); %#ok<AGROW>
    trialYmax(end+1,1) = max(y); %#ok<AGROW>
    valid(k) = true;
end

if ~any(valid)
    warning('No valid trials after preprocessing.');
    y_diverge = NaN;
    stats = struct();
    return;
end

% Define common Y overlap region to avoid extrapolation bias
yLow  = max(trialYmin);
yHigh = min(trialYmax);
if ~(isfinite(yLow) && isfinite(yHigh) && yHigh > yLow)
    warning('Insufficient overlap in Y across trials (yLow=%.3f, yHigh=%.3f).', yLow, yHigh);
    y_diverge = NaN;
    stats = struct();
    return;
end
yGrid = linspace(yLow, yHigh, NBins);

% Interpolate X at each yGrid for each trial, separated by group
X_left  = nan(NBins, sum(groupFlag==1 & valid));
X_right = nan(NBins, sum(groupFlag==0 & valid));

iL = 1; iR = 1;
for k = 1:numel(allTrials)
    if ~valid(k), continue; end
    tr = allTrials(k);
    [x,y] = getXY(tr);
    if isempty(x), continue; end
    % Interpolate X as function of Y along the grid
    xq = interp1(y, x, yGrid, 'linear', 'extrap');  % we ensured overlap, so little to no extrap
    if groupFlag(k)==1
        X_left(:,iL)  = xq; iL = iL+1;
    else
        X_right(:,iR) = xq; iR = iR+1;
    end
end

% Trim columns that remained NaN (in case)
X_left  = X_left(:, all(~isnan(X_left),1));
X_right = X_right(:,all(~isnan(X_right),1));

nL = size(X_left,2);
nR = size(X_right,2);

if nL < MinTrials || nR < MinTrials
    warning('Skipping session: not enough trials (Left=%d, Right=%d)', nL, nR);
    y_diverge = NaN;
    stats = struct('nLeft',nL,'nRight',nR);
    return;
end

% ------------------- Bin-wise statistics -------------------
p_raw     = nan(NBins,1);
dCohen    = nan(NBins,1);
meanLeft  = nan(NBins,1);
meanRight = nan(NBins,1);

for b = 1:NBins
    Lb = X_left(b, :);
    Rb = X_right(b,:);

    % Basic summaries
    meanLeft(b)  = mean(Lb, 'omitnan');
    meanRight(b) = mean(Rb, 'omitnan');

    switch testStr
        case "ranksum"
            try
                p_raw(b) = ranksum(Lb, Rb); % nonparametric MW-U
            catch
                % Fallback to ttest2 if Statistics Toolbox not available
                [~, p_raw(b)] = ttest2(Lb, Rb, 'Vartype','unequal');
            end
        otherwise
            [~, p_raw(b)] = ttest2(Lb, Rb, 'Vartype','unequal');
    end

    % Cohen's d (pooled SD)
    mL = mean(Lb,'omitnan'); mR = mean(Rb,'omitnan');
    sL = std(Lb,[],'omitnan'); sR = std(Rb,[],'omitnan');
    n_l = sum(~isnan(Lb)); n_r = sum(~isnan(Rb));
    sp = sqrt(((n_l-1)*sL^2 + (n_r-1)*sR^2) / max(n_l+n_r-2,1));
    if sp>0
        dCohen(b) = (mL - mR)/sp;
    end
end

% Benjamini-Hochberg FDR
p_fdr = bhFDR(p_raw);

% Find first Y with ≥ MinRun consecutive significant bins
sig = p_fdr < Alpha;
y_diverge = NaN;
if any(sig)
    % Find runs of true
    dSig = diff([false; sig; false]);
    runStarts = find(dSig==1);
    runEnds   = find(dSig==-1)-1;
    runLens   = runEnds - runStarts + 1;
    ok = find(runLens >= MinRun, 1, 'first');
    if ~isempty(ok)
        firstRunStart = runStarts(ok);
        y_diverge = yGrid(firstRunStart);
    end
end

% ------------------- Plot overlay -------------------
if showPlot && isfinite(y_diverge)
    yl = y_diverge;
    xl = [min([X_left(:); X_right(:)], [], 'omitnan'), max([X_left(:); X_right(:)], [], 'omitnan')];
    % Draw a horizontal guideline at the divergence Y across current axes
    hold on;
    yline(yl, '--', sprintf('Choice diverges @ Y=%.2f', yl), 'LabelVerticalAlignment','bottom');
end

% ------------------- Output stats -------------------
stats = struct('yGrid',yGrid,'p_raw',p_raw,'p_fdr',p_fdr, ...
               'dCohen',dCohen,'meanLeft',meanLeft,'meanRight',meanRight, ...
               'nLeft',nL,'nRight',nR,'Alpha',Alpha,'MinRun',MinRun);

end

% ---------- Helper: Benjamini-Hochberg FDR ----------
function p_adj = bhFDR(p)
    p = p(:);
    [ps, idx] = sort(p);
    m = numel(p);
    % Avoid zeros/NaNs issues
    ps(isnan(ps)) = 1;
    q = (m ./ (1:m)') .* ps;
    q = min(accumarray((m:-1:1)', q(end:-1:1), [], @min), [], 'all'); %#ok<NASGU> % not used
    % Classic BH
    bh = ps .* m ./ (1:m)';
    bh = cummin(bh(end:-1:1));
    bh = bh(end:-1:1);
    p_adj = nan(size(p));
    p_adj(idx) = min(bh, 1);
end
