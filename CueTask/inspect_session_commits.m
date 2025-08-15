function inspect_session_commits

basepath = pwd;
% Minimal per-session analysis:
% (1) Build baseline template from BC trials
% (2) Classify each trial and find first commit-Y
% (3) Make 3 quick figures

%% ---- params (keep simple) ----
NBins      = 120;
dprimeThr  = 1.0;     % separation needed to classify
dotSize    = 20;

clrBC = [0.20 0.60 0.80];  % BC correct
clrBE = [0.30 0.75 0.40];  % BC error
clrST = [0.85 0.40 0.40];  % Stim (all)

%% ---- load ----
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

%% ---- baseline template from BC only ----
L = idxBC(choice(idxBC)==1);
R = idxBC(choice(idxBC)==0);
assert(~isempty(L) && ~isempty(R), 'Need both Left and Right BC trials.');

bandX  = 4;                                        % ±5 px around center; tweak if needed
centerX = estimate_center_x(tracking, behavTrials, [L(:); R(:)]);  % session-specific

% Y grid spans only the center-arm samples (no side ports)
allYmin = inf; allYmax = -inf;
for tr = [L(:); R(:)]'
    [~,y] = xy_unique_center(tracking, behavTrials, tr, centerX, bandX);
    if ~isempty(y)
        allYmin = min(allYmin, min(y));
        allYmax = max(allYmax, max(y));
    end
end
assert(isfinite(allYmin) && isfinite(allYmax) && allYmax>allYmin, ...
    'Could not define Y-range for center arm.');
yG = linspace(allYmin, allYmax, NBins).';

interp_no_extrap = @(x,y,yq) interp1(y, x, yq, 'linear', NaN);

Xl = nan(NBins, numel(L));
Xr = nan(NBins, numel(R));
for i=1:numel(L)
    [x,y]   = xy_unique_center(tracking, behavTrials, L(i), centerX, bandX);
    if numel(y)>=2, Xl(:,i) = interp_no_extrap(x,y,yG); end
end
for i=1:numel(R)
    [x,y]   = xy_unique_center(tracking, behavTrials, R(i), centerX, bandX);
    if numel(y)>=2, Xr(:,i) = interp_no_extrap(x,y,yG); end
end

muL = mean(Xl,2,'omitnan');
muR = mean(Xr,2,'omitnan');
sdL = std(Xl,0,2,'omitnan'); sdL(~isfinite(sdL)|sdL==0) = 1e-6;
sdR = std(Xr,0,2,'omitnan'); sdR(~isfinite(sdR)|sdR==0) = 1e-6;

dprime = abs(muL - muR) ./ sqrt(0.5*(sdL.^2 + sdR.^2));

%% ---- commit detector (scan from low Y upward) ----
function ci = commit_idx(xq)
    % llr > 0 => looks Left; < 0 => looks Right
    llL = -0.5*((xq - muL)./sdL).^2 - log(sdL);
    llR = -0.5*((xq - muR)./sdR).^2 - log(sdR);
    llr = llL - llR;
    ok  = isfinite(llr) & isfinite(dprime) & (dprime >= dprimeThr);
    ci  = find(ok & (llr ~= 0), 1, 'first');   % earliest decisive bin
end

%% ---- per-trial commit points ----
conds = {idxBC, idxBE, idxST};
S = struct('trials',[],'commitY',[],'commitX',[]);
for c = 1:3
    ids = conds{c};
    commY = nan(size(ids));
    commX = nan(size(ids));
    for k = 1:numel(ids)
        tr = ids(k);
        [x,y] = xy_unique_center(tracking, behavTrials, tr, centerX, bandX);
        xq    = interp_no_extrap(x,y,yG);
        ci    = commit_idx(xq);
        if ~isempty(ci) && ci>=1 && ci<=numel(yG)
            commY(k) = yG(ci);
            commX(k) = xq(ci);
        end
    end
    S(c).trials  = ids;
    S(c).commitY = commY;
    S(c).commitX = commX;
end

%% ---- FIG 1: BC template over BC trajectories ----
figure; hold on;
for k = 1:numel(idxBC)
    tr  = idxBC(k);
    idx = InIntervals(tracking.timestamps,[behavTrials.timestamps(tr,1), behavTrials.choiceTS(tr)]);
    plot(tracking.position.x(idx), tracking.position.y(idx), '-', 'Color', [0.7 0.7 0.7 0.25]);
end
plot(muL, yG, 'b', 'LineWidth', 2);
plot(muR, yG, 'r', 'LineWidth', 2);
xlabel('X'); ylabel('Y'); title('Baseline template + BC trajectories');
axis equal; grid on;

%% ---- FIG 2: all trials with commit dots ----
figure; hold on;
colors = {clrBC, clrBE, clrST}; labels = {'BC','BE','Stim'};
for c = 1:3
    ids = conds{c};
    for k=1:numel(ids)
        tr  = ids(k);
        idx = InIntervals(tracking.timestamps,[behavTrials.timestamps(tr,1), behavTrials.choiceTS(tr)]);
        plot(tracking.position.x(idx), tracking.position.y(idx), '-', 'Color', [colors{c} 0.18]);
        if isfinite(S(c).commitY(k))
            plot(S(c).commitX(k), S(c).commitY(k), 'o', 'MarkerSize', 5, ...
                'MarkerFaceColor', colors{c}, 'MarkerEdgeColor','k');
        end
    end
end
xlabel('X'); ylabel('Y'); title('Trajectories with commit points');
axis equal; grid on;

%% ---- FIG 3: commit‑Y scatter by condition ----
figure; hold on;
for c=1:3
    yv = S(c).commitY;
    scatter(c + 0.12*(rand(size(yv))-0.5), yv, dotSize, 'filled', ...
            'MarkerFaceColor', colors{c}, 'MarkerEdgeAlpha',0.9);
    yline(mean(yv,'omitnan'), '--', 'Color', colors{c});
end
xlim([0.5 3.5]); xticks(1:3); xticklabels(labels);
ylabel('Commit Y'); title('Commit Y per condition');
grid on;

% quick printout (useful for your early vs late hypothesis)
fprintf('Mean commit Y  BC: %.2f | BE: %.2f | ST: %.2f\n', ...
    mean(S(1).commitY,'omitnan'), mean(S(2).commitY,'omitnan'), mean(S(3).commitY,'omitnan'));

end

%% ---- helper: sort/unique Y so interp1 is happy ----
function [x,y] = xy_unique(tracking, behavTrials, tr)
    idx = InIntervals(tracking.timestamps,[behavTrials.timestamps(tr,1), behavTrials.choiceTS(tr)]);
    x = tracking.position.x(idx); 
    y = tracking.position.y(idx);
    [y,ord] = sort(y(:)); 
    x = x(ord);
    [y,ui] = unique(y);
    x = x(ui);
end

function cx = estimate_center_x(tracking, behavTrials, trials)
% Robust center-arm X from the lower stem (before branching).
vals = [];
for tr = trials(:)'
    [x,y] = xy_unique(tracking, behavTrials, tr);
    if numel(y) < 5, continue; end
    ylo = quantile(y, 0.05);    % bottom 5–40% of Y = straight stem
    yhi = quantile(y, 0.40);
    keep = (y >= ylo) & (y <= yhi);
    if any(keep), vals(end+1) = median(x(keep)); end %#ok<AGROW>
end
cx = median(vals, 'omitnan');
end

function [x,y] = xy_unique_center(tracking, behavTrials, tr, centerX, bandX)
% Sort/unique by Y, then keep only samples near the stem center (|X-centerX|<=bandX).
[x,y] = xy_unique(tracking, behavTrials, tr);
keep  = abs(x - centerX) <= bandX;
x = x(keep); y = y(keep);
end