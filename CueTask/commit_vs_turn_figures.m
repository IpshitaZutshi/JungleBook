function OUT = commit_vs_turn_figures(basepath)
% See header in previous version.

if nargin==0 || isempty(basepath), basepath = pwd; end

%% ---------- Params ----------
NBins       = 120;
dprimeThr   = 1.0;
llrThr      = 0.75;
bandX       = 5;
dotSize     = 10;

% Geometric turn detector thresholds (applied to FULL trajectory)
offsetThr       = 2.0;     % |X - centerX| in pixels
slopeThr_dydy   = 0.25;    % // NEW: |dx/dy| (pixels of X per pixel of Y)
minRun          = 8;       % persistence in bins
yMinTurn = 70; 

% Plot colors
clrBC = [0.20 0.60 0.80]; clrBE = [0.30 0.75 0.40]; clrST = [0.85 0.40 0.40];
colors = {clrBC, clrBE, clrST}; labels = {'BC','BE','Stim'};

% Gates
yMinCommit   = 70;   % // NEW: only mark commits at or above this Y
turnYMaxFrac = 0.92; % // NEW: search for turns only in lower ~92% of y-range (pre-port)

interp_no_extrap = @(x,y,yq) interp1(y, x, yq, 'linear', NaN);

%% ---------- Load ----------
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

L = idxBC(choice(idxBC)==1);
R = idxBC(choice(idxBC)==0);
assert(~isempty(L) && ~isempty(R), 'Need both Left and Right BC trials.');

%% ---------- Center-arm estimate ----------
centerX = estimate_center_x(tracking, behavTrials, [L(:); R(:)]);

%% ---------- Build template (center-arm only) ----------
allYmin = inf; allYmax = -inf;
for tr = [L(:); R(:)]'
    [~,y] = xy_unique_center(tracking, behavTrials, tr, centerX, bandX);
    if ~isempty(y), allYmin = min(allYmin, min(y)); allYmax = max(allYmax, max(y)); end
end
assert(isfinite(allYmin)&&isfinite(allYmax)&&allYmax>allYmin, 'Cannot set Y grid.');
yG = linspace(allYmin, allYmax, NBins).';

% // NEW: cap for turn search (avoid port zone)
turnYMax = yG(max(1, round(turnYMaxFrac * numel(yG))));

% Interpolate BC center-arm trajectories to yG
Xl = nan(NBins, numel(L));  Xr = nan(NBins, numel(R));
for i=1:numel(L)
    [x,y] = xy_unique_center(tracking, behavTrials, L(i), centerX, bandX);
    if numel(y)>=2, Xl(:,i) = interp_no_extrap(x,y,yG); end
end
for i=1:numel(R)
    [x,y] = xy_unique_center(tracking, behavTrials, R(i), centerX, bandX);
    if numel(y)>=2, Xr(:,i) = interp_no_extrap(x,y,yG); end
end

% Template stats
muL = mean(Xl,2,'omitnan'); muR = mean(Xr,2,'omitnan');
sdL = std(Xl,0,2,'omitnan'); sdR = std(Xr,0,2,'omitnan');
sdL(~isfinite(sdL)|sdL==0) = 1e-6;  sdR(~isfinite(sdR)|sdR==0) = 1e-6;
dprime = abs(muL - muR) ./ sqrt(0.5*(sdL.^2 + sdR.^2));

%% ---------- Per-trial COMMIT & TURN ----------
conds = {idxBC, idxBE, idxST};
S = struct('trials',[],'commitY',[],'commitX',[],'turnY',[],'turnX',[]);
for c = 1:3
    ids   = conds{c};
    commY = nan(size(ids)); commX = nan(size(ids));
    turnY = nan(size(ids)); turnX = nan(size(ids));
    for k = 1:numel(ids)
        tr = ids(k);

        % --- COMMIT: center-arm space
        [xc, yc] = xy_unique_center(tracking, behavTrials, tr, centerX, bandX);
        if numel(yc) >= 2
            xq_center = interp_no_extrap(xc, yc, yG);
            ci = commit_idx(xq_center);                 % returns [] if none
            if ~isempty(ci)
                if yG(ci) >= yMinCommit                 % // NEW: enforce Y gate
                    commY(k) = yG(ci);
                    commX(k) = xq_center(ci);
                else
                    % try the next qualifying bin at/above yMinCommit
                    above = find((1:numel(yG))'>=ci & yG>=yMinCommit, 1, 'first');
                    if ~isempty(above) && isfinite(xq_center(above))
                        commY(k) = yG(above);
                        commX(k) = xq_center(above);
                    end
                end
            end
        end

        % --- TURN: FULL trajectory, search only up to turnYMax
        [xf, yf] = xy_unique(tracking, behavTrials, tr);
        if numel(yf) >= 2
            xq_full = interp_no_extrap(xf, yf, yG);
            [yt, xt] = turn_point(xq_full, yG);
            turnY(k) = yt;  turnX(k) = xt;
        end
    end
    S(c).trials  = ids;
    S(c).commitY = commY;  S(c).commitX = commX;
    S(c).turnY   = turnY;  S(c).turnX   = turnX;
end

%% ---------- Figures ----------
figure('Name','Commit vs Turn'); 
t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% Panel A
ax1 = nexttile; hold(ax1,'on');
for k=1:numel(idxBC)
    tr=idxBC(k);
    idx = InIntervals(tracking.timestamps,[behavTrials.timestamps(tr,1), behavTrials.choiceTS(tr)]);
    plot(ax1, tracking.position.x(idx), tracking.position.y(idx), '-', 'Color',[0 0 0 0.15]);
end
shadeMeanSD(ax1, muL, sdL, yG, [0.25 0.45 0.95]);
shadeMeanSD(ax1, muR, sdR, yG, [0.95 0.35 0.35]);
plot(ax1, muL, yG, 'b', 'LineWidth',2);
plot(ax1, muR, yG, 'r', 'LineWidth',2);
title(ax1,'BC trajectories + template (mean±SD)'); xlabel(ax1,'X'); ylabel(ax1,'Y');
axis(ax1,'equal'); grid(ax1,'on'); box(ax1,'on');

% Panel B
ax2 = nexttile; hold(ax2,'on');
for c=1:3
    ids=conds{c};
    for k=1:numel(ids)
        tr  = ids(k);
        idx = InIntervals(tracking.timestamps,[behavTrials.timestamps(tr,1), behavTrials.choiceTS(tr)]);
        plot(ax2, tracking.position.x(idx), tracking.position.y(idx), '-', 'Color',[colors{c} 0.18]);
        if isfinite(S(c).commitY(k))
            plot(ax2, S(c).commitX(k), S(c).commitY(k), 'o', 'MarkerSize',dotSize, ...
                 'MarkerFaceColor', colors{c}, 'MarkerEdgeColor','k');
        end
    end
end
shadeMeanSD(ax2, muL, sdL, yG, [0.25 0.45 0.95]);
shadeMeanSD(ax2, muR, sdR, yG, [0.95 0.35 0.35]);
plot(ax2, muL, yG, 'b', 'LineWidth',2);
plot(ax2, muR, yG, 'r', 'LineWidth',2);
title(ax2,'BC/BE/Stim + commit points (Y≥70)'); xlabel(ax2,'X'); ylabel(ax2,'Y');
axis(ax2,'equal'); grid(ax2,'on'); box(ax2,'on');

% Panel C
ax3 = nexttile; hold(ax3,'on');
for c=1:3
    ids=conds{c};
    for k=1:numel(ids)
        tr  = ids(k);
        idx = InIntervals(tracking.timestamps,[behavTrials.timestamps(tr,1), behavTrials.choiceTS(tr)]);
        plot(ax3, tracking.position.x(idx), tracking.position.y(idx), '-', 'Color',[colors{c} 0.18]);
        if isfinite(S(c).turnY(k))
            plot(ax3, S(c).turnX(k), S(c).turnY(k), 's', 'MarkerSize',dotSize, ...
                 'MarkerFaceColor', colors{c}, 'MarkerEdgeColor','k');
        end
    end
end
shadeMeanSD(ax3, muL, sdL, yG, [0.25 0.45 0.95]);
shadeMeanSD(ax3, muR, sdR, yG, [0.95 0.35 0.35]);
plot(ax3, muL, yG, 'b', 'LineWidth',2);
plot(ax3, muR, yG, 'r', 'LineWidth',2);
title(ax3,'BC/BE/Stim + turn points (pre‑port)'); xlabel(ax3,'X'); ylabel(ax3,'Y');
axis(ax3,'equal'); grid(ax3,'on'); box(ax3,'on');

%% ---------- Console summaries ----------
for c=1:3
    fprintf('%s  meanY(commit)=%.2f   meanY(turn)=%.2f   median ΔY(turn-commit)=%.2f  (n=%d)\n', ...
        labels{c}, mean(S(c).commitY,'omitnan'), mean(S(c).turnY,'omitnan'), ...
        median(S(c).turnY - S(c).commitY,'omitnan'), numel(S(c).trials));
end

%% ---------- Outputs ----------
OUT.S = S;
OUT.template = struct('muL',muL,'muR',muR,'sdL',sdL,'sdR',sdR, ...
                      'dprime',dprime,'yG',yG,'centerX',centerX);
OUT.labels = labels;

%% ---------- Nested ----------
    function ci = commit_idx(xq_center)
        llL = -0.5*((xq_center - muL)./sdL).^2 - log(sdL);
        llR = -0.5*((xq_center - muR)./sdR).^2 - log(sdR);
        llr = llL - llR;  % >0 Left, <0 Right

        startBin = find(isfinite(dprime) & dprime >= dprimeThr, 1, 'first');
        if isempty(startBin), ci = []; return; end

        ok = ((1:numel(llr))' >= startBin) & isfinite(llr) & (abs(llr) >= llrThr);
        % // NEW: also require Y >= yMinCommit directly here
        ok = ok & (yG >= yMinCommit);
        ci = find(ok, 1, 'first');
    end
end

%% ===================== Helpers (unchanged except turn_point) =====================
function cx = estimate_center_x(tracking, behavTrials, trials)
vals = [];
for tr = trials(:)'
    [x,y] = xy_unique(tracking, behavTrials, tr);
    if numel(y) < 10, continue; end
    ylo = quantile(y,0.05); yhi = quantile(y,0.40);
    keep = (y >= ylo) & (y <= yhi);
    if any(keep), vals(end+1) = median(x(keep)); end %#ok<AGROW>
end
cx = median(vals, 'omitnan');
end

function [x,y] = xy_unique_center(tracking, behavTrials, tr, centerX, bandX)
[x,y] = xy_unique(tracking, behavTrials, tr);
keep = abs(x - centerX) <= bandX;
x = x(keep); y = y(keep);
end

function [x,y] = xy_unique(tracking, behavTrials, tr)
idx = InIntervals(tracking.timestamps,[behavTrials.timestamps(tr,1), behavTrials.choiceTS(tr)]);
x = tracking.position.x(idx); 
y = tracking.position.y(idx);
[y,ord] = sort(y(:)); x = x(ord);
[y,ui] = unique(y);   x = x(ui);
end

function shadeMeanSD(ax, mu, sd, y, rgb)
if all(~isfinite(mu)), return; end
lo = mu - sd; hi = mu + sd;
keep = isfinite(lo) & isfinite(hi) & isfinite(y);
if ~any(keep), return; end
xBand = [lo(keep); flipud(hi(keep))];
yBand = [y(keep);  flipud(y(keep))];
ph = patch(ax, xBand, yBand, rgb, 'EdgeColor','none', 'FaceAlpha',0.15);
uistack(ph,'bottom');
end

function y0 = first_run_y(sig, yG, minRun)
d = diff([false; sig(:); false]);
starts = find(d==1); ends = find(d==-1)-1;
lens = ends - starts + 1;
k = find(lens>=minRun, 1, 'first');
if isempty(k), y0 = NaN; else, y0 = yG(starts(k)); end
end

function [yTurn, xTurn] = turn_point(xq, yG)
xq = fillmissing(xq,'linear','EndValues','nearest');
xq = smoothdata(xq,'movmedian',5);
xq = smoothdata(xq,'movmean',7);
dy   = median(diff(yG),'omitnan');
dxdy = gradient(xq) ./ dy;
dxdy_s = smoothdata(dxdy, 'movmean', 5);

thresh = 0.5 * max(abs(dxdy_s));  % 50% of max slope magnitude
sig = abs(dxdy_s) >= thresh;
% find earliest run of >= 3 bins above threshold
d = diff([false; sig(:); false]);
starts = find(d==1); ends = find(d==-1)-1;
lens = ends - starts + 1;
k = find(lens>=3, 1, 'first');
if isempty(k)
    yTurn = NaN; xTurn = NaN; return;
end
idx = starts(k);
yTurn = yG(idx);
xTurn = xq(idx);
end


