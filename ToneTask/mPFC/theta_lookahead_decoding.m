function out = theta_lookahead_decoding(sess, decodingPath, varargin)
% theta_lookahead_decoding
% Lookahead = decoded y-position minus true y-position.
% Aligns lookahead to theta peaks (from your saved thetaLFP) and
% builds peri-peak heatmaps/averages for baseline (stim==0) vs stim (stim==1).
%
% Required inputs
%   sess         : session folder containing *thetaLFP.mat, *Tracking.Behavior.mat,
%                  *TrialBehavior.Behavior.mat
%   decodingPath : NetCDF path with variables:
%                  y_position (posterior), time, y_position_value
%
% Options (Name,Value)
%   'Win'          : [-0.2 0.2]   peri-peak window (s)
%   'SpeedThresh'  : 2            forward-motion gate (cm/s)
%   'AlignTol'     : 0.05         |tracking t - posterior t| max (s)
%   'NumStd'       : 2            amplitude gate for peaks (see below)
%   'MinPeakDist'  : 0.06         minimum distance between peaks (s)
%   'UsePowerGate' : false        if true, gate peaks by lfp.thetapower instead of |filtered|
%   'SaveFig'      : ''           if nonempty, saves PNG+FIG with this basename
%
% Output struct:
%   out.timeAxis
%   out.base / out.stim: periErr, mean, sem, N
%   out.params

% -------------------- options --------------------
p = inputParser;
addParameter(p,'Win',[-0.2 0.2]);
addParameter(p,'SpeedThresh',2);
addParameter(p,'AlignTol',0.05);
addParameter(p,'NumStd',2);
addParameter(p,'MinPeakDist',0.06);
addParameter(p,'UsePowerGate',false);
addParameter(p,'SaveFig','');
parse(p,varargin{:});
opt = p.Results;

% -------------------- load session data --------------------
cwd0 = pwd; if exist(sess,'dir'); cd(sess); end

% theta LFP with fields: timestamps, filtered, thetapower, thetaphase (0..2π)
f = dir('*thetaLFP.mat');  assert(~isempty(f),'*.thetaLFP.mat not found.');
S = load(f(1).name);       lfp = fetchField(S,{'lfp','LFP','ThetaLFP'});

% tracking & trials
f = dir('*Tracking.Behavior.mat');      assert(~isempty(f),'*Tracking.Behavior.mat not found.');
S = load(f(1).name);                    tracking = fetchField(S,{'tracking','Tracking'});
f = dir('*TrialBehavior.Behavior.mat'); assert(~isempty(f),'*TrialBehavior.Behavior.mat not found.');
S = load(f(1).name);                    behavTrials = fetchField(S,{'behavTrials','Behavior'});
if isfield(behavTrials,'Timestamps') && ~isfield(behavTrials,'timestamps')
    behavTrials.timestamps = behavTrials.Timestamps;
end
assert(isfield(behavTrials,'timestamps') && isfield(behavTrials,'stim'),...
    'behavTrials.timestamps or behavTrials.stim missing');

ts_track = tracking.timestamps(:);
ypos     = tracking.position.y(:);
if isfield(tracking,'position') && isfield(tracking.position,'vy')
    vy = tracking.position.vy(:);
else
    vy = [0; diff(ypos)] ./ [1; diff(ts_track)]; % crude fallback (cm/s)
end
dt_track = median(diff(ts_track),'omitnan');

% -------------------- read posterior --------------------
posterior_pos = ncread(decodingPath,'y_position');         % (nPos x nTime)
post_time     = ncread(decodingPath,'time');               % (nTime)
post_pos      = ncread(decodingPath,'y_position_value');   % (nPos, cm)

% -------------------- compute lookahead at tracking timestamps --------------------
lookahead = nan(size(ts_track));
for i = 1:numel(ts_track)
    [d, tscur] = min(abs(post_time - ts_track(i)));
    if d <= opt.AlignTol
        curPost = posterior_pos(:,tscur);
        [~, idxMax] = max(curPost);
        dec_cm = post_pos(idxMax);
        lookahead(i) = dec_cm - ypos(i);  % + = decoded ahead of animal
    end
end

% -------------------- theta peaks from your saved fields --------------------
lfp_ts = lfp.timestamps(:);
xf     = double(lfp.filtered(:));     % bandpassed theta signal
ph     = double(lfp.thetaphase(:));   % 0..2π -> convert to -π..π
ph(ph>pi) = ph(ph>pi) - 2*pi;

% amplitude gate
if opt.UsePowerGate
    pow = double(lfp.thetapower(:));
    thr = mean(pow,'omitnan') + opt.NumStd*std(pow,[],'omitnan');
    ampGate = pow > thr;
else
    rmssig = abs(xf);
    thr = opt.NumStd*std(rmssig,'omitnan');
    ampGate = xf > thr;               % positive peaks only (as requested)
end

% rising zero-crossings of phase
phaseUp = [false; diff(ph>0)] > 0;

pkIdx = find(ampGate & phaseUp);

% enforce minimum inter-peak distance
if ~isempty(pkIdx)
    keep = true(size(pkIdx));
    lastT = -Inf;
    for k = 1:numel(pkIdx)
        t = lfp_ts(pkIdx(k));
        if (t - lastT) < opt.MinPeakDist
            keep(k) = false;
        else
            lastT = t;
        end
    end
    pkIdx = pkIdx(keep);
end
thetaTimes_all = lfp_ts(pkIdx);

% -------------------- assign peaks to baseline vs stim forward runs --------------------
basePeaks = peaksInTrials(thetaTimes_all, behavTrials.timestamps(behavTrials.stim==0,:), ts_track, vy, opt.SpeedThresh);
stimPeaks = peaksInTrials(thetaTimes_all, behavTrials.timestamps(behavTrials.stim==1,:), ts_track, vy, opt.SpeedThresh);

% -------------------- peri-event matrices --------------------
timeAxis = opt.Win(1):dt_track:opt.Win(2);
periBase = buildPeri(lookahead, ts_track, basePeaks, timeAxis);
periStim = buildPeri(lookahead, ts_track, stimPeaks, timeAxis);

% -------------------- plot --------------------
fh = figure('Color','w','Position',[180 120 1150 780]);

subplot(2,2,1);
imagesc(timeAxis, 1:size(periBase,1), periBase); axis xy;
colormap(gca,'magma'); colorbar; caxisForBoth(gca, periBase, periStim);
title(sprintf('Baseline (stim==0), N=%d peaks', size(periBase,1)));
xlabel('Time from theta peak (s)'); ylabel('Event #');

subplot(2,2,2);
imagesc(timeAxis, 1:size(periStim,1), periStim); axis xy;
colormap(gca,'magma'); colorbar; caxisForBoth(gca, periBase, periStim);
title(sprintf('mPFC silencing (stim==1), N=%d peaks', size(periStim,1)));
xlabel('Time from theta peak (s)'); ylabel('Event #');

subplot(2,2,3);
[mB,seB] = meanSem(periBase);
[mS,seS] = meanSem(periStim);
hold on;
plot(timeAxis, mB, 'k', 'LineWidth', 2);
plot(timeAxis, mB+seB, 'k--', timeAxis, mB-seB, 'k--');
plot(timeAxis, mS, 'g', 'LineWidth', 2);
plot(timeAxis, mS+seS, 'g--', timeAxis, mS-seS, 'g--');
yline(0,'k:');
xlabel('Time from theta peak (s)'); ylabel('Decode - true pos (cm)');
legend({'Baseline','Baseline \pm SEM','Silencing','Silencing \pm SEM'},'Location','best');
title('Theta lookahead (decoded - true)');

subplot(2,2,4);
% simple modulation score: post - pre
preWin  = timeAxis>=-0.08 & timeAxis<=-0.02;
postWin = timeAxis>=+0.02 & timeAxis<=+0.08;
msB = rowfun(@(x) nanmean(x(postWin)) - nanmean(x(preWin)), periBase);
msS = rowfun(@(x) nanmean(x(postWin)) - nanmean(x(preWin)), periStim);
hold on;
histogram(msB, 'Normalization','probability');
histogram(msS, 'Normalization','probability');
xlabel('Modulation score (cm)'); ylabel('Probability');
legend({'Baseline','Silencing'}); box off;
title('Per-peak modulation (post - pre)');

if ~isempty(opt.SaveFig)
    saveas(fh, [opt.SaveFig '.png']);
    savefig(fh, [opt.SaveFig '.fig']);
end

% -------------------- outputs --------------------
out = struct();
out.timeAxis = timeAxis;
out.base = packOut(periBase, mB, seB);
out.stim = packOut(periStim, mS, seS);
out.params = opt;

if exist('cwd0','var'); cd(cwd0); end
end

% =================== helpers ===================
function val = fetchField(S, names)
for k = 1:numel(names)
    if isfield(S,names{k}), val = S.(names{k}); return; end
end
error('Missing field. Tried: %s', strjoin(names,', '));
end

function peaks = peaksInTrials(pkTimes, intervals, ts_track, vy, speedThresh)
peaks = [];
if isempty(intervals) || isempty(pkTimes), return; end
for i = 1:size(intervals,1)
    t1 = intervals(i,1); t2 = intervals(i,2);
    sel = pkTimes>=t1 & pkTimes<=t2;
    if any(sel)
        tk = pkTimes(sel);
        % forward-motion gate at peak time
        vAtPk = interp1(ts_track, vy, tk, 'linear', 'extrap');
        tk = tk(vAtPk > speedThresh);
        peaks = [peaks; tk]; %#ok<AGROW>
    end
end
end

function peri = buildPeri(signal, sigT, eventTimes, tAxis)
peri = nan(numel(eventTimes), numel(tAxis));
if isempty(eventTimes), return; end
dt_median = median(diff(sigT),'omitnan');
for i = 1:numel(eventTimes)
    tgt = eventTimes(i) + tAxis;
    idx = interp1(sigT, 1:numel(sigT), tgt, 'nearest', 'extrap');
    dt  = abs(sigT(idx) - tgt);
    bad = dt > 1.5*dt_median;         % reject poor alignments
    row = signal(idx);
    row(bad) = NaN;
    peri(i,:) = row;
end
end

function [m,se] = meanSem(X)
m  = nanmean(X,1);
n  = sum(~isnan(X),1);
se = nanstd(X,[],1) ./ sqrt(max(n,1));
end

function S = packOut(peri, m, se)
S = struct('periErr',peri,'mean',m,'sem',se,'N',size(peri,1));
end

function r = rowfun(f, M)
if isempty(M), r = []; return; end
r = nan(size(M,1),1);
for i = 1:size(M,1), r(i) = f(M(i,:)); end
end

function caxisForBoth(ax, A, B)
lims = [nanmin([A(:);B(:)]), nanmax([A(:);B(:)])];
if all(isfinite(lims)), caxis(ax, lims); end
end
