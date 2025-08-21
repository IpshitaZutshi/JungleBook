function out = theta_lookahead_decoding(sess, decodingPath, varargin)
% theta_lookahead_decoding
% Lookahead = decoded y-position minus true y-position.
% Key change: compute lookahead on the DECODER timebase (post_time), not the
% tracking timebase, then align to theta peaks (from thetaLFP).
%
% Inputs
%   sess         : session folder with *thetaLFP.mat, *Tracking.Behavior.mat,
%                  *TrialBehavior.Behavior.mat
%   decodingPath : NetCDF with variables:
%                  y_position (posterior), time, y_position_value
%
% Options (Name,Value)
%   'Win'            : [-0.2 0.2] peri-peak window (s)
%   'SpeedThresh'    : 2          forward-motion gate (cm/s)
%   'NumStd'         : 2          amplitude gate for peaks
%   'MinPeakDist'    : 0.06       minimum distance between peaks (s)
%   'UsePowerGate'   : false      gate by thetapower instead of |filtered|
%   'UsePosteriorMean': false     if true, use E[pos] (mean), else MAP
%   'SaveFig'        : ''         basename to save PNG+FIG
%
% Output:
%   out.timeAxis
%   out.base / out.stim: periErr, mean, sem, N
%   out.params

% ---------- options ----------
p = inputParser;
addParameter(p,'Win',[-0.2 0.2]);
addParameter(p,'SpeedThresh',2);
addParameter(p,'NumStd',2);
addParameter(p,'MinPeakDist',0.06);
addParameter(p,'UsePowerGate',false);
addParameter(p,'UsePosteriorMean',false);
addParameter(p,'SaveFig','');
parse(p,varargin{:});
opt = p.Results;

% ---------- load session data ----------
cwd0 = pwd; if exist(sess,'dir'); cd(sess); end

% theta LFP (fields: timestamps, filtered, thetapower, thetaphase 0..2π)
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

% Tracking series
ts_track = tracking.timestamps(:);
y_track  = tracking.position.y(:);
if isfield(tracking,'position') && isfield(tracking.position,'vy')
    vy_track = tracking.position.vy(:);
else
    vy_track = [0; diff(y_track)]./[1; diff(ts_track)];
end

% ---------- read decoder ----------
posterior_pos = ncread(decodingPath,'y_position');         % (nPos x nTime)
post_time     = ncread(decodingPath,'time');               % (nTime)
post_pos      = ncread(decodingPath,'y_position_value');   % (nPos, cm)

% Decoder time step & peri window axis (fixed length)
dts = diff(post_time); dts = dts(isfinite(dts) & dts>0);
dt_post = median(dts,'omitnan');
nSamp   = max(5, round((opt.Win(2)-opt.Win(1))/dt_post) + 1);
timeAxis = linspace(opt.Win(1), opt.Win(2), nSamp);

% ---------- compute lookahead on decoder timebase ----------
% 1) True y at decoder times
y_at_post  = interp1(ts_track, y_track, post_time, 'linear', NaN);
vy_at_post = interp1(ts_track, vy_track, post_time, 'linear', NaN);

% 2) Decoded y at decoder times
if opt.UsePosteriorMean
    % Expected value of position (cm)
    normPost = posterior_pos ./ max(eps, sum(posterior_pos,1));
    dec_y = post_pos(:)' * normPost;    % 1 x nTime
    dec_y = dec_y(:);                   % column
else
    % MAP
    [~,idxMax] = max(posterior_pos, [], 1);
    dec_y = post_pos(idxMax)';          % column
end

% 3) Lookahead at decoder times
lookahead_post = dec_y - y_at_post;     % column
validDecode = isfinite(lookahead_post) & ~isnan(lookahead_post) & (sum(posterior_pos,1)'>0);
% keep only valid decoder samples
lookahead_post(~validDecode) = NaN;

% ---------- theta peaks (using your lfp fields) ----------
lfp_ts = lfp.timestamps(:);
xf     = double(lfp.filtered(:));
ph     = double(lfp.thetaphase(:)); ph(ph>pi) = ph(ph>pi) - 2*pi;

if opt.UsePowerGate
    pow = double(lfp.thetapower(:));
    thr = mean(pow,'omitnan') + opt.NumStd*std(pow,[],'omitnan');
    ampGate = pow > thr;
else
    rmssig = abs(xf);
    thr = opt.NumStd*std(rmssig,'omitnan');
    ampGate = xf > thr; % positive peaks
end
phaseUp = [false; diff(ph>0)] > 0;
pkIdx = find(ampGate & phaseUp);

% min inter-peak spacing
if ~isempty(pkIdx)
    keep = true(size(pkIdx)); lastT = -Inf;
    for k=1:numel(pkIdx)
        t=lfp_ts(pkIdx(k));
        if (t - lastT) < opt.MinPeakDist, keep(k)=false; else, lastT=t; end
    end
    pkIdx = pkIdx(keep);
end
thetaTimes_all = lfp_ts(pkIdx);

% ---------- select peaks within trials & forward motion ----------
% Gates evaluated at peak times using tracking speed
basePeaks = peaksInTrials(thetaTimes_all, behavTrials.timestamps(behavTrials.stim==0,:), ...
                          ts_track, vy_track, opt.SpeedThresh);
stimPeaks = peaksInTrials(thetaTimes_all, behavTrials.timestamps(behavTrials.stim==1,:), ...
                          ts_track, vy_track, opt.SpeedThresh);

% ---------- build peri-event matrices on decoder timebase ----------
periBase = buildPeri(lookahead_post, post_time, basePeaks, timeAxis);
periStim = buildPeri(lookahead_post, post_time, stimPeaks, timeAxis);

% Optional: also require forward motion **within** the window (decoder speed)
% (uncomment if desired)
% vperiBase = buildPeri(vy_at_post, post_time, basePeaks, timeAxis);
% vperiStim = buildPeri(vy_at_post, post_time, stimPeaks, timeAxis);
% spmaskB = vperiBase > opt.SpeedThresh; periBase(~spmaskB) = NaN;
% spmaskS = vperiStim > opt.SpeedThresh; periStim(~spmaskS) = NaN;

% ---------- plot ----------
fh = figure('Color','w','Position',[180 120 1150 780]);
subplot(2,2,1); heatmapSafe(timeAxis, periBase, 'Baseline (stim==0)');
subplot(2,2,2); heatmapSafe(timeAxis, periStim, 'mPFC silencing (stim==1)');

subplot(2,2,3);
[mB,seB] = meanSem(periBase); [mS,seS] = meanSem(periStim);
hold on;
plot(timeAxis, mB, 'k', 'LineWidth', 2);
plot(timeAxis, mB+seB, 'k--', timeAxis, mB-seB, 'k--');
plot(timeAxis, mS, 'g', 'LineWidth', 2);
plot(timeAxis, mS+seS, 'g--', timeAxis, mS-seS, 'g--');
yline(0,'k:'); xlabel('Time from theta peak (s)'); ylabel('Decode - true pos (cm)');
legend({'Baseline','Baseline \pm SEM','Silencing','Silencing \pm SEM'},'Location','best');
title('Theta lookahead (decoded - true)');

subplot(2,2,4);
preWin  = timeAxis>=-0.08 & timeAxis<=-0.02;
postWin = timeAxis>=+0.02 & timeAxis<=+0.08;
msB = rowfun(@(x) nanmean(x(postWin)) - nanmean(x(preWin)), periBase);
msS = rowfun(@(x) nanmean(x(postWin)) - nanmean(x(preWin)), periStim);
hold on;
if ~isempty(msB), histogram(msB,'Normalization','probability'); end
if ~isempty(msS), histogram(msS,'Normalization','probability'); end
xlabel('Modulation score (cm)'); ylabel('Probability');
legend({'Baseline','Silencing'}); box off;
title('Per-peak modulation (post - pre)');

if ~isempty(opt.SaveFig)
    saveas(fh, [opt.SaveFig '.png']);
    savefig(fh, [opt.SaveFig '.fig']);
end

% ---------- outputs ----------
out = struct();
out.timeAxis = timeAxis;
out.base = packOut(periBase, mB, seB);
out.stim = packOut(periStim, mS, seS);
out.params = opt;

if exist('cwd0','var'); cd(cwd0); end
end

% ================= helpers =================
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
        % forward-motion gate evaluated at peak time via tracking
        vAtPk = interp1(ts_track, vy, tk, 'linear', 'extrap');
        tk = tk(vAtPk > speedThresh);
        peaks = [peaks; tk]; %#ok<AGROW>
    end
end
end

function peri = buildPeri(signal, sigT, eventTimes, tAxis)
% Interpolate a signal defined on sigT onto windows centered at eventTimes.
signal = signal(:); sigT = sigT(:); tAxis = tAxis(:)';    % orient
[sigTuniq, iu] = unique(sigT,'stable'); signal = signal(iu);
peri = nan(numel(eventTimes), numel(tAxis));
if isempty(eventTimes), return; end
for i = 1:numel(eventTimes)
    tgt = eventTimes(i) + tAxis;                           % row
    peri(i,:) = interp1(sigTuniq, signal, tgt, 'linear', NaN);
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

function heatmapSafe(timeAxis, periMat, ttl)
if isempty(periMat)
    imagesc(timeAxis, 1, nan(1,numel(timeAxis))); axis xy;
    nstr = 'N=0 peaks';
else
    imagesc(timeAxis, 1:size(periMat,1), periMat); axis xy;
    nstr = sprintf('N=%d peaks', size(periMat,1));
end
colormap(gca,'magma'); colorbar;
vals = periMat(:); vals = vals(isfinite(vals));
if ~isempty(vals), caxis([min(vals) max(vals)]); end
xlabel('Time from theta peak (s)'); ylabel('Event #');
title(sprintf('%s, %s', ttl, nstr));
end
