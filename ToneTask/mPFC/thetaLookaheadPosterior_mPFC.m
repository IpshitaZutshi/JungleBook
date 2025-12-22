function thetaLookaheadPosterior_mPFC()

% Average posterior in RELATIVE coords over theta phase/time (2 cycles),
% split by stim vs baseline
%
% Outputs:
%   Figure 1a: mean±SEM decoded-true in 3 s pre-lick (baseline vs stim)
%   Figure 2b: Contourf of avg posterior(y_rel, phase) over 2 cycles (baseline vs stim)
%   Figure 2c: Quantification of lookahead
%
% Dependencies:
%   - catpad (to append and match matrix dimensions)


%% ----------------------- USER SETTINGS ----------------------- %%
basepath = pwd;
fileloc = strsplit(basepath,'\');
sessRel = fileloc{end};
mouseId = sessRel(1:4);

decodingPath     = strcat('Z:\Homes\zz737\ipshita_data\Auditory_Task\',mouseId,'\Final');
decodingName     = 'py_data/theta_decoding_lickLoc_y/up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].nc';
changePointName  = 'py_data/theta_decoding_lickLoc_y/change_point_posterior_up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].mat';

sessAbs = strcat('Z:\Homes\zutshi01\Recordings\Auditory_Task\mPFC\',mouseId,'\',sessRel);     % <-- your session folder
cd(sessAbs);

% Relative-position axis
yRelMin = -125;  % behind (past)
yRelMax =  125;  % ahead  (future)
nYrel   = 60;  % edges count; bins = nYrel-1

% Theta visualization bins
nPhaseBinsPerCycle = 1800;            % per theta cycle
sigmaSmooth = 1.5;                  % smoothing for plot

% Pre-lick window
preLickSec = 2.5;

velThresh = 5;

%% ----------------------- LOAD BEHAV / LFP ----------------------- %%
file = dir('*.Tracking.Behavior.mat');          load(file.name);  
file = dir('*.TrialBehavior.Behavior.mat');     load(file.name);  
file = dir('*.thetaLFP.mat');                   load(file.name);  

% Load decoding data (posterior over y-position)
[posterior_pos, post_time, post_pos] = loadDecodingData( ...
    decodingPath, sessRel, decodingName, changePointName);

%% ---- Session-level decoder QC (median absolute error) ----
[~, idxMax] = max(posterior_pos, [], 1);
decY = post_pos(idxMax).';   % MAP decoded position at each decoder time bin

trueY = interp1(tracking.timestamps(:), tracking.position.y(:), post_time(:), 'linear', nan);
v = interp1(tracking.timestamps(:), tracking.position.v(:), post_time(:), 'linear', nan);
move = v > 5;

err  = abs(decY(:) - trueY(:));
medAE = median(err(move), 'omitnan');

fprintf('Decoder QC: median |decoded-true| = %.2f (units of tracking.position.y)\n', medAE);

% Skip bad sessions (choose threshold in same units as position)
medAE_thresh = 10;  % e.g., 15 cm (EDIT)
if medAE > medAE_thresh
    warning('Skipping session: decoder too poor (medAE=%.2f > %.2f).', medAE, medAE_thresh);
    return
end

% Ensure column vectors
post_time = post_time(:);
post_pos  = post_pos(:);

% Only select tone trials 
trialNum = find(behavTrials.linTrial == 0);

%% ----------------------- ACCUMULATORS ----------------------- %%
% Pre-lick decoded-true traces (variable length handled by padding)
frameRate = 30;
nSamp = max(1, round(preLickSec/(1/frameRate)) + 1);
decodedDiff  = nan(numel(trialNum),nSamp);                              % [nTrials x nSamp]
stimTrial    = false(numel(trialNum),1);
avgVel       = nan(numel(trialNum),1);

% Posterior accumulators (relative position × 2 cycles)
yRelEdges = linspace(yRelMin, yRelMax, nYrel);
yRelC     = (yRelEdges(1:end-1)+yRelEdges(2:end))/2;

nX = 2*nPhaseBinsPerCycle;                      % two cycles stitched (0..4?)
basePostSum = zeros(numel(yRelC), nX);
stimPostSum = zeros(numel(yRelC), nX);
baseN = 0; stimN = 0;

%% ----------------------- MAIN LOOP ----------------------- %%
for ii = 1:numel(trialNum)

    % Skip trials where the first port was licked because there's not enough time
    if behavTrials.lickLoc(ii)==0
        continue
    end

    curTrial = trialNum(ii);
    tWin = behavTrials.timestamps(curTrial,:);

    % Shift start to first v>velThresh within trial
    idxTrialTracking = find(tracking.timestamps > tWin(1) & tracking.timestamps < tWin(2));
    idxVel = find(tracking.position.v(idxTrialTracking) > velThresh, 1, 'first');
    if ~isempty(idxVel)
        tWin(1) = tracking.timestamps(idxTrialTracking(1) + idxVel - 1);
    end

    stimTrial(ii) = logical(behavTrials.stim(curTrial));

    % ----- Build tracking time segment for velocity + pre-lick trace -----
    [posTrk, velTrk, tTrk] = extractTrackingSegment(tWin, tracking);
    if isempty(tTrk)
        continue
    end

    % Average velocity in pre-lick window
    idxV = tTrk >= (tWin(2)-preLickSec) & tTrk <= tWin(2);
    if any(idxV)
        avgVel(ii) = mean(velTrk(idxV), 'omitnan');
    end

    % ----- Decoder indices for this trial -----
    [idxDec1, idxDec2] = decoderIndicesForWindow(tWin,post_time);
    if isempty(idxDec1)
        continue
    end

    % Decoder time base for this trial
    tDec = post_time(idxDec1:idxDec2);

    % True position on decoder time base
    trueY = interp1(tTrk, posTrk, tDec, 'linear', nan);

    % Theta phase on decoder time base
    phiDec = interp1(lfp.timestamps(:), lfp.thetaphase(:), tDec, 'linear', nan);

    % Posterior for this trial window
    Pwin = posterior_pos(:, idxDec1:idxDec2);   % [nPos x nTdec]

    if isempty(Pwin) || all(~isfinite(trueY)) || all(~isfinite(phiDec))
        continue
    end

    %% ----- Pre-lick decoded-true trace (decoder -> interpolated to tracking) -----
    % Position on decoder bins
    [~, idxPosMax] = max(Pwin, [], 1);
    decY = post_pos(idxPosMax).';

    % Interpolate decoder onto tracking timestamps for this window
    decY_trk = interp1(tDec, decY(:), tTrk(:), 'linear', nan);

    % Compute pre-lick diff row (padded)
    decodedDiff(ii,:) = computePreWindowDiff(posTrk, decY_trk, tTrk, tWin(2), preLickSec, nSamp);

    % ----- Accumulate RELATIVE posterior over 2 theta cycles -----
    troughIdx = findThetaTroughs(phiDec);
    if numel(troughIdx) < 3
        continue
    end

    for k = 1:(numel(troughIdx)-2)

        i1 = troughIdx(k);
        i3 = troughIdx(k+2) - 1;     % end of 2 cycles
        if i3 <= i1, continue; end

        % slice 2-cycle posterior and true pos
        P2 = Pwin(:, i1:i3);          % [nPos x nT2]
        y2 = trueY(i1:i3);            % [nT2 x 1]
        if any(~isfinite(y2)) || isempty(P2), continue; end
        
                % --- speed threshold (cm/s) ---
        v2 = interp1(tracking.timestamps, tracking.position.v, ...
                     post_time(idxDec1+i1-1:idxDec1+i3-1), 'linear', nan);

        if mean(v2,'omitnan') < 5
            continue
        end            

        % Convert absolute posterior -> relative posterior, then resample to 2*nPhaseBins
        R2 = relPosterior2cycles(P2, post_pos, y2, yRelEdges, nPhaseBinsPerCycle);

        if stimTrial(ii)
            stimPostSum(:,:,end+1) = R2; 
            stimN = stimN+1;
        else
            basePostSum(:,:,end+1) = R2;
            baseN = baseN+1;
        end
    end
end

%% ----------------------- VELOCITY STATS (OPTIONAL) ----------------------- %%
baseVel = avgVel(~stimTrial & isfinite(avgVel));
stimVel = avgVel(stimTrial  & isfinite(avgVel));
fprintf('Avg velocity (%.1fs pre-lick): baseline n=%d, stim n=%d\n', preLickSec, mean(baseVel), mean(stimVel));
if exist('ranksum','file') == 2
    [p,~,stats] = ranksum(baseVel, stimVel);
    fprintf('ranksum p = %.3g, z = %.3f\n', p, stats.zval);
else
    fprintf('(ranksum not available: Statistics toolbox missing)\n');
end

%% ----------------------- FIGURE 1: PRE-LICK TRACE ----------------------- %%
nSamp = size(decodedDiff,2);
t = linspace(-preLickSec, 0, nSamp);

baseMat = decodedDiff(~stimTrial, :);
stimMat = decodedDiff(stimTrial,  :);

mBase = mean(baseMat, 1, 'omitnan');
mStim = mean(stimMat, 1, 'omitnan');

nBase = sum(isfinite(baseMat), 1);
nStim = sum(isfinite(stimMat), 1);

semBase = std(baseMat, 0, 1, 'omitnan') ./ sqrt(max(nBase,1));
semStim = std(stimMat, 0, 1, 'omitnan') ./ sqrt(max(nStim,1));

figure; hold on
shadedErrorBar(t, mBase, semBase);
shadedErrorBar(t, mStim, semStim,[0 0 1]);
xlabel(sprintf('Time (s) relative to lick/end'));
ylabel('Decoded - true position (MAP)');
legend({'Baseline','Stim'}, 'Location','best');
title(sprintf('Decoded difference (%.1f s pre-lick): mean \\pm SEM', preLickSec));
grid on


%% ----------------------- FIGURE 2: AVG POSTERIOR ----------------------- %%
baseAvg = mean(basePostSum,3);
stimAvg = mean(stimPostSum,3);

baseAvg = smooth2D(baseAvg, sigmaSmooth);
stimAvg = smooth2D(stimAvg, sigmaSmooth);

phaseAxis = linspace(0, 4*pi, nX);

%% ---------- Posterior mean offset per phase (E[y_rel] vs phase) ----------
% basePostSum / stimPostSum are [nY x nX x nWindows]
% yRelC is [nY x 1]

ycol = yRelC(:);                          % nY x 1

% Expected y_rel at each phase for each window: E[y_rel] = sum_y y * P(y|phase)
baseExp = squeeze( sum(basePostSum .* reshape(ycol,[],1,1), 1) );  % nX x nBaseWin
stimExp = squeeze( sum(stimPostSum .* reshape(ycol,[],1,1), 1) );  % nX x nStimWin

% Mean + SEM across windows (3rd-dim slices)
mBase = mean(baseExp, 2, 'omitnan');
mStim = mean(stimExp, 2, 'omitnan');

nBase = sum(isfinite(baseExp), 2);
nStim = sum(isfinite(stimExp), 2);

semBase = std(baseExp, 0, 2, 'omitnan') ./ sqrt(max(nBase,1));
semStim = std(stimExp, 0, 2, 'omitnan') ./ sqrt(max(nStim,1));

% Plot
figure; hold on
subplot(1,3,1)
shadedErrorBar(phaseAxis, mBase, semBase, [0 0 0]);   % baseline black
shadedErrorBar(phaseAxis, mStim, semStim, [1 0 0]);   % stim red
yline(0,'k:');
xticks([0 2*pi 4*pi]); xticklabels({'0','2\pi','4\pi'})
xlabel('Theta phase/time (2 cycles)');
ylabel('Posterior mean offset E[y_{rel}] (cm)');
title('Posterior mean lookahead/lookbehind vs theta phase (mean ± SEM across windows)');
legend({'Baseline','Stim'}, 'Location','best');
grid on

cl = [min([baseAvg(:); stimAvg(:)]), max([baseAvg(:); stimAvg(:)])];
levels = linspace(cl(1), cl(2), 20);

subplot(1,3,2)
contourf(phaseAxis, yRelC, baseAvg, levels, 'LineColor','none');
set(gca,'YDir','normal')
xticks([0 2*pi 4*pi]); xticklabels({'0','2\pi','4\pi'})
xlabel('Theta phase/time (2 cycles)'); ylabel('Relative position (y - current y)');
title(sprintf('Baseline (all cycles), N=%d', baseN));
colormap(jet); caxis(cl); colorbar
ylim([-20 20])

subplot(1,3,3)
contourf(phaseAxis, yRelC, stimAvg, levels, 'LineColor','none');
set(gca,'YDir','normal')
xticks([0 2*pi 4*pi]); xticklabels({'0','2\pi','4\pi'})
xlabel('Theta phase/time (2 cycles)'); ylabel('Relative position (y - current y)');
title(sprintf('Stim (all cycles), N=%d', stimN));
colormap(jet); caxis(cl); colorbar
ylim([-20 20])

end

%% ======================= HELPERS ======================= %%
function [posterior_pos, post_time, post_pos] = loadDecodingData(decodingPath, sessRel, decodingName, changePointName) %#ok<INUSD>
file_nc = fullfile(decodingPath, sessRel, decodingName);
posterior_pos = ncread(file_nc, 'y_position');
post_time     = ncread(file_nc, 'time');
post_pos      = ncread(file_nc, 'y_position_value');
end

function [pos, vel, t] = extractTrackingSegment(tWin, tracking)
[~, i1] = min(abs(tracking.timestamps - tWin(1)));
[~, i2] = min(abs(tracking.timestamps - tWin(2)));
i1 = max(1,i1); i2 = min(numel(tracking.timestamps),i2);
if i2 <= i1
    pos=[]; vel=[]; t=[];
    return
end
pos = tracking.position.y(i1:i2);
vel = tracking.position.v(i1:i2);
t   = tracking.timestamps(i1:i2);
end

function [idxDec1, idxDec2] = decoderIndicesForWindow(tWin,post_time)
[~, idxDec1] = min(abs(post_time - tWin(1)));
[~, idxDec2] = min(abs(post_time - tWin(2)));

% keep "half-second" guard behavior, because the decoder just skips time
% periods, so it might end up aligning with a different trial
if abs(post_time(idxDec1) - tWin(1)) > 0.5 && idxDec1 < numel(post_time)
    idxDec1 = idxDec1 + 1;
end
if abs(post_time(idxDec2) - tWin(2)) > 0.5 && idxDec2 > 1
    idxDec2 = idxDec2 - 1;
end

idxDec1 = max(1, idxDec1);
idxDec2 = min(numel(post_time), idxDec2);

if (idxDec2 <= idxDec1) || any(diff(post_time(idxDec1:idxDec2)) > 0.011)
    idxDec1 = [];
    idxDec2 = [];
end
end

function diffRow = computePreWindowDiff(posTrk, decTrk, tTrk, tEnd, winSec, nSamp)
posTrk = posTrk(:);
decTrk = decTrk(:);
tTrk   = tTrk(:);

diffSig = decTrk - posTrk;

t0 = tEnd - winSec;
t1 = tEnd;

idx = find(tTrk >= t0 & tTrk <= t1);

diffRow = nan(1, nSamp);
if isempty(idx), return; end

seg = diffSig(idx).';
if numel(seg) >= nSamp
    diffRow = seg(end-nSamp+1:end);
else
    diffRow(end-numel(seg)+1:end) = seg;
end
end

function troughIdx = findThetaTroughs(phaseVec)
% Troughs near 0 using circular distance; no toolbox required.
phaseVec = phaseVec(:);
phaseVec(~isfinite(phaseVec)) = nan;

d = abs(angle(exp(1i*(phaseVec-pi)))); % 0 at phase == -pi
th = 0.5;                                % radians tolerance
% unwrap for cycle labeling
good = isfinite(phaseVec);
phu = nan(size(phaseVec));
phu(good) = unwrap(phaseVec(good));

cyc = floor(phu/(2*pi));   % cycle labels

% candidates near trough
cand = (d < th) & isfinite(d) & isfinite(cyc);

troughIdx = [];
if ~any(cand), return; end

% pick ONE index per cycle: the minimum d within that cycle
uCyc = unique(cyc(cand));
troughIdx = nan(numel(uCyc),1);

for k = 1:numel(uCyc)
    ii = find(cand & (cyc == uCyc(k)));
    [~, m] = min(d(ii));
    troughIdx(k) = ii(m);
end

troughIdx = sort(troughIdx);
end

function R2 = relPosterior2cycles(P2, post_pos, trueY, yRelEdges, nPhaseBinsPerCycle)
% Convert ABS posterior -> REL posterior (y_rel bins) for each time bin,
% normalize, then resample to fixed length = 2 cycles for visualization.

post_pos = post_pos(:);
trueY    = trueY(:);

nT = size(P2,2);
nY = numel(yRelEdges)-1;

Rfine = zeros(nY, nT);

for t = 1:nT
    yrel = post_pos - trueY(t);
    p    = P2(:,t);
    if ~all(isfinite(yrel)) || ~any(isfinite(p)), continue; end

    bin = discretize(yrel, yRelEdges);
    good = isfinite(bin) & isfinite(p);

    if any(good)
        Rfine(:,t) = accumarray(bin(good), p(good), [nY 1], @sum, 0);
    end
end

% normalize each time bin to sum to 1
cs = sum(Rfine,1);
cs(cs==0) = 1;
Rfine = Rfine ./ cs;

% resample time to fixed 2-cycle length
x0 = linspace(0,1,nT);
x1 = linspace(0,1,2*nPhaseBinsPerCycle);

R2 = interp1(x0, Rfine.', x1, 'linear', 0).';  % [nY x 2*nBins]
end

function A = smooth2D(A, sigma)
r = ceil(3*sigma);
x = -r:r;
g = exp(-(x.^2)/(2*sigma^2));
g = g/sum(g);
A = conv2(conv2(A,g,'same'),g','same');
end

function shadedErrorBar(x,y,e,col)
if nargin<4, col=[0 0 0]; end
x=x(:)'; y=y(:)'; e=e(:)';
ok=isfinite(x)&isfinite(y)&isfinite(e);
x=x(ok); y=y(ok); e=e(ok);
hold on
patch([x fliplr(x)],[y-e fliplr(y+e)],col,...
      'FaceAlpha',0.3,'EdgeColor','none');
plot(x,y,'Color',col,'LineWidth',2);
end