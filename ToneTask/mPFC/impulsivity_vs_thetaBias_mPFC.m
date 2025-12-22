function impulsivity_vs_thetaBias_mPFC()
% impulsivity_vs_thetaBias_mPFC
% Standalone analysis:
%   impulsivity(trial) = behavTrials.toneGain(target) - behavTrials.lickLoc(chosen)
%   neuralBias(trial)  = mean over time window of [Ahead - Behind] using FULL posterior
%                        Ahead = sum_{yrel>0} P(yrel)*yrel
%                        Behind= sum_{yrel<0} P(yrel)*|yrel|
%                        Bias  = Ahead - Behind
% Then correlate impulsivity with neuralBias (robust fit if available).
%
% Notes:
% - Assumes tracking.position.y, tracking.position.v, behavTrials.timestamps exist.
% - Assumes decoding .nc contains y_position (posterior), time, y_position_value.
% - No theta phase needed here; this is time-to-lick aligned.
%
% Paste into a NEW .m file and run.

%% ---------------- USER SETTINGS ----------------
sessRel = 'IZ39_220707_sess17';
mouseId = sessRel(1:4);

decodingPath     = strcat('Z:\Homes\zz737\ipshita_data\Auditory_Task\',mouseId,'\Final');
decodingName     = 'py_data/theta_decoding_lickLoc_y/up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].nc';

sessAbs = strcat('Z:\Homes\zutshi01\Recordings\Auditory_Task\mPFC\',mouseId,'\',sessRel);     % <-- your session folder
cd(sessAbs);

% Trial selection
useToneTrialsOnly = true;     % uses behavTrials.linTrial==0 (like your older scripts)

% Time window (seconds) relative to lick/end time (tWin(2))
tWinRel = [-2, -0.1];       % "around ~1s before lick"
minWinFracGood = 0.5;         % require at least 50% of bins finite in window

% Optional speed control (cm/s)
useSpeedThresh = true;
speedThresh = 5;

% Optional: also report partial correlation controlling for speed
doPartialCorr = true;

%% ---------------- LOAD FILES ----------------
file = dir('*.Tracking.Behavior.mat');      load(file.name); %#ok<LOAD>
file = dir('*.TrialBehavior.Behavior.mat'); load(file.name); %#ok<LOAD>

% Decoding (posterior over y-position)
[posterior_pos, post_time, post_pos] = loadDecodingData(decodingPath, sessRel, decodingName);
post_time = post_time(:);
post_pos  = post_pos(:);

% Trial list
if useToneTrialsOnly && isfield(behavTrials,'linTrial')
    trialList = find(behavTrials.linTrial == 0);
else
    trialList = 1:size(behavTrials.timestamps,1);
end

nT = numel(trialList);

%% ---------------- PREALLOC ----------------
impulsivity = nan(nT,1);
neuralBias  = nan(nT,1);
meanSpeed   = nan(nT,1);
stimTrial   = false(nT,1);

%% ---------------- MAIN LOOP ----------------
for ii = 1:nT
    tr = trialList(ii);

    % trial window and lick time (end)
    tWin = behavTrials.timestamps(tr,:);
    tLick = tWin(2);

    % condition label
    if isfield(behavTrials,'stim')
        stimTrial(ii) = logical(behavTrials.stim(tr));
    end

    % ----- impulsivity metric -----
    impulsivity(ii) = computeImpulsivity(behavTrials, tr);

    % ----- decoder indices for analysis window -----
    t0 = tLick + tWinRel(1);
    t1 = tLick + tWinRel(2);
    if t1 <= t0, continue; end

    [iDec0, iDec1] = decoderIndicesForWindow(post_time, [t0 t1]);
    if isempty(iDec0), continue; end

    % Guard against big decoding gaps
    if any(diff(post_time(iDec0:iDec1)) > 0.011)
        continue
    end

    tDec = post_time(iDec0:iDec1);
    Pseg = posterior_pos(:, iDec0:iDec1);  % [nPos x nBins]
    if isempty(Pseg), continue; end

    % True position on decoder time base
    trueY = interp1(tracking.timestamps(:), tracking.position.y(:), tDec, 'linear', nan);

    % Speed on decoder time base
    vDec = interp1(tracking.timestamps(:), tracking.position.v(:), tDec, 'linear', nan);
    meanSpeed(ii) = mean(vDec, 'omitnan');

    if useSpeedThresh
        if ~isfinite(meanSpeed(ii)) || meanSpeed(ii) < speedThresh
            continue
        end
    end

    % Require enough finite samples
    good = isfinite(trueY) & all(isfinite(Pseg),1).';
    if mean(good) < minWinFracGood
        continue
    end

    % ----- neural bias in that window -----
    neuralBias(ii) = computePosteriorBias(Pseg, post_pos, trueY);
end

%% ---------------- CLEAN / REPORT ----------------
ok = isfinite(impulsivity) & isfinite(neuralBias);

fprintf('Trials total: %d\n', nT);
fprintf('Trials kept (impulsivity+neuralBias finite): %d\n', nnz(ok));
if useSpeedThresh
    fprintf('Speed threshold: %.2f cm/s applied\n', speedThresh);
end
fprintf('Analysis window: [%.2f, %.2f] s relative to lick\n', tWinRel(1), tWinRel(2));

%% ---------------- PLOTS: SCATTER + FIT ----------------
x = impulsivity(ok);
y = neuralBias(ok);
s = stimTrial(ok);
v = meanSpeed(ok);

figure; hold on
plot(x(~s), y(~s), 'k.', 'MarkerSize', 14);
plot(x(s),  y(s),  'r.', 'MarkerSize', 14);
xlabel('Impulsivity = toneGain(target) - lickLoc(chosen)');
ylabel('Neural bias (Ahead - Behind), cm·prob');
title('Impulsivity vs neural prospective bias (posterior-based)');
grid on
legend({'Baseline','Stim'}, 'Location','best');

%% ===== Separate fits by condition + Spearman per condition =====
isBase = ~s;   % s = stimTrial(ok) from your script
isStim =  s;

% --- Spearman correlations (per condition) ---
[rB,pB] = corr(x(isBase), y(isBase), 'Type','Spearman', 'Rows','complete');
[rS,pS] = corr(x(isStim), y(isStim), 'Type','Spearman', 'Rows','complete');
fprintf('Spearman baseline: rho = %.3f, p = %.3g (n=%d)\n', rB, pB, nnz(isBase));
fprintf('Spearman stim:     rho = %.3f, p = %.3g (n=%d)\n', rS, pS, nnz(isStim));

% --- Fit separate lines (robust if available) ---
xx = linspace(min(x), max(x), 200);

if exist('robustfit','file') == 2
    bB = robustfit(x(isBase), y(isBase));   % [intercept; slope]
    bS = robustfit(x(isStim), y(isStim));
    yB = bB(1) + bB(2)*xx;
    yS = bS(1) + bS(2)*xx;

    plot(xx, yB, 'k-', 'LineWidth', 2);
    plot(xx, yS, 'r-', 'LineWidth', 2);

    fprintf('Robust fit baseline: intercept = %.3f, slope = %.3f\n', bB(1), bB(2));
    fprintf('Robust fit stim:     intercept = %.3f, slope = %.3f\n', bS(1), bS(2));
else
    pBfit = polyfit(x(isBase), y(isBase), 1);
    pSfit = polyfit(x(isStim), y(isStim), 1);
    yB = polyval(pBfit, xx);
    yS = polyval(pSfit, xx);

    plot(xx, yB, 'k-', 'LineWidth', 2);
    plot(xx, yS, 'r-', 'LineWidth', 2);

    fprintf('LS fit baseline: intercept = %.3f, slope = %.3f\n', pBfit(2), pBfit(1));
    fprintf('LS fit stim:     intercept = %.3f, slope = %.3f\n', pSfit(2), pSfit(1));
end

legend({'Baseline','Stim','All-data fit','Base fit','Stim fit'}, 'Location','best');

% --- Test slope/intercept difference with interaction model (OLS) ---
% Model: y = b0 + b1*x + b2*stim + b3*(x*stim)
% b2 = intercept shift (stim vs base at x=0)
% b3 = slope difference
X = [ones(numel(x),1), x, double(isStim), x.*double(isStim)];
b = X \ y;
yhat = X*b;
res  = y - yhat;
n = numel(y);
p = size(X,2);
s2 = sum(res.^2) / (n - p);
covb = s2 * inv(X'*X);
se = sqrt(diag(covb));
tval = b ./ se;
pval = 2*tcdf(-abs(tval), n-p);

fprintf('\nInteraction model (OLS): y = b0 + b1*x + b2*stim + b3*x*stim\n');
fprintf('b0 (base intercept): %.3f (p=%.3g)\n', b(1), pval(1));
fprintf('b1 (base slope):     %.3f (p=%.3g)\n', b(2), pval(2));
fprintf('b2 (stim intercept shift at x=0): %.3f (p=%.3g)\n', b(3), pval(3));
fprintf('b3 (stim slope diff):             %.3f (p=%.3g)\n\n', b(4), pval(4));


% Robust fit if available
if exist('robustfit','file') == 2
    X = [ones(numel(x),1), x];
    b = robustfit(x, y);                  % b(1)=intercept, b(2)=slope
    xx = linspace(min(x), max(x), 200);
    yy = b(1) + b(2)*xx;
    plot(xx, yy, 'b-', 'LineWidth', 2);

    % p-value for slope (from robustfit stats)
    [~, stats] = robustfit(x, y);
    pSlope = stats.p(2);
    fprintf('robustfit slope = %.4g, p = %.3g\n', b(2), pSlope);
else
    % fallback: least squares
    p = polyfit(x, y, 1);
    xx = linspace(min(x), max(x), 200);
    yy = polyval(p, xx);
    plot(xx, yy, 'b-', 'LineWidth', 2);
    fprintf('polyfit slope = %.4g (no p-value; robustfit not available)\n', p(1));
end

% Optional: partial correlation controlling for speed
if doPartialCorr && exist('partialcorr','file') == 2
    [rpar, ppar] = partialcorr(x, y, v, 'Rows','complete', 'Type','Spearman');
    fprintf('Partial Spearman corr (control speed): r = %.3f, p = %.3g\n', rpar, ppar);
elseif doPartialCorr
    fprintf('partialcorr not available (Stats toolbox). Skipping partial correlation.\n');
end

end

%% ================= HELPERS =================

function [posterior_pos, post_time, post_pos] = loadDecodingData(decodingPath, sessRel, decodingName)
file_nc = fullfile(decodingPath, sessRel, decodingName);
posterior_pos = ncread(file_nc, 'y_position');
post_time     = ncread(file_nc, 'time');
post_pos      = ncread(file_nc, 'y_position_value');
end

function [i0, i1] = decoderIndicesForWindow(post_time, tWin)
% tWin = [t0 t1]
[~, i0] = min(abs(post_time - tWin(1)));
[~, i1] = min(abs(post_time - tWin(2)));
i0 = max(1, i0);
i1 = min(numel(post_time), i1);
if i1 <= i0
    i0 = []; i1 = [];
end
end

function imp = computeImpulsivity(behavTrials, tr)
% Impulsivity = toneGain(target port) - lickLoc(chosen port)
% We implement this robustly:
% - If toneGain is a vector: use toneGain(tr)
% - If toneGain is a matrix: try to index with a target-port field; else error
% - lickLoc assumed scalar per trial; if matrix, try chosen-port field similarly

% --- toneGain(target) ---
tg = behavTrials.toneGain;
if isvector(tg)
    toneGainTarget = tg(tr);
else
    % matrix: need target port index
    targetPort = pickPortIndex(behavTrials, tr, {'targetPort','tonePort','toneID','tone','target'});
    toneGainTarget = tg(tr, targetPort);
end

% --- lickLoc(chosen) ---
ll = behavTrials.lickLoc;
if isvector(ll)
    lickLocChosen = ll(tr);
else
    chosenPort = pickPortIndex(behavTrials, tr, {'chosenPort','lickPort','choice','response','chosen'});
    lickLocChosen = ll(tr, chosenPort);
end

imp = toneGainTarget - lickLocChosen;
end

function port = pickPortIndex(behavTrials, tr, candidates)
% Tries to find a per-trial port index field.
port = [];
for i = 1:numel(candidates)
    f = candidates{i};
    if isfield(behavTrials, f)
        val = behavTrials.(f);
        if isvector(val)
            port = val(tr);
        else
            % if it's a matrix with one column that matches trials, try take first col
            if size(val,1) >= tr
                port = val(tr,1);
            end
        end
        if ~isempty(port) && isfinite(port)
            port = round(port);
            return
        end
    end
end
error(['toneGain/lickLoc is a matrix but I could not find a port index field. ' ...
       'Add a field name to candidates in pickPortIndex() that matches your behavTrials structure.']);
end

function bias = computePosteriorBias(Pseg, post_pos, trueY)
% Computes mean bias across time bins:
% bias_t = sum_{yrel>0} p*yrel - sum_{yrel<0} p*|yrel|
post_pos = post_pos(:);
trueY    = trueY(:);

nBins = size(Pseg,2);
b = nan(nBins,1);

for t = 1:nBins
    yrel = post_pos - trueY(t);
    p = Pseg(:,t);

    if ~isfinite(trueY(t)) || ~any(isfinite(p)), continue; end

    ia = (yrel > 0) & isfinite(yrel) & isfinite(p);
    ib = (yrel < 0) & isfinite(yrel) & isfinite(p);

    ahead  = sum(p(ia) .* yrel(ia));
    behind = sum(p(ib) .* abs(yrel(ib)));

    b(t) = ahead - behind;
end

bias = mean(b, 'omitnan');
end
