function wrapper_strategy_metrics
%% ==================== USER SETTINGS ====================
direc = '\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\Cue task';

sessions = {'T20\Final\T20_241101_171908', ...
    'T20\Final\T20_241102_163911', ...
    'T20\Final\T20_241103_131201',...
    'T20\Final\T20_241111_154144',...
    'T22\Final\T22_241022_150631', ... 
    'T22\Final\T22_241029_145741',...
    'T22\Final\T22_241030_144753',...
    };

% Template/metric parameters
NBins      = 120;     % Y grid
MinRun     = 8;       % consecutive bins needed for a "commit"
llrThr     = 2.5;     % per-bin LLR threshold (use 2.5–4)
bandwidth  = 2;       % smoothing for raincloud violin
minTrialsTemplate = 20; % require this many baseline-correct trials to build template
yFloorUse = false; yFloor = 0; % optional: ignore bins below yFloor
postThr    = 0.85;    % posterior >= 0.80 (or <= 0.20) counts as “decided”
dprimeThr  = 0.6;     % template separability gate
stablePct  = 0.60;    % after commit, must stay same sign ≥80% of remaining bins
usePosteriorHysteresis = true;  % set false to use cumulative-LLR variant below

% Plot settings
clrBC = [0.20 0.60 0.80];  % Baseline Correct
clrBE = [0.30 0.75 0.40];  % Baseline Error
clrST = [0.85 0.40 0.40];  % Stim (correct+error)

%% ==================== ACCUMULATORS ACROSS SESSIONS ====================
commitY_BC = [];   % baseline-correct (for sanity check; leave-one-out)
commitY_BE = [];   % baseline-error
commitY_ST = [];   % stim (pooled correct+error)

indec_BC = []; crossings_BC = []; lateSlope_BC = [];
indec_BE = []; crossings_BE = []; lateSlope_BE = [];
indec_ST = []; crossings_ST = []; lateSlope_ST = [];

%% ==================== LOOP SESSIONS ====================
for s = 1:numel(sessions)
    basepath = fullfile(direc, sessions{s});
    noseFile = dir(fullfile(basepath, '*Tracking.Nose.mat'));
    behFile  = dir(fullfile(basepath, '*TrialBehavior.Events.mat'));
    if isempty(noseFile) || isempty(behFile)
        fprintf('[%d] Missing files: %s\n', s, basepath);
        continue;
    end

    load(fullfile(noseFile.folder, noseFile.name), 'tracking');
    load(fullfile(behFile.folder, behFile.name), 'behavTrials');

    % --------- SELECT TRIAL SETS ---------
    stim  = behavTrials.stim(:);
    corr  = behavTrials.correct(:);
    ch    = behavTrials.choice(:); % 1=L, 0=R

    idxBC = find(stim==0 & corr==1);
    idxBE = find(stim==0 & corr==0);
    idxST = find(stim==1); % stim trials, both correct & error

    if numel(idxBC) < minTrialsTemplate
        fprintf('[%d] Skip: baseline-correct too few (%d)\n', s, numel(idxBC));
        continue;
    end

    % --------- BUILD BASELINE TEMPLATE (from BC only) ---------
    [muL, muR, sdL, sdR, mu0, yG] = buildBaselineTemplate_session(tracking, behavTrials, idxBC, NBins);
    % Optional yFloor
    if yFloorUse
        useBins = yG >= yFloor;
    else
        useBins = true(size(yG));
    end
    muL=muL(useBins); muR=muR(useBins); sdL=sdL(useBins); sdR=sdR(useBins); mu0=mu0(useBins); yG=yG(useBins);

    % --------- METRICS PER CONDITION (single-trial friendly) ---------
    % Baseline-Correct (leave-one-out to avoid circularity if desired)
    [cY, indec, crossN, lateS] = metrics_for_trials(tracking, behavTrials, idxBC, ...
    muL,sdL,muR,sdR,mu0,yG, llrThr, MinRun, postThr, dprimeThr, stablePct, usePosteriorHysteresis);
    commitY_BC   = [commitY_BC;   cY(:)];
    indec_BC     = [indec_BC;     indec(:)];
    crossings_BC = [crossings_BC; crossN(:)];
    lateSlope_BC = [lateSlope_BC; lateS(:)];

    % Baseline-Error
    [cY, indec, crossN, lateS] = metrics_for_trials(tracking, behavTrials, idxBE, ...
        muL,sdL,muR,sdR,mu0,yG, llrThr, MinRun, postThr, dprimeThr, stablePct, usePosteriorHysteresis);
    commitY_BE   = [commitY_BE;   cY(:)];
    indec_BE     = [indec_BE;     indec(:)];
    crossings_BE = [crossings_BE; crossN(:)];
    lateSlope_BE = [lateSlope_BE; lateS(:)];

    % Stimulation (pooled)
    [cY, indec, crossN, lateS] = metrics_for_trials(tracking, behavTrials, idxST, ...
        muL,sdL,muR,sdR,mu0,yG, llrThr, MinRun, postThr, dprimeThr, stablePct, usePosteriorHysteresis);
    commitY_ST   = [commitY_ST;   cY(:)];
    indec_ST     = [indec_ST;     indec(:)];
    crossings_ST = [crossings_ST; crossN(:)];
    lateSlope_ST = [lateSlope_ST; lateS(:)];
end

%% ==================== Commit Y only + simple ANOVA ====================
figure('Color','w'); 
raincloud3({commitY_BC, commitY_BE, commitY_ST}, {'BC','BE','Stim'}, ...
           {clrBC, clrBE, clrST}, bandwidth);
ylabel('Commit Y (LLR)');
title('Commit Y across conditions');

%% ---------- One-way ANOVA + Tukey (ignores NaNs) ----------
y   = [commitY_BC(:); commitY_BE(:); commitY_ST(:)];
grp = [repmat("BC",  numel(commitY_BC),1); ...
       repmat("BE",  numel(commitY_BE),1); ...
       repmat("Stim",numel(commitY_ST),1)];
ok = ~isnan(y);
y = y(ok); grp = grp(ok);

[p, tbl, S] = anova1(y, grp, 'off');         % one-way ANOVA
F  = tbl{2,5}; df1 = tbl{2,3}; df2 = tbl{3,3};
fprintf('ANOVA: F(%d,%d) = %.3f, p = %.4g\n', df1, df2, F, p);

% Pairwise Tukey-Kramer posthoc
C = multcompare(S, 'CType','tukey-kramer', 'Display','off');
labels = {'BC','BE','Stim'};
pairwise = array2table(C, 'VariableNames', ...
    {'group1','group2','lowerCI','diff','upperCI','pAdj'});
pairwise.group1 = categorical(labels(pairwise.group1).');
pairwise.group2 = categorical(labels(pairwise.group2).');
disp('Pairwise (Tukey-Kramer), adjusted p-values:');
disp(pairwise(:,{'group1','group2','diff','pAdj'}));

end

%% ==================== HELPERS ====================

function [muL, muR, sdL, sdR, mu0, yGrid] = buildBaselineTemplate_session(tracking, behavTrials, idxBC, NBins)
% Build L/R templates from baseline-correct trials (idxBC) of THIS session.

% Split BC by choice
L = idxBC(behavTrials.choice(idxBC)==1);
R = idxBC(behavTrials.choice(idxBC)==0);

% Find common Y overlap
mins=[]; maxs=[];
for tr=[L(:);R(:)]'
    idx = InIntervals(tracking.timestamps,[behavTrials.timestamps(tr,1), behavTrials.choiceTS(tr)]);
    y   = tracking.position.y(idx);
    if numel(y)>=5, mins(end+1)=min(y); maxs(end+1)=max(y); end %#ok<AGROW>
end
yLow = max(mins); yHigh = min(maxs);
yGrid = linspace(yLow, yHigh, NBins);

% Interpolate X(Y)
Xl = nan(NBins,numel(L)); Xr = nan(NBins,numel(R));
for i=1:numel(L), [x,y]=xy(tracking,behavTrials,L(i)); Xl(:,i)=interp1(y,x,yGrid,'linear','extrap'); end
for i=1:numel(R), [x,y]=xy(tracking,behavTrials,R(i)); Xr(:,i)=interp1(y,x,yGrid,'linear','extrap'); end

% Means/SDs (smoothed, floored)
muL = movmean(mean(Xl,2,'omitnan'),5);
muR = movmean(mean(Xr,2,'omitnan'),5);
sdL = max(movmean(std(Xl,0,2,'omitnan'),5), 1e-6);
sdR = max(movmean(std(Xr,0,2,'omitnan'),5), 1e-6);
mu0 = (muL+muR)/2;
end

function [commitY, indec, crossN, lateSlope] = metrics_for_trials(tracking, behavTrials, idxList, ...
    muL,sdL,muR,sdR,mu0,yG, llrThr, MinRun, postThr, dprimeThr, stablePct, usePosteriorHysteresis)

% Outputs are one value per trial in idxList:
%   commitY   : earliest Y of sustained commitment (robust to noise)
%   indec     : fraction of Y with weak evidence (|LLR|<1) in separable bins
%   crossN    : number of evidence sign flips (only in separable bins)
%   lateSlope : slope of |LLR| over last 20% of Y

commitY   = nan(numel(idxList),1);
indec     = nan(numel(idxList),1);
crossN    = nan(numel(idxList),1);
lateSlope = nan(numel(idxList),1);

% Enforce column orientation
yG  = yG(:);
muL = muL(:); muR = muR(:);
sdL = sdL(:); sdR = sdR(:);
mu0 = mu0(:);

if isempty(idxList), return; end

for k = 1:numel(idxList)
    tr = idxList(k);

    % ---- interpolate this trial onto template grid ----
    [xq, yq] = xy(tracking, behavTrials, tr);
    xq = interp1(yq, xq, yG, 'linear', 'extrap');
    xq = xq(:);

    % (Optional) tiny smoothing to kill one-bin spikes
    if numel(xq) >= 3
        xq = movmedian(xq, 3);
    end

    % ---- per-bin LLR (Gaussian) ----
    % floor SDs already applied when building template; still guard divisions
    sdLsafe = max(sdL, 1e-6);
    sdRsafe = max(sdR, 1e-6);
    llL = -0.5*((xq - muL)./sdLsafe).^2 - log(sdLsafe);
    llR = -0.5*((xq - muR)./sdRsafe).^2 - log(sdRsafe);
    llr = llL - llR;
    llr(isnan(llr)) = 0;  % avoid NaNs propagating

    % ---- separability gate: d' >= threshold ----
    dprime  = abs(muL - muR) ./ sqrt(0.5*(sdLsafe.^2 + sdRsafe.^2));
    goodBin = dprime >= dprimeThr;

    % If nothing is separable, skip metrics for this trial
    if ~any(goodBin)
        commitY(k)   = NaN;
        indec(k)     = NaN;
        crossN(k)    = NaN;
        lateSlope(k) = NaN;
        continue;
    end

    % ---- robust commit detection ----
    if usePosteriorHysteresis
        % Posterior threshold with MinRun + stability after commit
        post = 1 ./ (1 + exp(-llr)); % P(Left | x,y)
        sig  = ((post >= postThr) | (post <= 1-postThr)) & goodBin;

        y0 = first_run_start_y(sig, yG, MinRun);    % earliest sustained bin
        commitY(k) = NaN;
        if isfinite(y0)
            b0  = find(yG==y0,1,'first');
            sgn = sign(llr(b0) + eps);
            stay = mean(sign(llr(b0:end) + eps) == sgn, 'omitnan'); % stability
            if stay >= stablePct
                commitY(k) = y0;
            else
                % try a later stable run
                sig(1:b0) = false;
                y1 = first_run_start_y(sig, yG, MinRun);
                if isfinite(y1)
                    b1   = find(yG==y1,1,'first');
                    sgn1 = sign(llr(b1) + eps);
                    stay1= mean(sign(llr(b1:end) + eps) == sgn1, 'omitnan');
                    if stay1 >= stablePct, commitY(k) = y1; end
                end
            end
        end
    else
        % Cumulative LLR with enter/leave hysteresis
        enter = 5; leave = 2;           % tweak if desired
        cum = cumsum(llr);
        sig = (abs(cum) >= enter) & goodBin;
        y0  = first_run_start_y(sig, yG, MinRun);
        commitY(k) = NaN;
        if isfinite(y0)
            b0 = find(yG==y0,1,'first');
            if all(abs(cum(b0:end)) >= leave)
                commitY(k) = y0;
            end
        end
    end

    % ---- other metrics ----
    indec(k) = mean( (abs(llr) < 1) & goodBin, 'omitnan');   % weak-evidence fraction in separable bins

    % count sign flips only in separable region
    llr_g = llr;
    llr_g(~goodBin) = NaN;
    crossN(k) = sum(diff(sign(llr_g + eps)) ~= 0, 'omitnan');

    % late slope over last 20% (if sufficient bins)
    last = yG >= quantile(yG, 0.8);
    if nnz(last) >= 2 && all(isfinite(llr(last)))
        p = polyfit(yG(last), abs(llr(last)), 1);
        lateSlope(k) = p(1);
    else
        lateSlope(k) = NaN;
    end
end
end

function [x,y] = xy(tracking,behavTrials,tr)
idx = InIntervals(tracking.timestamps,[behavTrials.timestamps(tr,1), behavTrials.choiceTS(tr)]);
x = tracking.position.x(idx); y = tracking.position.y(idx);
[y,ord] = sort(y(:)); x = x(ord); [y,ui] = unique(y); x = x(ui);
end

function y0 = first_run_start_y(sig,yGrid,MinRun)
d = diff([false; sig(:); false]); st=find(d==1); en=find(d==-1)-1;
len = en-st+1; k=find(len>=MinRun,1,'first');
y0 = NaN; if ~isempty(k), y0 = yGrid(st(k)); end
end

function raincloud3(groups, labels, colors, bw)
% groups: cell{1xK} of vectors
% labels: cellstr 1xK
% colors: cell{1xK} RGB
if nargin<4, bw=2; end
K = numel(groups);
hold on;
for k=1:K
    data = groups{k}; data = data(isfinite(data));
    if isempty(data), continue; end
    [f,xi] = ksdensity(data,'Bandwidth',bw);
    f = f / max(f) * 0.35; % width scale
    x0 = k; 
    fill(x0 + [f, -f(end:-1:1)], [xi, xi(end:-1:1)], colors{k}, ...
        'FaceAlpha', 0.25, 'EdgeColor','none');
    % jittered points
    xs = x0 + 0.14*(rand(size(data))-0.5);
    scatter(xs, data, 28, 'filled', 'MarkerFaceColor', colors{k}, 'MarkerFaceAlpha', 0.8);
    % mean line
    ym = mean(data); plot([x0-0.25, x0+0.25], [ym, ym], '--', 'Color', colors{k});
end
xlim([0.5 K+0.5]); xticks(1:K); xticklabels(labels); grid on; box on;
end
