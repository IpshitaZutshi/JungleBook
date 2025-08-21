%% Compare bias-arm turn Y across trials (BC, BE, Stim)
% Loops over given session folders and collects bias-arm turn Y for each trial.
% Does stats on TRIALS (pooled) rather than per-session means.

clear; clc;

% -------------------- Configure your sessions here --------------------
sessionDirs = {'C:\Users\Ipshita\NYU Langone Health Dropbox\Ipshita Zutshi\Cue Task\Data\T20\Final\T20_241101_171908',...
    'C:\Users\Ipshita\NYU Langone Health Dropbox\Ipshita Zutshi\Cue Task\Data\T20\Final\T20_241102_163911',...
    'C:\Users\Ipshita\NYU Langone Health Dropbox\Ipshita Zutshi\Cue Task\Data\T22\Final\T22_241022_150631',...
    'C:\Users\Ipshita\NYU Langone Health Dropbox\Ipshita Zutshi\Cue Task\Data\T22\Final\T22_241029_145741'};  % default: just the current folder

% ----------------------------------------------------------------------

% Parameters
yMinSearch = 100;   % skip base
bandX      = 6;     % center-arm band half-width
yRef       = 90;    % centerX reference Y
smoothWin  = 2;     % smoothing window for movmedian

% Colors / labels
labels = {'BC','BE','Stim'};
colors = {[0.20 0.60 0.80], [0.30 0.75 0.40], [0.85 0.40 0.40]};

% Trial-level accumulators
turnY_all = {[], [], []};  % BC, BE, Stim
trialGroups = [];  % 1=BC, 2=BE, 3=Stim
trialVals   = [];

for s = 1:numel(sessionDirs)
    basepath = sessionDirs{s};
    try
        % --- Load files ---
        noseFile = dir(fullfile(basepath, '*Tracking.Body.mat'));
        behFile  = dir(fullfile(basepath, '*TrialBehavior.Events.mat'));
        assert(~isempty(noseFile) && ~isempty(behFile), 'Missing files in %s', basepath);
        load(fullfile(noseFile.folder, noseFile.name), 'tracking');
        load(fullfile(behFile.folder, behFile.name), 'behavTrials');

        % convenience handles
        stim   = behavTrials.stim(:);
        choice = behavTrials.choice(:);   % 1 = left, 0 = right
        corr   = behavTrials.correct(:);

        idxBC = find(stim==0 & corr==1);
        idxBE = find(stim==0 & corr==0);
        idxST = find(stim==1);
        conds = {idxBC, idxBE, idxST};

        % --- centerX from BC at Y closest to 90 ---
        centerX_vals = nan(numel(idxBC),1);
        for i = 1:numel(idxBC)
            tr = idxBC(i);
            idx_bc = InIntervals(tracking.timestamps, ...
                [behavTrials.timestamps(tr,1)+2, behavTrials.choiceTS(tr)]);
            x = tracking.position.x(idx_bc);
            y = tracking.position.y(idx_bc);
            if ~isempty(x) && ~isempty(y) && any(isfinite(y))
                [~, j] = min(abs(y - yRef));
                centerX_vals(i) = x(j);
            end
        end
        centerX = mean(centerX_vals, 'omitnan');

        % --- bias from STIM trials ---
        nLeft_ST  = sum(choice(idxST) == 1);
        nRight_ST = sum(choice(idxST) == 0);
        biasDir = (nLeft_ST > nRight_ST);   % 1=left, 0=right

        % --- collect bias-arm trial Ys ---
        for c = 1:3
            ids = conds{c};
            ids = ids(choice(ids) == biasDir); % only bias-arm

            for k = 1:numel(ids)
                tr  = ids(k);
                idx = InIntervals(tracking.timestamps, ...
                    [behavTrials.timestamps(tr,1)+2, behavTrials.choiceTS(tr)]);
                x = tracking.position.x(idx);
                y = tracking.position.y(idx);
                t = tracking.timestamps(idx);

                if numel(y) < 3, continue; end
                validIdx = (y >= yMinSearch) & (abs(x - centerX) <= bandX);
                if nnz(validIdx) < 3, continue; end

                % smooth
                x_valid = smoothdata(x(validIdx), 'movmedian', smoothWin);
                y_valid = smoothdata(y(validIdx), 'movmedian', smoothWin);
                t_valid = t(validIdx);

                dx_dt = [diff(x_valid) ./ diff(t_valid); 0];
                dy_dt = [diff(y_valid) ./ diff(t_valid); 0];
                dx_dy = dx_dt ./ dy_dt;

                [~, idxMax] = max(abs(dx_dy));
                y_valid_orig = y(validIdx);
                yTurn = y_valid_orig(idxMax);

                % store
                turnY_all{c}(end+1) = yTurn; %#ok<SAGROW>
                trialGroups(end+1,1) = c; %#ok<SAGROW>
                trialVals(end+1,1)   = yTurn; %#ok<SAGROW>
            end
        end

        fprintf('[%d/%d] %s: collected %d/%d/%d bias trials\n', ...
            s, numel(sessionDirs), basepath, ...
            numel(turnY_all{1}), numel(turnY_all{2}), numel(turnY_all{3}));

    catch ME
        warning('Session failed: %s\n  -> %s', sessionDirs{s}, ME.message);
    end
end

%% -------------------- Stats --------------------
nCounts = cellfun(@numel, turnY_all); % number of trials in each condition

[pKW, tbl, statsKW] = kruskalwallis(trialVals, trialGroups, 'off');
chi2KW = tbl{2,5};
dfKW   = tbl{2,3};

posthocStr = '';
if pKW < 0.05
    c = multcompare(statsKW, 'CType','dunn-sidak', 'Display','off');
    for r = 1:size(c,1)
        posthocStr = sprintf('%s%s vs %s: p=%.4g\n', posthocStr, ...
            labels{c(r,1)}, labels{c(r,2)}, c(r,6));
    end
end

statsStr = sprintf(['Kruskal-Wallis: \\chi^2(%d) = %.3f, p = %.4g\n' ...
                    'n_{BC}=%d, n_{BE}=%d, n_{Stim}=%d\n%s'], ...
    dfKW, chi2KW, pKW, nCounts(1), nCounts(2), nCounts(3), posthocStr);

%% -------------------- Plot --------------------
figure('Name','Bias-arm turn Y (trial level)','Position',[100 100 800 500]);
set(gcf,'Renderer','painters')
hold on;
boxplot(trialVals, trialGroups, 'Labels', labels, 'Widths', 0.55, 'Symbol','');
for c = 1:3
    mask = (trialGroups==c);
    xj = c + 0.15*(rand(sum(mask),1)-0.5);
    scatter(xj, trialVals(mask), 20, 'MarkerFaceColor', colors{c}, ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end
ylabel('Turn Y (bias-arm trials)');
title('Bias-arm turn Y (trial-level, pooled across sessions)');
grid on;

% annotation of stats
yl = ylim;
text(2, yl(1) - 0.1*(yl(2)-yl(1)), statsStr, ...
    'Units','data','FontSize',10, 'VerticalAlignment','top');

ylim([yl(1) - 0.15*(yl(2)-yl(1)), yl(2)]);