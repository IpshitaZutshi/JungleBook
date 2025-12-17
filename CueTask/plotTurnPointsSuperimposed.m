%% Plot bias-arm turn points and compare Y-values across conditions
% Auto-run from folder

basepath = pwd;

% --- Load files ---
noseFile = dir(fullfile(basepath, '*Tracking.Body.mat'));
behFile  = dir(fullfile(basepath, '*TrialBehavior.Events.mat'));
assert(~isempty(noseFile) && ~isempty(behFile), 'Missing files in %s', basepath);
load(fullfile(noseFile.folder, noseFile.name), 'tracking');
load(fullfile(behFile.folder, behFile.name), 'behavTrials');

stim   = behavTrials.stim(:);
choice = behavTrials.choice(:); % 1 = left, 0 = right
corr   = behavTrials.correct(:);

idxBC = find(stim==0 & corr==1);
idxBE = find(stim==0 & corr==0);
idxST = find(stim==1);

conds  = {idxBC, idxBE, idxST};
colors = {[0.20 0.60 0.80], [0.30 0.75 0.40], [0.85 0.40 0.40]};
labels = {'BC','BE','Stim'};

% --- Parameters ---
yMinSearch = 100;   % skip base
bandX      = 4;    % center-arm band
bandXmin = 1;

% --- Estimate center arm X from BC trials ---
% --- Estimate center arm X as mean X at Y closest to 90 (BC trials) ---
yRef = 90;
centerX_vals = nan(numel(idxBC),1);

for i = 1:numel(idxBC)
    tr = idxBC(i);

    % use center-arm window: cue to choice
    idx_bc = InIntervals(tracking.timestamps, ...
        [behavTrials.timestamps(tr,1)+2, behavTrials.choiceTS(tr)]);

    x = tracking.position.x(idx_bc);
    y = tracking.position.y(idx_bc);

    if ~isempty(x) && ~isempty(y) && any(isfinite(y))
        [~, j] = min(abs(y - yRef));   % index of Y closest to 90
        centerX_vals(i) = x(j);
    end
end

centerX = mean(centerX_vals, 'omitnan');

% --- Determine bias from STIM trials ---
nLeft_ST  = sum(choice(idxST) == 1);
nRight_ST = sum(choice(idxST) == 0);
if nLeft_ST > nRight_ST
    biasDir = 1; % left
else
    biasDir = 0; % right
end
fprintf('Bias direction from STIM trials: %s\n', ternary(biasDir==1, 'Left', 'Right'));

% --- Collect turn Y values ---
turnY_all = cell(1,3);

figure;
set(gcf,'Renderer','painters')
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% === Subplot 1: All trials with bias vs non-bias turn points ===
ax1 = nexttile; hold(ax1,'on');

for c = 1:3
    ids = conds{c};
    for k = 1:numel(ids)
        tr  = ids(k);
        idx = InIntervals(tracking.timestamps, ...
            [behavTrials.timestamps(tr,1)+2, behavTrials.choiceTS(tr)]);
        x = tracking.position.x(idx);
        y = tracking.position.y(idx);
        t = tracking.timestamps(idx);

        if numel(y) >= 3
            validIdx = (y >= yMinSearch) & (abs(x - centerX) <= bandX) & (abs(x - centerX) >= bandXmin);
            if nnz(validIdx) >= 3
                % Smooth valid segment
                x_valid = smoothdata(x(validIdx), 'movmedian', 2);
                y_valid = smoothdata(y(validIdx), 'movmedian', 2);
                t_valid = t(validIdx);

                % Manual derivatives
                dx_dt = [diff(x_valid) ./ diff(t_valid); 0];
                dy_dt = [diff(y_valid) ./ diff(t_valid); 0];

                % dx/dy
                dx_dy = dx_dt ./ dy_dt;

                % Max |dx/dy|
                [~, idxMax] = max(abs(dx_dy));

                % Map to original coords
                x_valid_orig = x(validIdx);
                y_valid_orig = y(validIdx);
                xTurn = x_valid_orig(idxMax);
                yTurn = y_valid_orig(idxMax);

                % Store bias turn Y
                if choice(tr) == biasDir
                    turnY_all{c}(end+1) = yTurn;
                end

                % Plot trajectory faint
                plot(ax1, x, y, '-', 'Color', [colors{c} 0.1], 'LineWidth', 0.8);

                % Plot turn point
                if choice(tr) == biasDir
                    % Bias arm: bold square
                    plot(ax1, xTurn, yTurn, 's', 'MarkerSize', 6, ...
                        'MarkerFaceColor', colors{c}, 'MarkerEdgeColor', 'k', 'LineWidth', 1.0);
                else
                    % Non-bias arm: faint circle
                    scatter(ax1, xTurn, yTurn, 25, ...
                        'MarkerFaceColor', colors{c} * 0.6 + 0.4, ...
                        'MarkerEdgeColor', 'none', ...
                        'MarkerFaceAlpha', 0.8);
                end
            end
        end
    end
end
xlabel(ax1,'X'); ylabel(ax1,'Y'); axis(ax1,'equal');
ylim([100 118])
xlim([84 102])
title(ax1, sprintf('All Trials (Bias = %s)', ternary(biasDir==1, 'Left', 'Right')));
grid(ax1,'on');

% === Subplot 2: Boxplot + jitter for bias-arm turns ===
ax2 = nexttile; hold(ax2,'on');
allY = [];
allG = [];
for c = 1:3
    allY = [allY, turnY_all{c}];
    allG = [allG, repmat(c, 1, numel(turnY_all{c}))];
end
boxplot(allY, allG, 'Colors', cell2mat(colors'), 'Labels', labels, 'Widths', 0.5);
for c = 1:3
    xj = c + 0.15*(rand(1,numel(turnY_all{c})) - 0.5);
    scatter(xj, turnY_all{c}, 25, 'MarkerFaceColor', colors{c}, ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end
ylabel(ax2,'Turn Y (Bias Arm)');
title(ax2,'Bias Arm Turn Y by Condition');
grid(ax2,'on');

% --- Helper ---
function out = ternary(cond, valTrue, valFalse)
if cond
    out = valTrue;
else
    out = valFalse;
end
end
