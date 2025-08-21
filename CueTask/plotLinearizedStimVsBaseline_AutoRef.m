function plotLinearizedStimVsBaseline_AutoRef(tracking, behavTrials)

tracking.position.y = tracking.position.y - min(tracking.position.y);  % normalize Y

% Trial indices
bL = find(behavTrials.stim'==0 & behavTrials.choice == 1);
bR = find(behavTrials.stim'==0 & behavTrials.choice == 0);
sL = find(behavTrials.stim'==1 & behavTrials.choice == 1);
sR = find(behavTrials.stim'==1 & behavTrials.choice == 0);

% All baseline trials (for reference)
allB = [bL; bR];

% Define Y-bins
yBins = linspace(20, 94, 100)';
dy = mean(diff(yBins));

% Extract all baseline X positions into bins for reference
ref_X = nan(length(yBins), length(allB));

for i = 1:length(allB)
    idx = InIntervals(tracking.timestamps, [behavTrials.timestamps(allB(i),1), behavTrials.choiceTS(allB(i))]);
    x = tracking.position.x(idx);
    y = tracking.position.y(idx);
    valid = y >= yBins(1) & y <= yBins(end);
    x = x(valid); y = y(valid);
    for j = 1:length(yBins)
        binIdx = abs(y - yBins(j)) < dy;
        if any(binIdx)
            ref_X(j,i) = nanmedian(x(binIdx));
        end
    end
end

% Reference trajectory = median X at each Y bin
refX_interp = nanmedian(ref_X, 2);

% Helper to linearize a trial
linearize = @(x,y) arrayfun(@(yb) nanmean(x(abs(y - yb) < dy)), yBins);

% Initialize matrices
bL_X = nan(length(yBins), length(bL));
bR_X = nan(length(yBins), length(bR));
sL_X = nan(length(yBins), length(sL));
sR_X = nan(length(yBins), length(sR));

% Compute linearized X trajectories
for i = 1:length(bL)
    idx = InIntervals(tracking.timestamps, [behavTrials.timestamps(bL(i),1), behavTrials.choiceTS(bL(i))]);
    bL_X(:,i) = linearize(tracking.position.x(idx), tracking.position.y(idx));
end
for i = 1:length(bR)
    idx = InIntervals(tracking.timestamps, [behavTrials.timestamps(bR(i),1), behavTrials.choiceTS(bR(i))]);
    bR_X(:,i) = linearize(tracking.position.x(idx), tracking.position.y(idx));
end
for i = 1:length(sL)
    idx = InIntervals(tracking.timestamps, [behavTrials.timestamps(sL(i),1), behavTrials.choiceTS(sL(i))]);
    sL_X(:,i) = linearize(tracking.position.x(idx), tracking.position.y(idx));
end
for i = 1:length(sR)
    idx = InIntervals(tracking.timestamps, [behavTrials.timestamps(sR(i),1), behavTrials.choiceTS(sR(i))]);
    sR_X(:,i) = linearize(tracking.position.x(idx), tracking.position.y(idx));
end

%% Plotting

figure;

% LEFT TRIALS (stim vs baseline)
subplot(1,2,1); hold on;
plot(bL_X, yBins, 'Color', [0 0.6 0 0.1])
plot(sL_X, yBins, 'Color', [0 0 1 0.1])
plot(nanmedian(bL_X,2), yBins, 'g', 'LineWidth', 2)
plot(nanmedian(sL_X,2), yBins, 'b', 'LineWidth', 2)
xlabel('Linearized X'); ylabel('Y'); title('Left Trials: Baseline vs Stim');

% Divergence
p = arrayfun(@(j) ttest2(bL_X(j,:), sL_X(j,:)), 1:length(yBins));
divIdx = find(p < 0.05, 1, 'first');
if ~isempty(divIdx)
    scatter(mean([nanmedian(bL_X(divIdx,:)), nanmedian(sL_X(divIdx,:))]), yBins(divIdx), ...
        50, 'r', 'filled');
end

% RIGHT TRIALS (stim vs baseline)
subplot(1,2,2); hold on;
plot(bR_X, yBins, 'Color', [0 0.6 0 0.1])
plot(sR_X, yBins, 'Color', [0 0 1 0.1])
plot(nanmedian(bR_X,2), yBins, 'g', 'LineWidth', 2)
plot(nanmedian(sR_X,2), yBins, 'b', 'LineWidth', 2)
xlabel('Linearized X'); ylabel('Y'); title('Right Trials: Baseline vs Stim');

% Divergence
p = arrayfun(@(j) ttest2(bR_X(j,:), sR_X(j,:)), 1:length(yBins));
divIdx = find(p < 0.05, 1, 'first');
if ~isempty(divIdx)
    scatter(mean([nanmedian(bR_X(divIdx,:)), nanmedian(sR_X(divIdx,:))]), yBins(divIdx), ...
        50, 'r', 'filled');
end

end
