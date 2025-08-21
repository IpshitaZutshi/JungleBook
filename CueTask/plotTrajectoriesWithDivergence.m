function plotTrajectoriesWithDivergence(tracking, behavTrials)

% Preprocess Y to start from zero
tracking.position.y = tracking.position.y - min(tracking.position.y);

% Trial indices
bL = find(behavTrials.stim'==0 & behavTrials.choice == 1); % baseline left
bR = find(behavTrials.stim'==0 & behavTrials.choice == 0); % baseline right
sL = find(behavTrials.stim'==1 & behavTrials.choice == 1); % stim left
sR = find(behavTrials.stim'==1 & behavTrials.choice == 0); % stim right

% Binning along Y
yBins = 96:-1:20;
dy = 1;

% Initialize arrays
bL_X = nan(length(yBins), length(bL)); bR_X = nan(length(yBins), length(bR));
sL_X = nan(length(yBins), length(sL)); sR_X = nan(length(yBins), length(sR));

figure;
titles = {'Left Trials: Baseline vs Stim', 'Right Trials: Baseline vs Stim'};
xlimVals = [88 98];
ylimVals = [70 100];

for trialType = 1:2
    subplot(1,2,trialType); hold on

    if trialType == 1  % LEFT trials
        [bX, bColor, stimIdx] = deal(bL_X, [0 1 0], sL);
        baseIdx = bL;
    else               % RIGHT trials
        [bX, bColor, stimIdx] = deal(bR_X, [0 1 0], sR);
        baseIdx = bR;
    end

    % Plot baseline
    for i = 1:length(baseIdx)
        idx = InIntervals(tracking.timestamps, ...
              [behavTrials.timestamps(baseIdx(i),1), behavTrials.choiceTS(baseIdx(i))]);
        x = tracking.position.x(idx); y = tracking.position.y(idx);
        plot(x, y, 'Color', [bColor 0.2])  % low opacity green
        [bX(:,i)] = binXoverY(x, y, yBins, dy);
    end

    % Plot stim
    stimX = nan(length(yBins), length(stimIdx));
    for i = 1:length(stimIdx)
        idx = InIntervals(tracking.timestamps, ...
              [behavTrials.timestamps(stimIdx(i),1), behavTrials.choiceTS(stimIdx(i))]);
        x = tracking.position.x(idx); y = tracking.position.y(idx);
        plot(x, y, 'Color', [0 0 1 0.2])  % low opacity blue
        stimX(:,i) = binXoverY(x, y, yBins, dy);
    end

    % Mean trajectories
    meanBase = nanmedian(bX, 2);
    meanStim = nanmedian(stimX, 2);
    plot(meanBase, yBins, 'Color', [0.2 0.6 0.2], 'LineWidth', 2);  % dark green
    plot(meanStim, yBins, 'Color', [0 0 0.6], 'LineWidth', 2);      % dark blue

    % Statistical divergence
    p = nan(size(yBins));
    for j = 1:length(yBins)
        [~, p(j)] = ttest2(bX(j,:), stimX(j,:));
    end
    sigBin = find(p < 0.05, 1, 'first');
    if ~isempty(sigBin)
        scatter(mean([meanBase(sigBin), meanStim(sigBin)]), ...
                yBins(sigBin), 60, 'r', 'filled');
    end

    title(titles{trialType});
    xlim([90 97])
    xlabel('X'); ylabel('Y');
    xlim(xlimVals); ylim(ylimVals);
end

end

function xBinned = binXoverY(x, y, yBins, dy)
    xBinned = nan(length(yBins),1);
    
    validIdx = y >= yBins(end) & y <= yBins(1) & x >= 91 & x <= 96;
    x = x(validIdx);
    y = y(validIdx);
    
    for j = 1:length(yBins)
        binIdx = abs(y - yBins(j)) < dy;
        if sum(binIdx) >= 3  % Require at least 3 samples
            xBinned(j) = nanmedian(x(binIdx));
        end
    end
end