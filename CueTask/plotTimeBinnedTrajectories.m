function plotTimeBinnedTrajectories(tracking, behavTrials)

% Parameters
nBins = 50;  % Number of time bins to divide each trajectory
alpha = 0.05; % Significance level

% Trial indices
bL = find(behavTrials.stim' == 0 & behavTrials.choice == 1);
bR = find(behavTrials.stim' == 0 & behavTrials.choice == 0);
sL = find(behavTrials.stim' == 1 & behavTrials.choice == 1);
sR = find(behavTrials.stim' == 1 & behavTrials.choice == 0);

% Initialize
[bL_X, bL_Y] = getTimeBinnedXY(bL, tracking, behavTrials, nBins);
[bR_X, bR_Y] = getTimeBinnedXY(bR, tracking, behavTrials, nBins);
[sL_X, sL_Y] = getTimeBinnedXY(sL, tracking, behavTrials, nBins);
[sR_X, sR_Y] = getTimeBinnedXY(sR, tracking, behavTrials, nBins);

% Plot
figure;
titles = {'Left Trials: Baseline vs Stim', 'Right Trials: Baseline vs Stim'};
conds = {{bL_X, bL_Y, sL_X, sL_Y}, {bR_X, bR_Y, sR_X, sR_Y}};

for i = 1:2
    subplot(1,2,i); hold on
    BLX = conds{i}{1}; BLY = conds{i}{2};
    STX = conds{i}{3}; STY = conds{i}{4};

    % Plot raw trajectories
    for tr = 1:size(BLX,2)
        plot(BLX(:,tr), BLY(:,tr), 'Color', [0 1 0 0.2]);  % green
    end
    for tr = 1:size(STX,2)
        plot(STX(:,tr), STY(:,tr), 'Color', [0 0 1 0.2]);  % blue
    end

    % Plot average
    meanBLX = nanmean(BLX, 2);
    meanBLY = nanmean(BLY, 2);
    meanSTX = nanmean(STX, 2);
    meanSTY = nanmean(STY, 2);

    plot(meanBLX, meanBLY, 'Color', [0.2 0.6 0.2], 'LineWidth', 2);  % dark green
    plot(meanSTX, meanSTY, 'Color', [0 0 0.6], 'LineWidth', 2);      % dark blue

    % Divergence point
    p = nan(nBins,1);
    for t = 1:nBins
        [~, p(t)] = ttest2(BLX(t,:), STX(t,:), 'Alpha', alpha);
    end
    sigIdx = find(p < alpha, 1, 'first');
    if ~isempty(sigIdx)
        scatter(mean([meanBLX(sigIdx), meanSTX(sigIdx)]), ...
                mean([meanBLY(sigIdx), meanSTY(sigIdx)]), ...
                60, 'r', 'filled');
    end

    title(titles{i});
    xlabel('X'); ylabel('Y');
    xlim([min(tracking.position.x)-5, max(tracking.position.x)+5]);
    ylim([min(tracking.position.y)-5, max(tracking.position.y)+5]);
end

end

function [xBinned, yBinned] = getTimeBinnedXY(trialIdx, tracking, behavTrials, nBins)
    xBinned = nan(nBins, length(trialIdx));
    yBinned = nan(nBins, length(trialIdx));
    for i = 1:length(trialIdx)
        tStart = behavTrials.timestamps(trialIdx(i), 1);
        tEnd = behavTrials.choiceTS(trialIdx(i));
        idx = InIntervals(tracking.timestamps, [tStart tEnd]);

        if sum(idx) < nBins
            continue  % skip if too few points
        end

        x = tracking.position.x(idx);
        y = tracking.position.y(idx);

        % Bin in time
        xResampled = resample(x, nBins, length(x));
        yResampled = resample(y, nBins, length(y));

        xBinned(:,i) = xResampled(:);
        yBinned(:,i) = yResampled(:);
    end
end
