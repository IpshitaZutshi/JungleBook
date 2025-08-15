function plotTurnPointsSuperimposed(tracking, behavTrials)

% Define trial types
trialTypes = struct( ...
    'bL', find(behavTrials.stim'==0 & behavTrials.choice==1), ...
    'sL', find(behavTrials.stim'==1 & behavTrials.choice==1), ...
    'bR', find(behavTrials.stim'==0 & behavTrials.choice==0), ...
    'sR', find(behavTrials.stim'==1 & behavTrials.choice==0) ...
);

fields = fieldnames(trialTypes);
colors = struct('bL', [0.3 0.8 0.3], 'sL', [0 0.5 0], ...
                'bR', [0.3 0.3 1], 'sR', [0 0 0.5]);

figure;
for f = 1:length(fields)
    field = fields{f};
    trialIdx = trialTypes.(field);
    isRight = contains(field, 'R');

    subplot(2,2,f); hold on;
    title(sprintf('Turn points for %s (up to choice)', field), 'FontWeight', 'bold');
    xlabel('X'); ylabel('Y');

    turnXs = [];
    turnYs = [];

    for i = 1:length(trialIdx)
        tIdx = trialIdx(i);
        tStart = behavTrials.timestamps(tIdx,1);
        tEnd = behavTrials.choiceTS(tIdx);
        idx = InIntervals(tracking.timestamps, [tStart, tEnd]);
        x = tracking.position.x(idx);
        y = tracking.position.y(idx);

        % Remove nans
        valid = ~(isnan(x) | isnan(y));
        x = x(valid);
        y = y(valid);

        % Skip if not enough points
        if numel(x) < 5, continue; end

        % Compute angle between vectors (smoothed diff)
        dx = smooth(diff(x),5);
        dy = smooth(diff(y),5);
        angles = atan2d(dy, dx);
        dAngles = abs([0; diff(angles)]);

        % Only consider portion before committing to turn
        if isRight
            idxRestrict = find(x < 100 & y > 90);  % restrict to junction
        else
            idxRestrict = find(x > 90 & y > 90);   % restrict to junction
        end
        dAnglesRestrict = dAngles(idxRestrict);
        if isempty(dAnglesRestrict), continue; end

        [~, maxTurnIdxLocal] = max(dAnglesRestrict);
        maxTurnIdx = idxRestrict(maxTurnIdxLocal);

        % Store and plot
        plot(x, y, 'Color', [0.5 0.5 0.5 0.2]);
        plot(x(maxTurnIdx), y(maxTurnIdx), 'ro');
        turnXs(end+1) = x(maxTurnIdx);
        turnYs(end+1) = y(maxTurnIdx);
    end

    if ~isempty(turnXs)
        meanX = mean(turnXs);
        meanY = mean(turnYs);
        plot(meanX, meanY, 'ko', 'MarkerFaceColor', 'k');
        text(meanX+1, meanY, 'Mean turn', 'FontSize', 8);
    end

    xlim([min(tracking.position.x)-2, max(tracking.position.x)+2]);
    ylim([min(tracking.position.y)-2, max(tracking.position.y)+2]);
end

end
