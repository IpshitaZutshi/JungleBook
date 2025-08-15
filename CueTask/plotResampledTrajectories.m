function plotResampledTrajectories(tracking, behavTrials)

% Settings
nPoints = 100;    % Number of resampled points
opacity = 0.2;    % Line opacity for individual trials

% Get trial indices
bL = find(behavTrials.stim'==0 & behavTrials.choice == 1);  % Baseline Left
sL = find(behavTrials.stim'==1 & behavTrials.choice == 1);  % Stim Left
bR = find(behavTrials.stim'==0 & behavTrials.choice == 0);  % Baseline Right
sR = find(behavTrials.stim'==1 & behavTrials.choice == 0);  % Stim Right

groups = {bL, sL, bR, sR};
labels = {'Baseline Left', 'Stim Left', 'Baseline Right', 'Stim Right'};
colors = {[0 0.6 0], [0 1 0], [0 0 0.8], [0 0.5 1]};  % Dark green, light green, dark blue, light blue

figure; hold on;

for g = 1:4
    trials = groups{g};
    resampledX = nan(nPoints, length(trials));
    resampledY = nan(nPoints, length(trials));

    for i = 1:length(trials)
        idx = InIntervals(tracking.timestamps, ...
               [behavTrials.timestamps(trials(i),1), behavTrials.choiceTS(trials(i))]);
        x = tracking.position.x(idx);
        y = tracking.position.y(idx);

        if length(x) < 2, continue; end

        % Plot raw trajectory
        plot(x, y, 'Color', [colors{g}, opacity]);

        % Path length
        d = [0; cumsum(sqrt(diff(x).^2 + diff(y).^2))];
        [d, keepIdx] = unique(d);  % Ensure monotonic increasing
        x = x(keepIdx);
        y = y(keepIdx);
        if length(d) < 2, continue; end

        % Resample
        dResampled = linspace(0, d(end), nPoints);
        xq = interp1(d, x, dResampled, 'linear', 'extrap');
        yq = interp1(d, y, dResampled, 'linear', 'extrap');

        resampledX(:,i) = xq;
        resampledY(:,i) = yq;
    end

    % Plot mean trajectory
    plot(nanmean(resampledX,2), nanmean(resampledY,2), ...
         'Color', colors{g}, 'LineWidth', 2, 'DisplayName', labels{g});
end

% Formatting
xlabel('X'); ylabel('Y');
% xlim([89.5 97.5]);
% ylim([90 120]);
% xlim([89.5 97.5]);
% ylim([90 120]);
title('Trajectories: Baseline vs Stim, Left vs Right');
axis square;

end
