
basepath = '\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\Cue task\T20\Final\T20_241117_132754';

noseFile = dir(fullfile(basepath, '*Tracking.Nose.mat'));
behFile  = dir(fullfile(basepath, '*TrialBehavior.Events.mat'));
if isempty(noseFile) || isempty(behFile)
    warning('Missing tracking or behavior file in %s', basepath);
    return
end

load(fullfile(noseFile.folder, noseFile.name), 'tracking');
load(fullfile(behFile.folder, behFile.name), 'behavTrials');

fig1 = figure('Color','w');
ax1 =  subplot(1,2,1);
hold on
plotResampledTrajectories(tracking, behavTrials)
[y_diverge, stats] = findChoiceDivergenceY(tracking, behavTrials, 'TrialFilter', behavTrials.correct==1);

%% Batch Y-divergence computation across sessions
direc = '\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\Cue task';

sessions = {
    'T20\Final\T20_241030_153043', ...
    'T20\Final\T20_241031_173852', ...
    'T20\Final\T20_241101_171908', ...
    'T20\Final\T20_241102_163911', ...
    'T20\Final\T20_241103_131201', ...
    'T20\Final\T20_241107_145457', ...
    'T20\Final\T20_241111_154144', ...
    'T20\Final\T20_241114_133417', ...
    'T20\Final\T20_241115_172023', ...
    'T20\Final\T20_241117_132754', ...
    'T20\Final\T20_241118_144118', ...
    'T20\Final\T20_241119_124533', ...
    'T20\Final\T20_241213_152020', ...
    'T22\Final\T22_241022_150631', ... 
    'T22\Final\T22_241023_175232', ...
    'T22\Final\T22_241024_153843', ...
    'T22\Final\T22_241029_145741', ...
    'T22\Final\T22_241030_144753', ...
    'T22\Final\T22_241102_173932', ...
    'T22\Final\T22_241103_123921', ...
    'T22\Final\T22_241104_150209', ...
    'T22\Final\T22_241113_150043', ...
    'T22\Final\T22_241114_143845', ...
    'T22\Final\T22_241119_144043', ...
    'T22\Final\T22_241209_173821', ...
    'T22\Final\T22_241209_182341', ... 
    'T22\Final\T22_241209_190750', ...
    'T22\Final\T22_241212_154105', ...
    'T22\Final\T22_241213_160938'
};

stim = 0;           % Compare Left vs Right in baseline trials
NBins = 120;        % Number of Y bins
Alpha = 0.001;       % FDR significance level
MinRun = 6;         % Consecutive bins required
testStr = 'ranksum';

yVals = nan(numel(sessions),1);

for s = 1:numel(sessions)
    basepath = fullfile(direc, sessions{s});
    noseFile = dir(fullfile(basepath, '*Tracking.Nose.mat'));
    behFile  = dir(fullfile(basepath, '*TrialBehavior.Events.mat'));
    if isempty(noseFile) || isempty(behFile)
        continue;
    end

    load(fullfile(noseFile.folder, noseFile.name), 'tracking');
    load(fullfile(behFile.folder, behFile.name), 'behavTrials');

    % Correct trials only
    [y_corr, ~] = findChoiceDivergenceY(tracking, behavTrials, ...
        'NBins', NBins, 'Alpha', Alpha, 'MinRun', MinRun, ...
        'Test', testStr, 'ShowPlot', false, ...
        'TrialFilter', behavTrials.correct==1 & behavTrials.stim'==0);

    % Error trials only
    [y_err, ~] = findChoiceDivergenceY(tracking, behavTrials, ...
        'NBins', NBins, 'Alpha', Alpha, 'MinRun', MinRun, ...
        'Test', testStr, 'ShowPlot', false, ...
        'TrialFilter', behavTrials.correct==0 & behavTrials.stim'==0);

    yVals_correct(s) = y_corr;
    yVals_error(s)   = y_err;
end

% --- Plot results ---
% --- Raincloud plot for Y-divergence values ---
valid = isfinite(yVals);
data = yVals(valid);

% --- Side-by-side raincloud plot for correct vs error ---
subplot(1,2,2); cla; hold on;

valid_corr = isfinite(yVals_correct);
valid_err  = isfinite(yVals_error);

data_corr = yVals_correct(valid_corr);
data_err  = yVals_error(valid_err);

% Parameters
violinColors = {[0.2 0.6 0.8], [0.8 0.4 0.4]}; % [Correct, Error]
jitterAmount = 0.15;
bandwidth = 2;  % adjust for smoothness

% ---- Correct half violin ----
[fC, xiC] = ksdensity(data_corr, 'Bandwidth', bandwidth);
fC = fC / max(fC) * 0.3;
fill(1 + [fC, -fC(end:-1:1)], [xiC, xiC(end:-1:1)], violinColors{1}, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Scatter (Correct)
xScatterC = 1 + jitterAmount * (rand(size(data_corr))-0.5);
scatter(xScatterC, data_corr, 60, 'filled', ...
    'MarkerFaceColor', violinColors{1}, 'MarkerFaceAlpha', 0.8);

% Mean line (Correct)
yMeanC = mean(data_corr);
plot([0.7, 1.3], [yMeanC, yMeanC], '--', 'Color', violinColors{1});
text(1.35, yMeanC, sprintf('%.2f', yMeanC), 'VerticalAlignment', 'bottom');

% ---- Error half violin ----
if sum(~isnan(yVals_error))>2
    [fE, xiE] = ksdensity(data_err, 'Bandwidth', bandwidth);
    fE = fE / max(fE) * 0.3;
    fill(2 + [fE, -fE(end:-1:1)], [xiE, xiE(end:-1:1)], violinColors{2}, ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    % Scatter (Error)
    xScatterE = 2 + jitterAmount * (rand(size(data_err))-0.5);
    scatter(xScatterE, data_err, 60, 'filled', ...
        'MarkerFaceColor', violinColors{2}, 'MarkerFaceAlpha', 0.8);
    
    % Mean line (Error)
    yMeanE = mean(data_err);
    plot([1.7, 2.3], [yMeanE, yMeanE], '--', 'Color', violinColors{2});
    text(2.35, yMeanE, sprintf('%.2f', yMeanE), 'VerticalAlignment', 'bottom');
end

% ---- Formatting ----
xlim([0.5 2.5]);
xticks([1 2]);
xticklabels({'Correct', 'Error'});
ylim([20 120]);
ylabel('Earliest significant Y');
title('Choice divergence points (Raincloud plot)');
box on; grid on;