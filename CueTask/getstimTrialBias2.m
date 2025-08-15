function stimBiasSummary = getStimTrialBiasWithBiasStrength

expPath = 'Z:\Homes\zutshi01\Recordings\Cue task\';
sessName = {'T2\T2__220824_130617','T6\T6_221001_173735','T6\T6_221021_170656','T8\T8_230316_164421','T9\T9_230217_190516',...
    'T11\T11_230723_142338','T20\Final\T20_241102_163911','T22\Final\T22_241022_150631',...
    'T2\T2__220826_131832','T9\T9_230302_173934','T11\T11_230729_160436','T22\T22_241106_141549'};

baselineChoicePct = [];
stimChoicePct = [];
baselineBiasStrength = [];
stimBiasStrength = [];
biasShift = [];

for ii = 1:length(sessName)
    cd(fullfile(expPath, sessName{ii}));

    file = dir('*.TrialBehavior.Events.mat');
    load(file.name);

    % Right cue (1 = right, 0 = left)
    baselineCue = behavTrials.cue(behavTrials.stim == 0);
    stimCue = behavTrials.cue(behavTrials.stim == 1);

    % Right choices (1 = right, 0 = left)
    baselineChoice = behavTrials.choice(behavTrials.stim == 0);
    stimChoice = behavTrials.choice(behavTrials.stim == 1);

    % Percent right choices
    pctRightBaselineCue = sum(baselineCue) / length(baselineCue);
    pctRightStimCue = sum(stimCue) / length(stimCue);
    
    % Percent right choices
    pctRightBaseline = sum(baselineChoice) / length(baselineChoice);
    pctRightStim = sum(stimChoice) / length(stimChoice);

    % Store values
    baselineChoicePct(ii) = pctRightBaseline;
    stimChoicePct(ii) = pctRightStim;
    baselineCuePct(ii) = pctRightBaselineCue;
    stimCuePct(ii) = pctRightStimCue;
end


%% Updated scatter plot: Bias Strength (color-coded by target region)

% Define mPFC vs hippocampus indices
hippocampusIdx = 1:8;
mPFCIdx = 9:12;

figure;
set(gcf,'Renderer','painters')
hold on;

subplot(1,2,1)
hold on
% Hippocampus sessions = blue
scatter(baselineChoicePct(hippocampusIdx), stimChoicePct(hippocampusIdx), ...
    50, 'b', 'filled');
% mPFC sessions = red
scatter(baselineChoicePct(mPFCIdx), stimChoicePct(mPFCIdx), ...
    50, 'r', 'filled');

% Unity line
plot([0.5 0.5], [0 1], 'k--');
plot([0 1], [0.5 0.5], 'k--');

xlabel('Baseline right-ward');
ylabel('Stim right-ward');
legend('Hippocampus', 'mPFC', 'Location', 'NorthWest');
axis square;

subplot(1,2,2)
hold on
% Hippocampus sessions = blue
scatter(baselineCuePct(hippocampusIdx), stimCuePct(hippocampusIdx), ...
    50, 'b', 'filled');
% mPFC sessions = red
scatter(baselineCuePct(mPFCIdx), stimCuePct(mPFCIdx), ...
    50, 'r', 'filled');

% Unity line
plot([0.5 0.5], [0 1], 'k--');
plot([0 1], [0.5 0.5], 'k--');

xlabel('Baseline right-ward');
ylabel('Stim right-ward');
legend('Hippocampus', 'mPFC', 'Location', 'NorthWest');
axis square;

% %% 🧪 Optional: Statistical test
% [p, h, stats] = signrank(biasShift);
% fprintf('Wilcoxon signed-rank test p = %.4f\n', p);

end
