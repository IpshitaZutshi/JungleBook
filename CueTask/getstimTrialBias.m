function getstimTrialBias

expPath = 'Z:\Homes\zutshi01\Recordings\Cue task\';
sessName = {'T2\T2__220824_130617','T6\T6_221001_173735','T6\T6_221021_170656','T8\T8_230316_164421','T9\T9_230217_190516',...
    'T11\T11_230723_142338','T20\Final\T20_241102_163911','T22\Final\T22_241022_150631'};

for ii = 1:length(sessName)
    cd(strcat(expPath,sessName{ii}))

    %% Within each folder, load the behavTrials files
    file = dir('*.TrialBehavior.Events.mat');
    load(file.name);

    baselineCue = behavTrials.cue(behavTrials.stim==0);
    stimCue = behavTrials.cue(behavTrials.stim==1);

    baselineCueMice(ii) = sum(baselineCue)./length(baselineCue);
    stimCueMice(ii) = sum(stimCue)./length(stimCue);

    baselineChoice = behavTrials.choice(behavTrials.stim==0);
    stimChoice = behavTrials.choice(behavTrials.stim==1);

    baselineChoiceMice(ii) = sum(baselineChoice)./length(baselineChoice);
    stimChoiceMice(ii) = sum(stimChoice)./length(stimChoice);

end

end