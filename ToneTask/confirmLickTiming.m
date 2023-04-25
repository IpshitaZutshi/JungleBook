function confirmLickTiming(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);

parse(p,varargin{:});
basepath = p.Results.basepath;

if ~isempty(dir([basepath filesep '*TrialBehavior.Behavior.mat'])) 
    disp('Behavior already detected! Loading file.');
    file = dir([basepath filesep '*TrialBehavior.Behavior.mat']);
    load(file(1).name);
end

file = dir([basepath filesep '*MergePoints.events.mat']);
load(file(1).name);

lfpT = bz_GetLFP(75,'noprompts',true);

% generate spectogram
params.Fs = lfpT.samplingRate; params.fpass = [1 100]; params.tapers = [3 5]; params.pad = 1;
[S,t,f] = mtspecgramc(single(lfpT.data),[1 0.5],params);
S = log10(S); % in Db
S_det= bsxfun(@minus,S,polyval(polyfit(f,mean(S,1),2),f)); 

figure
subplot(2,1,1)
imagesc(t,f,S_det',[-1.5 1.5]);
xlim([t(1) t(end)])

subplot(2,1,2)
plot(t,nanmean(S_det,2))
xlim([t(1) t(end)])

power60 = nanmean(S_det,2);
power60High = power60>1;
%find the first instances when power becomes high
powIdx = find(diff(power60High)==1);

hold on
scatter(t(powIdx+1),power60(powIdx+1))

timestampsHigh = t(powIdx+1);

% Find the psth of these timestamps with licks
[stccg, t] = CCG({behavTrials.timestamps(:,2) timestampsHigh'},[],'binSize',0.01,'duration',1);
psth = stccg(:,2,1);
figure
plot(t,psth);
end