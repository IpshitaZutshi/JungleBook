
function [hippocampalLayers] = getHippocampalLayers_IZ(varargin)
% Identify hippocampal layers based in spectral hallmarks
% 
% INPUTS
% <optional>
% 'basepath'          Default pwd
% 'lfp'               lfp        a buzcode-formatted lfp structure (use bz_GetLFP)
%                       needs fields: lfp.data, lfp.timestamps, lfp.samplingRate.
%                       If empty or no exist, look for lfp in basePath folder
% 'saveSummary'       Default true
% 'saveMat'           Detault true
% 
% OUTPUT
% hippocampalLayers   channelinfo structure with best channel for stratum
%                       oriens, piramidale, radiatum and slm, and an all
%                       field.
%
% LAYER DEFINITION
% Pyramidal layer is the highest ripple SNR channel.
% Oriens is the highest theta power channel above pyramidal channel.
% Slm is the highest theta power channel below pyramidal channel.
% Radiatum is channel with highest current sink during SPW-ripples.
%
% Manu Valero 2020
% WORK IN PROGRESS!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse options
p = inputParser;
addParameter(p,'basepath',pwd,@isstruct);
addParameter(p,'lfp',[],@isstruct);
addParameter(p,'saveSummary',true,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',true,@islogical);

parse(p,varargin{:})
basepath = p.Results.basepath;
lfp = p.Results.lfp;
saveMat = p.Results.saveMat;
saveSummary = p.Results.saveSummary;
force = p.Results.force;

%% Deal with inputs
cd(basepath);

targetFile = dir('*.hippocampalLayers.channelinfo.mat');
if ~isempty(targetFile) && ~force
    disp('Hippocampal layers already computed! Loading file.');
    load(targetFile.name);
    return
end

sessionInfo = bz_getSessionInfo;
file = dir([basepath filesep '*.region.mat']);
load(file.name);
file = dir(('*.session.mat'));
load(file.name);

% For now, keep predetermined PyrCh
if strcmp(session.animal.name,'IZ11') || strcmp(session.animal.name,'IZ15') || strcmp(session.animal.name,'IZ23') ||...
        strcmp(session.animal.name,'IZ29') || strcmp(session.animal.name,'IZ30') || strcmp(session.animal.name,'IZ32')
    overwritePyrCh = 0;
else
    overwritePyrCh = 1;
end

pyrCh = region.CA1sp; % Predefined manually
% Only select channels from the shank that includes the pyramidal channel
for ch = 1:size(sessionInfo.AnatGrps,2)
    if ismember(pyrCh, sessionInfo.AnatGrps(ch).Channels)
        channel_order = sessionInfo.AnatGrps(ch).Channels;
    end
end

% Silly exception for these animals because shank was broken
if strcmp(session.animal.name,'IZ11') || strcmp(session.animal.name,'IZ15') 
    channel_order = channel_order(1:11);
end

% If there's a reference channel, define
if isfield(session.channelTags,'RippleNoise')
    refChannel = session.channelTags.RippleNoise.channels-1;
else
    refChannel = [];
end

%% Compute channel features
lfp = bz_GetLFP('all','noPrompts', true);
% Correct noise and interpolate broken channels
lfp = bz_interpolateLFP(lfp,'refChan',refChannel);

% For slm, calculate theta
channels.pyramidal = pyrCh;
channels.channelOrder = channel_order;
powerProfile_theta = bz_PowerSpectrumProfile([5 12],'channels',channel_order,'lfp',lfp,'showfig',false,'saveMat',false); % [0:63]
powerProfile_hfo = bz_PowerSpectrumProfile([120 240],'channels',channel_order,'lfp',lfp,'showfig',false,'saveMat',false); % [0:63]
powerProfile_highpass = bz_PowerSpectrumProfile([500 1200],'channels',channel_order,'lfp',lfp,'showfig',false,'saveMat',false);

%Save the calculated profiles
channels.theta = powerProfile_theta;
channels.hfo = powerProfile_hfo;
channels.highpass = powerProfile_highpass;

%slm is channel of max theta power
channels.slm = powerProfile_theta.channels(find(powerProfile_theta.mean == max(powerProfile_theta.mean)));

%hfo near pyr
pyrChIdx = find(channel_order==channels.pyramidal);
hfo_pyr = abs(powerProfile_hfo.mean(pyrChIdx-min(7,(pyrChIdx-1)):pyrChIdx+1));
if overwritePyrCh
    %channels.pyramidal = powerProfile_theta.channels(find(hfo_theta == min(hfo_theta)));
    channels.pyramidal = powerProfile_hfo.channels(find(powerProfile_hfo.mean == max(hfo_pyr)));
end

channelsAbovePyr = channel_order(1:find(channel_order==channels.pyramidal)-1);
if isempty(channelsAbovePyr)|| length(channelsAbovePyr) == 1
    channels.oriens = NaN;
else
    [~,idxs] = ismember(channelsAbovePyr,powerProfile_theta.channels);
    channels.oriens = powerProfile_theta.channels(find(powerProfile_theta.mean == max(powerProfile_theta.mean(idxs))));
end

% For s.r., calculate sharp wave profile
ripples = bz_FindRipples(pwd,channels.pyramidal,'noise',refChannel);
twin = 0.1;

csdLFP = lfp;
[~,idxs] = ismember(channel_order,lfp.channels);
csdLFP.data = lfp.data(:,idxs);
csdLFP.channels = lfp.channels(:,idxs);
[evCsd,lfpAvg] = bz_eventCSD(csdLFP, ripples.peaks,'twin',[twin twin],'plotLFP',false,'plotCSD',false);
csdRippleProfile = [0 mean(evCsd.data(find(evCsd.timestamps > -10 & evCsd.timestamps < 10),:)) 0];
csdRippleProfile(1:find(channel_order==channels.pyramidal)) = 0;
csdRippleProfile(find(channel_order==channels.slm):end) = 0;
channels.radiatum = channel_order(find(csdRippleProfile == min(csdRippleProfile)));
channels.all = [channels.oriens; channels.pyramidal; channels.radiatum; channels.slm];

%% Summary plot
figure
set(gcf,'renderer','painters')
subplot(1,2,1)
hold on
dev1 = powerProfile_theta.mean - powerProfile_theta.std;
dev2 = powerProfile_theta.mean + powerProfile_theta.std;   
nC = 1:1:length(powerProfile_theta.channels);
hold on
fill([dev1 flip(dev2)],[nC flip(nC)],[.8 .2 .2],'FaceAlpha',.2,'EdgeColor','none');
p1 = plot(powerProfile_theta.mean,nC,'color',[.8 .2 .2]);

dev1 = powerProfile_hfo.mean - powerProfile_hfo.std;
dev2 = powerProfile_hfo.mean + powerProfile_hfo.std;
fill([dev1 flip(dev2)]+ 20,[nC flip(nC)],[.2 .2 .2],'FaceAlpha',.2,'EdgeColor','none');
p2 = plot(powerProfile_hfo.mean + 20,nC,'color',[.2 .2 .2]);
ylim([min(nC) max(nC)]);

dev1 = powerProfile_highpass.mean - powerProfile_highpass.std;
dev2 = powerProfile_highpass.mean + powerProfile_highpass.std;
fill([dev1 flip(dev2)]+ 40,[nC flip(nC)],[.2 .2 .6],'FaceAlpha',.2,'EdgeColor','none');
p3 = plot(powerProfile_highpass.mean + 40,nC,'color',[.2 .2 .6]);
ylim([min(nC) max(nC)]);

ax = axis;

plot(ax(1:2),[find(channel_order==channels.pyramidal) find(channel_order==channels.pyramidal)],'color',[.8 .2 1]);
text(ax(2),find(channel_order==channels.pyramidal),'Pyr','color',[.8 .2 1]);

plot(ax(1:2),[find(channel_order==channels.oriens) find(channel_order==channels.oriens)],'color',[.2 .2 1]);
text(ax(2),find(channel_order==channels.oriens),'Or','color',[.2 .2 1]);

plot(ax(1:2),[find(channel_order==channels.radiatum) find(channel_order==channels.radiatum)],'color',[.5 .5 .1]);
text(ax(2),find(channel_order==channels.radiatum),'Rad','color',[.5 .5 .1]);

plot(ax(1:2),[find(channel_order==channels.slm) find(channel_order==channels.slm)],'color',[.1 .8 .1]);
text(ax(2),find(channel_order==channels.slm),'Slm','color',[.1 .8 .1]);

legend([p1 p2 p3], '[5-12Hz]', '[120-240Hz]','[500-1000Hz]','Location','northwest');
set(gca,'TickDir','out'); set(gca,'YDir','rev'); xlabel('dB'); ylabel('Channels');

subplot(1,2,2)
contourf(evCsd.timestamps,(nC(2:end-1)),evCsd.data',40,'LineColor','none');hold on;
box off; colormap(jet); caxis([-max(abs(evCsd.data(:))) max(abs(evCsd.data(:)))]);
hold on
for kk = 1:size(lfpAvg.data,2)
    plot(lfpAvg.timestamps,(lfpAvg.data(:,kk)/1000)+kk-1,'k')
end

xs = [evCsd.timestamps(1) evCsd.timestamps(end)];
plot(xs, [find(channel_order==channels.pyramidal) find(channel_order==channels.pyramidal)],'color',[.8 .2 1]);
text(xs(2), find(channel_order==channels.pyramidal),'Pyr','color',[.8 .2 1]);

plot(xs,[find(channel_order==channels.oriens) find(channel_order==channels.oriens)],'color',[.2 .2 1]);
text(xs(2),find(channel_order==channels.oriens),'Or','color',[.2 .2 1]);

plot(xs,[find(channel_order==channels.radiatum) find(channel_order==channels.radiatum)],'color',[.5 .5 .1]);
text(xs(2),find(channel_order==channels.radiatum), 'Rad','color',[.5 .5 .1]);

plot(xs,[find(channel_order==channels.slm) find(channel_order==channels.slm)],'color',[.1 .8 .1]);
text(xs(2),find(channel_order==channels.slm),'Slm','color',[.1 .8 .1]);
ylim([min(nC) max(nC)]);
set(gca,'TickDir','out','YDir','reverse'); ylabel('Channels'); xlabel('Time [s]');

if saveSummary
    mkdir('SummaryFigures'); % create folder
    saveas(gcf,['SummaryFigures\hippocampalLayers.png']);
    saveas(gcf,['SummaryFigures\hippocampalLayers.eps'],'epsc');
end

hippocampalLayers = channels;
%hippocampalLayers.channels = lfp.channels;

if saveMat
    filename = split(pwd,filesep); filename = filename{end};
    save([filename '.hippocampalLayers.channelinfo.mat'],'hippocampalLayers');
end

end
