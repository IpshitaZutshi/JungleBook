function [DS1,DS2] = autoDSdetectionIZ(varargin)
% The previous return values [DS2amplitudes,ml_chan,h_chan] are found in
%   the recent version as DS2.amplitudes, DS2.detectorinfo.h_channel and 
%   .ml_channel
p = inputParser;
addParameter(p,'hilus_chan',[],@isnumeric);
addParameter(p,'ml_chan',[],@isnumeric);

parse(p,varargin{:})

hilus_chan = p.Results.hilus_chan;
ml_chan = p.Results.ml_chan;

basepath = pwd;
basename = bz_BasenameFromBasepath(basepath);
outputfolder = 'Analysis/LFP';

if ~isfolder(outputfolder)
    mkdir(outputfolder);
end

%% Load from file if run already
% if isfile([basename '.DS1.events.mat']) && isfile([basename '.DS2.events.mat'])
%     load([basename '.DS1.events.mat'])
%     load([basename '.DS2.events.mat'])
%     fprintf('%s dentate spikes already detected, loading from file\n',basepath);
%     return
% end

%% Exclude badchannels
% Note: .badchannels are 1-INDEXED. SUBTRACT 1 when transfering to bz_fkts.
selected = true(1,size(bz_GetLFP('all','basepath',basepath,'intervals',[1 2]).data,2)); % Load a snipet to determine channel count

[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true); 

if isfile(fullfile(basepath,filesep,[basename '.session.mat']))
    load(fullfile(basepath,filesep,[basename '.session.mat']));
    if isfield(session.channelTags,'Bad')
        badchannels = session.channelTags.Bad.channels;
        fprintf('Bad channels loaded from %s.session.mat\n',basename);
    end
elseif isfile(fullfile(basepath,filesep,[basename '.sessionInfo.mat']))
    load(fullfile(basepath,filesep,[basename '.sessionInfo.mat']));
    if isfield(sessionInfo,'badchannels')
        badchannels = sessionInfo.badchannels;
        fprintf('Bad channels loaded from %s.sessionInfo.mat\n',basename);
    else
        badchannels = [];
        fprintf('No badchannels found, leaving empty\n');
    end
else
    if ~exist('badchannels','var')
        badchannels = [];
        fprintf('No badchannels found, leaving empty\n');
    end
end

%% Get regions
if isfile(fullfile(basepath,filesep,'region_template.mat'))
    load(fullfile(basepath,filesep,'region_template.mat'),'region_template')
    isinDG = strcmp(region_template,'DG');
    fprintf('Loaded region template from %s\n',basepath);
elseif isfile(fullfile(basepath,filesep,'..',filesep,'region_template.mat'))
    load(fullfile(basepath,filesep,'..',filesep,'region_template.mat'));
    isinDG = strcmp(region_template,'DG');
    fprintf('Loaded region template from top folder.\n');
else
    isinDG = [];
    fprintf('No region labels found, all in.\n');
end

selected(badchannels) = false;
if ~isempty(isinDG)
    selected = find(and(selected,isinDG));
else
    selected = find(selected);
end

%% Determine the best channels  
if isfile('chanMap.mat')
    load('chanMap.mat', 'xcoords')
    load('chanMap.mat', 'ycoords')
else
    chanMap = createChannelMap(basepath,basename,'staggered');
    xcoords = chanMap.xcoords;
    ycoords = chanMap.ycoords;
end

file = dir(('*.hippocampalLayersCSD.channelinfo.mat'));
load(file.name);  
pyrCh = hippocampalLayers.pyramidal; 

% Only select channels from the shank that includes the pyramidal channel
for ch = 1:size(sessionInfo.AnatGrps,2)
    if ismember(pyrCh, sessionInfo.AnatGrps(ch).Channels)
        channel_order = sessionInfo.AnatGrps(ch).Channels;
    end
end

ml_chan = hippocampalLayers.ml;
idx_ml = find(channel_order==ml_chan);
if length(channel_order) == 16
    if isempty(idx_ml)
        idx_ml = find(channel_order==hippocampalLayers.slm);
        idx_ml = idx_ml+1;
        ml_chan = channel_order(idx_ml);
    end
    hilus_chan = channel_order(idx_ml+2);
elseif length(channel_order) == 64
    if idx_ml<48
        hilus_chan = channel_order(idx_ml+12);
    else
        hilus_chan = channel_order(end);
    end
end

if isempty(hilus_chan) || isempty(ml_chan)
    
    lfp = bz_GetLFP(channel_order,'basepath',basepath,'intervals',[60 660]); % Get 10 minutes of LFP    
    lfp = bz_interpolateLFP(lfp);
    gammalfp = bz_Filter(lfp,'passband',double([30 120]));
    gammapow(:,1) = 1:size(gammalfp.amp,2);
    gammapow(:,2) = nanmean(gammalfp.amp);
    
    figure; title('Best Gamma Channel');
    plot(gammapow(:,2),gammapow(:,1))
    set(gca,'YDir','reverse')
    %scatter(xcoords,ycoords,[],gammapow(:,2),'filled');
    
    [~,idx]=sort(gammapow(:,2),'descend');
    gammapow = gammapow(idx,:);
    gammachan = gammapow(1,1); %Use highest gamma power channel
    hold on; line([min(gammapow(:,2)) max(gammapow(:,2))],[gammachan gammachan],'Color','k');

    %% Find raw events to do stratification
    trace = bz_GetLFP(gammachan,'basepath',basepath,'intervals',[60 660]);%[60 inf]);
    trace = bz_Filter(trace,'passband',[10 inf]); % Filter out low frequency
    evts = find(diff(abs(double(trace.data))>...
        5*nanstd(double(trace.data)))==1); % Threshold crossings
    mamps = median(double(lfp.data(evts,:)),1);
    
    [~,hidx] = max(mamps);
    [~,midx] = min(mamps);
    fprintf('Selected %i as hilus and %i as molecular-layer channel (0-Index).\n',hilus_chan,ml_chan);
else
    fprintf('%i was given as hilus- and %i was given as molecular-layer channel (0-Index).\n',hilus_chan,ml_chan);
end

clear amps mamps

%% Yuta's DS detection
%LImit detection to non-rem periods

[DS2] = DetectDSpikes_v4(basename,hilus_chan,ml_chan,[0 Inf]);

%% Generate and save DS2 amplitude maps
% fprintf('%s DS detection complete, getting DS2 amplitudes\n',basepath);
% lfp = bz_GetLFP(channel_order,'basename',basename);
% 
% DS2amplitudes = median(double(lfp.data(round(DS2.peaks*1250),:)),1)*.00038; % Convert to mV.
% 
% DS2mapFig = figure; hold on; title('DS2 amplitudes')
% scatter(xcoords(1:length(channel_order)),ycoords(1:length(channel_order)),[],DS2amplitudes,'filled');
% scatter(xcoords(hilus_chan+1),ycoords(hilus_chan+1),'r');
% scatter(xcoords(ml_chan+1),ycoords(ml_chan+1),'r');
% colorbar
% 
% save([basename '.DS2map.mat'],'xcoords','ycoords','DS2amplitudes');
% saveas(DS2mapFig,fullfile(outputfolder,filesep,'DS2map.fig'));
% saveas(DS2mapFig,fullfile(outputfolder,filesep,'DS2map.png'));
end
