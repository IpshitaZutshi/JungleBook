function [ripples] = rippleDetectorAndrea(varargin)

%% Default values
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'rippleChannels',[],@isnumeric);
addParameter(p,'threshold',0.6,@isnumeric);
addParameter(p,'passband',[130 200],@isnumeric);
addParameter(p,'thresholds',[2 5],@isnumeric)
addParameter(p,'stdev',[],@isnumeric);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'padding',0.1,@isnumeric);
addParameter(p,'minDuration',0.020*1250,@isnumeric)

parse(p,varargin{:})
basepath = p.Results.basepath;
rippleChannels = p.Results.rippleChannels;
passband = p.Results.passband;
threshold = p.Results.threshold;
saveMat = p.Results.saveMat;
sd = p.Results.stdev;
lowThresholdFactor = p.Results.thresholds(1);
highThresholdFactor = p.Results.thresholds(2);
padding = p.Results.padding;
minDuration = p.Results.minDuration;

%% Load Session Metadata and several variables if not provided
cd(basepath);
session = sessionTemplate(basepath,'showGUI',false);
if isempty(rippleChannels)
    if ~isempty(dir([session.general.name,'.hippocampalLayers.channelinfo.mat']))
        file = dir([session.general.name,'.hippocampalLayers.channelinfo.mat']);
        load(file.name);
    end
    pyrCh = find(hippocampalLayers.channelOrder==hippocampalLayers.pyramidal);
    rippleChannels = hippocampalLayers.channelOrder(pyrCh-3:pyrCh+4);
end

%Load lfp
lfp = bz_GetLFP(hippocampalLayers.channelOrder,'noPrompts', true);
% Correct noise and interpolate broken channels
lfp = bz_interpolateLFP(lfp);
data = double(lfp.data(:,pyrCh-3:pyrCh+4));

%% Calculate SPW-R probability using CNN
ripples.swrProb = detect_ripples_cnn(data, lfp.samplingRate);
timestamps = get_intervals(ripples.swrProb,'threshold',threshold);

%% Calculate peak timestamps and amplitude, exclude noise ripples
signal = bz_Filter(double(lfp.data(:,pyrCh)),'filter','butter','passband',passband,'order', 3);
% Parameters
windowLength = lfp.samplingRate/lfp.samplingRate*11;

% Square and normalize signal
squaredSignal = signal.^2;
% squaredSignal = abs(opsignal);
window = ones(windowLength,1)/windowLength;
keep = [];
keep = logical(keep); 

[normalizedSquaredSignal,sd] = unity(Filter0(window,sum(squaredSignal,2)),sd,keep);
thresholded = normalizedSquaredSignal > lowThresholdFactor;

%for each interval, calculate peak timestamp
curRip = 0;
for ii = 1:length(timestamps)
    % add padding
    int = (timestamps(ii,1)-(padding*1250)):(timestamps(ii,2)+(padding*1250));
    [maxValue,peakSig] = max(normalizedSquaredSignal(int));
    start = find(diff(thresholded(int))>0,1,'first');
    stop = find(diff(thresholded(int))<0,1,'last');
    if (maxValue < lowThresholdFactor) || isempty(start) || isempty(stop) || (stop<start) ||(stop-start)<minDuration
        continue
    else
        curRip = curRip+1;    
        ripples.peaks(curRip) = lfp.timestamps(int(1)+ peakSig);
        ripples.timestamps(curRip,1) = lfp.timestamps(int(1) + start);
        ripples.timestamps(curRip,2) = lfp.timestamps(int(1) + stop);    
    end    
end

%% Ripple Stats
if size(ripples.timestamps,1)>2
    %lfp = bz_GetLFP('all');
    filtered = bz_Filter(lfp,'channels',hippocampalLayers.pyramidal,'filter','butter','passband',passband,'order',3);
    ripples.detectorinfo.detectionparms.frequency = lfp.samplingRate;
    [maps,data,stats] = bz_RippleStats(filtered.data,filtered.timestamps,ripples);
    ripples.maps = maps;
    ripples.data = data;
    ripples.stats = stats;
    %bz_PlotRippleStats(ripples.maps, ripples.data, ripples.stats);
end

if saveMat
    disp('Saving Ripples Results...');
    save([session.general.name , '.ripplesCNN.events.mat'],'ripples');    
end

function y = Filter0(b,x)

if size(x,1) == 1
	x = x(:);
end

if mod(length(b),2)~=1
	error('filter order should be odd');
end

shift = (length(b)-1)/2;

[y0 z] = filter(b,1,x);

y = [y0(shift+1:end,:) ; z(1:shift,:)];

function [U,stdA] = unity(A,sd,restrict)

if ~isempty(restrict),
	meanA = mean(A(restrict));
	stdA = std(A(restrict));
else
	meanA = mean(A);
	stdA = std(A);
end
if ~isempty(sd),
	stdA = sd;
end

U = (A - meanA)/stdA;
