function [ripple_power] = rippleBandPower(varargin)


%% Modified from bz_FindRipples

% USAGE
% call the function the same way you would call bz_findRipples
% 
% pyrCh = 121; 
% noiseCh = 111;
% [ripple_power] = rippleBandPower(pwd,pyrCh,'noise',noiseCh,'savemat',true,'durations',[30 100],'passband',[130 200]);
% 
% 
% Lucy Anderson 2025
% 
% 
%
%    Ripples are detected using the normalized squared signal (NSS) by
%    thresholding the baseline, merging neighboring events, thresholding
%    the peaks, and discarding events with excessive duration.
%    Thresholds are computed as multiples of the standard deviation of
%    the NSS. Alternatively, one can use explicit values, typically obtained
%    from a previous call.  The estimated EMG can be used as an additional
%    exclusion criteria
%
% INPUTS - note these are NOT name-value pairs... just raw values
%    lfp            unfiltered LFP (one channel) to use
%	 timestamps	    timestamps to match filtered variable
%    <options>      optional list of property-value pairs (see tables below)
%
%    OR
%
%    basepath       path to a single session to run findRipples on
%    channel      	Ripple channel to use for detection (0-indexed, a la neuroscope)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'thresholds'  thresholds for ripple beginning/end and peak, in multiples
%                   of the stdev (default = [2 5]); must be integer values
%     'durations'   min inter-ripple interval and max ripple duration, in ms
%                   (default = [30 100]). 
%     'minDuration' min ripple duration. Keeping this input nomenclature for backwards
%                   compatibility
%     'restrict'    interval used to compute normalization (default = all)
%     'frequency'   sampling rate (in Hz) (default = 1250Hz)
%     'stdev'       reuse previously computed stdev
%     'show'        plot results (default = 'off')
%     'noise'       noisy unfiltered channel used to exclude ripple-
%                   like noise (events also present on this channel are
%                   discarded)
%     'passband'    N x 2 matrix of frequencies to filter for ripple detection 
%                   (default = [130 200])
%     'EMGThresh'   0-1 threshold of EMG to exclude noise
%     'saveMat'     logical (default=false) to save in buzcode format
%     'plotType'   1=original version (several plots); 2=only raw lfp
%    =========================================================================
%
% OUTPUT
%
%    ripples        buzcode format .event. struct with the following fields
%                   .timestamps        Nx2 matrix of start/stop times for
%                                      each ripple
%                   .detectorName      string ID for detector function used
%                   .peaks             Nx1 matrix of peak power timestamps 
%                   .stdev             standard dev used as threshold
%                   .noise             candidate ripples that were
%                                      identified as noise and removed
%                   .peakNormedPower   Nx1 matrix of peak power values
%                   .detectorParams    struct with input parameters given
%                                      to the detector

warning('this function is under development and may not work... yet')

% Default values
p = inputParser;
addParameter(p,'thresholds',[2 5],@isnumeric)
addParameter(p,'durations',[30 100],@isnumeric)
addParameter(p,'restrict',[],@isnumeric)
addParameter(p,'frequency',1250,@isnumeric)
addParameter(p,'stdev',[],@isnumeric)
addParameter(p,'show','off',@isstr)
addParameter(p,'noise',[],@ismatrix)
addParameter(p,'passband',[130 200],@isnumeric)
addParameter(p,'EMGThresh',.9,@isnumeric);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'minDuration',20,@isnumeric)
addParameter(p,'plotType',2,@isnumeric)

if isstr(varargin{1})  % if first arg is basepath
    addRequired(p,'basepath',@isstr)
    addRequired(p,'channel',@isnumeric)    
    parse(p,varargin{:})
    basename = bz_BasenameFromBasepath(p.Results.basepath);
    basepath = p.Results.basepath;
    passband = p.Results.passband;
    EMGThresh = p.Results.EMGThresh;
    lfp = bz_GetLFP(p.Results.channel,'basepath',p.Results.basepath,'basename',basename);%currently cannot take path inputs
    signal = bz_Filter(double(lfp.data),'filter','butter','passband',passband,'order', 3);
    timestamps = lfp.timestamps;
elseif isnumeric(varargin{1}) % if first arg is filtered LFP
    addRequired(p,'lfp',@isnumeric)
    addRequired(p,'timestamps',@isnumeric)
    parse(p,varargin{:})
    passband = p.Results.passband;
    EMGThresh = p.Results.EMGThresh;
    signal = bz_Filter(double(p.Results.lfp),'filter','butter','passband',passband,'order', 3);
    timestamps = p.Results.timestamps;
    basepath = pwd;
    basename = bz_BasenameFromBasepath(basepath);
end

% assign parameters (either defaults or given)
frequency = p.Results.frequency;
show = p.Results.show;
restrict = p.Results.restrict;
sd = p.Results.stdev;
noise = p.Results.noise;
lowThresholdFactor = p.Results.thresholds(1);
highThresholdFactor = p.Results.thresholds(2);
minInterRippleInterval = p.Results.durations(1);
maxRippleDuration = p.Results.durations(2);
minRippleDuration = p.Results.minDuration;
plotType = p.Results.plotType;

%% filter and calculate noise


% Parameters
windowLength = frequency/frequency*11;

% Square and normalize signal
squaredSignal = signal.^2;
% squaredSignal = abs(opsignal);
window = ones(windowLength,1)/windowLength;
keep = [];
if ~isempty(restrict)
    for i=1:size(restrict,1)
        keep = InIntervals(timestamps,restrict);
    end
end
keep = logical(keep); 

[normalizedSquaredSignal,sd] = unity(Filter0(window,sum(squaredSignal,2)),sd,keep);


ripple_power.normed_power_trace = normalizedSquaredSignal;
ripple_power.timestamps = timestamps; 


%Save
if p.Results.saveMat
    save(fullfile(basepath, [basename '.rippleBandPower.mat']),'ripple_power')
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

