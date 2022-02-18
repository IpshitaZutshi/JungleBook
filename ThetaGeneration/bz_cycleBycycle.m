
function [shape_features, data, lfp_lowpass] = bz_cycleBycycle (varargin)

% Calculates lfp frequency, amplitude and waveform assymetry using a cycle
% by cycle analysis - inspired by Cole & Voytek 2019 (https://www.physiology.org/doi/pdf/10.1152/jn.00273.2019)
%
% USAGE
%   shape_features = bz_cycleBycycle (<options>)
%
%
% INPUT
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'lowpass'         frequency for lowpass filter. [0 X], default [0.1 25]
%     'passband'        pass frequency range. [X Y], default [5 10] for theta
%     'samplingRate'    lfp sampling rate. default 1250  
%     'channels'        if input is a buzcode lfp structure with field
%                       samples.channels, will only filter the selected
%                       channels
%     'overwrite'       If precomputed, overwrite existing file. true or false 
%     'powerThresh'     integer power threshold to use as cut off, 
%                       measured in standard deviations (default = 2)
%    =========================================================================

% OUTPUT:
%    shape_features    a struct with each datapoint being one identified cycle.
%                              Structure subfields:
%                             - sample_peak : sample of 'x' at which the peak occurs
%                             - sample_zerox_decay : sample of the decaying zerocrossing
%                             - sample_zerox_rise : sample of the rising zerocrossing
%                             - sample_last_trough : sample of the last trough
%                             - sample_next_trough : sample of the next trough
%                             - period : period of the cycle
%                             - time_decay : time between peak and next trough
%                             - time_rise : time between peak and previous trough
%                             - time_peak : time between rise and decay zerocrosses
%                             - time_trough : duration of previous trough estimated by zerocrossings
%                             - volt_decay : voltage change between peak and next trough
%                             - volt_rise : voltage change between peak and previous trough
%                             - volt_amp : average of rise and decay voltage
%                             - volt_peak : voltage at the peak
%                             - volt_trough : voltage at the last trough
%                             - time_rdsym : fraction of cycle in the rise period
%                             - time_ptsym : fraction of cycle in the peak period
%                             - band_amp : average analytic amplitude of the oscillation
%                               computed using narrowband filtering and the Hilbert
%                               transform. Filter length is 3 cycles of the low
%                               cutoff frequency. Average taken across all time points
%                               in the cycle.
%                             - is_burst : True if the cycle is part of a detected oscillatory burst

% Ipshita Zutshi, 7/19

%% Parse inputs
sess = bz_getSessionInfo(pwd,'noPrompts',true);      

p = inputParser;
addParameter(p,'lowpass',[0.1 25],@isnumeric);
addParameter(p,'passband',[5 12],@isnumeric);
addParameter(p,'samplingRate',1250,@isnumeric);
addParameter(p,'channels',0:(sess.nChannels-1),@isvector);
addParameter(p,'overwrite',false,@islogical);
addParameter(p,'powerThresh',2,@isnumeric)

parse(p,varargin{:});
lowpass = p.Results.lowpass;
passband = p.Results.passband;
samplingRate = p.Results.samplingRate;
channels = p.Results.channels;
overwrite = p.Results.overwrite;
powerThresh = p.Results.powerThresh;

% Check cutoff frequency input is valid
if length(passband) ~= 2
    disp('Two cutoff frequencies required for bandpass and bandstop filters')
    return;
end
if passband(1) >= passband(2)
    disp('Second cutoff frequency must be greater than first for bandpass and bandstop filters')
    return;
end

% if exist([sess.FileName '.Bycycle.mat'],'file') && ~overwrite
%     disp('Cycle by cycle analysis already performed! Loading file.');
%     load([sess.FileName '.Bycycle.mat']);
%     return
% end

if isempty(dir('*.lfp'))
    try 
        bz_LFPfromDat(pwd,'outFs',1250); % generating lfp
    catch
        disp('Problems with bz_LFPfromDat, resampling...');
        ResampleBinary(strcat(sessionInfo.session.name,'.dat'),strcat(sessionInfo.session.name,'.lfp'),...
            sessionInfo.nChannels,1,sessionInfo.rates.wideband/sessionInfo.rates.lfp);
    end
end

lfp = bz_GetLFP(channels,'noPrompts', true);

% Go back from 0 indexing for channels
channels = channels+1;

%lfp input
if isstruct(lfp)
    data = lfp.data;
    timestamps = lfp.timestamps;
    samplingRate = lfp.samplingRate;
elseif iscell(lfp) %for multiple trials
    celllengths = cellfun(@length,lfp);
    data = vertcat(lfp{:});
elseif isnumeric(lfp)
    data = lfp;
    timestamps = [1:length(lfp)]'./samplingRate;
end


%% First low-pass filter to remove high-frequency artifacts

Wn = [lowpass(1)/(samplingRate/2) lowpass(2)/(samplingRate/2)];
[b,a] = butter(3,Wn);
for ii = 1:size(data,2)
    lfp_lowpass(:,ii) = filtfilt(b,a,double(data(:,ii)));
end

%% Next, filter in the desired frequency range

Wn = [passband(1)/(samplingRate/2) lowpass(2)/(samplingRate/2)];
[b,a] = butter(3,Wn);
for ii = 1:length(channels)
    lfp_narrow(:,ii) = filtfilt(b,a,double(data(:,ii)));
end

%% Find peaks, troughs and zero-crossings

for ii = 1:length(channels)

    %% Find rising and falling zerocrossings (narrowband)    % Rising zerocrossings
    pos = diff(lfp_narrow(:,ii) < 0); 
    zeroriseN = find(pos==-1);

    % Falling zerocrossings
    zerofallN = find(pos==1);

    %Compute number of peaks and troughs
    if zeroriseN(end) > zerofallN(end)
        P = length(zeroriseN)-1;
        T = length(zerofallN);
    else
        P = length(zeroriseN);
        T = length(zerofallN)-1;
    end
    
%% Calculate peaks (using lowpass data)
    for pp = 1:P
        % find sample range between most recent zero rise and the next zero fall
        startzerorise = zeroriseN(pp);
        endzerorise = zerofallN(zerofallN>startzerorise);
        % find peak
        [~, idx] = max(lfp_lowpass(startzerorise:endzerorise(1),ii));
        Ps(pp) = idx+startzerorise;
    end

%% Calculate troughs (using lowpass data)
    for tt = 1:T
        % find sample range between most recent zero rise and the next zero fall
        startzerofall = zerofallN(tt);
        endzerofall = zeroriseN(zeroriseN>startzerofall);
        % find trough
        [~, idx] = min(lfp_lowpass(startzerofall:endzerofall,ii));
        Ts(tt) = idx+startzerofall;
    end    

%% Find zerocrossings within each cycle after peaks and troughs are identified
    [zR, zD] = find_zerox(lfp_lowpass(:,ii), Ps, Ts);
    
%% Compute features
    shape_features = compute_features(lfp_lowpass(:,ii),lfp_narrow(:,ii), samplingRate,passband,Ps,Ts,zR,zD,ii,'powerThresh',powerThresh); 
end
   % save([sess.FileName '.Bycycle.mat'],'shape_features');
end

function [zR, zD] = find_zerox(x,Ps,Ts)

%     A rising zerocrossing occurs when the voltage crosses
%     midway between the trough voltage and subsequent peak voltage.
%     A decay zerocrossing is defined similarly.
%     If this voltage is crossed at multiple times, the temporal median is taken
%     as the zerocrossing.

%     Parameters
%     ----------
%     x : array-like 1d
%         voltage time series
%     Ps : array 1d
%         time points of oscillatory peaks
%     Ts : arrays 1d
%         time points of osillatory troughs
% 
%     Returns
%     -------
%     zeroxR : array-like 1d
%         indices at which oscillatory rising zerocrossings occur
%     zeroxD : array-like 1d
%         indices at which oscillatory decaying zerocrossings occur


% Calculate the number of rises and decays
if Ps(1) < Ts(1)
    N_rises = length(Ps) - 1;
    N_decays = length(Ts);
    idx_bias = 0;
else
    N_rises = length(Ps);
    N_decays = length(Ts) - 1;
    idx_bias = 1;
end

% Find zerocrossings for rise
zR(1:N_rises) = 0;
for ii = 1:N_rises
    x_temp = x(Ts(ii):Ps(ii + 1 - idx_bias) + 1);
    if isempty(x_temp)
        x_temp = 0;
    end
    x_temp  = x_temp - ((x_temp(1) + x_temp(end))/2);
    
% If data is all 0s, just set the zerocrossing to be halfway between.    
    if sum(x_temp)==0
        zR(ii) = Ts(ii)+(length(x_temp)/2);
% If rise is actually decay, just set the zerocrossing to be halfway between.        
    elseif x_temp(1) > x_temp(end)
        zR(ii) = Ts(ii)+(length(x_temp)/2);
    else
         pos = diff(x_temp < 0); 
        zeroriseN = find(pos==-1);
        zR(ii) = Ts(ii)+(median(zeroriseN));
    end
end

% Find zerocrossings for decay
zD(1:N_decays) = 0;
for ii = 1:N_decays
    x_temp = x(Ps(ii):Ts(ii + idx_bias) + 1);
    x_temp  = x_temp - ((x_temp(1) + x_temp(end))/2);
    if isempty(x_temp)
        x_temp = 0;
    end    
% If data is all 0s, just set the zerocrossing to be halfway between.    
    if sum(x_temp)==0
        zD(ii) = Ps(ii)+(length(x_temp)/2);
% If rise is actually decay, just set the zerocrossing to be halfway between.        
    elseif x_temp(1) > x_temp(end)
        zD(ii) = Ps(ii)+(length(x_temp)/2);
    else
        pos = diff(x_temp < 0); 
        zerofallN = find(pos==1);
        zD(ii) = Ps(ii)+(median(zerofallN));
    end
end

end
