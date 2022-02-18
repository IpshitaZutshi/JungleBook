function shape_features = compute_features (x,x_narrow, Fs, f_range, Ps, Ts, zR, zD, chNum, varargin)

% Quantify the shape of oscillatory waveforms on a cycle-by-cycle basis
% from Cole & Voytek 2019 (https://www.physiology.org/doi/pdf/10.1152/jn.00273.2019)
%
% USAGE
%   shape_features = compute_features (signal, Fs, f_range, varargin)
%
%
% INPUT
% 
%     x               1d array
%                     voltage time series
%     x_narrow        1d array. lfp filtered within f_range
%     Fs              float
%                     sampling rate (Hz)
%     f_range         frequency range for narrowband signal of interest (Hz)
%     Ps              Array of detected peaks  
%     Ts              Array of detected troughs
%     zR              Array of detected rising zero-crossings
%     zD              Array of detected falling zero-crossings
%     chNum           Current channel number in a multi-channel probe, 1 if using a single time-series.  
%     <options>       optional list of property-value pairs (see table below)
%        
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%        'powerThresh'               integer power threshold to use as cut off, 
%                                    measured in standard deviations (default = 2)
%    =========================================================================

% OUTPUT:
%    shape_features            a struct with each datapoint being one identified cycle.
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
p = inputParser;
addParameter(p,'powerThresh',2,@isnumeric)

parse(p,varargin{:});
powerThresh = p.Results.powerThresh;

% First double check lengths
lenPT = min(length(Ts),length(Ps));
lenZ = min(length(zR),length(zD));
finLen = min(lenZ,lenPT);

Ts = Ts(1:finLen);
Ps = Ps(1:finLen);
zR = zR(1:finLen);
zD = zD(1:finLen);

% Next check that the first cycle is aligned correctly
% i.e, Ts(1)<zR(1)<Ps(1)<zD(1)
if Ts(1) > Ps(1)
    if zR(1) >= zD(1)
        Ps = Ps(2:end);
        Ts = Ts(1:(finLen-1));
        zR = zR(1:(finLen-1));
        zD = zD(2:end);
    elseif zR(1) < zD(1)
        Ps = Ps(2:end);
        Ts = Ts(1:(finLen-1));
        zR = zR(1:(finLen-1));
        zD = zD(1:(finLen-1));
    end
elseif Ts(1) == Ps(1)
        Ps = Ps(2:end);
        Ts = Ts(2:end);
        zR = zR(2:end);
        zD = zD(2:end);
end
   
shape_features.sample_peak(:,chNum) = Ps;
shape_features.sample_zerox_decay(:,chNum) = zD;
shape_features.sample_zerox_rise(:,chNum) = zR;
shape_features.sample_last_trough(:,chNum) = Ts(1:(end-1));
shape_features.sample_next_trough(:,chNum) = Ts(2:end);
    
%% Compute duration of period
shape_features.period(:,chNum) = shape_features.sample_next_trough - shape_features.sample_last_trough;       

%% Compute duration of peak
shape_features.time_peak(:,chNum) = zD(1:(end-1)) - zR(1:(end-1));

%% Compute duration of last trough
shape_features.time_trough(:,chNum) = zR(2:end) - zD(1:(end-1));

%% Determine extrema voltage
shape_features.volt_peak(:,chNum) = x(Ps(1:end));
shape_features.volt_trough(:,chNum) = x(Ts(1:end));

%% Determine rise and decay characteristics
shape_features.time_decay(:,chNum) = -(Ps(1:(end-1)) - Ts(2:end));
shape_features.time_rise(:,chNum) = (Ps(1:(end-1))) -(Ts(1:(end-1)));

shape_features.volt_decay(:,chNum) = (x(Ps(1:(end-1))) - x(Ts(2:end)));
shape_features.volt_rise(:,chNum) = abs(x(Ts(1:(end-1))) - x(Ps(1:(end-1))));
shape_features.volt_amp(:,chNum) = (shape_features.volt_decay + shape_features.volt_rise) / 2;

%% Compute rise-decay symmetry features
shape_features.time_rdsym(:,chNum) = shape_features.time_rise ./ shape_features.period;

%% Compute peak-trough symmetry features
shape_features.time_ptsym(:,chNum) = shape_features.time_peak ./ (shape_features.time_peak + shape_features.time_trough);

%% Compute average oscillatory amplitude estimate during cycle
amp = fastrms(x_narrow,ceil(Fs./f_range(1)));  
for ii = 1:(length(shape_features.sample_peak(:,chNum))-1)
    shape_features.band_amp(ii,chNum) = mean(amp(Ts(ii):Ts(ii+1)));    
end

%% update intervals to remove sub-threshold power periods
%amp_delta = amp_by_time(x,Fs,[1 4]);
power = amp;

%remove periods of high noise with points above 99% cutoff as 0
power_sorted = sort(amp, 'ascend');
cutoff99 = power_sorted(ceil(0.99*length(power_sorted)));
power(power>cutoff99) = 0;

disp('finding intervals below power threshold...')
thresh = mean(power) + std(power)*powerThresh;
minWidth = (Fs./f_range(2)) * 2; % set the minimum width to two cycles
intervals = [1 length(power)];

below=find(power<thresh);
if max(diff(diff(below))) == 0
    below_thresh = [below(1) below(end)];
elseif length(below)>0
    ends=find(diff(below)~=1);
    ends(end+1)=length(below);
    ends=sort(ends);
    lengths=diff(ends);
    stops=below(ends);
    starts=lengths;
    starts = [1; starts];
    below_thresh(:,2)=stops;
    below_thresh(:,1)=stops-starts;
else
    below_thresh=[];
end

% now merge interval sets from input and power threshold
intervals = SubtractIntervals(intervals,below_thresh);  % subtract out low power intervals

intervals = intervals(diff(intervals')>minWidth,:); % only keep min width epochs

positive(1:length(amp))= 0;
for ii = 1:size(intervals,1)
    positive(intervals(ii,1):intervals(ii,2)) = 1;
end
shape_features.is_burst(:,chNum) = positive;

% %% Define whether or not each cycle is part of a burst
% if strcmp(burst_detection_method,'cycles') == 1
% %To do        df = detect_bursts_cycles(df, x, **burst_detection_kwargs)
% elseif strcmp(burst_detection_method,'amp') == 1
%     shape_features = detect_bursts_df_amp(shape_features, x, Fs, f_range, chNum);
% else
%     disp('Invalid entry for "burst_detection_method"')
%     return;
% end
end

function amp = amp_by_time(x,Fs,F_range)

% Calculate the amplitude time series
[b, a] = butter(3,[F_range(1)/(Fs/2) F_range(2)/(Fs/2)],'bandpass');
filt = filtfilt(b,a,double(x));
amp = fastrms(filt,ceil(Fs./F_range(1)));  % approximate power is frequency band
% hilb = hilbert(filt);
% lfpphase = mod(angle(hilb),2*pi);
end

function shape_features = detect_bursts_df_amp(shape_features, x, Fs, f_range, chNum)

%  Determine which cycles in a signal are part of an oscillatory
%  burst using an amplitude thresholding approach
thresh_hi = 1;
thresh_lo = 0.5;
N_cycles_min = 6;

%% Compute amplitude time series
x_magnitude = amp_by_time(x, Fs, f_range);

%% Rescale magnitude by median
x_magnitude = x_magnitude / median(x_magnitude);

%% Identify time periods of oscillation using the 2 thresholds
x_magnitude(1) = 0;
x_magnitude(end) = 0;
idx_over_hi = find(x_magnitude >= thresh_hi);

positive(1:length(x_magnitude)) = 0;
positive(idx_over_hi) = 1;

for ii = 1:length(idx_over_hi)
    j_down = idx_over_hi(ii)-1;
    if positive(j_down) == 0
        j_down_done = 0;
        while j_down_done == 0
            if x_magnitude(j_down) >=thresh_lo
                positive(j_down) = 1;
                j_down = j_down-1;
                if j_down<0
                    j_down_done = 1;
                end
            else
                j_down_done = 1;
            end
        end
    end
    j_up = idx_over_hi(ii)+1;
    
    if positive(j_up) == 0
        j_up_done = 0;
        while j_up_done == 0
            if x_magnitude(j_up) >= thresh_lo
                positive(j_up) = 1;
                j_up = j_up+1;
                if j_up >= length(x_magnitude)
                    j_up_done = 1;
                end
            else
                j_up_done = 1;
            end
        end
    end
end

%% Remove short time periods of oscillation
min_period_length = ceil(N_cycles_min * Fs / f_range(1));
shape_features.is_burst(:,chNum) = rmv_short_periods(positive, min_period_length);

end

function isosc = rmv_short_periods(x, N)

if sum(x) == 0
    isosc = x;
    return 
end

osc_changes = diff(x);
osc_starts = find(osc_changes == 1);
osc_ends = find(osc_changes == -1);

if isempty(osc_starts)
    osc_starts = 0;
end
if isempty(osc_ends)
    osc_ends = length(osc_changes);
end

if osc_ends(1) < osc_starts(1)
    osc_starts = [0 osc_starts];
end
if osc_ends(end) < osc_starts(end)
    osc_ends = [osc_ends length(osc_changes)];
end

osc_length = osc_ends - osc_starts;
osc_starts_long = osc_starts(osc_length >= N);
osc_ends_long = osc_ends(osc_length >= N);

isosc(1:length(x)) = 0;
for osc  = 1:length(osc_starts_long)
    isosc(osc_starts_long(osc):osc_ends_long(osc)) = 1;
end
   
end
