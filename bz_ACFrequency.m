function [ACData] = bz_ACFrequency(spikes,varargin)
% USAGE
%[ACData] = bz_ACFrequency(varargin)
% 
% INPUTS
% spikes        -spike time cellinfo struct
%
% lfp           -lfp struct with a single channel from bz_GetLFP()
%
% passband      -frequency range for AC modulation [lowHz highHz] form
%
% intervals     -(optional) may specify timespans over which to calculate 
%               autocorrelation.  Formats accepted: tstoolbox intervalSet
%               or a 2column matrix of [starts stops] in seconds
%
% samplingRate  -specifies lfp sampling frequency default=1250
%
%
% plotting      -logical if you want to plot, false if not, default=true
%
% saveMat       -logical to save cellinfo .mat file with results, default=true
%
%
% OUTPUTS
%
% cor           - mxn matrix of the autocorrelation
%               for each cell. m, number of cells, n,range of lags
%
% lag           - AC lags
%
% spectrum      - mxn matrix of the spectral power of the autocorrelation
%               for each cell. m, number of cells, n,range of frequencies
%
% frequency     - Range of frequency bins
%
% lfpfreq       - mx1 array, If lfp is provided, peak lfp theta frequency whenever
%                 cell spikes. m, number of cells
%
% unitfreq      - mx1 array, Peak unit AC frequency
%
% unitpower     - mx1 array. Power at peak unit AC frequency.
%
%
% Calculates peak frequency of the autocorrelation of spikes using the
% chronux toolbox
% If provided, also calculates the peak LFP theta frequency when the cell
% spikes.
%
% Ipshita Zutshi 2021

%% defaults
p = inputParser;
% addRequired(p,'spikes',[],@isstruct);
addParameter(p,'lfp',[],@isstruct)
addParameter(p,'passband',[6 14],@isnumeric)
addParameter(p,'intervals',[0 inf],@isnumeric)
addParameter(p,'samplingRate',1250,@isnumeric)
addParameter(p,'tbin',0.005,@isnumeric)
addParameter(p,'forceReload',true,@islogical)
addParameter(p,'saveMat',false,@islogical)

parse(p,varargin{:})

% spikes = p.Results.spikes;
lfp = p.Results.lfp;
passband = p.Results.passband;
intervals = p.Results.intervals; % interval(s) over which to calculate
samplingRate = p.Results.samplingRate; % sampling rate of continuous signal (LFP)
tbin = p.Results.tbin; % time resolution for AC
saveMat = p.Results.saveMat;

if exist([spikes.sessionName '.ACData.cellinfo.mat'],'file') && forceReload == false
    disp('loading AC Data from cellinfo file..')
    load([spikes.sessionName '.ACData.cellinfo.mat'])
    
else

    %% Get peak frequency for every time point in LFP
    if ~isempty(lfp)
        [wave,f_lfp,t_lfp]=getWavelet(double(lfp.data(:,1)),samplingRate,passband(1),passband(2),40,0);
        [~,mIdx]=max(wave);%get index max power for each timepoint
        freq = f_lfp(mIdx);
        t_lfp = sort(t_lfp);
    end

    %% Extract spikes in the specified intervals

    for a = 1:length(spikes.times)

        bools = InIntervals(spikes.times{a},intervals);
        s =spikes.times{a}(bools);
        if length(s)<5
%             ACData.spectrum(a,:) = zeros(1,655); 
%             ACData.cor(a,:) = zeros(1,201);
%             ACData.lag(1,:) =  zeros(1,201);            
%             ACData.frequency(1,:) = zeros(1,655);  
            ACData.unitfreq(a) = nan;  
            ACData.unitpower(a) = nan;  
            if ~isempty(lfp)
                ACData.lfpfreq(a) = nan;  
            end
        else
        %% calculate crosscorr    
            [cor, lag] = CrossCorr(s, s);
            [cor, lag, smooth] = SmoothCor(cor, lag, tbin); % solve for smoothed signal, and chop off center peak and negative half of autocorrelation
            [S_units,f_units] = getChronuxSpectrum(cor, tbin);

            S_units = S_units(f_units>passband(1)&f_units<passband(2));
            f_units = f_units(f_units>passband(1)&f_units<passband(2));
            [maxP,freqP] = max(S_units);
            
            ACData.cor(a,:) = cor;
            ACData.lag(1,:) = lag;
            ACData.spectrum(a,:) = S_units;
            ACData.frequency(1,:) = f_units;
            ACData.unitfreq(a) = f_units(freqP);  
            ACData.unitpower(a) = maxP;
            
         %% calculate lfp frequency   
            if ~isempty(lfp)           
                %Find lfp power at the spiketimes
                s = sort(s);
                [~, ~, bins] = histcounts(t_lfp,s);
                [~,y] = unique(bins);
                closestIndex = y(2:end);
                lfpFreq = freq(closestIndex);
                ACData.lfpfreq(a) = nanmean(lfpFreq);
            end

        end
    end

    ACData.UID = spikes.UID;
    ACData.sessionName = spikes.sessionName;

    if saveMat
     save([spikes.sessionName '.ACData.cellinfo.mat'],'ACData');
    end
end
end
 
 function [cor, lag, smooth] = SmoothCor(cor, lag, t_bin)
    
    std_smooth_kernel = .005;
    
    kernel = pdf('normal', -std_smooth_kernel*10/t_bin:std_smooth_kernel*10/t_bin, 0, std_smooth_kernel/t_bin); % convolve with gaussian

    if isempty(cor), smooth = []; return; % if no signals, then send it back
    
    else
        
        smooth = zeros(size(cor));
        
        for i = 1:size(cor,2)  

            smooth(:,i) = ndnanfilter(cor(:,i), kernel(:), []);
        
        end
        
        lag = lag(round(end/2):end,:);
        
        cor = cor(round(end/2):end,:);
        
        smooth = smooth(round(end/2):end, :);

    end 
    
end

function [S, f] = getChronuxSpectrum(cor, t_bin)

    params.tapers = [1.25, 1];
    params.fpass = [];
    params.Fs = 1/t_bin;
    params.pad = 6;  

    [S,f]=mtspectrumc(cor,params); % compute spectrum of cor

end
