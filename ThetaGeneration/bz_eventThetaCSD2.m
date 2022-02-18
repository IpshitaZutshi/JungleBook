
function [ csd, lfpAvg ] = bz_eventThetaCSD2 (lfp, events, varargin)

% [ CSD ] = bz_eventThetaCSD (lfp, events, varargin)
% Calculates event-triggered (i.e. opto stim) CSD map from a linear array of LFPs

% INPUT
%    lfp            a buzcode structure with fields lfp.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%                   -lfp can also be a [t x 1] timeseries signal. in which
%                   case you need to input 'samplingRate'
%   events          events timestamps (in sec)

%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%       channels    vector with channels to inlcude. If empty take all (default)
%       twin        time window around events to calculate average. Default [0.1 0.1]
%       spat_sm     degree of spatial smoothing. Default = 0.
%       temp_sm     degree of temporal smoothing. Default = 0.
%       plotCSD     true/false. Default true.
%       plotLFP     true/false. Default true.
%       NumStd      number of standard deviations above filtered rms signal to
%                   detect theta peaks, default 1.5
%       fRange      Frequency band to filter theta, default [6 10]
%       thetaCh     Specify channel to align peaks from other channels
%    =========================================================================

% OUTPUT:
%    CSD            a buzcode structure with fields csd.data,
%                                                   csd.timestamps
%                                                   csd.samplingRate
%                                                   csd.channels 
%                                                   csd.params
%    lfpAvg         a buzcode structure with fields lpf.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%                                                   lfp.channels 
%                                                   lfp.params
%
%

%% Parse inputs

p = inputParser;
addParameter(p,'channels',1:size(lfp.data,2),@isvector);
addParameter(p,'samplingRate',1250,@isnumeric);
addParameter(p,'twin',[0.1 0.1],@isnumeric);
addParameter(p,'spat_sm',11,@isnumeric);
addParameter(p,'temp_sm',11,@isnumeric);
addParameter(p,'doDetrend',false,@islogical);
addParameter(p,'plotCSD',true,@islogical);
addParameter(p,'plotLFP',true,@islogical);
addParameter(p,'NumStd',1.5,@isnumeric);
addParameter(p,'fRange',[6 12],@isnumeric);
addParameter(p,'thetaCh',ceil(size(lfp.data,2)/2),@isnumeric);
addParameter(p,'invert',false,@islogical);


parse(p,varargin{:});
channels = p.Results.channels;
samplingRate = p.Results.samplingRate;
spat_sm = p.Results.spat_sm;
temp_sm = p.Results.temp_sm;
doDetrend = p.Results.doDetrend;
plotCSD = p.Results.plotCSD;
plotLFP = p.Results.plotLFP;
NumStd = p.Results.NumStd;
fRange = p.Results.fRange;
thetaCh = p.Results.thetaCh;
invert = p.Results.invert;

%lfp input
if isstruct(lfp)
    data = lfp.data;
    if invert
        data = -1*data;
    end
    timestamps = lfp.timestamps;
    samplingRate = lfp.samplingRate;
    channels = lfp.channels;
elseif iscell(lfp) %for multiple trials
    celllengths = cellfun(@length,lfp);
    data = vertcat(lfp{:});
elseif isnumeric(lfp)
    data = lfp;
    timestamps = [1:length(lfp)]'./samplingRate;
end

twin = p.Results.twin*samplingRate;
events = round(events*samplingRate);

%% Filter the LFP
for ii = 1:length(channels)
    data_filt(:,ii) = bandpassFilterCSD(double(data(:,ii)),samplingRate,fRange);
    data_phase(:,ii) = instPhaseCSD(double(data_filt(:,ii)));
    if channels(ii) == thetaCh
        %Take the absolute of the filtered signal and calculate the standard
        %deviation
        rmssig  = abs(data_filt(:,ii));
        stdsig = std(rmssig);
        %Find the index of peaks that are more than NumStd standard deviations 
        peakIdx = find(data_filt(:,ii)>(NumStd*stdsig) & [0; diff(data_phase(:,ii)>0)]);             
    end 
end

%events = events((events + (5*1250) <= size(data,1)) & (events - (5*1250) > 0));

% Now, only select peaks within event timestamps
idx = [];
for ii = 1:length(peakIdx)
    for jj = 1:size(events,2)
        if peakIdx(ii) >= (events(1,jj)) && peakIdx(ii) <= (events(2,jj))
            idx = [idx; peakIdx(ii)];
        end
    end
end

%% Conpute event-triggered LFP average
% events = events((events + twin(2) <= size(data,1)) & (events - twin(1) > 0));
lfp = nan(twin(1)+twin(2)+1,length(channels),length(idx));
% 
for e = 1:length(idx)
    lfp(:,:,e) = data(idx(e)-twin(1):idx(e)+twin(2),1:length(channels));
end
 
lfp_avg = nanmean(lfp,3)*-1;
num = length(idx);

%% Compute CSD

% detrend
if doDetrend
   lfp_avg = detrend(lfp_avg')';
end
    
% temporal smoothing
if temp_sm > 0
   for ch = 1:size(lfp_avg,2) 
       lfp_avg(:,ch) = smooth(lfp_avg(:,ch),temp_sm,'sgolay');
   end
end

% spatial smoothing
if spat_sm > 0
   for t = 1:size(lfp_avg,1) 
       lfp_avg(t,:) = smooth(lfp_avg(t,:),spat_sm,'lowess');
   end
end

% calculate CSD 
CSD = diff(lfp_avg,2,2);

% generate output structure
csd.data = CSD;
csd.timestamps = -twin(1):twin(2);
csd.samplingRate = samplingRate;
csd.channels = channels; 
csd.params.spat_sm = spat_sm;
csd.params.temp_sm = temp_sm;
csd.params.detrend = doDetrend;
csd.num = num;

lfpAvg.data = lfp_avg;
lfpAvg.timestamps = -twin(1):twin(2);
lfpAvg.samplingRate = samplingRate;
lfpAvg.channels = channels; 
lfpAvg.params.spat_sm = spat_sm;
lfpAvg.params.temp_sm = temp_sm;
lfpAvg.params.detrend = doDetrend;

%% Plot

if plotLFP
    
    taxis = (-(twin(1)/samplingRate):(1/samplingRate):(twin(2)/samplingRate))*1e3;
    cmax = max(max(CSD)); 
    figure;
    subplot(1,2,1);
    contourf(taxis,1:size(CSD,2),CSD',40,'LineColor','none');hold on;
    colormap jet; caxis([-cmax cmax]);
    set(gca,'YDir','reverse');xlabel('time (s)');ylabel('channel');title('CSD'); 
    plot([0 0],[1 size(CSD,2)],'--k');hold on;
    
    subplot(1,2,2);
    for ch=1:size(lfp_avg,2)
        offset = 400*(ch-1);
        sh_tmp = 1e0*(lfp_avg(:,ch)) + offset;
        plot(taxis,sh_tmp,'k','LineWidth',1.5); hold on;
        clear sh_tmp
    end
    set(gca,'YDir','reverse','YTickLabel',[]);ylim([-1000 offset+1000]);xlim([taxis(1) taxis(end)]);
    xlabel('time (ms)');ylabel('channel');title('LFP');   
    plot([0 0],ylim,'--r');hold on;

       
elseif plotCSD  
    
     cmax = max(max(CSD)); 
   
     figure;
     contourf(taxis,1:size(CSD,2),CSD',40,'LineColor','none');hold on;
     colormap jet; caxis([-cmax cmax]);
     set(gca,'YDir','reverse','YTickLabel',[]);ylim([-1000 offset+1000]);xlim([taxis(1) taxis(end)]);
     plot([0 0],[1 size(CSD,2)],'--k');hold on;
   
end

end
