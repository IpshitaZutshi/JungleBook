
function [CCGmean,t_lag,CCGstd,CCGevents,skippedevents ] = EventVsContinousCCG2(data,data_t,events,lagwin,varargin)
% Data should be continuous time series associated with time stamps data_t
% events should be in seconds, as should lagwin
% lagwin = 2 element vector specifying pre and post time surrounding event (make both positive!!) 

% To do
% Save raster of surrounding events?
% Need way to say want a specific timestamps so can average across videos of similar frame rate / ie impose sampling interval 
% Add memory mapping option to IsolateEpochs2
% Eventually write so uses frametimes themselves (interp to pre and post
% values) instead of av diff between frame times - problem: will result in different length segments 

% REQUIRED INPUT
% =========================================================================
%     Properties        Values
% -------------------------------------------------------------------------
% 'data'                continuous time series data (vector)
% 'data_t'              timestamps for data
% 'events'              event times used as 'triggers', in seconds
% 'lagwin'              2 element vector specifying pre an post time
%                       surrounding event (make both positive)

%
%
% OPTIONAL PROPERTY-VALUE INPUT PAIRS
% =========================================================================
%     Properties        Values
% -------------------------------------------------------------------------
% 'cutOff'             Exclude intervals that include data values below or
%                      greater than given standard deviation, used to exclude outliers 
% 'units'               Specify if units are 'samples' or 'time'. Time in
%                       seconds is the default. Use units if wnt to keep size of CCGs
%                       consistent and sampling rate is the same across
%                       vids (true if upsampled).

% RS 2018

%% Parse Inputs

% Defaults
p = inputParser;
addParameter(p,'cutOff',[])
addParameter(p,'units',[])

% User defined
parse(p,varargin{:});
cutOff = p.Results.cutOff;
units = p.Results.units;

%%
%Find frames closest to event of interest 
event_ints = interp1(data_t,data_t,events,'nearest'); % redundant, but keep so stand-alone fx 

%estimate dt, number of samples around events to grab - assumes relatively consistent sampling
if isempty(units)
    data_dt = median(diff(data_t)); %average dt
    t_lag_dt = round(lagwin./data_dt); % estimate lag window in units samples ...
    t_lag = [-t_lag_dt(1):t_lag_dt(2)].*data_dt;
elseif strcmp(units,'samples')
    data_dt = median(diff(data_t)); %average dt
    t_lag_dt = lagwin;
    t_lag = [-t_lag_dt(1):t_lag_dt(2)].*data_dt;
else 
    disp('Input units must be *samples* or *time* or empty') 
end
    
%Get the time points closest to events, and samples in window around them
[eventinds] = find(ismember(data_t,event_ints)); %get INDEX where interp events are member of data timestamps
if size(eventinds,2)>1
    eventinds = eventinds';
end
events_dt(:,1) = bsxfun(@(A,B) A+B,eventinds,t_lag_dt(1).*-1); %events_dt must be column; add/subtract t_lag_dt (works if symmetric)
events_dt(:,2) = bsxfun(@(A,B) A+B,eventinds,t_lag_dt(2)); % in dt units (actual indices into data instead of time)

skippedevents = (isnan(event_ints)); %?

[CCGevents,~] = IsolateEpochs3(data,events_dt,0,1,'includeNaN'); %data ints; written as if sf = 1; wonky work around, fix later
%[CCGevents_time,~] = IsolateEpochs3(data_t,events_dt,0,1,'includeNaN'); %time ints; written as if sf = 1; wonky work around, fix later

CCGevents = single(cat(2,CCGevents{:})); % might be too large 
%CCGevents_time = single(cat(2,CCGevents_time{:})); %%%% ADD IN TIME CHECK
%FOR OUTLIER TRIALS 

%%%%%%%%%% COME BACK TO THIS STD cut-off (exclude noise in data - useful specifically for wf)
if ~isempty(cutOff)
    dataSTD = nanstd(data);
    dataMean = nanmean(data);
    upperCutOff = dataMean + dataSTD*cutOff;
    lowerCutOff = dataMean - dataSTD*cutOff;
    [r,c]=find(CCGevents> upperCutOff | CCGevents < lowerCutOff);
    outlierTrials = unique(c);
    inds=zeros(1,size(CCGevents,2));
    inds(outlierTrials) = 1;
    CCGevents = CCGevents(:,~inds);
end

% Final output
CCGmean = nanmean(CCGevents,2);
CCGstd = nanstd(CCGevents,[],2);

% Find ripples/data in intervals 
% Check manu's code to see how he does this...? 
% Check my code to see how i did this? lollll where to find? 
end

% imagesc(t_lag,1:179,CCGevents');axis xy
% hold on; plot([0 0],get(gca,'ylim'),'w','linewidth',1.2)
