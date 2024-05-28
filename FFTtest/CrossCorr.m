function [cor, lags] = CrossCorr(ts1,ts2,binsize)
% [acor, lag] = CMBHOME.Spike.CrossCorr(ts1)
% [xcor, lag] = CMBHOME.Spike.CrossCorr(ts1, ts2, varargin)
%
% Calculates cross correlation (or auto correlation) for vectors ts1 and
% ts2 of spike times. ts1 and ts2 are not binned, or binary, spike trains but
% rather the time that spikes occured.
%
% If ts1 is the only argument, the autocorrelation is produced (also if
% ts1==ts2). In this case, all zeros in the set of latencies will be removed
% before calculating the xcorr histogram.
%
% ARGUMENTS
%   ts1         vector of spike times (seconds)
%   ts2         (optional) vector of spike times
%
% PROPERTIES
%   'binsize'       (default=.01 seconds) binsize for cross correlation
%   'norm'          (default='prob') determines normalization of cross
%                   correlogram. Can be 'prob', 'count', 'unbiased'. To use
%                   the unbiased option, 'epoch' property must be passed.
%   'lag'           (default [-.4 4] seconds) included to speed up algorithm.
%                   this defines the upper and lower limit in the peristumulus
%                   time histogram
%   'suppress_plot' (default 1) 1 with plot, 0 without plot
%   'epoch'         2 element vector indicating start and stop times of recording window
%
% RETURNS
%   cor         a col. vector of cross correlation values corresponding to 'lag'
%   lag         a col. vector of center-aligned lag values (seconds)
%
% alex wiltschko and andrew bogaard


ac = 0;

if isequal(ts1,ts2), ac = 1; end % is autocorr

db = nan(length(ts1), 3);

psth = [];
s1 = 1;
lag = [-0.5 0.5];
spkind = 1;

while spkind <= length(ts1)
   
    s = s1;
    
    while ts2(s) < ts1(spkind)+lag(1) && s < length(ts2)
    
        s = s+1;
    
    end
    
    s1 = s;
    
    f = s;
    
    while ts2(f) <= ts1(spkind)+lag(2) && f < length(ts2)
        
        f = f+1;
        
    end
    
    if ts2(s)<=ts1(spkind)+lag(2)
        db(spkind, :) = [s f-1 ts1(spkind)];
    end
    
    spkind = spkind+1;
    
end

dspk = diff(db(:,1:2), 1, 2);

for i = 0:max(dspk)
    
    where = dspk>=i;
    
    tf = db(where,1)+i;
    
    psth = [psth; ts2(tf)-db(where,3)];

end

if ac, psth(psth==0) = []; end % remove zeros in autocorrelation

lags = lag(1)+binsize/2:binsize:lag(2)+binsize/2;
cor = hist(psth, lags);
 
lags = lags(:);
cor = cor(:);

