function cscData = formatcscCSD(cscData,varargin)
% FORMATCSC   Format timestamp and sample data and convert to appropriate
% units
%
% CSCDATA = FORMATCSC(CSCDATA) accepts a scalar structure CSCDATA as
% returned by READCSCFILE or a scalar element of the structure array
% returned by READCSCLIST and performs the following formatting operations:
%
%       1) If FS (sampling frequency) exceeds 2034 Hz, TIMESTAMPS and SAMPLE
%       are subsampled to 2034Hz
%
%       2) TIMESTAMPS are extrapolated for all samples and converted into
%       seconds
%
%       3) All R SAMPLE records are concatenated (reshaped) into a single
%       column vector and converted from bits to millivolts.
%
%       4) Fields FS and CHANNEL are converted from R vectors to scalars.
%
%       5) The units of TIMESTAMPS and SAMPLE are saved in new fields
%       TUNITS and SUNITS, respectively.
%
% CSCDATA = FORMATCSC(...,'SUNITS',SUNITNAME) allows specification of the
% sample units. Valid values of SUNITNAME include:
%
%       - volts:            'V','Volts'
%       - millivolts:       'mV','millivolts'
%       - microvolts:       'uV','microvolts'
%       - bits/raw:         'bits','raw'
%
% CSCDATA = FORMATCSC(...,'TUNITS',TUNITNAME) allows specification of the
% timestamp units. Valid values of TUNITNAME include:
%
%       - seconds:          's','sec','seconds','second'
%       - milliseonds:      'ms','msec','milliseconds','millisecond'
%       - microseconds:     'us','usec','microsec','microseconds','microsecond'
%       - raw:              'raw'  (Currently microseconds, November 2014)
%
% See also: READCSCFILE, READCSCLIST, GETCSCDATA

if ~isempty(varargin)
    [~,props] = parseparams(varargin);
    assert(isequal(varargin,props) && mod(length(props),2)==0,...
        'Conversion units must be entered as name value pairs ');
    assert(length(props)<=4,'Too many input arguments');
    
    validTUnits = {...
        's','sec','seconds','second',...
        'ms','msec','milliseconds','millisecond',...
        'us','usec','microsec','microseconds','microsecond','raw'};
    validSUnits = {'v','volts','mv','millivolts','uv','microvolts','raw'};
    
    P = inputParser;
    addParamValue(P,'sUnits','mv',@(u)any(validatestring(lower(u),validSUnits)));
    addParamValue(P,'tUnits','s',@(u)any(validatestring(lower(u),validTUnits)));
    parse(P,varargin{:});
    sUnits = P.Results.sUnits;
    tUnits = P.Results.tUnits;         
else
    sUnits = 'mV';
    tUnits = 's';
end
% Convert sample to sUnits and convert to vector
switch lower(sUnits)
    case {'v','volts'}
        sampConversion = cscData.bits2volts;
        cscData.sUnits = 'V';
    case {'mv','millivolts'}
        sampConversion = cscData.bits2volts*1000;
        cscData.sUnits = 'mV';
    case {'uv','microvolts'}
        sampConversion = cscData.bits2volts*1000000;
        cscData.sUnits = 'uV';
    case {'raw','bits'}
        sampConversion = 1;
        cscData.sUnits = 'bits';
end
cscData.sample = reshape(cscData.sample*sampConversion,numel(cscData.sample),1);

% Number of raw time units per second
% %%%%%%%%%% CONSTANT, NOT A PARAMETER %%%%%%%%%%
RAW_UNITS_PER_SECOND = 1000000; % Microseconds (November 2014)
% MUST BE UPDATED IF AND ONLY IF RAW UNITS CHANGE
switch lower(tUnits)
    case {'s','sec','seconds','second'}
        tConversion = 1/RAW_UNITS_PER_SECOND;
        cscData.tUnits = 's';
    case {'ms','msec','milliseconds','millisecond'}
        tConversion = 1000/RAW_UNITS_PER_SECOND;
        cscData.tUnits = 'ms';
    case {'us','usec','microsec','microseconds','microsecond'}
        tConversion = 1000000/RAW_UNITS_PER_SECOND;
        cscData.tUnits = 'us';
    case {'raw'} % Equivalent to microseconds (November 2014)
        tConversion = 1;
        cscData.tUnits = 'us'; % November 2014
end

% Find all timestamps instead of every 512th and convert to tUnits
allts = arrayfun(@(ts2,ts1,n)linspace(ts1,ts2,(n+1)),...
    [cscData.timestamps(2:end) 2*cscData.timestamps(end)-cscData.timestamps(end-1)],...
    cscData.timestamps,512*ones(1,length(cscData.timestamps)),'UniformOutput',0);
cscData.timestamps = cell2mat(cellfun(@(x)x(1:end-1),allts,'UniformOutput',0))'*tConversion;

% If data wasn't subsampled it needs to be
if cscData.Fs(1)>2034
    rateCorrection = round(cscData.Fs(1)/2034);
    cscData.sample = cscData.sample(1:rateCorrection:end);
    cscData.timestamps = cscData.timestamps(1:rateCorrection:end);
    cscData.Fs = cscData.Fs/rateCorrection;
end
cscData.Fs = median(cscData.Fs);

cscData.ch = cscData.ch(1);

cscData = rmfield(cscData,'nValidSamp');