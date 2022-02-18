function [phasePrecession] = bz_phasePrecessionAvg(positions,spikes,lfp,varargin)

% USAGE
% [data, stats] = bz_phasePrecessionAvg(positions,spikes,lfp,varargin)
% Calculates averaged phase precession
%
% INPUTS
%
%   spikes    - buzcode format .cellinfo. struct with the following fields
%               .times 
%   positions - [t x y ] or [t x] position matrix or
%               cell with several of these matrices (for different conditions)
%   lfp       - a buzcode structure with fields lfp.data,
%                                               lfp.timestamps
%                                               lfp.samplingRate
%      or
%   behavior  - buzcode format behavior struct - NOT YET IMPLEMENTED
%   <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties        Values
%    -------------------------------------------------------------------------
%     'maxGap'          time gaps between successive position samples exceeding
%                       this threshold (e.g. undetects) will not be interpolated
%                       (default = 100 ms)
%     'boundaries'      onset and offset for single-lap phase precession can be
%                       determined either automatically based on spike count
%                       ('count', default) or using explicit firing field
%                       boundaries ([Xstart Xstop], where each X is in [0..1])
%     'slope'           search start value for slope (default = 0)
%     'passband'        Default = [5 12]  
%     'orderKalmanVel'	order of Kalman Velocity Filter (default 2)
%     'speedThresh'		speed threshold to compute firing rate
%     'saveMat'         true  
%    =========================================================================
%
%  OUTPUT
%
%    data.x                position samples
%
%    data.position.ok      spikes occurring at valid coordinates (logical)
%    data.position.t       spike times (only for valid coordinates)
%    data.position.x       x coordinate for each spike
%    data.position.phase   spike phase for each spike (in radians)
%    data.position.lap     lap number for each spike
%
%    data.rate.t           spike times
%    data.rate.r           spike rate for each spike
%    data.rate.phase       spike phase for each spike (in radians)
%    data.rate.lap         lap number for each spike
%
%    Additional statistics computed using phase vs position data:
%
%    stats.slope           phase precession slope (via circular regression)
%    stats.intercept       phase precession intercept (via circular regression)
%    stats.r2              coefficient of determination for circular regression
%    stats.p               p-value for circular regression
%    stats.x               center x for early, middle and late subfields
%    stats.mean            mean phase for early, middle and late subfields
%    stats.var             phase variance for early, middle and late subfields
%    stats.std             phase standard deviation for early, middle and late subfields
%    stats.conf            95% confidence intervals
%    stats.all             cell array containing all phases for each subfield
%                          (useful for population analyses)
%
%    stats.lap.slope       phase precopenession slope (via circular regression)
%    stats.lap.intercept   phase precession intercept (via circular regression)
%    stats.lap.r2          coefficient of determination for circular regression
%    stats.lap.p           p-value for circular regression
%
%    Additional statistics computed using phase vs rate data:
%
%    stats.rate.mean       mean phase for each spike rate
%    stats.rate.var        phase variance for each spike rate
%    stats.rate.conf       phase 95% confidence intervals for each spike rate
%
%    stats.boundaries      field boundaries

%% parse inputs
p=inputParser;
addParameter(p,'maxGap',0.1,@isnumeric);
addParameter(p,'boundaries','count');
addParameter(p,'slopebound',[-1.5 1.5],@isnumeric);
addParameter(p,'passband',[6 12]);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'orderKalmanVel',2,@isnumeric);
addParameter(p,'speedThresh',0.1,@isnumeric);
addParameter(p,'curate_PFs',true,@islogical);

parse(p,varargin{:});
maxGap = p.Results.maxGap;
slopebound = p.Results.slopebound;
boundaries = p.Results.boundaries;
passband = p.Results.passband;
saveMat = p.Results.saveMat;
order = p.Results.orderKalmanVel;
speedThresh = p.Results.speedThresh;
curate_PFs = p.Results.curate_PFs;

% number of conditions
if iscell(positions)
    conditions = length(positions); 
elseif isvector(positions)
    conditions = 1;
end

%%% TODO: conditions label

%% Filter lfp using hilbert transform

[b,a] = butter(3,[passband(1)/(lfp.samplingRate/2) passband(2)/(lfp.samplingRate/2)],'bandpass'); % order 3
filt = FiltFiltM(b,a,double(lfp.data(:,1)));
hilb = hilbert(filt);
lfpphase = mod(angle(hilb),2*pi);
lfpData = [lfp.timestamps lfpphase];

%% Calculate
%Erase positions below speed threshold
for iCond = 1:size(positions,2)
    % Compute speed
    post = positions{iCond}(:,1);
    % - 1D 
    if size(positions{iCond},2)==2
        posx = positions{iCond}(:,2);
        [~,~,~,vx,vy,~,~] = KalmanVel(posx,posx*0,post,order);
    elseif size(positions{iCond},2)==3
        posx = positions{iCond}(:,2);
        posy = positions{iCond}(:,3);
        [~,~,~,vx,vy,~,~] = KalmanVel(posx,posy,post,order);
    else
        warning('This is not a linear nor a 2D space!');
    end
    % Absolute speed
    v = sqrt(vx.^2+vy.^2);
    
    % Compute timestamps where speed is under threshold
    positions{iCond}(:,2) = posx;
    if size(positions{iCond},2) > 2
        positions{iCond}(:,3) = posy;
    end
    positions{iCond}(v<speedThresh,:) = [];
    
end

% Curate place field boundaries
if curate_PFs
    boundaries = bz_curatePlaceFieldStats(boundaries, positions, spikes, 'input_mapstats', true);
end

tic
for unit = 1:length(spikes.times)
    for c = 1:conditions   
        timenow = toc;
        fprintf('%.2d:%.2d  Unit %d, condition %d\n', floor(timenow/60), round(mod(timenow, 60)), unit, c)
        boundaryData = boundaries{unit}{c}.fieldX./size(boundaries{unit}{c}.field,1);
        if isnan(boundaryData(1))
           phasePrecession.data{unit}{c} = nan;
           phasePrecession.stats{unit}{c} = nan;
           continue
        end
        for bd = 1:size(boundaryData,1)
            if (boundaryData(bd,2)-boundaryData(bd,1))>0
            % Extend boundaries by 1 bin on either side to account for
            % extra spikes. 
                if boundaryData(bd,1) ~= 0.02
                    boundaryData(bd,1) = boundaryData(bd,1)-0.02;
                end
                if boundaryData(bd,2) ~= 1
                    boundaryData(bd,2) = boundaryData(bd,2)+0.02;
                end
                [phasePrecession.data{unit}{bd,c}, phasePrecession.stats{unit}{bd,c}] = PhasePrecessionSimplified(positions{c},spikes.times{unit},lfpData,boundaryData(bd,:),'maxGap',maxGap,'abound',slopebound);     
            %    [phasePrecession.data{unit}{bd,c}, phasePrecession.stats{unit}{bd,c}] = PhasePrecession_IZ(positions{c},spikes.times{unit},lfpData,'boundaries',boundaryData(bd,:));  
            end
        end
    end
end


%% restructure into cell info data type

phasePrecession.UID = spikes.UID;
phasePrecession.sessionName = spikes.sessionName;
try
    phasePrecession.region = spikes.region; 
catch
   %warning('spikes.region is missing') 
end

if saveMat
   save([phasePrecession.sessionName '.phasePrecessionAvg.cellinfo.mat'],'phasePrecession'); 
end

end

