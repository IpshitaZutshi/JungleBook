function [data,stats] = PhasePrecessionSimplified(positions,spikes,phases,boundaries,varargin)


%% Defaults and Parms
p = inputParser;
addParameter(p,'abound',[2 2]);
addParameter(p,'da',0.1);
addParameter(p,'maxGap',0.1,@isnumeric);

parse(p,varargin{:});

abound = p.Results.abound;
da = p.Results.da;
maxGap = p.Results.maxGap;


%% Default values
data.x = positions;
data.lindat = [];
data.phidat = [];
data.position.t = [];
data.position.x = [];
data.position.phase = [];
data.rate.t = [];
data.rate.r = [];
data.rate.phase = [];

stats.slope = nan;
stats.intercept = nan;
stats.r2 = nan;
stats.p = nan;

%% Boundary definitions
stats.boundaries = boundaries;

% Firing curve boundaries
fieldStart = boundaries(1);
fieldStop = boundaries(2);
fieldSize = abs(fieldStart-fieldStop);


%% Compute spike phases
if isempty(spikes), return; end
spikePhases = Interpolate(phases,spikes,'trim','off','type','circular');
if isempty(spikePhases), return; end


%% Interpolate positions at spike times
if ~isempty(positions)
	% Make sure positions are normalized
	if max(positions(:,2)) > 1 || min(positions(:,2)) < 0
		positions(:,2) = ZeroToOne(positions(:,2));
		warning('Parameter ''positions'' should contain values in [0 1]. The data will now be transformed accordingly.');
	end
	[x,ignored] = Interpolate(positions,spikes,'trim','off','maxGap',maxGap);
	data.position.ok = ~ignored;
	data.position.t = spikes(~ignored);
	data.position.x = x(:,2);
	data.position.phase = spikePhases(~ignored,2);
    
    % Get spikes and phases within the field
	x = data.position.x;    
    ok = ~isnan(x) & x>=fieldStart & x <=fieldStop;
    lindat = (x(ok)-fieldStart)./fieldSize;
    data.lindat = lindat;
    phidat = data.position.phase(ok); %In radians
    phasedeg = rad2deg(phidat);
    phasedeg = phasedeg-180;
    phidat = deg2rad(phasedeg);
    data.phidat = phidat;
end

%% Compute unwrapped spike phases
unwrapped = [phases(:,1) unwrap(phases(:,2))];
[spikeUnwrappedPhases,~] = Interpolate(unwrapped,spikes,'trim','off');
% Get beginning and end unwrapped phases of two cycles surrounding each spike
startUnwrappedPhases = spikeUnwrappedPhases(:,2)-2*pi;
stopUnwrappedPhases = spikeUnwrappedPhases(:,2)+2*pi;

% Compute instantaneous rate at spike times
data.rate.r = CountInIntervals(spikeUnwrappedPhases(:,2),[startUnwrappedPhases stopUnwrappedPhases]);
data.rate.t = spikes;
data.rate.phase = spikePhases(:,2);


%% Perform the actual computation
if ~isempty(positions)
    
    [aopt, phi0, ~, R, p] = rccc(lindat, phidat, abound, da); 
	stats.slope = aopt;
	stats.intercept = phi0;
	stats.r2 = R;
	stats.p = p;
end

