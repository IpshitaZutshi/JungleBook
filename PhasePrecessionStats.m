function  [aopt phi0 rho R p lindat phidat] = PhasePrecessionStats(Position, Phase, PositionField)     

% Position - of each spike, provided in bins
% Phase - of each spike, provided in phase 0-360

if ismember([1 116],PositionField) == [1,1], PositionField = [PositionField(PositionField>55) PositionField(PositionField<55)+116]; end
Position_norm = (Position-min(PositionField)) / (max(PositionField)-min(PositionField));
% lindat - location of each spike (normalized 0-1)
lindat = 1-Position_norm;

Phase = Phase-180;
%Let's normalize phase
Phase = degtorad(Phase);
% phidat - phase of each spike
phidat = Phase;

%Take out NaNs
lindat = lindat(~isnan(phidat));
phidat = phidat(~isnan(phidat));

% abound - slope boundaries.  Does not allow alogorithm to fit lines with slopes outside of this range
abound = [-2 2]; % [-1.6 1.6] is suggested by Christian, [-2 2] was suggested by chris
%abound = [-1.6 1.6];

da = 0.1;

[aopt phi0 rho R p] = rccc(lindat, phidat, abound, da);


% output values
% aopt = slope value
% phi0 = y-intercept
