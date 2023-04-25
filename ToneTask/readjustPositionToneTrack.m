function readjustPositionToneTrack(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);

parse(p,varargin{:});
basepath = p.Results.basepath;

%% Look for the tracking and behavior files
if ~isempty(dir([basepath filesep '*.Tracking.Behavior.mat'])) 
    disp('Loading tracking');
    file = dir([basepath filesep '*.Tracking.Behavior.mat']);
    load(file(1).name);
end

if ~isempty(dir([basepath filesep '*TrialBehavior.Behavior.mat'])) 
    disp('Behavior already detected! Loading file.');
    file = dir([basepath filesep '*TrialBehavior.Behavior.mat']);
    load(file(1).name);
end
[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);
%% Now the track is 122 cm, each reward port is separated by 20.32 cm
%Find current range -
start = min(tracking.position.y);
tracking.position.y  = tracking.position.y-start;
if start>0
    tracking.position.y = (tracking.position.y*(122/(120-start)));
end

save([basepath filesep sessionInfo.FileName '.Tracking.Behavior.mat'],'tracking');
end