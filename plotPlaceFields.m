function plotPlaceFields(varargin)


%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'dimension',2,@isnumeric)
parse(p,varargin{:});
basepath = p.Results.basepath;
dimension = p.Results.dimension;

%% Deal with inputs
if ~isempty(dir([basepath filesep '*.SessionPulses.Events.mat']))
    disp('Session pulses already detected! Loading file.');
    file = dir([basepath filesep '*.SessionPulses.Events.mat']);
    load(file.name);
end

if ~isempty(dir([basepath filesep '*.Behavior.mat'])) 
    disp('Behavior already detected! Loading file.');
    file = dir([basepath filesep '*.Behavior.mat']);
    load(file(1).name);
end

if ~isempty(dir([basepath filesep '*.spikes.cellinfo.mat']))
    disp('Spikes already detected! Loading file.');
    file = dir([basepath filesep '*.spikes.cellinfo.mat']);
    load(file.name);
end

%% Find subfolder recordings
flds = fieldnames(sessionPulses);

for ii = 1:size(flds,1)
    cd([basepath filesep flds{ii}]);
    fprintf('Computing place fields in %s folder \n',flds{ii});
    
    %Get OFF timestamps
    Offidx = find(sessionPulses.(flds{ii}).stim==0);
    trialsIdx = (behavior.masks.trials == Offidx);
    trialsIdx = sum(trialsIdx,2);
    
    %Get ON timestamps
    Onidx = find(sessionPulses.(flds{ii}).stim==1);
    trialsIdxOn = (behavior.masks.trials == Onidx);
    trialsIdxOn = sum(trialsIdxOn,2);    
    
    %OFF trials
    xOFF_left = behavior.position.x(behavior.masks.arm==0 & trialsIdx & behavior.masks.recording==ii);
    xOFF_right = behavior.position.x(behavior.masks.arm==1 & trialsIdx & behavior.masks.recording==ii);
    
    yOFF_left = behavior.position.y(behavior.masks.arm==0 & trialsIdx & behavior.masks.recording==ii);
    yOFF_right = behavior.position.y(behavior.masks.arm==1 & trialsIdx & behavior.masks.recording==ii);

    linOFF_left = behavior.position.lin(behavior.masks.arm==0 & trialsIdx & behavior.masks.recording==ii);
    linOFF_right = behavior.position.lin(behavior.masks.arm==1 & trialsIdx & behavior.masks.recording==ii);
    
    tsOFF_left = behavior.timestamps(behavior.masks.arm==0 & trialsIdx & behavior.masks.recording==ii);
    tsOFF_right = behavior.timestamps(behavior.masks.arm==1 & trialsIdx & behavior.masks.recording==ii);
    
    %ON trials
    xON_left = behavior.position.x(behavior.masks.arm==0 & trialsIdxOn & behavior.masks.recording==ii);
    xON_right = behavior.position.x(behavior.masks.arm==1 & trialsIdxOn & behavior.masks.recording==ii);
    
    yON_left = behavior.position.y(behavior.masks.arm==0 & trialsIdxOn & behavior.masks.recording==ii);
    yON_right = behavior.position.y(behavior.masks.arm==1 & trialsIdxOn & behavior.masks.recording==ii);

    linON_left = behavior.position.lin(behavior.masks.arm==0 & trialsIdxOn & behavior.masks.recording==ii);
    linON_right = behavior.position.lin(behavior.masks.arm==1 & trialsIdxOn & behavior.masks.recording==ii);
    
    tsON_left = behavior.timestamps(behavior.masks.arm==0 & trialsIdxOn & behavior.masks.recording==ii);
    tsON_right = behavior.timestamps(behavior.masks.arm==1 & trialsIdxOn & behavior.masks.recording==ii);

    %Populate position matrix
    if dimension == 1 % If linear maps
        positions{1} = [tsOFF_left linOFF_left];
        positions{2} = [tsON_left linON_left];        
        positions{3} = [tsOFF_right linOFF_right];
        positions{4} = [tsON_right linON_right];
        
    elseif dimension ==2
        positions{1} = [tsOFF_left xOFF_left yOFF_left];
        positions{2} = [tsON_left xON_left yON_left];
        positions{3} = [tsOFF_right xOFF_right yOFF_right];
        positions{4} = [tsON_right xON_right yON_right];
    end
    
    firingMaps.(flds{ii}) = bz_firingMapAvg_IZ(positions,spikes,'plotFig',true,'saveMat',true);

end

end