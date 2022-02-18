function getPhasePrecession(varargin)


%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
parse(p,varargin{:});
basepath = p.Results.basepath;

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

if ~isempty(dir([basepath filesep '*.SessionArmChoice.Events.mat'])) 
    disp('Arm Choice already detected! Loading file.');
    file = dir([basepath filesep '*.SessionArmChoice.Events.mat']);
    load(file(1).name);
end

file = dir(('*.hippocampalLayers.channelinfo.mat'));
load(file.name);  

file = dir(('*.session.mat'));
load(file.name); 

%% Get lfp
pyrCh = hippocampalLayers.oriens; 
if ~isfield(session.channelTags,'RippleNoise')
    lfp = bz_GetLFP(pyrCh,'noPrompts', true);
else
    refChannel = session.channelTags.RippleNoise.channels-1;
    lfp = bz_GetLFP([pyrCh refChannel],'noPrompts', true);
    lfp = bz_interpolateLFP(lfp,'refChan',refChannel);
end

%% Find subfolder recordings
flds = fieldnames(sessionPulses);

for ii = 1:size(flds,1)
    cd([basepath filesep flds{ii}]);
    
    if ~isempty('*.placeFields.cellinfo.mat')
        fprintf('Loading place fields & computing phase precession in %s folder \n',flds{ii});
        file = dir('*.placeFields.cellinfo.mat');
        load(file(1).name);
    else
        fprintf('First calculate place field stats!');
        return
    end

    %Get OFF correct timestamps
     Offidx = find(sessionPulses.(flds{ii}).stim==0 & sessionArmChoice.(flds{ii}).choice(1:length(sessionPulses.(flds{ii}).stim))'==1);
    trialsIdx = (behavior.masks.trials == Offidx);
    trialsIdx = sum(trialsIdx,2);
    
    %Get ON correct timestamps
    Onidx = find(sessionPulses.(flds{ii}).stim==1 & sessionArmChoice.(flds{ii}).choice(1:length(sessionPulses.(flds{ii}).stim))'==1);
    trialsIdxOn = (behavior.masks.trials == Onidx);
    trialsIdxOn = sum(trialsIdxOn,2);    
    
    %OFF trials
    linOFF_left = behavior.position.lin(behavior.masks.arm==0 & trialsIdx & behavior.masks.recording==ii);
    linOFF_right = behavior.position.lin(behavior.masks.arm==1 & trialsIdx & behavior.masks.recording==ii);
    
    tsOFF_left = behavior.timestamps(behavior.masks.arm==0 & trialsIdx & behavior.masks.recording==ii);
    tsOFF_right = behavior.timestamps(behavior.masks.arm==1 & trialsIdx & behavior.masks.recording==ii);
    
    %ON trials
    linON_left = behavior.position.lin(behavior.masks.arm==0 & trialsIdxOn & behavior.masks.recording==ii);
    linON_right = behavior.position.lin(behavior.masks.arm==1 & trialsIdxOn & behavior.masks.recording==ii);
    
    tsON_left = behavior.timestamps(behavior.masks.arm==0 & trialsIdxOn & behavior.masks.recording==ii);
    tsON_right = behavior.timestamps(behavior.masks.arm==1 & trialsIdxOn & behavior.masks.recording==ii);

    %Populate position matrix
    positions{1} = [tsOFF_left linOFF_left];
    positions{2} = [tsON_left linON_left];        
    positions{3} = [tsOFF_right linOFF_right];
    positions{4} = [tsON_right linON_right];

    %Correct positions for linearized maze
    for kk = 1:length(positions)
        linMap = positions{kk}(:,2);
        linMapDouble = linMap+170;
        idxMap = (250-linMapDouble)>0; %min([abs(linMap-120),abs(linMapDouble-120)],[],2);
        linMapNew = linMap;
        linMapNew(idxMap) =  linMapDouble(idxMap);           
        linMapNew = linMapNew-75;
        %remove sections when the mouse went backwards
        a = diff(linMapNew);
        idxA= find(a>150);
        idxB= find(a<-150);
        idxtoKeep = ones(length(linMapNew),1);
        %Remove the sections between idxA and idxB
        for pp = 1:length(idxA)
            temp = idxB((idxB-idxA(pp))>0);
            idxEnd = min(temp);
            idxtoKeep(idxA(pp):idxEnd) = 0;
        end
        idxtoKeep = logical(idxtoKeep);
        positions{kk} = [positions{kk}(idxtoKeep,1) linMapNew(idxtoKeep)];
    end
    
    phasePrecession.(flds{ii}) = bz_phasePrecessionAvg(positions,spikes,lfp,'boundaries',placeFieldStats.mapStats,'saveMat',true);

end

end