function getPlaceFields(varargin)


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

if ~isempty(dir([basepath filesep '*.SessionArmChoice.Events.mat'])) 
    disp('Arm Choice already detected! Loading file.');
    file = dir([basepath filesep '*.SessionArmChoice.Events.mat']);
    load(file(1).name);
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
    
    %Get OFF correct timestamps
    Offidx = find(sessionPulses.(flds{ii}).stim==0 & sessionArmChoice.(flds{ii}).choice(1:length(sessionPulses.(flds{ii}).stim))'==1);
    trialsIdx = (behavior.masks.trials == Offidx);
    trialsIdx = sum(trialsIdx,2);
 
    %Get OFF incorrect timestamps
    OffidxIN = find(sessionPulses.(flds{ii}).stim==0 & sessionArmChoice.(flds{ii}).choice(1:length(sessionPulses.(flds{ii}).stim))'==0);
    if ~isempty(OffidxIN)
        trialsIdxIN = (behavior.masks.trials == OffidxIN);
        trialsIdxIN = sum(trialsIdxIN,2);
    else
        trialsIdxIN = zeros(length(behavior.masks.trials),1);
    end
    
    %Get ON correct timestamps
    Onidx = find(sessionPulses.(flds{ii}).stim==1 & sessionArmChoice.(flds{ii}).choice(1:length(sessionPulses.(flds{ii}).stim))'==1);
    trialsIdxOn = (behavior.masks.trials == Onidx);
    trialsIdxOn = sum(trialsIdxOn,2);    

    %Get ON incorrect timestamps
    OnidxIN = find(sessionPulses.(flds{ii}).stim==1 & sessionArmChoice.(flds{ii}).choice(1:length(sessionPulses.(flds{ii}).stim))'==0);
    if ~isempty(OnidxIN)
        trialsIdxOnIN = (behavior.masks.trials == OnidxIN);
        trialsIdxOnIN = sum(trialsIdxOnIN,2);  
    else
        trialsIdxOnIN = zeros(length(behavior.masks.trials),1);
    end
    
    %OFF trials
    xOFF_left = behavior.position.x(behavior.masks.arm==0 & trialsIdx & behavior.masks.recording==ii);
    xOFF_right = behavior.position.x(behavior.masks.arm==1 & trialsIdx & behavior.masks.recording==ii);
    
    yOFF_left = behavior.position.y(behavior.masks.arm==0 & trialsIdx & behavior.masks.recording==ii);
    yOFF_right = behavior.position.y(behavior.masks.arm==1 & trialsIdx & behavior.masks.recording==ii);

    linOFF_left = behavior.position.lin(behavior.masks.arm==0 & trialsIdx & behavior.masks.recording==ii);
    linOFF_right = behavior.position.lin(behavior.masks.arm==1 & trialsIdx & behavior.masks.recording==ii);
    
    tsOFF_left = behavior.timestamps(behavior.masks.arm==0 & trialsIdx & behavior.masks.recording==ii);
    tsOFF_right = behavior.timestamps(behavior.masks.arm==1 & trialsIdx & behavior.masks.recording==ii);
    
    %OFF incorrect trials
    xOFF_leftIN = behavior.position.x(behavior.masks.arm==0 & trialsIdxIN & behavior.masks.recording==ii);
    xOFF_rightIN = behavior.position.x(behavior.masks.arm==1 & trialsIdxIN & behavior.masks.recording==ii);
    
    yOFF_leftIN = behavior.position.y(behavior.masks.arm==0 & trialsIdxIN & behavior.masks.recording==ii);
    yOFF_rightIN = behavior.position.y(behavior.masks.arm==1 & trialsIdxIN & behavior.masks.recording==ii);

    linOFF_leftIN = behavior.position.lin(behavior.masks.arm==0 & trialsIdxIN & behavior.masks.recording==ii);
    linOFF_rightIN = behavior.position.lin(behavior.masks.arm==1 & trialsIdxIN & behavior.masks.recording==ii);
    
    tsOFF_leftIN = behavior.timestamps(behavior.masks.arm==0 & trialsIdxIN & behavior.masks.recording==ii);
    tsOFF_rightIN = behavior.timestamps(behavior.masks.arm==1 & trialsIdxIN & behavior.masks.recording==ii);
    
    
    %ON trials
    xON_left = behavior.position.x(behavior.masks.arm==0 & trialsIdxOn & behavior.masks.recording==ii);
    xON_right = behavior.position.x(behavior.masks.arm==1 & trialsIdxOn & behavior.masks.recording==ii);
    
    yON_left = behavior.position.y(behavior.masks.arm==0 & trialsIdxOn & behavior.masks.recording==ii);
    yON_right = behavior.position.y(behavior.masks.arm==1 & trialsIdxOn & behavior.masks.recording==ii);

    linON_left = behavior.position.lin(behavior.masks.arm==0 & trialsIdxOn & behavior.masks.recording==ii);
    linON_right = behavior.position.lin(behavior.masks.arm==1 & trialsIdxOn & behavior.masks.recording==ii);
    
    tsON_left = behavior.timestamps(behavior.masks.arm==0 & trialsIdxOn & behavior.masks.recording==ii);
    tsON_right = behavior.timestamps(behavior.masks.arm==1 & trialsIdxOn & behavior.masks.recording==ii);

    %ON incorrect trials
    xON_leftIN = behavior.position.x(behavior.masks.arm==0 & trialsIdxOnIN & behavior.masks.recording==ii);
    xON_rightIN = behavior.position.x(behavior.masks.arm==1 & trialsIdxOnIN & behavior.masks.recording==ii);
    
    yON_leftIN = behavior.position.y(behavior.masks.arm==0 & trialsIdxOnIN & behavior.masks.recording==ii);
    yON_rightIN = behavior.position.y(behavior.masks.arm==1 & trialsIdxOnIN & behavior.masks.recording==ii);

    linON_leftIN = behavior.position.lin(behavior.masks.arm==0 & trialsIdxOnIN & behavior.masks.recording==ii);
    linON_rightIN = behavior.position.lin(behavior.masks.arm==1 & trialsIdxOnIN & behavior.masks.recording==ii);
    
    tsON_leftIN = behavior.timestamps(behavior.masks.arm==0 & trialsIdxOnIN & behavior.masks.recording==ii);
    tsON_rightIN = behavior.timestamps(behavior.masks.arm==1 & trialsIdxOnIN & behavior.masks.recording==ii);


    %Populate position matrix
    if dimension == 1 % If linear maps
        positions{1} = [tsOFF_left linOFF_left];
        positions{2} = [tsON_left linON_left];        
        positions{3} = [tsOFF_right linOFF_right];
        positions{4} = [tsON_right linON_right];
        positions{5} = [tsOFF_leftIN linOFF_leftIN];
        positions{6} = [tsON_leftIN linON_leftIN];        
        positions{7} = [tsOFF_rightIN linOFF_rightIN];
        positions{8} = [tsON_rightIN linON_rightIN];
        
    elseif dimension ==2
        positions{1} = [tsOFF_left xOFF_left yOFF_left];
        positions{2} = [tsON_left xON_left yON_left];
        positions{3} = [tsOFF_right xOFF_right yOFF_right];
        positions{4} = [tsON_right xON_right yON_right];
        positions{5} = [tsOFF_leftIN xOFF_leftIN yOFF_leftIN];
        positions{6} = [tsON_leftIN xON_leftIN yON_leftIN];
        positions{7} = [tsOFF_rightIN xOFF_rightIN yOFF_rightIN];
        positions{8} = [tsON_rightIN xON_rightIN yON_rightIN];
    end
    
    %Correct positions for linearized maze
    if dimension == 1
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
    end

    firingMaps.(flds{ii}) = bz_firingMapAvg_IZ(positions,spikes,'plotFig',true,'saveMat',true);

end

end