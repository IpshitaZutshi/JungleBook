function SessTriallByTrialMaps(varargin)

%%Firingmap must already be calculated, merely adds a field in that
%%structure!
%% Defaults and Parms
p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'dimension',1,@isnumeric)
addParameter(p,'downsample',false,@islogical)
parse(p,varargin{:});

expPath = p.Results.expPath;
dimension = p.Results.dimension;
downsample = p.Results.downsample;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');

%% Start collecting data
for ss = 1:size(allSess,1)
    cd(strcat(allSess(ss).folder,'\',allSess(ss).name));

    %% Deal with inputs
    file = dir('*.SessionPulses.Events.mat');
    load(file.name);
    
    file = dir('*.SessionArmChoice.Events.mat') ;
    load(file(1).name);

    file = dir('*.Behavior.mat'); 
    load(file(1).name);
    
    if ~downsample
        file = dir('*.spikes.cellinfo.mat');
        load(file.name);
    end

    %% Find subfolder recordings
    flds = fieldnames(sessionPulses);

    for ii = 1:size(flds,1)
        cd(strcat(allSess(ss).folder,'\',allSess(ss).name,'\',flds{ii}));
        fprintf('Updating place fields in %s folder \n',flds{ii});

        if downsample
            load('DSspikes.cellinfo.mat')
            file = dir('*_DS.firingMapsAvg.cellinfo.mat');
            disp('Loading firing map');
            load(file.name);           
        else
            file = dir('*.firingMapsAvg.cellinfo.mat');
            disp('Loading firing map');
            load(file.name);
        end
       
        %Populate position matrix
        if dimension == 1 % If linear maps
            for pf = 1:size(sessionPulses.(flds{ii}).stim,2)
                positions{pf} = [behavior.timestamps(behavior.masks.trials==pf & behavior.masks.recording==ii) ...
                    behavior.position.lin(behavior.masks.trials==pf & behavior.masks.recording==ii)];
            end     
        end

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

        if ~downsample
            trialMaps = bz_firingMapAvg_IZ(positions,spikes,'minTime',0,'plotFig',false,'saveMat',false);        
            firingMaps.trialMaps = trialMaps.rateMaps;
            firingMaps.stim = sessionPulses.(flds{ii}).stim;
            firingMaps.choice = sessionArmChoice.(flds{ii}).choice;
            firingMaps.arm = sessionArmChoice.(flds{ii}).visitedArm;
            save([firingMaps.sessionName '.firingMapsAvg.cellinfo.mat'],'firingMaps'); 
        else
            trialMaps = bz_firingMapAvg_IZ(positions,DSspikes,'minTime',0,'plotFig',false,'saveMat',false);        
            firingMaps.trialMaps = trialMaps.rateMaps;
            firingMaps.stim = sessionPulses.(flds{ii}).stim;
            firingMaps.choice = sessionArmChoice.(flds{ii}).choice;
            firingMaps.arm = sessionArmChoice.(flds{ii}).visitedArm;
            save([firingMaps.sessionName '_DS.firingMapsAvg.cellinfo.mat'],'firingMaps');          
            clear DSspikes
        end
    end
end
end