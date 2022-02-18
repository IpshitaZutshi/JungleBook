function SessTrialsPlaceFields(varargin)

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
        %Get OFF correct timestamps
        Offidx = find(sessionPulses.(flds{ii}).stim==0 & sessionArmChoice.(flds{ii}).choice(1:length(sessionPulses.(flds{ii}).stim))'==1);
        trialsIdx_first = (behavior.masks.trials == Offidx(1:floor(length(Offidx)/2)));
        trialsIdx_first = sum(trialsIdx_first,2);
        trialsIdx_last = (behavior.masks.trials == Offidx(floor(length(Offidx)/2)+1:end));
        trialsIdx_last = sum(trialsIdx_last,2);
        
        %Get OFF incorrect timestamps
        OffidxIN = find(sessionPulses.(flds{ii}).stim==0 & sessionArmChoice.(flds{ii}).choice(1:length(sessionPulses.(flds{ii}).stim))'==0);
        if ~isempty(OffidxIN)
            trialsIdxIN_first = (behavior.masks.trials == OffidxIN(1:floor(length(OffidxIN)/2)));
            trialsIdxIN_first = sum(trialsIdxIN_first,2);
            trialsIdxIN_last = (behavior.masks.trials == OffidxIN(floor(length(OffidxIN)/2)+1:end));
            trialsIdxIN_last = sum(trialsIdxIN_last,2);
        else
            trialsIdxIN_first = zeros(length(behavior.masks.trials),1);
            trialsIdxIN_last = zeros(length(behavior.masks.trials),1);
        end

        %Get ON correct timestamps
        Onidx = find(sessionPulses.(flds{ii}).stim==1 & sessionArmChoice.(flds{ii}).choice(1:length(sessionPulses.(flds{ii}).stim))'==1);
        trialsIdxOn_first = (behavior.masks.trials == Onidx(1:floor(length(Onidx)/2)));
        trialsIdxOn_first = sum(trialsIdxOn_first,2);
        trialsIdxOn_last = (behavior.masks.trials == Onidx(floor(length(Onidx)/2)+1:end));
        trialsIdxOn_last = sum(trialsIdxOn_last,2);  

        %Get ON incorrect timestamps
        OnidxIN = find(sessionPulses.(flds{ii}).stim==1 & sessionArmChoice.(flds{ii}).choice(1:length(sessionPulses.(flds{ii}).stim))'==0);
        if ~isempty(OnidxIN)
            trialsIdxOnIN_first = (behavior.masks.trials == OnidxIN(1:floor(length(OnidxIN)/2)));
            trialsIdxOnIN_first = sum(trialsIdxOnIN_first,2);
            trialsIdxOnIN_last = (behavior.masks.trials == OnidxIN(floor(length(OnidxIN)/2)+1:end));
            trialsIdxOnIN_last = sum(trialsIdxOnIN_last,2);
        else
            trialsIdxOnIN_first = zeros(length(behavior.masks.trials),1);
            trialsIdxOnIN_last = zeros(length(behavior.masks.trials),1);
        end

        %OFF correct
        linOFF_left_first = behavior.position.lin(behavior.masks.arm==0 & trialsIdx_first & behavior.masks.recording==ii);
        linOFF_left_last = behavior.position.lin(behavior.masks.arm==0 & trialsIdx_last & behavior.masks.recording==ii);
        linOFF_right_first = behavior.position.lin(behavior.masks.arm==1 & trialsIdx_first & behavior.masks.recording==ii);
        linOFF_right_last = behavior.position.lin(behavior.masks.arm==1 & trialsIdx_last & behavior.masks.recording==ii);

        tsOFF_left_first = behavior.timestamps(behavior.masks.arm==0 & trialsIdx_first & behavior.masks.recording==ii);
        tsOFF_left_last = behavior.timestamps(behavior.masks.arm==0 & trialsIdx_last & behavior.masks.recording==ii);
        tsOFF_right_first = behavior.timestamps(behavior.masks.arm==1 & trialsIdx_first & behavior.masks.recording==ii);
        tsOFF_right_last = behavior.timestamps(behavior.masks.arm==1 & trialsIdx_last & behavior.masks.recording==ii);

        %OFF incorrect trials
        linOFF_leftIN_first = behavior.position.lin(behavior.masks.arm==0 & trialsIdxIN_first & behavior.masks.recording==ii);
        linOFF_leftIN_last = behavior.position.lin(behavior.masks.arm==0 & trialsIdxIN_last & behavior.masks.recording==ii);
        linOFF_rightIN_first = behavior.position.lin(behavior.masks.arm==1 & trialsIdxIN_first & behavior.masks.recording==ii);
        linOFF_rightIN_last = behavior.position.lin(behavior.masks.arm==1 & trialsIdxIN_last & behavior.masks.recording==ii);

        tsOFF_leftIN_first = behavior.timestamps(behavior.masks.arm==0 & trialsIdxIN_first & behavior.masks.recording==ii);
        tsOFF_leftIN_last = behavior.timestamps(behavior.masks.arm==0 & trialsIdxIN_last & behavior.masks.recording==ii);
        tsOFF_rightIN_first = behavior.timestamps(behavior.masks.arm==1 & trialsIdxIN_first & behavior.masks.recording==ii);
        tsOFF_rightIN_last = behavior.timestamps(behavior.masks.arm==1 & trialsIdxIN_last & behavior.masks.recording==ii);

        %ON correct    
        linON_left_first = behavior.position.lin(behavior.masks.arm==0 & trialsIdxOn_first & behavior.masks.recording==ii);
        linON_left_last = behavior.position.lin(behavior.masks.arm==0 & trialsIdxOn_last & behavior.masks.recording==ii);
        linON_right_first = behavior.position.lin(behavior.masks.arm==1 & trialsIdxOn_first & behavior.masks.recording==ii);
        linON_right_last = behavior.position.lin(behavior.masks.arm==1 & trialsIdxOn_last & behavior.masks.recording==ii);

        tsON_left_first = behavior.timestamps(behavior.masks.arm==0 & trialsIdxOn_first & behavior.masks.recording==ii);
        tsON_left_last = behavior.timestamps(behavior.masks.arm==0 & trialsIdxOn_last & behavior.masks.recording==ii);
        tsON_right_first = behavior.timestamps(behavior.masks.arm==1 & trialsIdxOn_first & behavior.masks.recording==ii);
        tsON_right_last = behavior.timestamps(behavior.masks.arm==1 & trialsIdxOn_last & behavior.masks.recording==ii);

        %ON incorrect trials
        linON_leftIN_first = behavior.position.lin(behavior.masks.arm==0 & trialsIdxOnIN_first & behavior.masks.recording==ii);
        linON_leftIN_last = behavior.position.lin(behavior.masks.arm==0 & trialsIdxOnIN_last & behavior.masks.recording==ii);
        linON_rightIN_first = behavior.position.lin(behavior.masks.arm==1 & trialsIdxOnIN_first & behavior.masks.recording==ii);
        linON_rightIN_last = behavior.position.lin(behavior.masks.arm==1 & trialsIdxOnIN_last & behavior.masks.recording==ii);

        tsON_leftIN_first = behavior.timestamps(behavior.masks.arm==0 & trialsIdxOnIN_first & behavior.masks.recording==ii);
        tsON_leftIN_last = behavior.timestamps(behavior.masks.arm==0 & trialsIdxOnIN_last & behavior.masks.recording==ii);
        tsON_rightIN_first = behavior.timestamps(behavior.masks.arm==1 & trialsIdxOnIN_first & behavior.masks.recording==ii);
        tsON_rightIN_last = behavior.timestamps(behavior.masks.arm==1 & trialsIdxOnIN_last & behavior.masks.recording==ii);
        


        %Populate position matrix
        if dimension == 1 % If linear maps
            positions{1} = [tsOFF_left_first linOFF_left_first];
            positions{2} = [tsOFF_left_last linOFF_left_last];       
            positions{3} = [tsON_left_first linON_left_first];
            positions{4} = [tsON_left_last linON_left_last];
            positions{5} = [tsOFF_right_first linOFF_right_first];
            positions{6} = [tsOFF_right_last linOFF_right_last];       
            positions{7} = [tsON_right_first linON_right_first];
            positions{8} = [tsON_right_last linON_right_last];
            positions{9} = [tsOFF_leftIN_first linOFF_leftIN_first];
            positions{10} = [tsOFF_leftIN_last linOFF_leftIN_last];       
            positions{11} = [tsON_leftIN_first linON_leftIN_first];
            positions{12} = [tsON_leftIN_last linON_leftIN_last];
            positions{13} = [tsOFF_rightIN_first linOFF_rightIN_first];
            positions{14} = [tsOFF_rightIN_last linOFF_rightIN_last];       
            positions{15} = [tsON_rightIN_first linON_rightIN_first];
            positions{16} = [tsON_rightIN_last linON_rightIN_last];
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
            sessCorrMaps = bz_firingMapAvg_IZ(positions,spikes,'plotFig',false,'saveMat',false);        
            firingMaps.sessCorrMaps = sessCorrMaps.rateMaps;
            save([firingMaps.sessionName '.firingMapsAvg.cellinfo.mat'],'firingMaps');             
        else
            sessCorrMaps = bz_firingMapAvg_IZ(positions,DSspikes,'plotFig',false,'saveMat',false);        
            firingMaps.sessCorrMaps = sessCorrMaps.rateMaps;
            save([firingMaps.sessionName '_DS.firingMapsAvg.cellinfo.mat'],'firingMaps');   
            clear DSspikes
        end
    end
end
end