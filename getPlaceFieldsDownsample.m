function getPlaceFieldsDownsample(varargin)


%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'dimension',1,@isnumeric)
addParameter(p,'isCA3',false,@islogical)
addParameter(p,'saveDSspikes',true,@islogical)
parse(p,varargin{:});

basepath = p.Results.basepath;
dimension = p.Results.dimension;
isCA3 = p.Results.isCA3;
saveDSspikes = p.Results.saveDSspikes;

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

for ii = 1:1:size(flds,1)
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
    linOFF_left = behavior.position.lin(behavior.masks.arm==0 & trialsIdx & behavior.masks.recording==ii);
    linOFF_right = behavior.position.lin(behavior.masks.arm==1 & trialsIdx & behavior.masks.recording==ii);
    
    tsOFF_left = behavior.timestamps(behavior.masks.arm==0 & trialsIdx & behavior.masks.recording==ii);
    tsOFF_right = behavior.timestamps(behavior.masks.arm==1 & trialsIdx & behavior.masks.recording==ii);
    
    %OFF incorrect trials
    linOFF_leftIN = behavior.position.lin(behavior.masks.arm==0 & trialsIdxIN & behavior.masks.recording==ii);
    linOFF_rightIN = behavior.position.lin(behavior.masks.arm==1 & trialsIdxIN & behavior.masks.recording==ii);
    
    tsOFF_leftIN = behavior.timestamps(behavior.masks.arm==0 & trialsIdxIN & behavior.masks.recording==ii);
    tsOFF_rightIN = behavior.timestamps(behavior.masks.arm==1 & trialsIdxIN & behavior.masks.recording==ii);
    
    
    %ON trials
    linON_left = behavior.position.lin(behavior.masks.arm==0 & trialsIdxOn & behavior.masks.recording==ii);
    linON_right = behavior.position.lin(behavior.masks.arm==1 & trialsIdxOn & behavior.masks.recording==ii);
    
    tsON_left = behavior.timestamps(behavior.masks.arm==0 & trialsIdxOn & behavior.masks.recording==ii);
    tsON_right = behavior.timestamps(behavior.masks.arm==1 & trialsIdxOn & behavior.masks.recording==ii);

    %ON incorrect trials
    linON_leftIN = behavior.position.lin(behavior.masks.arm==0 & trialsIdxOnIN & behavior.masks.recording==ii);
    linON_rightIN = behavior.position.lin(behavior.masks.arm==1 & trialsIdxOnIN & behavior.masks.recording==ii);
    
    tsON_leftIN = behavior.timestamps(behavior.masks.arm==0 & trialsIdxOnIN & behavior.masks.recording==ii);
    tsON_rightIN = behavior.timestamps(behavior.masks.arm==1 & trialsIdxOnIN & behavior.masks.recording==ii);


    %Populate position matrix
    if dimension == 1 % If linear maps
        positions{ii,1} = [tsOFF_left linOFF_left];
        positions{ii,2} = [tsON_left linON_left];        
        positions{ii,3} = [tsOFF_right linOFF_right];
        positions{ii,4} = [tsON_right linON_right];
        positions{ii,5} = [tsOFF_leftIN linOFF_leftIN];
        positions{ii,6} = [tsON_leftIN linON_leftIN];        
        positions{ii,7} = [tsOFF_rightIN linOFF_rightIN];
        positions{ii,8} = [tsON_rightIN linON_rightIN];
    end
    
    %Correct positions for linearized maze
    if dimension == 1
        for kk = 1:size(positions,2)
            linMap = positions{ii,kk}(:,2);
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
            positions{ii,kk} = [positions{ii,kk}(idxtoKeep,1) linMapNew(idxtoKeep)];
        end
    end
    
    % assign position to spikes   
    for posID = 1:4
       	numSpikesSide = [];
        numSpikesCenter = []; 
        for cellnum = 1:size(spikes.times,2)
            spikePos{ii,posID}{cellnum}.pos = [];
            spikePos{ii,posID}{cellnum}.idx = [];
            %Find position at the spiketimes
            s = sort(spikes.times{cellnum});
            tBehav = positions{ii,posID}(:,1); 
            for tSp = 1:length(s)
                [minT,idx] = min(abs(tBehav-s(tSp)));
                if minT<0.05
                    spikePos{ii,posID}{cellnum}.pos  = [spikePos{ii,posID}{cellnum}.pos; positions{ii,posID}(idx,2)];
                    spikePos{ii,posID}{cellnum}.idx  = [spikePos{ii,posID}{cellnum}.idx; tSp];
                end
            end
            numSpikesSide(cellnum) = sum(spikePos{ii,posID}{cellnum}.pos<110);
            numSpikesCenter(cellnum) = sum(spikePos{ii,posID}{cellnum}.pos>=110);
        end                   
        spkCount{ii}.Side(:,posID) = numSpikesSide;
        spkCount{ii}.Center(:,posID) = numSpikesCenter;
    end
                
%             [~, ~, bins1] = histcounts(tBehav,s);
%             [~,y] = unique(bins1);
%             closestIndex = y(2:end);
%             spikesDS{ii,posID}{cellnum} = positions{ii,posID}(closestIndex,2);

end

for ii = 1:1:size(flds,1)
    
    cd([basepath filesep flds{ii}]);
    DSspikes = spikes;
    if ~isCA3              
        %ratematch return arms separately
        minRate_left = min(spkCount{ii}.Side(:,1:2),[],2);
        minRate_right = min(spkCount{ii}.Side(:,3:4),[],2);    
        %ratematch central arm combined for left & right
        spksCenter = min([spkCount{ii}.Center(:,1)+spkCount{ii}.Center(:,3) spkCount{ii}.Center(:,2)+spkCount{ii}.Center(:,4)],[],2);                 
    else
        
        % if CA3, just need to calculate minSpikes across morning and
        % afternoon sessions;
        minRate_left = min([spkCount{1}.Side(:,1) spkCount{1}.Side(:,2) spkCount{2}.Side(:,1) spkCount{2}.Side(:,2)],[],2);
        minRate_right = min([spkCount{1}.Side(:,3) spkCount{1}.Side(:,4) spkCount{2}.Side(:,3) spkCount{2}.Side(:,4)],[],2);   
        %ratematch central arm combined for left & right
        spksCenter = min([spkCount{1}.Center(:,1)+spkCount{1}.Center(:,3) spkCount{1}.Center(:,2)+spkCount{1}.Center(:,4) ...
            spkCount{2}.Center(:,1)+spkCount{2}.Center(:,3) spkCount{2}.Center(:,2)+spkCount{2}.Center(:,4)],[],2);                 
    end
    
    for cellnum = 1:size(spikes.times,2)
        spikeTimes = [];
        for posID = 1:4
            idxRet = find(spikePos{ii,posID}{cellnum}.pos<110);
            if ~isempty(idxRet)
                if posID <3
                    y = randsample(idxRet,minRate_left(cellnum));
                else
                    y = randsample(idxRet,minRate_right(cellnum));
                end
                spikeTimes = [spikeTimes; spikes.times{cellnum}(spikePos{ii,posID}{cellnum}.idx(y))];
            end                        
        end
        for posID = 1:2
            spkPos = [spikePos{ii,posID}{cellnum}.pos;spikePos{ii,posID+2}{cellnum}.pos];
            spkIdx = [spikePos{ii,posID}{cellnum}.idx;spikePos{ii,posID+2}{cellnum}.idx];
            idxCenter = find(spkPos>=110);
            if ~isempty(idxCenter)
                if length(idxCenter)> spksCenter(cellnum)
                    y = randsample(idxCenter,spksCenter(cellnum));
                    spikeTimes = [spikeTimes; spikes.times{cellnum}(spkIdx(y))];
                else
                    spikeTimes = [spikeTimes; spikes.times{cellnum}(spkIdx(idxCenter))];
                end      
            end
        end
        DSspikes.times{cellnum} = sort(spikeTimes);
    end
    
    for kk = 1:8
       DSposn{kk} = positions{ii,kk};
    end
    [firingMaps.(flds{ii})]= bz_firingMapAvg_IZ(DSposn,DSspikes,'downsample',true,'plotFig',false,'saveMat',true);
    
    if saveDSspikes
        save(['DSspikes.cellinfo.mat'],'DSspikes'); 
    end
    clear DSspikes
end

end