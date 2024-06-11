function thetaComp = sessionThetaCompression(varargin)

% This function extracts tone cells and place cells during forward tone
% trials. It calculates the CCG in theta timescale and correlates it to the
% distance between field field peaks. 
% For tone cells, it does this in the tone peaks

warning off
%% defaults
p = inputParser;
addParameter(p,'expPath',pwd,@isfolder)
addParameter(p,'plotfig',true,@islogical)

parse(p,varargin{:})
expPath = p.Results.expPath;
plotfig = p.Results.plotfig;

%% Load relevant stuff
cd(expPath)
file = dir(['*.spikes.cellinfo.mat']);
load(file.name);
file = dir(['*.spikeData.cellinfo.mat']);
load(file.name);
file = dir(['*.Tracking.Behavior.mat']);
load(file(1).name);
file = dir(['*TrialBehavior.Behavior.mat']);
load(file.name);
file = dir(['*cell_metrics.cellinfo.mat']);
load(file.name);
file = dir(['*.rateMapsAvgnotLog.cellinfo.mat']);
load(file(1).name);
file = dir(['*.rateMapsTrial.cellinfo.mat']);
trialMap = load(file(1).name);

%% Label cells that are either place or tone
placeCellLog = zeros(spikes.numcells,1);
toneCellLog =  zeros(spikes.numcells,1);
fieldLoc = nan(spikes.numcells,1);
tonefieldLoc = nan(spikes.numcells,1);
fieldEdges= nan(spikes.numcells,2);
tonefieldEdges= nan(spikes.numcells,2);

for ii = 1:spikes.numcells

    if strcmp(cell_metrics.putativeCellType{ii},'Pyramidal Cell') ~= 1 
        continue
    end

    %% Classify cells as place/tone cells and extract the field location
    dataMat = [];
    dataMatTone = [];
    for jj = 2:7    
        dataMat = [dataMat;firingMaps.forward.rateMaps{ii}{jj}];
        a = fillmissing(firingMaps.tone.rateMaps{ii}{jj},'linear'); % tone
        dataMatTone = [dataMatTone;a];
    end           
    spaceMap = nanmean(dataMat,1);
    toneMap = nanmean(dataMatTone,1);        
    [~,idxMax] = max(toneMap); 

    corrSpace = [];
    corrTone = []; 
    for pp = 1:6
       for jj = (pp+1):6     
           a = corrcoef(dataMat(pp,:),dataMat(jj,:),'rows','complete');
               corrSpace = [corrSpace a(1,2)];
           a = corrcoef(dataMatTone(pp,:),dataMatTone(jj,:),'rows','complete');         
           corrTone = [corrTone a(1,2)];
       end
    end     
    spaceCorr = nanmean(corrSpace);  
    toneCorr = nanmean(corrTone);  

    % Detect fields
    Field_Info_space = detectFields(spaceMap);
    if ~isempty(Field_Info_space) 
        [~,idxMax] = max(Field_Info_space(:,1)); % Take the location that has the highest rate
    end

    if ~isempty(Field_Info_space) && (spaceCorr > 0.1) && Field_Info_space(idxMax,4)>5
        placeCellLog(ii) = 1;                
        fieldLoc(ii) = Field_Info_space(idxMax,4)*125/50;
        fieldEdges(ii,:) = Field_Info_space(idxMax,2:3)*125/50;
    end

    Field_Info_tone = detectFields(toneMap,'maxFieldSize',40);
    if ~isempty(Field_Info_tone) 
        [~,idxMax] = max(Field_Info_tone(:,1)); % Take the location that has the highest rate
    end

    if ~isempty(Field_Info_tone) && (toneCorr > 0.1) && Field_Info_tone(idxMax,4)>20
        toneCellLog(ii) = 1;                
        tonefieldLoc(ii) = Field_Info_tone(idxMax,4)*125/50;
        tonefieldEdges(ii,:) = Field_Info_tone(idxMax,2:3)*125/50;
    end   
end

if ~isfield(behavTrials,'probe')
    behavTrials.probe(1:length(behavTrials.lickLoc),1) = 0;
end

if ~isfield(behavTrials,'stim')
    behavTrials.stim(1:length(behavTrials.lickLoc),1) = 0;
end


%% Calculate peak-peak estimation between all place cells 
tSp = [];
toneSp = [];
for a = 1:length(spikes.times)  
    [idx_ts] = InIntervals(spikes.times{a},[tracking.timestamps(1) tracking.timestamps(end)]); 
    tsBehav = spikes.times{a}(idx_ts);
    idx = behavTrials.linTrial == 0 & behavTrials.probe == 0 & behavTrials.stim == 0;
    intervals = behavTrials.timestamps(idx,:);

    %% Convert spike position to tone field position
    toneSpos = spikeData.pos{a};
    gain = [122/9, 122/32 122/55.53 122/79.62 122/102.79 122/122];
    for tt = 1:6
        idx = behavTrials.linTrial==0 & behavTrials.lickLoc == (tt-1);
        ts_spike = tracking.timestamps(spikeData.posIdx{a});
        [idxPos] = InIntervals(ts_spike,behavTrials.timestamps(idx,:));   
        toneSpos(idxPos) = toneSpos(idxPos)*gain(tt);
    end
    toneSpos(toneSpos>125) = nan;
    toneSpos = reshape(toneSpos,length(toneSpos),1);

    % Only take spikes that occur in the forward run
    bools = InIntervals(tsBehav,intervals);
    posData = spikeData.pos{a};
    posData = reshape(posData,length(posData),1);
    if placeCellLog(a)        
        tSp{a} = tsBehav(bools & posData>fieldEdges(a,1) & posData<fieldEdges(a,2));
    else
        tSp{a} = [];
    end
    if toneCellLog(a)
        toneSp{a} = tsBehav(bools & toneSpos>tonefieldEdges(a,1) & toneSpos<tonefieldEdges(a,2));
    else
        toneSp{a} = [];
    end

end

binSize = 0.005;
if sum(placeCellLog)>1
    [cor, lag] = CCG(tSp,[],'binSize',binSize,'duration',0.5);   
else
    cor = [];
end

if sum(toneCellLog)>1
    [corTone, lag] = CCG(toneSp,[],'binSize',binSize,'duration',0.5);   
else
    cor = [];
end

%% Get pairwise ccg information
pairIDs = [];
place_offset = [];
ccgs_time_offset = [];

pairIDs_tone = [];
tone_offset = [];
tone_ccgs_time_offset = [];


%% Place cell CCG pairs
for i = 1:(size(cor,3)-1)
    for j = (i+1):size(cor,3)  
        %% If not enough spikes for the ccg
        if sum(cor(:,i,j)) < 50
            continue
        end
        %% If both cells are place cells - but only take cells that are overlapping
        if placeCellLog(i)==1 && placeCellLog(j)==1
    
            if fieldEdges(i,1) <= fieldEdges(j,1)
                if fieldEdges(i,2)< (fieldEdges(j,1)-2)
                    continue
                end
            else 
                if fieldEdges(j,2)< (fieldEdges(i,1)-2)
                    continue
                end
            end

            [cor2, lag2, smooth] = SmoothCor(cor(:,i,j), lag, binSize);
            filtCor = bandpass(smooth,[6 12],1/binSize);
            
            pairIDs = [pairIDs;i j];
    
            % Get theta-cycle time offset
            %[~,locs] = findpeaks(filtCor,'MinPeakHeight',mean(abs(filtCor)));
            idxTheta = lag2>-0.1 & lag2<0.1;
            lags = lag2(idxTheta);
            [peakCor,maxTheta] = max(smooth(idxTheta));

            if peakCor>(mean(smooth(idxTheta))+2*std(smooth(idxTheta)))
                ccgs_time_offset = [ccgs_time_offset lags(maxTheta)];
            else
                ccgs_time_offset = [ccgs_time_offset nan];
            end
            % if isempty(locs)~=1
            %     [~,J]=sort(abs(lag2(locs)),'ascend');
            %     ccgs_time_offset = [ccgs_time_offset lag2(locs(J(1)))];
            % else
            %     ccgs_time_offset = [ccgs_time_offset nan];
            % end
            
            place_offset = [place_offset fieldLoc(j)-fieldLoc(i)];

        end
    end
end

%% Tone cell CCG pairs
for i = 1:(size(cor,3)-1)
    for j = (i+1):size(cor,3)  
        %% If not enough spikes for the ccg
        if sum(corTone(:,i,j)) <50
            continue
        end
        %% If both cells are place cells - but only take cells that are overlapping
        if toneCellLog(i)==1 && toneCellLog(j)==1

            if tonefieldEdges(i,1) <= tonefieldEdges(j,1)
                if tonefieldEdges(i,2)< (tonefieldEdges(j,1)-2)
                    continue
                end
            else 
                if tonefieldEdges(j,2)< (tonefieldEdges(i,1)-2)
                    continue
                end
            end

            [cor2, lag2, smooth] = SmoothCor(corTone(:,i,j), lag, binSize);
            filtCor = bandpass(smooth,[6 12],1/binSize);
            
            pairIDs_tone = [pairIDs_tone;i j];
    
            idxTheta = lag2>-0.1 & lag2<0.1;
            lags = lag2(idxTheta);
            [peakCor,maxTheta] = max(smooth(idxTheta));

            if peakCor>(mean(smooth(idxTheta))+2*std(smooth(idxTheta)))
                tone_ccgs_time_offset = [tone_ccgs_time_offset lags(maxTheta)];
            else
                tone_ccgs_time_offset = [tone_ccgs_time_offset nan];
            end

            % % Get theta-cycle time offset
            % [~,locs] = findpeaks(filtCor,'MinPeakHeight',mean(abs(filtCor)));
            % if isempty(locs)~=1
            %     [~,J]=sort(abs(lag2(locs)),'ascend');
            %     tone_ccgs_time_offset = [tone_ccgs_time_offset lag2(locs(J(1)))];
            % else
            %     tone_ccgs_time_offset = [tone_ccgs_time_offset nan];
            % end
            
            tone_offset = [tone_offset tonefieldLoc(j)-tonefieldLoc(i)];

        end
    end
end

%% Store important variables
thetaComp.place_offset = place_offset;
thetaComp.ccgs_time_offset = ccgs_time_offset;
thetaComp.pairIDs = pairIDs;
thetaComp.tone_offset = tone_offset;
thetaComp.tone_ccgs_time_offset = tone_ccgs_time_offset;
thetaComp.pairIDs_tone = pairIDs_tone;

if plotfig
   figure
   subplot(1,2,1)
   scatter(place_offset,ccgs_time_offset,'.')

   subplot(1,2,2)
   scatter(tone_offset, tone_ccgs_time_offset,'.' )
end

end