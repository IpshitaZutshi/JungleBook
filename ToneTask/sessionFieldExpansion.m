function fieldData = sessionFieldExpansion(varargin)

% This function extracts tone cells and place cells during forward tone
% trials. It calculates the CCG in "place space", and in "theta timescale"
% and compares across each port.

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
file = dir(['*.rateMapsAvg.cellinfo.mat']);
load(file(1).name);
file = dir(['*.rateMapsTrial.cellinfo.mat']);
trialMap = load(file(1).name);

%% Label cells that are either place or tone
placeCellLog = zeros(spikes.numcells,1);

for ii = 1:spikes.numcells

    if strcmp(cell_metrics.putativeCellType{ii},'Pyramidal Cell') == 1 
        cellType(ii) = 1;
    else
        cellType(ii) = 0;
    end

    %% Get the average place map across all tone trials
    idxTrials = find(behavTrials.linTrial==0);
    
    dataMat = [];
    for kk = 1:(length(idxTrials)-1)
        dataMat(kk,:) = trialMap.firingMaps.forward.rateMaps{ii}{idxTrials(kk)};
    end       
    C = corrcoef(dataMat','Rows','complete');
    B = tril(C,-1);
    D = B(abs(B)>0);
    spaceMap = nanmean(dataMat,1);
    spaceCorr(ii) = nanmean(D);

    %% Detect fields
    Field_Info_space = detectFields(spaceMap,'minFieldSize',4,'maxFieldSize',35);
    if ~isempty(Field_Info_space) 
        [~,idxMax] = max(Field_Info_space(:,1)); % Take the location that has the biggest field
    end

    if ~isempty(Field_Info_space) && (spaceCorr(ii) > 0) && cellType(ii) == 1 && Field_Info_space(idxMax,4)>3
        placeCellLog(ii) = 1;                
        fieldLoc(ii) = Field_Info_space(idxMax,4);      
        % find the center of mass
        extractedField = spaceMap(Field_Info_space(idxMax,2):Field_Info_space(idxMax,3));
        last_nonNaN_index = find(~isnan(extractedField), 1, 'last');
        extractedField = extractedField(1:last_nonNaN_index);
        avgCOM(ii) = (((1:length(extractedField))*extractedField')./sum(extractedField))+(Field_Info_space(idxMax,2)-1);
    end    

    % Also detect peak firing locations for each trial type
    if placeCellLog(ii) == 1 
        fieldEdges(ii,:) = Field_Info_space(idxMax,2:3);
        for jj = 1:6 
            field_Info = detectFields(firingMaps.forward.rateMaps{ii}{jj+1},'minFieldSize',2,'maxFieldSize',35,'maxRate',5);
            % Find the field location closest to the average location,
            % which means, the peak should be within XX bins on each side -
            idx1 = [];
            if ~isempty(field_Info)
                peak1 = field_Info(:,4);
                idx1 = find(peak1>(fieldLoc(ii)-15) & peak1<(fieldLoc(ii)+15));
            end

            % % As long as the peak firing rate > 5, find the max
            % [maxFR,idx] = max(firingMaps.forward.rateMaps{ii}{jj+1});
 
            if ~isempty(idx1)     
                
                fieldStart(ii,jj) = field_Info(idx1(1),2);
                fieldEnd(ii,jj) = field_Info(idx1(1),3);    
                % find the center of mass
                extractedField = firingMaps.forward.rateMaps{ii}{jj+1}(fieldStart(ii,jj):fieldEnd(ii,jj));
                last_nonNaN_index = find(~isnan(extractedField), 1, 'last');
                extractedField = extractedField(1:last_nonNaN_index);

                peakSpace(ii,jj) = peak1(idx1(1));
                
                fieldStartRate(ii,jj) = firingMaps.forward.rateMaps{ii}{jj+1}(fieldStart(ii,jj));

                
                fieldEndRate(ii,jj) = extractedField(end);

                
                trialCOM(ii,jj) = (((1:length(extractedField))*extractedField')./sum(extractedField))+(fieldStart(ii,jj)-1);                
            else
                peakSpace(ii,jj) = nan;
                fieldStart(ii,jj) = nan;
                fieldEnd(ii,jj) = nan;
                trialCOM(ii,jj) = nan;
                fieldEndRate(ii,jj) = nan;
                fieldStartRate(ii,jj) = nan;
            end
        end
    else
        peakSpace(ii,1:6) = nan;
        fieldStart(ii,1:6) = nan;
        fieldEnd(ii,1:6) = nan;
        fieldLoc(ii) = nan;
        fieldEdges(ii,1:2) = nan;
        trialCOM(ii,1:6) = nan;
        fieldEndRate(ii,1:6) = nan;
        fieldStartRate(ii,1:6) = nan;   
        avgCOM(ii) = nan;
    end

end
% 
%% Store important variables
%idx = find(~isnan(fieldLoc));
fieldData.placefield = peakSpace*(125/50);
fieldData.avgCOM = avgCOM*(125/50);
fieldData.avgfield = fieldLoc*(125/50);
fieldData.avgfieldEdges = fieldEdges*(125/50);
fieldData.fieldStart = fieldStart*(125/50);
fieldData.fieldEnd = fieldEnd*(125/50);
fieldData.trialCOM = trialCOM*(125/50);
fieldData.fieldEndRate = fieldEndRate;
fieldData.fieldStartRate = fieldStartRate;   


if plotfig
   figure
   set(gcf,'Renderer','painters')
   set(gcf,'Color','w')
   subplot(2,4,1)
   groupStats({fieldData.placefield(:,1),fieldData.placefield(:,2),fieldData.placefield(:,3),...
       fieldData.placefield(:,4),fieldData.placefield(:,5),fieldData.placefield(:,6)},1:6,'inAxis',true)
   title('Place field peak') 

   subplot(2,4,2)
   subMat  = (fieldData.placefield-fieldData.avgfield');
   groupStats({subMat(:,1),subMat(:,2),subMat(:,3),subMat(:,4),subMat(:,5),subMat(:,6)},1:6,'inAxis',true)
   title('Peak loc - avg loc') 

   subplot(2,4,3)
   groupStats({fieldData.fieldStart(:,1),fieldData.fieldStart(:,2),fieldData.fieldStart(:,3),...
       fieldData.fieldStart(:,4),fieldData.fieldStart(:,5),fieldData.fieldStart(:,6)},1:6,'inAxis',true)
   title('Place field start') 

   subplot(2,4,4)
   groupStats({fieldData.fieldEnd(:,1),fieldData.fieldEnd(:,2),fieldData.fieldEnd(:,3),...
       fieldData.fieldEnd(:,4),fieldData.fieldEnd(:,5),fieldData.fieldEnd(:,6)},1:6,'inAxis',true)
   title('Place field end') 

   subplot(2,4,5)
   groupStats({fieldData.trialCOM(:,1),fieldData.trialCOM(:,2),fieldData.trialCOM(:,3),...
       fieldData.trialCOM(:,4),fieldData.trialCOM(:,5),fieldData.trialCOM(:,6)},1:6,'inAxis',true)
   title('Trial COM') 

   subplot(2,4,6)
   subMat  = (fieldData.trialCOM-fieldData.avgCOM');
   groupStats({subMat(:,1),subMat(:,2),subMat(:,3),subMat(:,4),subMat(:,5),subMat(:,6)},1:6,'inAxis',true)
   title('Trial COM - avg COM') 

   subplot(2,4,7)
   groupStats({fieldData.fieldStartRate(:,1),fieldData.fieldStartRate(:,2),fieldData.fieldStartRate(:,3),...
       fieldData.fieldStartRate(:,4),fieldData.fieldStartRate(:,5),fieldData.fieldStartRate(:,6)},1:6,'inAxis',true)
   title('Place field start rate') 

   subplot(2,4,8)
   groupStats({fieldData.fieldEndRate(:,1),fieldData.fieldEndRate(:,2),fieldData.fieldEndRate(:,3),...
       fieldData.fieldEndRate(:,4),fieldData.fieldEndRate(:,5),fieldData.fieldEndRate(:,6)},1:6,'inAxis',true)
   title('Place field end rate')    
end

end