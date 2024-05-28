function thetaComp = sessionThetaCompression(varargin)

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

%% Label cells that are either place or tone
placeCellLog = zeros(spikes.numcells,1);
toneCellLog = zeros(spikes.numcells,1);

for ii = 1:spikes.numcells

    if strcmp(cell_metrics.putativeCellType{ii},'Pyramidal Cell') == 1 
        cellType(ii) = 1;
    else
        cellType(ii) = 0;
    end

    %% Check if its a stable place cell or tone cell
    dataMat = [];
    dataMatTone = [];
    info = [];
    for jj = 2:7    
        dataMat = [dataMat;firingMaps.forward.rateMaps{ii}{jj}]; % Place
        info = [info;firingMaps.forward.information{ii}{jj}];
        a = fillmissing(firingMaps.tone.rateMaps{ii}{jj},'linear'); % tone
        dataMatTone = [dataMatTone;a];
    end           
    spaceMap = nanmean(dataMat,1);
    toneMap = nanmean(dataMatTone,1);        
    spatialInfo(ii) = nanmean(info);

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
    toneCorr(ii) = nanmean(corrTone);  
    spaceCorr(ii) = nanmean(corrSpace);

    %% Detect fields
    Field_Info_tone = detectFields(toneMap);
    Field_Info_space = detectFields(spaceMap);

    if ~isempty(Field_Info_tone) && (toneCorr(ii) > 0.1)
        toneCellLog(ii) = 1;       
    end

    if ~isempty(Field_Info_space) && (spaceCorr(ii) > 0.1) && (spatialInfo(ii) > 0.3)
        placeCellLog(ii) = 1;
    end    

    % Also detect peak firing locations for each trial type
    if placeCellLog(ii) == 1 % If it is a place cell, find the location
        % Take the location that has the biggest field
        [~,idxMax] = max(Field_Info_space(:,1));
        fieldLoc = Field_Info_space(idxMax,4);
        fieldEdges(ii,:) = Field_Info_space(idxMax,2:3);
        for jj = 1:6 
            field_Info = detectFields(firingMaps.forward.rateMaps{ii}{jj+1});
            % Find the field location closest to the average location,
            % which means, the peak should be within XX bins on each side -
            % 10?
            idx1 = [];
            if ~isempty(field_Info)
                peak1 = field_Info(:,4);
                idx1 = find(peak1>(fieldLoc-10) & peak1<(fieldLoc+10));
            end
         
            if ~isempty(idx1)
                peakSpace(ii,jj) = peak1(idx1(1));
                fieldStart(ii,jj) = field_Info(idx1(1),2);
                fieldEnd(ii,jj) = field_Info(idx1(1),3);
            else
                peakSpace(ii,jj) = nan;
                fieldStart(ii,jj) = nan;
                fieldEnd(ii,jj) = nan;
            end
        end
    else
        peakSpace(ii,1:6) = nan;
        fieldStart(ii,1:6) = nan;
        fieldEnd(ii,1:6) = nan;
    end

end
% 
% %% First, CCGs between place cell pairs
if ~isfield(behavTrials,'probe')
    behavTrials.probe(1:length(behavTrials.lickLoc)) = 0;
end

% Calculate peak-peak estimation between all cells within a trial Type
for ii = 1:2
    tSp = [];
    for a = 1:length(spikes.times)  
        [idx_ts] = InIntervals(spikes.times{a},[tracking.timestamps(1) tracking.timestamps(end)]); 
        tsBehav = spikes.times{a}(idx);

        if ii == 1
            idx = (behavTrials.lickLoc ==3 | behavTrials.lickLoc ==2) & behavTrials.linTrial == 0 & behavTrials.probe == 0;
        elseif ii == 2
            idx = behavTrials.lickLoc >=5 & behavTrials.linTrial == 0 & behavTrials.probe == 0;
        end

        intervals = behavTrials.timestamps(idx,:);
        bools = InIntervals(tsBehav,intervals);
        if placeFieldLog(a)
            % Only take spikes that occur within the place field
            tSp{a} = tsBehav(bools & spikeData.pos{a}>fieldEdges(a,1) & spikeData.pos{a}>fieldEdges(a,2));
        else
            tSp{a} = [];
        end
    end

    % Longer timescale difference between cells
    [ccg2,t2] = CCG(tSp,[],'binSize',0.1,'duration',2);
    t2 = t2*1000;

    % Theta-cycle time difference between place fields
    [ccg,t] = CCG(tSp,[],'binSize',0.002,'duration',0.4);    
    t = t * 1000; % from sec to msec

    %% Initialize    
    pairIDs{ii} = NaN(nchoosek(length(spikes.times),2),2);

    ccgs_place_offset{ii} = NaN(nchoosek(length(spikes.times),2),1);
    place_offset{ii} = NaN(nchoosek(length(spikes.times),2),1);
    ccgs_time_offset{ii} = NaN(nchoosek(length(spikes.times),2),1);
    ccgs_placefield{ii} = NaN(nchoosek(length(spikes.times),2),1);
    ccgs_tonefield{ii} = NaN(nchoosek(length(spikes.times),2),1);

    %% Get pairwise ccg information
    kk = 1;
    for i = 1:(size(ccg,3)-1)
        for j = (i+1):size(ccg,3)   
            [ccg_trace2,~,smooth2] = SmoothCor(ccg2(:,i,j),t2,0.1); % slow timescale 
            [ccg_trace,~,smooth] = SmoothCor(ccg(:,i,j),t,0.002); % fast timescale 

            pairIDs{ii}(kk,:) = [i j];

            % Get place field time offset
            [~,indx] = max(ccg_trace2);
            if max(ccg_trace2)> (mean(smooth2)+(0.5*std(smooth2)))
                ccgs_place_offset{ii}(kk) = t2(indx);
            end

            % Get place field offset
            place_offset{ii}(kk) = peakSpace(j,ii)-peakSpace(i,ii);

            % Get theta-cycle time offset
            [~,locs] = findpeaks(smooth,'MinPeakHeight',mean(smooth),'MinPeakDistance',50);
            if isempty(locs)~=1
                [~,J]=sort(abs(t(locs)),'ascend');
                ccgs_time_offset{ii}(kk)=t(locs(J(1)));
            end   

            if (cellType(i)==1 && placeCellLog(i)==1 && cellType(j)==1 && placeCellLog(j)==1)
                ccgs_placefield{ii}(kk) = 1;
            end

            if (cellType(i)==1 && toneCellLog(i)==1 && cellType(j)==1 && toneCellLog(j)==1)
                ccgs_tonefield{ii}(kk) = 1;
            end

            kk = kk+1;
        end
    end
end

%% Store important variables

thetaComp.placefield = ccgs_placefield;
thetaComp.tonefield = ccgs_tonefield;
thetaComp.thetaDT = ccgs_time_offset;
thetaComp.placeDT = ccgs_place_offset;
thetaComp.placeDist = place_offset;

if plotfig
    % selectidx = ccgs_placefield{1}==1 & ccgs_placefield{2}==1 & ccgs_time_offset{1} <250 & ccgs_time_offset{1} >-250 ...
    %     & ccgs_time_offset{2} <250 & ccgs_time_offset{2} >-250 & ...
    %     ccgs_place_offset{1} <800 & ccgs_place_offset{1} >-800 ...
    %     & ccgs_place_offset{2} <800 & ccgs_place_offset{2} >-800;

    selectidx1 = (abs(place_offset{1})<10) & abs(ccgs_time_offset{1}) <250 & abs(ccgs_time_offset{2}) < 250;
    selectidx2 = (abs(place_offset{2})<10) & abs(ccgs_time_offset{1}) <250 & abs(ccgs_time_offset{2}) < 250;

    figure
    subplot(2,2,1)
    %scatter(ccgs_time_offset{1}(selectidx),ccgs_place_offset{1}(selectidx),'.')
    scatter(place_offset{1}(selectidx1),ccgs_time_offset{1}(selectidx1),'.')
    hold on
    lsline
    %[r,p]= corrcoef(ccgs_time_offset{1}(selectidx),ccgs_place_offset{1}(selectidx));
    [r,p]= corrcoef(place_offset{1}(selectidx1),ccgs_time_offset{1}(selectidx1),'Rows','complete');
    xlabel('delta Position bins')
    ylabel('Theta CCG(ms)')
    title(strcat('Short (1-3), r=', num2str(r(1,2)),', p= ',num2str(p(1,2))));
    
    subplot(2,2,2)
    scatter(place_offset{2}(selectidx2),ccgs_time_offset{2}(selectidx2),'.')
    %scatter(ccgs_time_offset{2}(selectidx),ccgs_place_offset{2}(selectidx),'.')
    hold on
    lsline
    [r,p]= corrcoef(place_offset{2}(selectidx2),ccgs_time_offset{2}(selectidx2),'Rows','complete');
    %[r,p]= corrcoef(ccgs_time_offset{2}(selectidx),ccgs_place_offset{2}(selectidx));
    % xlim([-250 250])
    % ylim([-800 800])
    xlabel('delta Position bins')
    ylabel('Theta CCG(ms)')
    title(strcat('Long (5-6), r=', num2str(r(1,2)),', p= ',num2str(p(1,2))));
    
    subplot(2,2,3)
    %groupStats([{ccgs_time_offset{1}} {ccgs_time_offset{2}}],[],'inAxis',true,'repeatedMeasures',true,'plotType','BoxLinesSEM');
    groupStats([{ccgs_time_offset{1}(selectidx1)} {ccgs_time_offset{2}(selectidx2)}],[],'inAxis',true);
    title('Theta timescale')
    
    subplot(2,2,4)
    %groupStats([{ccgs_place_offset{1}(selectidx)} {ccgs_place_offset{2}(selectidx)}],[],'inAxis',true,'repeatedMeasures',true,'plotType','BoxLinesSEM');
    groupStats([{place_offset{1}} {place_offset{2}}],[],'inAxis',true);
    title('Place field distance')
end

end