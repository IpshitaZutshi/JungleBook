function precessData = sessionPhasePrecession(varargin)

% This function extracts tone cells and place cells during forward tone
% trials. It calculates the phase precession slop for each trial type

warning off
%% defaults
p = inputParser;
addParameter(p,'expPath',pwd,@isfolder)
addParameter(p,'plotfig',true,@islogical)
addParameter(p,'tonecell',false,@islogical)
addParameter(p,'lfpChan',63)
addParameter(p,'speedThresh',5)

parse(p,varargin{:})
expPath = p.Results.expPath;
plotfig = p.Results.plotfig;
tonecell = p.Results.tonecell;
lfpChan = p.Results.lfpChan;
speedThresh = p.Results.speedThresh;

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
ss = 0;

for ii = 1:spikes.numcells

    if ~strcmp(cell_metrics.putativeCellType{ii},'Pyramidal Cell') == 1         
        continue
    end

    tSp = spikeData.pos{ii};
    tSp = reshape(tSp,[length(tSp),1]);
    vSp = tracking.position.vy(spikeData.posIdx{ii});
    pSp = rad2deg(spikeData.phase{ii}); 

    %% Get the average map across all tone trials
    idxTrials = find(behavTrials.linTrial==0);
    
    dataMatTone = [];
    dataMat = [];
    for kk = 1:(length(idxTrials)-1)
        dataMat(kk,:) = trialMap.firingMaps.forward.rateMaps{ii}{idxTrials(kk)};
        a = fillmissing(trialMap.firingMaps.tone.rateMaps{ii}{idxTrials(kk)},'linear');
        dataMatTone(kk,:) = a;
    end       
    
    spaceMap = nanmean(dataMat,1);
    toneMap = nanmean(dataMatTone,1);       
    
    C = corrcoef(dataMat','Rows','complete');
    B = tril(C,-1);
    D = B(abs(B)>0);
    spaceCorr = nanmean(D);

    C = corrcoef(dataMatTone','Rows','complete');
    B = tril(C,-1);
    D = B(abs(B)>0);
    toneCorr = nanmean(D);

    %% Detect fields
    Field_Info_space = detectFields(spaceMap,'minFieldSize',4,'maxFieldSize',35);
    if ~isempty(Field_Info_space) 
        [~,idxMaxSpace] = max(Field_Info_space(:,1)); % Take the location that has the biggest field
    end

    Field_Info_tone = detectFields(toneMap,'minFieldSize',4,'maxFieldSize',35);
    if ~isempty(Field_Info_tone) 
        [~,idxMaxTone] = max(Field_Info_tone(:,1)); % Take the location that has the biggest field
    end

    if ~tonecell %% Only look at place cells
        logicalIdx  = ~isempty(Field_Info_space) && spaceCorr > 0.1 && Field_Info_space(idxMaxSpace,4)>3;
    else %% Only look at tone cells
        logicalIdx  = ~isempty(Field_Info_tone) && toneCorr > 0.1 && Field_Info_tone(idxMaxTone,4)>25;
    end

    if logicalIdx % If the cell fulfils the criteria then loop through to get trial type data
        ss = ss+1;
        cellnum(ss) = ii; 

        if ~tonecell
            fieldLoc = Field_Info_space(idxMaxSpace,4); 
        else
            fieldLoc = Field_Info_tone(idxMaxTone,4); 
        end
        [idx] = InIntervals(spikes.times{ii},[tracking.timestamps(1) tracking.timestamps(end)]); 
        tsBehav1 = spikes.times{ii}(idx);        

        for trialType = 2:6 %% Loop through trial types
            if ~isfield(behavTrials,'probe')
                behavTrials.probe(1:length(behavTrials.linTrial),1) = 0;
            end
            if ~isfield(behavTrials,'stim')
                behavTrials.stim(1:length(behavTrials.linTrial),1) = 0;
            end

            idxTrial = behavTrials.lickLoc == (trialType - 1) & behavTrials.linTrial == 0 & behavTrials.correct == 1 & behavTrials.probe == 0 & behavTrials.stim == 0;
            intervals = behavTrials.timestamps(idxTrial,:) - 0.033;


            field_Info = detectFields(firingMaps.forward.rateMaps{ii}{1 + trialType}, 'minFieldSize', 2);
            if ~isempty(field_Info)

                % Only select relevant field
                peak1 = field_Info(:,4);
                if ~tonecell
                    idx1 = find(peak1>(fieldLoc-15) & peak1<(fieldLoc+15));
                else % just pick the field with the higher rate
                    [~,idx1] = max(field_Info(:,1));
                end

                if length(idx1)>1 %% case where it detects multiple overlapping fields
                    [~,pickidx] = max(field_Info(idx1,1));
                    idx1 = idx1(pickidx);
                end
                % Select spikes within trial intervals and field boundaries
                bools1 = InIntervals(tsBehav1, intervals);

                if ~isempty(idx1)   
                    fieldStart = field_Info(idx1, 2) * 125 / 50;
                    fieldEnd = field_Info(idx1, 3) * 125 / 50;
                    fieldSize = fieldEnd - fieldStart;

                    % Select spikes within place fields
                    bools_1 = bools1 & tSp > fieldStart & tSp < fieldEnd & vSp> speedThresh; % Skip early spikes that are not accurate 

                    % Position of those spikes
                    lindat = (tSp(bools_1) - fieldStart) ./ fieldSize;

                    % Also estimate phase of those spikes
                    phidat = pSp(bools_1);
                    phasedeg = phidat - 180;
                    phidat = deg2rad(phasedeg);

                    % Calculate circular-linear regression
                    [slope(ss,trialType-1), offset(ss,trialType-1), ~, corrR(ss,trialType-1), sigP(ss,trialType-1)] = rccc(lindat, phidat, [-2 2], 0.05); 
                else
                    slope(ss,trialType-1) = nan;
                    offset(ss,trialType-1) = nan;
                    corrR(ss,trialType-1) = nan;
                    sigP(ss,trialType-1) = nan;
                end
            else
                slope(ss,trialType-1) = nan;
                offset(ss,trialType-1) = nan;
                corrR(ss,trialType-1) = nan;
                sigP(ss,trialType-1) = nan;
            end
        end
    end
end

%% Store important variables
if ss>0
    precessData.slope = slope;
    precessData.offset = offset;
    precessData.corrT = corrR;
    precessData.sigP = sigP;
    precessData.cellnum = cellnum;
else
    precessData.slope = [];
    precessData.offset = [];
    precessData.corrT = [];
    precessData.sigP = [];
    precessData.cellnum = [];
end

if plotfig
   figure
   set(gcf,'Renderer','painters')
   set(gcf,'Color','w')
   col = [241/255 114/255 42/255;...
        247/255 149/255 33/255;...
        249/255 197/255 81/255;...
        143/255 189/255 107/255;...
        87/255 116/255 144/255];

   for ii = 1:5
        subplot(2,5,ii)
        histogram(precessData.slope(:,ii),-2:0.2:2,'FaceColor',[0.7 0.7 0.7])
        hold on
        line([nanmedian(precessData.slope(:,ii)) nanmedian(precessData.slope(:,ii))],[0 5],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
        histogram(precessData.slope(precessData.sigP(:,ii)<0.05,ii),-2:0.2:2,'FaceColor',col(ii,:))
        line([nanmedian(precessData.slope(precessData.sigP(:,ii)<0.05,ii)) nanmedian(precessData.slope(precessData.sigP(:,ii)<0.05,ii))],[0 5],'Color',col(ii,:),'LineWidth',1.5)
        idxSig{ii} = precessData.sigP(:,ii)<0.05;
   end

   subplot(2,5,6)
   groupStats({precessData.slope(:,1),precessData.slope(:,2),precessData.slope(:,3),...
       precessData.slope(:,4),precessData.slope(:,5)},1:5,'inAxis',true,'color',col)
   title('Slope all') 

   subplot(2,5,7)
   groupStats({precessData.offset(:,1),precessData.offset(:,2),precessData.offset(:,3),...
       precessData.offset(:,4),precessData.offset(:,5)},1:5,'inAxis',true,'color',col)
   title('Intercept all') 
   
   subplot(2,5,8)
   groupStats({precessData.slope(idxSig{1},1),precessData.slope(idxSig{2},2),precessData.slope(idxSig{3},3),...
       precessData.slope(idxSig{4},4),precessData.slope(idxSig{5},5)},1:5,'inAxis',true,'color',col)
   title('Slope sig') 
   
   subplot(2,5,9)
   groupStats({precessData.offset(idxSig{1},1),precessData.offset(idxSig{2},2),precessData.offset(idxSig{3},3),...
       precessData.offset(idxSig{4},4),precessData.offset(idxSig{5},5)},1:5,'inAxis',true,'color',col)
   title('Intercept sig') 
end

end