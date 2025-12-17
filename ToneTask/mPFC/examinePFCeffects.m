function examinePFCeffects

sess = {'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   
    'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...
    'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',...
    'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16'}; 
    
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';
var = 'stim'; % <<<< IMPORTANT: must match behavTrials field used for silencing/stim

AllcellType = [];
AllspatialCorr = [];
AlltoneCorr = [];
AllspatialCorrStim = [];
AlltoneCorrStim = [];
AllspaceMap = [];
AlltoneMap = [];
AllspaceMapStim = [];
AlltoneMapStim = [];
AllspaceField = [];
AlltoneField = [];
AllspaceFieldStim = [];
AlltoneFieldStim = [];
AllspatialInfo = [];
AlltoneInfo = [];
AllspatialInfoStim = [];
AlltoneInfoStim = [];

% ===== NEW: fair ceiling + observed + shuffle (baseline vs stim) =====
AllspaceCeilBase = [];
AllspaceCeilStim = [];
AlltoneCeilBase  = [];
AlltoneCeilStim  = [];

AllspaceBaseStimMatched = [];
AlltoneBaseStimMatched  = [];

AllspaceBaseStimShuf = [];
AlltoneBaseStimShuf  = [];

nSplit   = 20;   % split-half repeats for ceiling
nShuffle = 200;  % shifts for shuffle null
minTrialsPerCond = 4; % require at least this many trials within a toneGain to compute split-half reliably

for ii = 1:length(sess)
    %% Load files
    cd(strcat(expPath,sess{ii}))    
    file = dir(['*.rateMapsAvgnotLog.cellinfo.mat']);
    load(file(1).name); %#ok<LOAD>
    file = dir(['*.rateMapsTrial.cellinfo.mat']);
    trialMaps = load(file(1).name);    
    file = dir(['*.Tracking.Behavior.mat']);
    load(file(1).name); %#ok<LOAD>
    file = dir(['*TrialBehavior.Behavior.mat']);
    load(file(1).name); %#ok<LOAD>
    [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts', true); %#ok<NASGU>
    file = dir(['*cell_metrics.cellinfo.mat']);
    load(file(1).name); %#ok<LOAD>
     
    cellType(1:length(cell_metrics.UID),1) = nan;
    spaceMap(1:length(cell_metrics.UID),1:50) = nan;
    toneMap(1:length(cell_metrics.UID),1:50) = nan;
    spaceMapStim(1:length(cell_metrics.UID),1:50) = nan;
    toneMapStim(1:length(cell_metrics.UID),1:50) = nan;
    spatialCorr(1:length(cell_metrics.UID),1) = nan;
    toneCorr(1:length(cell_metrics.UID),1) = nan;    
    spatialCorrStim(1:length(cell_metrics.UID),1) = nan;
    toneCorrStim(1:length(cell_metrics.UID),1) = nan;    
    spaceField(1:length(cell_metrics.UID),1) = nan;
    toneField(1:length(cell_metrics.UID),1) = nan;   
    spaceFieldStim(1:length(cell_metrics.UID),1) = nan;
    toneFieldStim(1:length(cell_metrics.UID),1) = nan;       
    spatialInfo(1:length(cell_metrics.UID),1) = nan;
    toneInfo(1:length(cell_metrics.UID),1) = nan; 
    spatialInfoStim(1:length(cell_metrics.UID),1) = nan;
    toneInfoStim(1:length(cell_metrics.UID),1) = nan; 

    % ===== NEW per-session arrays =====
    spaceCeilBase(1:length(cell_metrics.UID),1) = nan;
    spaceCeilStim(1:length(cell_metrics.UID),1) = nan;
    toneCeilBase(1:length(cell_metrics.UID),1)  = nan;
    toneCeilStim(1:length(cell_metrics.UID),1)  = nan;

    spaceBaseStimMatched(1:length(cell_metrics.UID),1) = nan;
    toneBaseStimMatched(1:length(cell_metrics.UID),1)  = nan;

    spaceBaseStimShuf(1:length(cell_metrics.UID),1) = nan;
    toneBaseStimShuf(1:length(cell_metrics.UID),1)  = nan;

    % Trial count alignment: your idx code elsewhere uses (1:end-1)
    nT = length(behavTrials.linTrial) - 1;

    % Precompute baseline/stim trial indices by toneGain (0..5)
    baseIdx = cell(6,1);
    stimIdx = cell(6,1);
    for g = 0:5
        baseIdx{g+1} = find(behavTrials.toneGain(1:nT) == g & behavTrials.correct(1:nT) == 1 & ...
                            behavTrials.linTrial(1:nT) == 0 & behavTrials.(var)(1:nT) == 0);
        stimIdx{g+1} = find(behavTrials.toneGain(1:nT) == g & behavTrials.correct(1:nT) == 1 & ...
                            behavTrials.linTrial(1:nT) == 0 & behavTrials.(var)(1:nT) == 1);
    end

    for kk=1:length(cell_metrics.UID)
    %% Extract pyramidal cells and calculate rate map correlations and information        
        if strcmp(cell_metrics.putativeCellType{kk},'Pyramidal Cell') == 1 
            cellType(kk,1) = 1;
        else
            cellType(kk,1) = 0;
        end

        % Original: trial-type averaged maps (across 2:7 baseline, and 14:19 stim)
        dataMat =[];
        dataMatTone = [];
        dataMatStim =[];
        dataMatToneStim = [];   
        
        info = [];
        infoTone = [];
        infoStim = [];
        infoToneStim = [];
        
        for jj = 2:7    
            dataMat = [dataMat;firingMaps.forward.rateMaps{kk}{jj}];
            info = [info;firingMaps.forward.information{kk}{jj}];
            
            a = fillmissing(firingMaps.tone.rateMaps{kk}{jj},'linear');
            dataMatTone = [dataMatTone;a];
            infoTone = [infoTone;firingMaps.tone.information{kk}{jj}];
            
            dataMatStim = [dataMatStim;firingMaps.forward.rateMaps{kk}{jj+12}];
            infoStim = [infoStim;firingMaps.forward.information{kk}{jj+12}];
            
            a = fillmissing(firingMaps.tone.rateMaps{kk}{jj+12},'linear');
            dataMatToneStim = [dataMatToneStim;a];  
            infoToneStim = [infoToneStim;firingMaps.tone.information{kk}{jj+12}];
        end
        
        spaceMap(kk,:) = nanmean(dataMat,1);
        toneMap(kk,:)  = nanmean(dataMatTone,1);        
        
        spaceMapStim(kk,:) = nanmean(dataMatStim,1);
        toneMapStim(kk,:)  = nanmean(dataMatToneStim,1);                
        
        % Original “within condition” across trial types (not a ceiling, but keep it)
        corrTone = []; 
        corrSpace = [];
        corrToneStim = []; 
        corrSpaceStim = [];
        
        for pp = 1:6
           for jj = (pp+1):6
               a = corrcoef(dataMat(pp,:),dataMat(jj,:),'rows','complete');
               corrSpace = [corrSpace a(1,2)];               
               a = corrcoef(dataMatTone(pp,:),dataMatTone(jj,:),'rows','complete');         
               corrTone = [corrTone a(1,2)];

               a = corrcoef(dataMatStim(pp,:),dataMatStim(jj,:),'rows','complete');
               corrSpaceStim = [corrSpaceStim a(1,2)];               
               a = corrcoef(dataMatToneStim(pp,:),dataMatToneStim(jj,:),'rows','complete');         
               corrToneStim = [corrToneStim a(1,2)];               
           end
        end     
        toneCorr(kk,1) = nanmean(corrTone);  
        spatialCorr(kk,1) = nanmean(corrSpace);
        
        toneCorrStim(kk,1) = nanmean(corrToneStim);  
        spatialCorrStim(kk,1) = nanmean(corrSpaceStim);
        
        spatialInfo(kk,1) = nanmean(info);
        toneInfo(kk,1) = nanmean(infoTone);           
        
        spatialInfoStim(kk,1) = nanmean(infoStim);
        toneInfoStim(kk,1) = nanmean(infoToneStim);           

        % =====================================================================
        % NEW: FAIR ceiling + matched baseline↔stim + shuffle (using trialMaps)
        % trialMaps.*.rateMaps{kk}{trial} is a per-trial map (1x50)
        % We use the trial indices from baseIdx/stimIdx to select relevant trials.
        % =====================================================================

        % Collect per-toneGain quantities, then average across gains
        ceilB_sp = nan(6,1); ceilS_sp = nan(6,1);
        ceilB_to = nan(6,1); ceilS_to = nan(6,1);
        obs_sp   = nan(6,1); obs_to   = nan(6,1);
        shuf_sp  = nan(6,1); shuf_to  = nan(6,1);

        for g = 1:6
            trB = baseIdx{g};
            trS = stimIdx{g};

            % Need enough trials to estimate split-half ceilings
            if numel(trB) >= minTrialsPerCond
                XB = stackTrialMaps(trialMaps.firingMaps.forward.rateMaps{kk}, trB);
                if ~isempty(XB)
                    ceilB_sp(g) = splitHalfCeil(XB, nSplit);
                end
                XBt = stackTrialMaps(trialMaps.firingMaps.tone.rateMaps{kk}, trB);
                if ~isempty(XBt)
                    ceilB_to(g) = splitHalfCeil(XBt, nSplit);
                end
            end

            if numel(trS) >= minTrialsPerCond
                XS = stackTrialMaps(trialMaps.firingMaps.forward.rateMaps{kk}, trS);
                if ~isempty(XS)
                    ceilS_sp(g) = splitHalfCeil(XS, nSplit);
                end
                XSt = stackTrialMaps(trialMaps.firingMaps.tone.rateMaps{kk}, trS);
                if ~isempty(XSt)
                    ceilS_to(g) = splitHalfCeil(XSt, nSplit);
                end
            end

            % Matched baseline↔stim similarity (mean map per toneGain)
            if ~isempty(trB) && ~isempty(trS)
                XB = stackTrialMaps(trialMaps.firingMaps.forward.rateMaps{kk}, trB);
                XS = stackTrialMaps(trialMaps.firingMaps.forward.rateMaps{kk}, trS);
                if ~isempty(XB) && ~isempty(XS)
                    mb = nanmean(XB,1);
                    ms = nanmean(XS,1);
                    a = corrcoef(mb, ms, 'rows','complete');
                    if size(a,1)==2, obs_sp(g) = a(1,2); end

                    % Shuffle null: circularly shift stim mean map
                    tmp = nan(nShuffle,1);
                    for ss = 1:nShuffle
                        sh = randi(50);
                        ms_sh = circshift(ms, [0 sh]);
                        aa = corrcoef(mb, ms_sh, 'rows','complete');
                        if size(aa,1)==2, tmp(ss) = aa(1,2); end
                    end
                    shuf_sp(g) = nanmean(tmp);
                end

                XBt = stackTrialMaps(trialMaps.firingMaps.tone.rateMaps{kk}, trB);
                XSt = stackTrialMaps(trialMaps.firingMaps.tone.rateMaps{kk}, trS);
                if ~isempty(XBt) && ~isempty(XSt)
                    mb = nanmean(XBt,1);
                    ms = nanmean(XSt,1);
                    a = corrcoef(mb, ms, 'rows','complete');
                    if size(a,1)==2, obs_to(g) = a(1,2); end

                    tmp = nan(nShuffle,1);
                    for ss = 1:nShuffle
                        sh = randi(50);
                        ms_sh = circshift(ms, [0 sh]);
                        aa = corrcoef(mb, ms_sh, 'rows','complete');
                        if size(aa,1)==2, tmp(ss) = aa(1,2); end
                    end
                    shuf_to(g) = nanmean(tmp);
                end
            end
        end

        spaceCeilBase(kk) = nanmean(ceilB_sp);
        spaceCeilStim(kk) = nanmean(ceilS_sp);
        toneCeilBase(kk)  = nanmean(ceilB_to);
        toneCeilStim(kk)  = nanmean(ceilS_to);

        spaceBaseStimMatched(kk) = nanmean(obs_sp);
        toneBaseStimMatched(kk)  = nanmean(obs_to);

        spaceBaseStimShuf(kk) = nanmean(shuf_sp);
        toneBaseStimShuf(kk)  = nanmean(shuf_to);

        %% Detect fields (unchanged)
        Field_Info = detectFields(spaceMap(kk,:));
        if isempty(Field_Info), spaceField(kk,1) = 0; else, spaceField(kk,1) = 1; end
        
        Field_Info = detectFields(toneMap(kk,:));
        if isempty(Field_Info), toneField(kk,1) = 0; else, toneField(kk,1) = 1; end      
        
        Field_Info = detectFields(spaceMapStim(kk,:));
        if isempty(Field_Info), spaceFieldStim(kk,1) = 0; else, spaceFieldStim(kk,1) = 1; end
        
        Field_Info = detectFields(toneMapStim(kk,:));
        if isempty(Field_Info), toneFieldStim(kk,1) = 0; else, toneFieldStim(kk,1) = 1; end           

    end
    
    % Append originals
    AllcellType = [AllcellType;cellType];
    AllspatialCorr = [AllspatialCorr;spatialCorr];
    AlltoneCorr = [AlltoneCorr;toneCorr];
    AllspatialCorrStim = [AllspatialCorrStim;spatialCorrStim];
    AlltoneCorrStim = [AlltoneCorrStim;toneCorrStim];
    AllspaceMap = [AllspaceMap;spaceMap];
    AlltoneMap = [AlltoneMap;toneMap];
    AllspaceMapStim = [AllspaceMapStim;spaceMapStim];
    AlltoneMapStim = [AlltoneMapStim;toneMapStim];
    AllspaceField = [AllspaceField;spaceField];
    AlltoneField = [AlltoneField;toneField];     
    AllspaceFieldStim = [AllspaceFieldStim;spaceFieldStim];
    AlltoneFieldStim = [AlltoneFieldStim;toneFieldStim];    
    AllspatialInfo = [AllspatialInfo;spatialInfo];
    AlltoneInfo = [AlltoneInfo;toneInfo];
    AllspatialInfoStim = [AllspatialInfoStim;spatialInfoStim];
    AlltoneInfoStim = [AlltoneInfoStim;toneInfoStim];

    % Append NEW fair metrics
    AllspaceCeilBase = [AllspaceCeilBase; spaceCeilBase];
    AllspaceCeilStim = [AllspaceCeilStim; spaceCeilStim];
    AlltoneCeilBase  = [AlltoneCeilBase;  toneCeilBase];
    AlltoneCeilStim  = [AlltoneCeilStim;  toneCeilStim];

    AllspaceBaseStimMatched = [AllspaceBaseStimMatched; spaceBaseStimMatched];
    AlltoneBaseStimMatched  = [AlltoneBaseStimMatched;  toneBaseStimMatched];

    AllspaceBaseStimShuf = [AllspaceBaseStimShuf; spaceBaseStimShuf];
    AlltoneBaseStimShuf  = [AlltoneBaseStimShuf;  toneBaseStimShuf];
    
    clear cellType spatialCorr toneCorr spatialCorrStim toneCorrStim spaceMap toneMap spaceMapStim toneMapStim ...
        spaceField toneField spaceFieldStim toneFieldStim spatialInfo toneInfo spatialInfoStim toneInfoStim ...
        spaceCeilBase spaceCeilStim toneCeilBase toneCeilStim spaceBaseStimMatched toneBaseStimMatched ...
        spaceBaseStimShuf toneBaseStimShuf
    
end

% =========================
% Your original plotting (unchanged up to here)
% =========================
figure
set(gcf,'Color','w')
set(gcf,'renderer','painters')    
set(gcf,'Position',[408 38 1107 946])

YlGnBu=cbrewer('seq', 'YlGnBu', 11);
fractTune = [];
linPos = linspace(1,122,50);
linTone = linspace(2000,22000,50);

activeCells = AllcellType & (nanmean(AllspaceMap,2)>0.2);
idxMaps{1} = AllcellType & ((AllspaceField & AllspatialCorr>0.1));
idxMaps{2} = AllcellType & ((AlltoneField & AlltoneCorr>0.1));
idxMaps{3} = AllcellType & AllspaceFieldStim & AllspatialCorrStim>0.1;
idxMaps{4} = AllcellType & AlltoneFieldStim & AlltoneCorrStim>0.1;
fractTune = round(([sum(idxMaps{1}) sum(idxMaps{3}) sum(idxMaps{2}) sum(idxMaps{4})]/sum(activeCells))*100); 

subplot(5,4,1)
bar(fractTune)
text(1:length(fractTune),fractTune,num2str(fractTune'),'vert','bottom','horiz','center'); 
box off

clear idxMaps
idxMaps{1} = AllcellType & ((AllspaceField & AllspatialCorr>0.1) | (AllspaceFieldStim & AllspatialCorrStim>0.1));
idxMaps{2} = AllcellType & ((AlltoneField & AlltoneCorr>0.1) | (AlltoneFieldStim & AlltoneCorrStim>0.1));

for ii = 1:2
   
   if ii == 1
       infoData{1} = AllspatialInfo(idxMaps{ii});
       infoData{2} = AllspatialInfoStim(idxMaps{ii});
       subplot(5,4,2)
       stats.spatialInfo = groupStats(infoData,[],'inAxis',true,'labelSummary',false,'repeatedMeasures',true);       
       title('Spatial Information')
   elseif ii == 2
       infoData{1} = AlltoneInfo(idxMaps{ii});
       infoData{2} = AlltoneInfoStim(idxMaps{ii});
       subplot(5,4,3)
       stats.toneInfo = groupStats(infoData,[],'inAxis',true,'labelSummary',false,'repeatedMeasures',true);         
       title('Tone Information')
   end

   selectedtoneMap = AlltoneMap(idxMaps{ii},:);
   selectedspaceMap = AllspaceMap(idxMaps{ii},:);
   selectedtoneMapStim = AlltoneMapStim(idxMaps{ii},:);
   selectedspaceMapStim = AllspaceMapStim(idxMaps{ii},:);
   
   [~,idxSpace] = max(selectedspaceMap,[],2);        
   [~,idxTone]  = max(selectedtoneMap,[],2);
   [~,idxSpaceStim] = max(selectedspaceMapStim,[],2);        
   [~,idxToneStim]  = max(selectedtoneMapStim,[],2);   
   
   normtoneMap = (selectedtoneMap-nanmean(selectedtoneMap,2))./nanstd(selectedtoneMap,[],2);
   normspaceMap = (selectedspaceMap-nanmean(selectedspaceMap,2))./nanstd(selectedspaceMap,[],2);
   normtoneMapStim = (selectedtoneMapStim-nanmean(selectedtoneMapStim,2))./nanstd(selectedtoneMapStim,[],2);
   normspaceMapStim = (selectedspaceMapStim-nanmean(selectedspaceMapStim,2))./nanstd(selectedspaceMapStim,[],2);
    
   if ii == 1            
       [~,sortidx] = sort(idxSpace,'ascend');
       [~,sortidxstim] = sort(idxSpaceStim,'ascend');
   else            
       [~,sortidx] = sort(idxTone,'ascend'); 
       [~,sortidxstim] = sort(idxToneStim,'ascend'); 
   end
    
   colormap(YlGnBu)
   subplot(5,4,8*(ii-1)+5)
   imagesc(linPos, 1:length(sortidx),normspaceMap(sortidx,:))
   ylabel('Cell ID'); xlabel('Position'); title('Task spatial'); caxis([-1 4])
    
   subplot(5,4,8*(ii-1)+6)
   imagesc(linPos, 1:length(sortidx),normspaceMapStim(sortidx,:))
   ylabel('Cell ID'); xlabel('Position'); title('Task spatial stim'); caxis([-1 4])
    
   subplot(5,4,8*(ii-1)+7)
   imagesc(linTone, 1:length(sortidx),normtoneMap(sortidx,:))
   ylabel('Cell ID'); xlabel('Frequency'); title('Task auditory'); caxis([-1 4])
    
   subplot(5,4,8*(ii-1)+8)
   imagesc(linTone, 1:length(sortidx),normtoneMapStim(sortidx,:))
   ylabel('Cell ID'); xlabel('Frequency'); title('Task auditory stim'); caxis([-1 4])

   subplot(5,4,8*(ii-1)+9)
   imagesc(linPos, 1:length(sortidxstim),normspaceMap(sortidxstim,:))
   ylabel('Cell ID'); xlabel('Position'); title('Task spatial'); caxis([-1 4])
    
   subplot(5,4,8*(ii-1)+10)
   imagesc(linPos, 1:length(sortidxstim),normspaceMapStim(sortidxstim,:))
   ylabel('Cell ID'); xlabel('Position'); title('Task spatial stim'); caxis([-1 4])
    
   subplot(5,4,8*(ii-1)+11)
   imagesc(linTone, 1:length(sortidxstim),normtoneMap(sortidxstim,:))
   ylabel('Cell ID'); xlabel('Frequency'); title('Task auditory'); caxis([-1 4])
    
   subplot(5,4,8*(ii-1)+12)
   imagesc(linTone, 1:length(sortidxstim),normtoneMapStim(sortidxstim,:))
   ylabel('Cell ID'); xlabel('Frequency'); title('Task auditory stim'); caxis([-1 4])
end

% ===== NEW: quick fair comparison panels =====
pyr = AllcellType == 1;

figure
subplot(1,2,1); hold on; box off
d{1} = AllspaceCeilBase(pyr);
d{2} = AllspaceCeilStim(pyr);
d{3} = AllspaceBaseStimMatched(pyr);
d{4} = AllspaceBaseStimShuf(pyr);
groupStats(d,[],'inAxis',true,'labelSummary',false);
title('Spatial (FAIR): ceilBase / ceilStim / base-stim / shuffle')
ylabel('Correlation')

subplot(1,2,2); hold on; box off
d{1} = AlltoneCeilBase(pyr);
d{2} = AlltoneCeilStim(pyr);
d{3} = AlltoneBaseStimMatched(pyr);
d{4} = AlltoneBaseStimShuf(pyr);
groupStats(d,[],'inAxis',true,'labelSummary',false);
title('Tone (FAIR): ceilBase / ceilStim / base-stim / shuffle')
ylabel('Correlation')

end

% =========================
% Helper: stack per-trial maps into matrix (nTrials x 50)
% rateMapsPerTrial is a cell array where rateMapsPerTrial{t} is 1x50 (or 50x1)
% =========================
function X = stackTrialMaps(rateMapsPerTrial, trIdx)
    X = [];
    if isempty(trIdx), return; end
    keep = false(size(trIdx));
    tmp = cell(numel(trIdx),1);
    for i = 1:numel(trIdx)
        t = trIdx(i);
        if t > numel(rateMapsPerTrial) || isempty(rateMapsPerTrial{t})
            continue
        end
        v = rateMapsPerTrial{t};
        v = v(:)'; % force row
        if numel(v) < 50
            v(50) = nan; % pad if needed
        end
        tmp{i} = v(1:50);
        keep(i) = true;
    end
    tmp = tmp(keep);
    if isempty(tmp), return; end
    X = cell2mat(tmp);
end

% =========================
% Helper: split-half ceiling correlation for trial x bin matrix
% =========================
function c = splitHalfCeil(X, nSplit)
    n = size(X,1);
    if n < 4
        c = nan;
        return
    end
    cc = nan(nSplit,1);
    for r = 1:nSplit
        rp = randperm(n);
        h1 = rp(1:floor(n/2));
        h2 = rp(floor(n/2)+1:end);
        m1 = nanmean(X(h1,:),1);
        m2 = nanmean(X(h2,:),1);
        a = corrcoef(m1, m2, 'rows','complete');
        if size(a,1)==2
            cc(r) = a(1,2);
        end
    end
    c = nanmean(cc);
end
