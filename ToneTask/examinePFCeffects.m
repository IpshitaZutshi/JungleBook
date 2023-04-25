function examinePFCeffects

sess = {'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   
    'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...
    'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',...
    'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16'}; 
    
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';

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

for ii = 1:length(sess)
    %% Load files
    cd(strcat(expPath,sess{ii}))    
    file = dir(['*.rateMapsAvg.cellinfo.mat']);
    load(file(1).name);
    file = dir(['*.rateMapsTrial.cellinfo.mat']);
    trialMaps = load(file(1).name);    
    file = dir(['*.Tracking.Behavior.mat']);
    load(file(1).name);
    file = dir(['*TrialBehavior.Behavior.mat']);
    load(file(1).name);
    [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts', true);
    file = dir(['*cell_metrics.cellinfo.mat']);
    load(file(1).name);
     
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
    
    for kk=1:length(cell_metrics.UID)
    %% Extract pyramidal cells and calculate rate map correlations and information        
        if strcmp(cell_metrics.putativeCellType{kk},'Pyramidal Cell') == 1 
            cellType(kk,1) = 1;
        else
            cellType(kk,1) = 0;
        end
        % Average spatial and tone correlation
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
        toneMap(kk,:) = nanmean(dataMatTone,1);        
        
        spaceMapStim(kk,:) = nanmean(dataMatStim,1);
        toneMapStim(kk,:) = nanmean(dataMatToneStim,1);                
        
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
        
        %% Detect fields
        Field_Info = detectFields(spaceMap(kk,:));
        if isempty(Field_Info)
            spaceField(kk,1) = 0;
        else 
            spaceField(kk,1) = 1;
        end
        
        Field_Info = detectFields(toneMap(kk,:));
        if isempty(Field_Info)
            toneField(kk,1) = 0;
        else 
            toneField(kk,1) = 1;
        end      
        
        Field_Info = detectFields(spaceMapStim(kk,:));
        if isempty(Field_Info)
            spaceFieldStim(kk,1) = 0;
        else 
            spaceFieldStim(kk,1) = 1;
        end
        
        Field_Info = detectFields(toneMapStim(kk,:));
        if isempty(Field_Info)
            toneFieldStim(kk,1) = 0;
        else 
            toneFieldStim(kk,1) = 1;
        end           
            
    end
    
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
    
    clear cellType spatialCorr toneCorr spatialCorrStim toneCorrStim spaceMap toneMap spaceMapStim toneMapStim ...
        spaceField toneField spaceFieldStim toneFieldStim spatialInfo toneInfo spatialInfoStim toneInfoStim
    
end

figure
set(gcf,'Color','w')
set(gcf,'renderer','painters')    
set(gcf,'Position',[408 38 1107 946])

YlGnBu=cbrewer('seq', 'YlGnBu', 11);
fractTune = [];
linPos = linspace(1,122,50);
linTone = linspace(2000,22000,50);

activeCells = AllcellType & (nanmean(AllspaceMap,2)>0.2);
idxMaps{1} = AllcellType & ((AllspaceField & AllspatialCorr>0.1)); % Idx for spatial cells
idxMaps{2} = AllcellType & ((AlltoneField & AlltoneCorr>0.1));  % Idx for tone cells 
idxMaps{3} = AllcellType & AllspaceFieldStim & AllspatialCorrStim>0.1; % Idx for spatial cells
idxMaps{4} = AllcellType & AlltoneFieldStim & AlltoneCorrStim>0.1;  % Idx for tone cells 
fractTune = round(([sum(idxMaps{1}) sum(idxMaps{3}) sum(idxMaps{2}) sum(idxMaps{4})]/sum(activeCells))*100); 

subplot(5,4,1)
bar(fractTune)
text(1:length(fractTune),fractTune,num2str(fractTune'),'vert','bottom','horiz','center'); 
box off

clear idxMaps
idxMaps{1} = AllcellType & ((AllspaceField & AllspatialCorr>0.1) | (AllspaceFieldStim & AllspatialCorrStim>0.1)); % Idx for spatial cells
idxMaps{2} = AllcellType & ((AlltoneField & AlltoneCorr>0.1) | (AlltoneFieldStim & AlltoneCorrStim>0.1));  % Idx for tone cells 
% idxMaps{3} = AllcellType & AllspaceFieldStim & AllspatialCorrStim>0.1; % Idx for spatial cells
% idxMaps{4} = AllcellType & AlltoneFieldStim & AlltoneCorrStim>0.1;  % Idx for tone cells 


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
   
   [maxSpace,idxSpace] = max(selectedspaceMap,[],2);        
   [maxTone,idxTone] = max(selectedtoneMap,[],2);
   
   [maxSpaceStim,idxSpaceStim] = max(selectedspaceMapStim,[],2);        
   [maxToneStim,idxToneStim] = max(selectedtoneMapStim,[],2);   
   
    %Z score normalize
    normtoneMap = (selectedtoneMap-nanmean(selectedtoneMap,2))./nanstd(selectedtoneMap,[],2);
    normspaceMap = (selectedspaceMap-nanmean(selectedspaceMap,2))./nanstd(selectedspaceMap,[],2);
    normtoneMapStim = (selectedtoneMapStim-nanmean(selectedtoneMapStim,2))./nanstd(selectedtoneMapStim,[],2);
    normspaceMapStim = (selectedspaceMapStim-nanmean(selectedspaceMapStim,2))./nanstd(selectedspaceMapStim,[],2);
    
    sortidx = [];
    sortidxstim = [];
    if ii == 1            
        [~,sortidx] = sort(idxSpace,'ascend');
        [~,sortidxstim] = sort(idxSpaceStim,'ascend');
    elseif ii == 2            
        [~,sortidx] = sort(idxTone,'ascend'); 
        [~,sortidxstim] = sort(idxToneStim,'ascend'); 
    end
    
    colormap(YlGnBu)
    subplot(5,4,8*(ii-1)+5)
    imagesc(linPos, 1:length(sortidx),normspaceMap(sortidx,:))
    ylabel('Cell ID')
    xlabel('Frequency')  
    title('Task spatial')
    caxis([-1 4])
    
    subplot(5,4,8*(ii-1)+6)
    imagesc(linPos, 1:length(sortidx),normspaceMapStim(sortidx,:))
    ylabel('Cell ID')
    xlabel('Position')
    caxis([-1 4])
    title('Task spatial stim')  
    
    subplot(5,4,8*(ii-1)+7)
    imagesc(linTone, 1:length(sortidx),normtoneMap(sortidx,:))
    ylabel('Cell ID')
    xlabel('Position')
    caxis([-1 4])
    title('Task auditory')  
    
    subplot(5,4,8*(ii-1)+8)
    imagesc(linTone, 1:length(sortidx),normtoneMapStim(sortidx,:))
    ylabel('Cell ID')
    xlabel('Position')
    caxis([-1 4])
    title('Task auditory stim')  
   
    subplot(5,4,8*(ii-1)+9)
    imagesc(linPos, 1:length(sortidxstim),normspaceMap(sortidxstim,:))
    ylabel('Cell ID')
    xlabel('Frequency')  
    title('Task spatial')
    caxis([-1 4])
    
    subplot(5,4,8*(ii-1)+10)
    imagesc(linPos, 1:length(sortidxstim),normspaceMapStim(sortidxstim,:))
    ylabel('Cell ID')
    xlabel('Position')
    caxis([-1 4])
    title('Task spatial stim')  
    
    subplot(5,4,8*(ii-1)+11)
    imagesc(linTone, 1:length(sortidxstim),normtoneMap(sortidxstim,:))
    ylabel('Cell ID')
    xlabel('Position')
    caxis([-1 4])
    title('Task auditory')  
    
    subplot(5,4,8*(ii-1)+12)
    imagesc(linTone, 1:length(sortidxstim),normtoneMapStim(sortidxstim,:))
    ylabel('Cell ID')
    xlabel('Position')
    caxis([-1 4])
    title('Task auditory stim')      
end

saveas(gcf,strcat(expPath,'Compiled\pFCMaps.png'));
saveas(gcf,strcat(expPath,'Compiled\pFCMaps.eps'),'epsc');
saveas(gcf,strcat(expPath,'Compiled\pFCMaps.fig'));
end

function Field_Info = detectFields(SmoothedFiringRate)
    minFieldSize = 5;
    maxFieldSize = 40;
    % Pad on each end with zeros for edge effects
    SmoothedFiringRate = [0 0 SmoothedFiringRate 0 0];
    [peakValues, peakLocations] = findpeaks(SmoothedFiringRate, 'minpeakheight',5, 'minpeakdistance', 10);
    Field_Info = [];
    for j = 1:length(peakLocations)
        FieldPeak = peakLocations(j);
        % FieldPeak must be 5 Hz or more
        if peakValues(j) < 5, continue, end
        LookForward = FieldPeak+1:length(SmoothedFiringRate);
        LookBack = 1:FieldPeak-1;
        PercentPeakRate_Forward = SmoothedFiringRate(LookForward)./peakValues(j);
        PercentPeakRate_Back = SmoothedFiringRate(LookBack)./peakValues(j);
        tempInd1 = find(PercentPeakRate_Forward < .2);
        if isempty(tempInd1), continue, end
        FieldEnd = FieldPeak+tempInd1(1); % this is the first bin forward of the animal that has a FR less than 20% of the peak rate
        tempInd2 = find(PercentPeakRate_Back < .2);
        if isempty(tempInd2)
            FieldStart = 1;
        else
            FieldStart = tempInd2(end); % this is the first bin forward of the animal that has a FR less than 20% of the peak rate
        end
        % Field must be more than 10 cm and less than 80cm
        if FieldEnd>FieldStart && FieldEnd-FieldStart < minFieldSize, continue, end
        if FieldEnd>FieldStart && FieldEnd-FieldStart > maxFieldSize, continue, end        
        if FieldEnd>50
            FieldEnd = 50;
        else
            FieldEnd = FieldEnd-2;
        end
        if FieldStart<3
            FieldStart = 1;
        else
            FieldStart = FieldStart-2;
        end
        Field_Info = [Field_Info;FieldStart, FieldEnd, FieldPeak];
    end

end