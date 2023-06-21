function compareActionNoToneTrials

%sess= {'IZ44\Final\IZ44_220827_sess4','IZ39\Final\IZ39_220629_sess12'};
%
sess={'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
    'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
    'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   
    'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...
    'IZ40\Final\IZ40_220708_sess17','IZ40\Final\IZ40_220714_sess18',...
    'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
    'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
    'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
    'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    
    'IZ44\Final\IZ44_220827_sess4', 'IZ44\Final\IZ44_220828_sess5',...
    'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...
    'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
    'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',...
    }; 

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';

AllrateMapLin = [];
AllrateMapLin1 = [];
AllrateMapLin2 = [];
AllrateMapLinEnd = [];
AllrateMapSpace = [];
AllrateMaptone = [];
AlllinField = [];
AlllinFieldEnd = [];

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
     
    linMapInit = []; 
    linMapInit1 = [];
    linMapInit2 = [];
    linMapEnd =[];
    spaceMap1 = [];
    toneMap1 = [];
    linField = [];
    linFieldEnd = [];
    
    for kk=1:length(cell_metrics.UID)
    %% Extract pyramidal cells and calculate rate map correlations and information        
        if strcmp(cell_metrics.putativeCellType{kk},'Pyramidal Cell') == 1 
            cellType = 1;
        else
            cellType = 0;
        end
        % Average spatial and tone correlation
        dataMat =[];
        dataMatTone = [];
        for jj = 2:7    
            dataMat = [dataMat;firingMaps.forward.rateMaps{kk}{jj}];
            a = fillmissing(firingMaps.tone.rateMaps{kk}{jj},'linear');
            dataMatTone = [dataMatTone;a];
        end         
        toneMap = nanmean(dataMatTone,1);        
        corrTone = []; 
        corrSpace = [];
        
        for pp = 1:6
           for jj = (pp+1):6
                a = corrcoef(dataMat(pp,:),dataMat(jj,:),'rows','complete');
               corrSpace = [corrSpace a(1,2)];               
               a = corrcoef(dataMatTone(pp,:),dataMatTone(jj,:),'rows','complete');         
               corrTone = [corrTone a(1,2)];
           end
        end     
        toneCorr = nanmean(corrTone);  
        spaceCorr = nanmean(corrSpace);
        
        %% Detect fields
        Field_Info = detectFields(toneMap);
        if isempty(Field_Info)
            toneField = 0;
        else 
            toneField = 1;
        end                
        
        [~,idxMax] = max(toneMap);
        % If its a trajectory cell
        if cellType == 1 && (toneField == 1) && (toneCorr > 0.1) && idxMax>15
            
            spaceMap1 = [spaceMap1; nanmean(dataMat,1)];        
            toneMap1 = [toneMap1; nanmean(dataMatTone,1)];    
            %its a trajectory cell!
            %Look at the linear track trials split up by lickLoc            
            linMapInit = [linMapInit;firingMaps.forward.rateMaps{kk}{1}];
            linMapInit1 = [linMapInit1;firingMaps.forward.rateMaps{kk}{26}];
            linMapInit2 = [linMapInit2;firingMaps.forward.rateMaps{kk}{27}];
            
            if length(firingMaps.forward.rateMaps{kk})==28
                linMapEnd = [linMapEnd;firingMaps.forward.rateMaps{kk}{28}];
            else
                temp(1:size(dataMat,2)) = nan;
                linMapEnd = [linMapEnd; temp];
            end
            
            Field_Info = detectFields(firingMaps.forward.rateMaps{kk}{1});
            if isempty(Field_Info)
                linField = [linField;0];
            else 
                linField = [linField;1];
            end     
            
            if length(firingMaps.forward.rateMaps{kk})==28
                Field_Info = detectFields(firingMaps.forward.rateMaps{kk}{28});
                if isempty(Field_Info)
                    linFieldEnd = [linFieldEnd;0];
                else 
                    linFieldEnd = [linFieldEnd;1];
                end   
            else
                 linFieldEnd = [linFieldEnd;nan];   
            end
        end
            
    end
    
    AllrateMapLin = [AllrateMapLin;linMapInit];
    AllrateMapLin1 = [AllrateMapLin1;linMapInit1];
    AllrateMapLin2 = [AllrateMapLin2;linMapInit2];
    AllrateMapLinEnd = [AllrateMapLinEnd;linMapEnd];
    AllrateMapSpace = [AllrateMapSpace;spaceMap1];
    AllrateMaptone = [AllrateMaptone;toneMap1];
    AlllinField = [AlllinField;linField];
    AlllinFieldEnd = [AlllinFieldEnd;linFieldEnd];
    
    clear linMapInit linMapEnd spaceMap1 toneMap1 linField linFieldEnd linMapInit1 linMapInit2
end


normlinMap = (AllrateMapLin-nanmean(AllrateMapLin,2))./nanstd(AllrateMapLin,[],2);
normlinMap1 = (AllrateMapLin1-nanmean(AllrateMapLin1,2))./nanstd(AllrateMapLin1,[],2);
normlinMap2 = (AllrateMapLin2-nanmean(AllrateMapLin2,2))./nanstd(AllrateMapLin2,[],2);
normlinMapEnd = (AllrateMapLinEnd-nanmean(AllrateMapLinEnd,2))./nanstd(AllrateMapLinEnd,[],2);
normspaceMap = (AllrateMapSpace-nanmean(AllrateMapSpace,2))./nanstd(AllrateMapSpace,[],2);
normtoneMap = (AllrateMaptone-nanmean(AllrateMaptone,2))./nanstd(AllrateMaptone,[],2);
        
figure
set(gcf,'Renderer','painters')
set(gcf,'Color','w')
linPos = linspace(1,122,50);
YlGnBu=cbrewer('seq', 'YlGnBu', 11);
colormap(YlGnBu)

subplot(3,3,1)
[maxLin,idxLin] = max(AllrateMapLin,[],2);
[~,sortidx] = sort(idxLin,'ascend');
imagesc(linPos, 1:size(normlinMap,1),normlinMap(sortidx,:))
caxis([-1 3])
title('All lin trials')

subplot(3,3,2)
imagesc(linPos, 1:size(normlinMap1,1),normlinMap1(sortidx,:))
caxis([-1 3])
title('Lin 1st half')

subplot(3,3,3)
imagesc(linPos, 1:size(normlinMap2,1),normlinMap2(sortidx,:))
caxis([-1 3])
title('Lin 2nd half')

subplot(3,3,4)
imagesc(linPos, 1:size(normlinMapEnd,1),normlinMapEnd(sortidx,:))
caxis([-1 3])
title('Lin End')


subplot(3,3,5)
imagesc(linPos, 1:size(normspaceMap,1),normspaceMap(sortidx,:))
caxis([-1 3])
title('Space')

subplot(3,3,6)
imagesc(linPos, 1:size(normtoneMap,1),normtoneMap(sortidx,:))
caxis([-1 3])
title('Tone')

subplot(3,3,7)
bar([sum(AlllinField)./length(AlllinField) sum(AlllinFieldEnd)./length(AlllinFieldEnd)])

selecttoneMap = normtoneMap(AlllinField==1,:);
selectspaceMap = normspaceMap(AlllinField==1,:);
selectlinMap = normlinMap(AlllinField==1,:);
selectlinMap1 = normlinMap1(AlllinField==1,:);
selectlinMap2 = normlinMap2(AlllinField==1,:);
[maxLin,idxLin] = max(selectlinMap,[],2);
[~,sortidx] = sort(idxLin,'ascend');

subplot(3,3,8)
imagesc(linPos, 1:size(selectlinMap1,1),selectlinMap1(sortidx,:))
caxis([-1 3])

subplot(3,3,9)
imagesc(linPos, 1:size(selectlinMap2,1),selectlinMap2(sortidx,:))
caxis([-1 3])

% figure
% histogram(idxLin)

saveas(gcf,strcat(expPath,'Compiled\ToneCellsLinearTrials.png'));
saveas(gcf,strcat(expPath,'Compiled\ToneCellsLinearTrials.eps'),'epsc');
saveas(gcf,strcat(expPath,'Compiled\ToneCellsLinearTrials.fig'));

end

function Field_Info = detectFields(SmoothedFiringRate)
    minFieldSize = 10;
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