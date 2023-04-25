function examineTrajectoryCells

%sess= {'IZ44\Final\IZ44_220827_sess4','IZ39\Final\IZ39_220629_sess12'};
%
sess= {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220629_sess12',...
    'IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220714_sess18',...
    'IZ39\Final\IZ39_220702_sess14',...
    'IZ40\Final\IZ40_220707_sess16','IZ40\Final\IZ40_220714_sess18'...
    'IZ43\Final\IZ43_220828_sess4','IZ43\Final\IZ43_220919_sess14',...
    'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220901_sess8',...
    'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11',...
    'IZ44\Final\IZ44_220827_sess4','IZ44\Final\IZ44_220830_sess7',...
    'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   
    'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...
    'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',...
    'IZ40\Final\IZ40_220705_sess15'}; 


expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';

for cc = 1:6
    AllrateMapLin{cc} = [];
    AllrateMapSpace{cc} = [];
end
supDeep = [];

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
        if cellType == 1 && (toneField == 1) && (toneCorr > 0.1) && idxMax>40
            %its a trajectory cell!
            %Look at the linear track trials split up by lickLoc            
            for cc = 1:6
                idxLick = find(behavTrials.lickLoc(1:(end-1))==(cc-1) & behavTrials.linTrial(1:(end-1))==1);
                rateMap(1:50) = nan;
                for ll = 1:length(idxLick)
                    rateMap = [rateMap;trialMaps.firingMaps.forward.rateMaps{kk,1}{1,idxLick(ll)}];
                end
                AllrateMapLin{cc} = [AllrateMapLin{cc}; nanmean(rateMap,1)];
                AllrateMapSpace{cc} = [AllrateMapSpace{cc};firingMaps.forward.rateMaps{kk,1}{1,cc+1}];
                clear rateMap
            end                        
        end
            
    end
end

figure
set(gcf,'Renderer','painters')
set(gcf,'Color','w')
linPos = linspace(1,122,50);
YlGnBu=cbrewer('seq', 'YlGnBu', 11);
colormap(YlGnBu)
rateRange = [1 6;1 15;1 24;1 34;1 43;1 49];
for cc = 1:6
    normLinMap = (AllrateMapLin{cc}-nanmean(AllrateMapLin{cc},2))./nanstd(AllrateMapLin{cc},[],2);
    normSpaceMap = (AllrateMapSpace{cc}-nanmean(AllrateMapSpace{cc},2))./nanstd(AllrateMapSpace{cc},[],2);
    subplot(3,6,cc)
    idxLin = find(~isnan(AllrateMapLin{cc}(:,1)));
    linMap = normLinMap(idxLin,:);
%     [~,sortidx] = sort(idxLin,'ascend');    
%     [~,idxLin] = max(linMap,[],2);
%     [~,sortidx] = sort(idxLin,'ascend');      
%    h = imagesc(linMap(sortidx,:));
    h = imagesc(linMap);
    set(h, 'AlphaData', ~isnan(linMap))
    caxis([-1 4])
    subplot(3,6,cc+12)
    plot(linPos,nanmean(normLinMap,1),'LineWidth',1.4,'Color',[0.5 0.5 0.5])
    xlim([1 122])
    hold on    
    subplot(3,6,cc+6)
    [~,idxLin] = max(AllrateMapSpace{cc},[],2);
    [~,sortidx] = sort(idxLin,'ascend');       
    h = imagesc(normSpaceMap(sortidx,:));
    set(h, 'AlphaData', ~isnan(normSpaceMap(sortidx,:)))
    caxis([-1 4])    
    subplot(3,6,cc+12)
    plot(linPos,nanmean(normSpaceMap,1),'LineWidth',1.4,'Color',[187/243 86/243 149/243])
    xlim([1 122])
end

saveas(gcf,strcat(expPath,'Compiled\LinearTrials_Tone.png'));
saveas(gcf,strcat(expPath,'Compiled\LinearTrials_Tone.eps'),'epsc');
saveas(gcf,strcat(expPath,'Compiled\LinearTrials_Tone.fig'));

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