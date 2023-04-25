function examineTrajectoryDeepSuperficial

%sess= {'IZ44\Final\IZ44_220827_sess4','IZ39\Final\IZ39_220629_sess12'};
%
sess= {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220629_sess12',...
    'IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220714_sess18','IZ39\Final\IZ39_220702_sess14'...
    'IZ40\Final\IZ40_220707_sess16','IZ40\Final\IZ40_220714_sess18'...
    'IZ43\Final\IZ43_220828_sess4','IZ43\Final\IZ43_220919_sess14',...
    'IZ44\Final\IZ44_220827_sess4','IZ44\Final\IZ44_220830_sess7',...
    'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   
    'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...
    'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',...
    'IZ40\Final\IZ40_220705_sess15'}; 

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';

trajectCell = [];
spaceCell = [];
deepSup  = [];

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
    if ~isempty(dir('*.deepSuperficialfromRipple.channelinfo.mat'))        
        file = dir(['*.deepSuperficialfromRipple.channelinfo.mat']);
        load(file(1).name);
    else 
        continue
    end
     
    
    for kk=1:length(cell_metrics.UID)
    %% Extract pyramidal cells and calculate rate map correlations and information        
        if strcmp(cell_metrics.putativeCellType{kk},'Pyramidal Cell') == 1 
            cellType = 1;
        else
            cellType = 0;
        end

        if strcmp(deepSuperficialfromRipple.channelClass{cell_metrics.maxWaveformCh1(kk)},'Deep') == 1
            deepSup  = [deepSup 1];
        elseif strcmp(deepSuperficialfromRipple.channelClass{cell_metrics.maxWaveformCh1(kk)},'Superficial') == 1
            deepSup  = [deepSup 0];
        else
            deepSup  = [deepSup nan];
        end        
        % Average spatial and tone correlation
        dataMat =[];
        dataMatTone = [];
        for jj = 2:7    
            dataMat = [dataMat;firingMaps.forward.rateMaps{kk}{jj}];
            a = fillmissing(firingMaps.tone.rateMaps{kk}{jj},'linear');
            dataMatTone = [dataMatTone;a];
        end         
        
        spaceMap = nanmean(dataMat,1);     
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
 
        Field_Info = detectFields(spaceMap);
        if isempty(Field_Info)
            spaceField = 0;
        else 
            spaceField = 1;
        end 
        
        % If its a trajectory cell
        if (cellType == 1) && (toneField == 1) && toneCorr>0.1
            trajectCell = [trajectCell 1];
        else
            trajectCell = [trajectCell 0];
        end
 
        % If its a place cell
        if (cellType == 1) && (spaceField == 1) && spaceCorr>0.1
            spaceCell = [spaceCell 1];
        else
            spaceCell = [spaceCell 0];
        end
        
            
    end
end

figure
set(gcf,'Renderer','painters')
set(gcf,'Color','w')

deepSupSpace = deepSup(logical(spaceCell));
deepSupTone = deepSup(logical(trajectCell));

sum(deepSup)/length(deepSup)

prop = [sum(deepSupSpace)/length(deepSupSpace) (length(deepSupSpace)-sum(deepSupSpace))/length(deepSupSpace);...
    sum(deepSupTone)/length(deepSupTone) (length(deepSupTone)-sum(deepSupTone))/length(deepSupTone)]

% saveas(gcf,strcat(expPath,'Compiled\LinearTrials_Tone.png'));
% saveas(gcf,strcat(expPath,'Compiled\LinearTrials_Tone.eps'),'epsc');
% saveas(gcf,strcat(expPath,'Compiled\LinearTrials_Tone.fig'));

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