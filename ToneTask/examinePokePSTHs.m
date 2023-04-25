function examinePokePSTHs

% Look at 2 types of pokes - 
%   1) Return run, in any of the middle ports while facing forward
%   2) Return run, in any of the middle ports while facing backwards

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

for tt = 1:3
    psthReward{tt} = [];
end
sumPoke1 = [];
sumPoke2 = [];
sumLicks = [];

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
    file = dir(['*spikes.cellinfo.mat']);
    load(file(1).name);       
    [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts', true);
    file = dir(['*cell_metrics.cellinfo.mat']);
    load(file(1).name);
    file = dir(['*.DigitalIn.Events.mat']);
    load(file(1).name);
    
    if size(digitalIn.ints{4},1)==2
        d = [digitalIn.ints{4} digitalIn.ints{5} digitalIn.ints{6} ...
            digitalIn.ints{7} digitalIn.ints{8}];        
    else
        d = [digitalIn.ints{4}' digitalIn.ints{5}' digitalIn.ints{6}'...
            digitalIn.ints{7}' digitalIn.ints{8}'];        
    end
    
    intsPeriods =[];
    intsPeriods(1,1) = d(1,1); % find stimulation intervals
    intPeaks =find(abs(diff(d(1,:)))>0.75);
    for jj = 1:length(intPeaks)
        intsPeriods(jj,2) = d(2,intPeaks(jj));
        intsPeriods(jj+1,1) = d(1,intPeaks(jj)+1);
    end
    intsPeriods(end,2) = d(2,end); 
    tsLick2 = intsPeriods;    
        
    for kk=1:length(cell_metrics.UID)
    %% Extract pyramidal cells and calculate rate map correlations and information 
        
        if strcmp(cell_metrics.putativeCellType{kk},'Pyramidal Cell') == 1 
            cellType = 1;
        else
            cellType = 0;
        end
        
        dataMatTone = [];
        for jj = 2:7    
            a = fillmissing(firingMaps.tone.rateMaps{kk}{jj},'linear');
            dataMatTone = [dataMatTone;a];
        end                 
        toneMap = nanmean(dataMatTone,1);        
        
        corrTone = []; 
        for pp = 1:6
           for jj = (pp+1):6           
               a = corrcoef(dataMatTone(pp,:),dataMatTone(jj,:),'rows','complete');         
               corrTone = [corrTone a(1,2)];
           end
        end     
        toneCorr = nanmean(corrTone);  
        
        %% Detect fields
        Field_Info = detectFields(toneMap);
        if isempty(Field_Info)
            toneField = 0;
        else 
            toneField = 1;
        end                
        
        % If its a trajectory cell, calculate its PSTH in response to each
        % trial type, correct, incorrect, return
        [~,idxMax] = max(toneMap);
        if cellType == 1 && (toneField == 1) && (toneCorr > 0.1) %&& idxMax>40
            sumLicks = [sumLicks; size(tsLick2,1)];
            if size(tsLick2,1)>140
                disp('wut')
            end
            for tt = 1:3

                if tt ==1 
                    st = behavTrials.timestamps(behavTrials.linTrial==0,2);
                elseif tt == 2 %Middle ports, linear trials
                    intPer = InIntervals(tsLick2(:,1),behavTrials.timestamps(behavTrials.linTrial==1,:));
                    st = tsLick2(intPer,1);   
                    if isempty(st)
                        sumPoke1 = [sumPoke1;0];
                    else
                        sumPoke1 = [sumPoke1;length(st)];
                    end                        
                elseif tt == 3 %Middle ports,return run
                    idxLin = behavTrials.linTrial==0;
                    tsFwd = behavTrials.timestamps(idxLin,1);
                    tsRet = behavTrials.timestamps(idxLin,2);
                    ints = [tsRet(1:(end-1)) tsFwd(2:end)];
                    intPer = InIntervals(tsLick2(:,1),ints);  
                    st = tsLick2(intPer,1);    
                    if isempty(st)
                        sumPoke2 = [sumPoke2;0];
                    else
                        sumPoke2 = [sumPoke2;length(st)];
                    end                       
                end
                    
                if ~isempty(st)
                    [stccg, tPSTH] = CCG({spikes.times{kk} st},[],'binSize',0.01,'duration',1,'norm','rate');
                    psthReward{tt} = [psthReward{tt}; stccg(:,2,1)'];               
                else
                    fillArr(1,1:21) = nan;
                    psthReward{tt} = [psthReward{tt}; fillArr];               
                end
            end

        end
      
    end
end


figure
set(gcf,'Renderer','painters')
set(gcf,'Color','w')

YlGnBu=cbrewer('seq', 'YlGnBu', 11);
colormap(YlGnBu)

idxT = tPSTH<0 & tPSTH>=-0.2;
avgRate = nanmean(psthReward{1}(:,idxT),2);
for tt = 1:3
    newpsth{tt} = psthReward{tt}(avgRate>1,:);
end
idxT = tPSTH<0 & tPSTH>=-0.5;
[~,idxMax2] = max(newpsth{1}(:,idxT),[],2);
[~,idxMax] = sort(idxMax2);
for tt = 1:3
    subplot(2,3,tt)
    %temp = zscore(newpsth{tt},[],2);
    %h = imagesc(tPSTH, 1:size(newpsth{tt},1),temp(idxMax,:));
    %set(h, 'AlphaData', ~isnan(temp))
    imagesc(tPSTH, 1:size(newpsth{tt},1),newpsth{tt});
    caxis([0 10])
    %ylim([125 350])
    colorbar
        
    subplot(2,3,tt+3)
    col = [0.5 0.5 0.5];
    meanpsth = nanmean(newpsth{tt},1);
    stdpsth = nanstd(newpsth{tt},1)./sqrt(size(newpsth{tt},1));         
    hold on
    fill([tPSTH; flipud(tPSTH)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],col,'linestyle','none','FaceAlpha',0.5);                    
    hold on
    hi = line(tPSTH,meanpsth,'LineWidth',1.5,'Color',col);
    line([0 0],[1 11],'Color','r')
    ylim([1 11])    
end

saveas(gcf,strcat(expPath,'Compiled\portPokePSTHs.png'));
saveas(gcf,strcat(expPath,'Compiled\portPokePSTHs.eps'),'epsc');
saveas(gcf,strcat(expPath,'Compiled\portPokePSTHs.fig'));

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