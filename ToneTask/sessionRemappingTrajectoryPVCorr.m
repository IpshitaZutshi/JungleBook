function sessionRemappingTrajectoryPVCorr

sess = {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
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
plotfig = 0;

YlGnBu=cbrewer('seq', 'YlGnBu', 11);

AllCellsCorr = [];
AllCellsSimPre = [];
AllCellsSim = [];

AllCellsRateMapLin1 = [];
AllCellsRateMapLin2 = [];
AllCellsRateMapTone = [];

for ii = 1:length(sess)
    
    %% Load files
    cd(strcat(expPath,sess{ii}))    
    file = dir(['*.rateMapsAvg.cellinfo.mat']);
    load(file(1).name);
    file = dir(['*.rateMapsTrial.cellinfo.mat']);
    firingMapsTrial = load(file(1).name);    
    file = dir(['*.Tracking.Behavior.mat']);
    load(file(1).name);
    file = dir(['*TrialBehavior.Behavior.mat']);
    load(file(1).name);
    [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts', true);
    file = dir(['*cell_metrics.cellinfo.mat']);
    load(file(1).name);    
    
    %Skip sessions if less than 5 post trials
    linIdx = find(behavTrials.linTrial==1);
    jumpLin = find(diff(linIdx)>1);  

    if isempty(jumpLin)
        lin_end = [];
    else
        lin_end = linIdx(jumpLin+1:end);
    end

    if isempty(lin_end)||length(lin_end)<=5
        continue
    end
    
    % Cluster trajectories and extract similar trials
    [trialNum, trialLabel, trajectLabel] = clusterTrajectories;
    
    % Check if the post-tone trials have more similarity to the no tone or
    % tone    
    typePre = mode(trajectLabel(trialLabel==1));
    typeTone = mode(trajectLabel(trialLabel==2));
    typePost = mode(trajectLabel(trialLabel==3));

    if typePost == typeTone
        sim = 1;
    else
        sim = 0 ;
    end
    
    if typePre == typeTone
        simPre  = 1;
    else
        simPre  = 0;
    end
    
    for kk=1:length(cell_metrics.UID)
        for zz = 1:6
            Maps{zz} = [];
        end
        %% First, only extract cells with a field in either of the three
        %trial types
        linMapInit = firingMaps.forward.rateMaps{kk}{1};
        linMapEnd = firingMaps.forward.rateMaps{kk}{28};
        toneMap = firingMaps.forward.rateMaps{kk}{7};
        
        %% Detect fields
        Field1 = detectFields(linMapInit);
        Field2 = detectFields(toneMap);
        Field3 = detectFields(linMapEnd);
        
        if isempty(Field1) && isempty(Field2) && isempty(Field3)
            continue
        end
        

        corrMap = [];
        simil  = [];
        similPre = [];
        corr1 = corrcoef(linMapInit',linMapEnd','rows','complete');
        corrMap(1) = corr1(1,2); %
        corr1 = corrcoef(linMapInit',toneMap','rows','complete');
        corrMap(2) = corr1(1,2); %        
        corr1 = corrcoef(toneMap',linMapEnd','rows','complete');
        corrMap(3) = corr1(1,2); %        
            
        % Skip sessions that have less than 5 post-tone trial linear trials of
        % a particular type
        if sum(trialLabel==3 & trajectLabel==typePost)<5
            AllCellsCorr = [AllCellsCorr; corrMap nan nan nan];
            AllCellsSim = [AllCellsSim;sim];
            AllCellsSimPre = [AllCellsSimPre;simPre];
            AllCellsRateMapLin1 = [AllCellsRateMapLin1;linMapInit];
            AllCellsRateMapLin2 = [AllCellsRateMapLin2;linMapEnd];
            AllCellsRateMapTone = [AllCellsRateMapTone;toneMap];
            continue;
        end
            
        %% Now extract firing map of each category according to trajectories
        
        for zz = 1:3
            switch zz
                case 1
                    idx = trialLabel == 1 & trajectLabel == typePre;
                case 2
                    idx = trialLabel == 3 & trajectLabel == typePost;                     
                case 3
                    idx = trialLabel == 2 & trajectLabel == typeTone;
                                 
            end
            
            trialIdx  = trialNum(idx);
            for tt = 1:length(trialIdx)
                Maps{zz} = [Maps{zz};firingMapsTrial.firingMaps.forward.rateMaps{kk}{trialIdx(tt)}];
            end
            
            if isempty(Maps{zz})
                avgMaps(zz,1:50) = nan;
            else
                avgMaps(zz,:) = nanmean(Maps{zz},1);
            end

        end   
           
        AllCellsRateMapLin1 = [AllCellsRateMapLin1;avgMaps(1,:)];
        AllCellsRateMapLin2 = [AllCellsRateMapLin2;avgMaps(2,:)];
        AllCellsRateMapTone = [AllCellsRateMapTone;avgMaps(3,:)];
                    
        corrMatrix = corrcoef(avgMaps','Rows','pairwise');
        C = tril(corrMatrix,-1);
        C = C(isnan(C) | C~=0);
        
        corrMap = [corrMap C'];
        simil = [simil sim];
        similPre = [similPre simPre];
        
        if plotfig
            saveLoc = strcat(expPath,sess{ii},'\TrajectMaps');
            if ~isfolder('TrajectMaps')
                mkdir('TrajectMaps')
            end  
            
            figure
            subplot(6,2,1)
            imagesc(linMapInit)
            colormap(YlGnBu)
            title(strcat('Max:',num2str(max(linMapInit)),' Hz'))
            
            subplot(6,2,5)
            imagesc(toneMap)
            title(strcat('Max:',num2str(max(toneMap)),' Hz'))
            
            subplot(6,2,9)
            imagesc(linMapEnd)
            title(strcat('Max:',num2str(max(linMapEnd)),' Hz'))
            
            for zz = 1:6
                subplot(6,2,2*(zz-1)+2)
                imagesc(Maps{zz})
                title(strcat('Max:',num2str(max(max(Maps{zz}))),' Hz'))
            end
            
            subplot(6,2,11)
            imagesc(avgMaps)
            title('Average trajectory maps')
            
            saveas(gcf,[saveLoc,filesep ,'Cell_', num2str(kk),'.png'],'png');
            saveas(gcf,[saveLoc,filesep ,'Cell_', num2str(kk),'.fig'],'fig');
            saveas(gcf,[saveLoc,filesep ,'Cell_', num2str(kk),'.eps'],'epsc');
            close all
        end
            
        AllCellsCorr = [AllCellsCorr;corrMap];
        AllCellsSim = [AllCellsSim;simil];
        AllCellsSimPre = [AllCellsSimPre;similPre];
    end
end

for ii = 1:3
    data1{ii} = AllCellsCorr(:,ii); 
end
dat = AllCellsCorr(:,2);
data4{1} = dat(AllCellsSimPre==0);
data4{2} = dat(AllCellsSimPre==1);

for ii = 4:6
    dat = AllCellsCorr(:,ii);
    data2{ii-3} = dat(AllCellsSim == 1 & AllCellsSimPre==0);
    data3{ii-3} = dat(AllCellsSim == 0 & AllCellsSimPre==0);
end

figure
set(gcf,'Renderer','painters')
subplot(3,4,1)
stats{1} = groupStats(data1,{},'inAxis',true,'repeatedMeasures',true);
subplot(3,4,2)
stats{2} = groupStats(data2,{},'inAxis',true,'repeatedMeasures',true);
subplot(3,4,3)
stats{3} = groupStats(data3,{},'inAxis',true,'repeatedMeasures',true);
subplot(3,4,4)
stats{4} = groupStats(data4,{},'inAxis',true);


YlGnBu=cbrewer('seq', 'YlGnBu', 11);
for kk = 1:2
    if kk == 1
        idx = (AllCellsSim == 1 & AllCellsSimPre==0);
    else
        idx = (AllCellsSim == 0 & AllCellsSimPre==0);
    end
    selectedlinMap = AllCellsRateMapLin1(idx,:);
    selectedlinMapEnd = AllCellsRateMapLin2(idx,:);
    selectedtoneMap = AllCellsRateMapTone(idx,:);

    %Deal with nans
    selectedlinMap = selectedlinMap(:,1:48);
    selectedlinMapEnd = selectedlinMapEnd(:,1:48);
    selectedtoneMap = selectedtoneMap(:,1:48);
    
    avglinMap = nanmean(selectedlinMap,2);
    avglinMapEnd = nanmean(selectedlinMapEnd,2);
    avgtoneMap = nanmean(selectedtoneMap,2);
    
    idxtoUseforcorr = ~isnan(avglinMap)&~isnan(avglinMapEnd)&~isnan(avgtoneMap) & avglinMap>0 &avglinMapEnd>0 &avgtoneMap>0;
    selectedlinMap = selectedlinMap(idxtoUseforcorr,:);
    selectedlinMapEnd = selectedlinMapEnd(idxtoUseforcorr,:);
    selectedtoneMap = selectedtoneMap(idxtoUseforcorr,:);   
    
    [maxLin,~] = max(selectedlinMap,[],2);
    [maxLinEnd,idxLin] = max(selectedlinMapEnd,[],2);            
    [maxTone,idxTone] = max(selectedtoneMap,[],2);
    %Z score normalize
    normlinMap = (selectedlinMap);%-nanmean(selectedlinMap,2))./nanstd(selectedlinMap,[],2);
    normlinMapEnd = (selectedlinMapEnd);%-nanmean(selectedlinMapEnd,2))./nanstd(selectedlinMapEnd,[],2);
    normtoneMap = (selectedtoneMap);%-nanmean(selectedtoneMap,2))./nanstd(selectedtoneMap,[],2);   
    
    [~,sortidx] = sort(idxLin,'ascend');

    linPos = linspace(1,122,50);
    
    colormap(YlGnBu)
    
    subplot(3,4,4*(kk-1)+4+1)
    imagesc(linPos, 1:length(sortidx),normlinMap(sortidx,:))
    ylabel('Cell ID')
    xlabel('Position')
    caxis([0 20])

    subplot(3,4,4*(kk-1)+4+2)
    imagesc(linPos, 1:length(sortidx),normtoneMap(sortidx,:))
    ylabel('Cell ID')
    xlabel('Position')
    caxis([-1 4])

    subplot(3,4,4*(kk-1)+4+3)
    imagesc(linPos, 1:length(sortidx),normlinMapEnd(sortidx,:))
    ylabel('Cell ID')
    xlabel('Position')
    caxis([-1 4])
   
    subplot(3,4,4*(kk-1)+4+4)
    a = corr(selectedlinMap,selectedlinMapEnd);
    Rvec= diag(a);
    plot(linPos(1:48),Rvec,'Color',[0.5 0.5 0.5],'LineWidth',1.5)
    
    hold on
    a = corr(selectedtoneMap,selectedlinMapEnd);
    Rvec= diag(a);
    plot(linPos(1:48),Rvec,'Color','b','LineWidth',1.5)
    
    a = corr(selectedtoneMap,selectedlinMap);
    Rvec= diag(a);
    plot(linPos(1:48),Rvec,'Color','r','LineWidth',1.5)
    ylim([0 1])
    
end



saveas(gcf,strcat(expPath,'Compiled\TrajectoryRemappingMaps.png'));
saveas(gcf,strcat(expPath,'Compiled\TrajectoryRemappingMaps.eps'),'epsc');
saveas(gcf,strcat(expPath,'Compiled\TrajectoryRemappingMaps.fig'));
save(strcat(expPath,'Compiled\TrajectoryRemappingMaps.mat'),'stats')   

end

function Field_Info = detectFields(SmoothedFiringRate)
    minFieldSize = 10;
    maxFieldSize = 50;
    % Pad on each end with zeros for edge effects
    SmoothedFiringRate = [0 0 SmoothedFiringRate 0 0];
    [peakValues, peakLocations] = findpeaks(SmoothedFiringRate, 'minpeakheight',5, 'minpeakdistance', 10);
    Field_Info = [];
    for j = 1:length(peakLocations)
        FieldPeak = peakLocations(j);
        % FieldPeak must be 8 Hz or more
        if peakValues(j) < 10, continue, end
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