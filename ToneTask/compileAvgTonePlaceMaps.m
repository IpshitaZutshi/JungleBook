function compileAvgTonePlaceMaps

% sess= {'IZ41\Final\IZ41_220624_sess5','IZ41\Final\IZ41_220629_sess7','IZ41\Final\IZ41_220701_sess9',...
%     'IZ41\Final\IZ41_220704_sess11','IZ41\Final\IZ41_220708_sess13','IZ41\Final\IZ41_220714_sess14',...

%    'IZ46\IZ46_230406_sess9','IZ46\IZ46_230420_sess19','IZ46\IZ46_230424_sess21',...
sess= {'IZ45\IZ45_230414_sess15',...
    'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
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


AllcellType = [];
AllsessType = [];
AlllinCorr = [];
AllspatialCorr = [];
AlltoneCorr = [];
AlllinMapInit = [];
AlllinMapEnd = [];
AllspaceMap = [];

AllretMapLinInit = [];
AllretMapCorr = [];
AllretMapIncorr = [];

AlltoneMap = [];
AlllinInfo = [];
AllspatialInfo = [];
AlltoneInfo = [];
AlllinpeakRate = [];
AllspacepeakRate = [];
AlltonepeakRate = [];
AlllinField = [];
AllspaceField = [];
AlltoneField = [];
AllretFieldlin = [];
AllretFieldCorrect = [];
AllsessID = [];
AllcellID = [];

AlltoneNoToneCorr = [];
AlllinlinEndCorr = [];
AlltonelinEndCorr = [];
AllretlinCorrCorr = [];
AllretCorrInCorrCorr = [];

for ii = 1:length(sess)
    %% Load files
    cd(strcat(expPath,sess{ii}))    
    file = dir(['*.rateMapsAvg.cellinfo.mat']);
    load(file(1).name);
    file = dir(['*.Tracking.Behavior.mat']);
    load(file(1).name);
    file = dir(['*TrialBehavior.Behavior.mat']);
    load(file(1).name);
    [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts', true);
    file = dir(['*cell_metrics.cellinfo.mat']);
    load(file(1).name);
     
    cellType(1:length(cell_metrics.UID),1) = nan;
    sessType(1:length(cell_metrics.UID),1) = nan;
    
    linCorr(1:length(cell_metrics.UID),1) = nan;
    spatialCorr(1:length(cell_metrics.UID),1) = nan;
    toneCorr(1:length(cell_metrics.UID),1) = nan;
    
    linMapInit(1:length(cell_metrics.UID),1:50) = nan;
    linMapEnd(1:length(cell_metrics.UID),1:50) = nan;
    spaceMap(1:length(cell_metrics.UID),1:50) = nan;
    toneMap(1:length(cell_metrics.UID),1:50) = nan;   
    
    retMapLinInit(1:length(cell_metrics.UID),1:50) = nan;
    retMapCorr(1:length(cell_metrics.UID),1:50) = nan;
    retMapIncorr(1:length(cell_metrics.UID),1:50) = nan;
    
    linInfo(1:length(cell_metrics.UID),1) = nan;
    spatialInfo(1:length(cell_metrics.UID),1) = nan;
    toneInfo(1:length(cell_metrics.UID),1) = nan;    
    
    linpeakRate(1:length(cell_metrics.UID),1) = nan;
    spacepeakRate(1:length(cell_metrics.UID),1) = nan;
    tonepeakRate(1:length(cell_metrics.UID),1) = nan;        

    linField(1:length(cell_metrics.UID),1) = nan;
    spaceField(1:length(cell_metrics.UID),1) = nan;
    toneField(1:length(cell_metrics.UID),1) = nan; 
    retFieldlin(1:length(cell_metrics.UID),1) = nan;
    retFieldCorrect(1:length(cell_metrics.UID),1) = nan;

    sessID(1:length(cell_metrics.UID),1) = nan;
    cellID(1:length(cell_metrics.UID),1) = nan;  
    
    toneNoToneCorr(1:length(cell_metrics.UID),1) = nan;
    linlinEndCorr(1:length(cell_metrics.UID),1) = nan;
    tonelinEndCorr(1:length(cell_metrics.UID),1) = nan;
    retlinCorrCorr(1:length(cell_metrics.UID),1) = nan;
    retCorrInCorrCorr(1:length(cell_metrics.UID),1) = nan;
    
    for kk=1:length(cell_metrics.UID)
        sessID(kk,1) = ii;
        cellID(kk,1) = kk;
    %% Extract relevant rate maps
        linMapInit(kk,:) = firingMaps.forward.rateMaps{kk}{1};
        retMapLinInit(kk,:) = firingMaps.reverse.rateMaps{kk}{1};
        retMapCorr(kk,:) = firingMaps.reverse.rateMaps{kk}{2};
        retMapIncorr(kk,:) = firingMaps.reverse.rateMaps{kk}{3};
        
        if length(firingMaps.forward.rateMaps{kk})==28
          linMapEnd(kk,:) = firingMaps.forward.rateMaps{kk}{28};
        end
    %% Extract pyramidal cells and calculate rate map correlations and information        
        if strcmp(cell_metrics.putativeCellType{kk},'Pyramidal Cell') == 1 
            cellType(kk,1) = 1;
        else
            cellType(kk,1) = 0;
        end
        if strcmp(sess{ii}(1:4),'IZ41')==1 || strcmp(sess{ii}(1:4),'IZ36')==1 || strcmp(sess{ii}(1:4),'IZ45')==1 || strcmp(sess{ii}(1:4),'IZ46')==1
            sessType(kk,1) = 0;
        else
            sessType(kk,1) = 1;
        end
        % Average linear track correlation
        corrLin = corrcoef(firingMaps.forward.rateMaps{kk}{26}',firingMaps.forward.rateMaps{kk}{27}','rows','complete');
        linCorr(kk,1) = corrLin(1,2);

        % Average linear track correlation (beginning and end)
        if length(firingMaps.forward.rateMaps{kk})==28
            corrLin = corrcoef(firingMaps.forward.rateMaps{kk}{1}',firingMaps.forward.rateMaps{kk}{28}','rows','complete');
            linlinEndCorr(kk,1) = corrLin(1,2);
        else
            linlinEndCorr(kk,1) = nan;
        end

        % Average linear track correlation (no tone and tone)
        corrLin = corrcoef(firingMaps.forward.rateMaps{kk}{1}',firingMaps.forward.rateMaps{kk}{7}','rows','complete');
        toneNoToneCorr(kk,1) = corrLin(1,2);

        % Average linear track correlation (tone and no tone end)
        if length(firingMaps.forward.rateMaps{kk})==28
            corrLin = corrcoef(firingMaps.forward.rateMaps{kk}{7}',firingMaps.forward.rateMaps{kk}{28}','rows','complete');
            tonelinEndCorr(kk,1) = corrLin(1,2);
        else
            tonelinEndCorr(kk,1) = nan;
        end                

        % Average return track correlation (no tone and correct tone)
        corrLin = corrcoef(firingMaps.reverse.rateMaps{kk}{1}',firingMaps.reverse.rateMaps{kk}{2}','rows','complete');
        retlinCorrCorr(kk,1) = corrLin(1,2);
        
        % Average return track correlation (incorrect vs correct tone)
        corrLin = corrcoef(firingMaps.reverse.rateMaps{kk}{2}',firingMaps.reverse.rateMaps{kk}{3}','rows','complete');
        retCorrInCorrCorr(kk,1) = corrLin(1,2);
        
        % Average spatial and tone correlation
        dataMat = [];
        dataMatTone = [];
        info = [];
        infoTone = [];
        for jj = 2:7           
            dataMat = [dataMat;firingMaps.forward.rateMaps{kk}{jj}];
            info = [info;firingMaps.forward.information{kk}{jj}];
            a = fillmissing(firingMaps.tone.rateMaps{kk}{jj},'linear');
            dataMatTone = [dataMatTone;a];
            infoTone = [infoTone;firingMaps.tone.information{kk}{jj}];
        end 
        
        spaceMap(kk,:) = nanmean(dataMat,1);        
        toneMap(kk,:) = nanmean(dataMatTone,1);        
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
        spatialCorr(kk,1) = nanmean(corrSpace);
        toneCorr(kk,1) = nanmean(corrTone);     
        
        linInfo(kk,1) = firingMaps.forward.information{kk}{1};
        spatialInfo(kk,1) = nanmean(info);
        toneInfo(kk,1) = nanmean(infoTone);   
        
        linpeakRate(kk,1) = nanmax(linMapInit(kk,:));
        spacepeakRate(kk,1) = nanmax(spaceMap(kk,:));
        tonepeakRate(kk,1) = nanmax(toneMap(kk,:));
        
        %% Detect fields
        Field_Info = detectFields(linMapInit(kk,:));
        if isempty(Field_Info)
            linField(kk,1) = 0;
        else 
            linField(kk,1) = 1;
        end

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
        
        Field_Info = detectFields(retMapLinInit(kk,:));
        if isempty(Field_Info)
            retFieldlin(kk,1) = 0;
        else 
            retFieldlin(kk,1) = 1;
        end  
        
        Field_Info = detectFields(retMapCorr(kk,:));
        if isempty(Field_Info)
            retFieldCorrect(kk,1) = 0;
        else 
            retFieldCorrect(kk,1) = 1;
        end         

    end
    AllcellType = [AllcellType;cellType];
    AllsessType = [AllsessType;sessType];
    AlllinCorr = [AlllinCorr;linCorr];
    AllspatialCorr = [AllspatialCorr;spatialCorr];
    AlltoneCorr = [AlltoneCorr;toneCorr];
    AlllinMapInit = [AlllinMapInit;linMapInit];
    AlllinMapEnd = [AlllinMapEnd;linMapEnd];
    AllspaceMap = [AllspaceMap;spaceMap];
    AlltoneMap = [AlltoneMap;toneMap];
    AlllinInfo = [AlllinInfo;linInfo];
    AllspatialInfo = [AllspatialInfo;spatialInfo];
    AlltoneInfo = [AlltoneInfo;toneInfo];
    AlllinpeakRate = [AlllinpeakRate;linpeakRate];
    AllspacepeakRate = [AllspacepeakRate;spacepeakRate];
    AlltonepeakRate = [AlltonepeakRate;tonepeakRate];
    AlllinField = [AlllinField;linField];
    AllspaceField = [AllspaceField;spaceField];
    AlltoneField = [AlltoneField;toneField];    
    AllsessID = [AllsessID;sessID];
    AllcellID = [AllcellID;cellID];    
    AlltoneNoToneCorr = [AlltoneNoToneCorr;toneNoToneCorr];
    AlllinlinEndCorr = [AlllinlinEndCorr;linlinEndCorr];
    AlltonelinEndCorr = [AlltonelinEndCorr; tonelinEndCorr];
    
    AllretMapLinInit = [AllretMapLinInit;retMapLinInit];
    AllretMapCorr = [AllretMapCorr;retMapCorr];
    AllretMapIncorr = [AllretMapIncorr;retMapIncorr];
    AllretFieldlin = [AllretFieldlin;retFieldlin];
    AllretFieldCorrect = [AllretFieldCorrect;retFieldCorrect];
    AllretlinCorrCorr = [AllretlinCorrCorr;retlinCorrCorr];
    AllretCorrInCorrCorr = [AllretCorrInCorrCorr;retCorrInCorrCorr];
       
    clear cellType sessType linCorr spatialCorr toneCorr linMapInit linMapEnd spaceMap toneMap toneNoToneCorr linlinEndCorr linInfo ...
        spatialInfo toneInfo linpeakRate tonelinEndCorr spacepeakRate tonepeakRate linField spaceField toneField sessID cellID ...
        retMapLinInit retMapCorr retMapIncorr retFieldlin retFieldCorrect retlinCorrCorr retCorrInCorrCorr
end

YlGnBu=cbrewer('seq', 'YlGnBu', 11);
fractTune = [];
linPos = linspace(1,122,50);
linTone = linspace(2000,22000,50);
tag = {'Control','ACgN'};

for ss = 1:2
    idxSess = AllsessType==(ss-1) & AllcellType == 1;
    idxActive = AllsessType==(ss-1) & AllcellType == 1 & nanmean(AllspaceMap,2)>0.2;
    idxMaps{1} = idxSess & AlllinField & AlllinCorr>0.1;
    idxMaps{2} = idxSess & AllspaceField & AllspatialCorr>0.1;
    idxMaps{3} = idxSess & AlltoneField & AlltoneCorr>0.1;    
    figure
    set(gcf,'Color','w')
    set(gcf,'renderer','painters')    
    set(gcf,'Position',[2384 318 617 547])
    for ii = 1:3
        fractTune(ss,ii) = sum(idxMaps{ii})/sum(idxActive);
        selectedlinMap = AlllinMapInit(idxMaps{ii},:);
        selectedlinMapEnd = AlllinMapEnd(idxMaps{ii},:);
        selectedtoneMap = AlltoneMap(idxMaps{ii},:);
        selectedspaceMap = AllspaceMap(idxMaps{ii},:);
        selectedsessID = AllsessID(idxMaps{ii});
        selectedcellID = AllcellID(idxMaps{ii});
        
        [maxLin,idxLin] = max(selectedlinMap,[],2);
        [maxLinEnd,~] = max(selectedlinMapEnd,[],2);
        [maxSpace,idxSpace] = max(selectedspaceMap,[],2);        
        [maxTone,idxTone] = max(selectedtoneMap,[],2);
        
        %Z score normalize
        normlinMap = (selectedlinMap-nanmean(selectedlinMap,2))./nanstd(selectedlinMap,[],2);
        normlinMapEnd = (selectedlinMapEnd-nanmean(selectedlinMapEnd,2))./nanstd(selectedlinMapEnd,[],2);
        normtoneMap = (selectedtoneMap-nanmean(selectedtoneMap,2))./nanstd(selectedtoneMap,[],2);
        normspaceMap = (selectedspaceMap-nanmean(selectedspaceMap,2))./nanstd(selectedspaceMap,[],2);
        
        if ii ==1            
            [~,sortidx] = sort(idxLin,'ascend');
        elseif ii == 2            
            [~,sortidx] = sort(idxSpace,'ascend');
        elseif ii == 3            
            [~,sortidx] = sort(idxTone,'ascend');
        end
 
        colormap(YlGnBu)
        subplot(3,4,4*(ii-1)+1)
        imagesc(linPos, 1:length(sortidx),normlinMap(sortidx,:))
        ylabel('Cell ID')
        xlabel('Position')
        caxis([-1 4])
        title('Linear early')
        
        subplot(3,4,4*(ii-1)+2)
        imagesc(linPos, 1:length(sortidx),normlinMapEnd(sortidx,:))
        caxis([-1 4])
        ylabel('Cell ID')
        xlabel('Position')
        title('Linear late')
        
        subplot(3,4,4*(ii-1)+3)
        imagesc(linPos, 1:length(sortidx),normspaceMap(sortidx,:))
        caxis([-1 4])
        ylabel('Cell ID')
        xlabel('Position')
        title('Task spatial')
        
        subplot(3,4,4*(ii-1)+4)
        imagesc(linTone, 1:length(sortidx),normtoneMap(sortidx,:))      
        caxis([-1 4])
        ylabel('Cell ID')
        xlabel('Frequency')  
        title('Task auditory')
    end    
    saveas(gcf,strcat(expPath,'Compiled\populationMaps',tag{ss},'.png'));
    saveas(gcf,strcat(expPath,'Compiled\populationMaps',tag{ss},'eps'),'epsc');
    saveas(gcf,strcat(expPath,'Compiled\populationMaps',tag{ss},'fig'));
end

figure
set(gcf,'Color','w')
set(gcf,'renderer','painters')
set(gcf,'Position',[2384 318 617 547])
subplot(2,2,1)
idxSess = AllsessType==0 & AllcellType == 1 & (AlltoneField==1) & (AlltoneCorr >0.1);
scatter(AllspatialCorr(idxSess),AlltoneCorr(idxSess),10,[187/243 86/243 149/243],'filled')
hold on
idxSess = AllsessType==0 & AllcellType == 1 & (AllspaceField==1) & (AllspatialCorr >0.1);
scatter(AllspatialCorr(idxSess),AlltoneCorr(idxSess),10,[0.5 0.5 0.5],'filled')
ylim([-0.2 1])
xlim([-0.2 1])
refline(1)
title('Control')
xlabel('Space correlation')
ylabel('Tone correlation')

subplot(2,2,2)
idxSess = AllsessType==1 & AllcellType == 1 & (AllspaceField==1);% & (AllspatialCorr >0.1);
scatter(AllspatialCorr(idxSess),AlltoneCorr(idxSess),10,[0.5 0.5 0.5],'filled')
hold on
idxSess = AllsessType==1 & AllcellType == 1 & (AlltoneField==1);% & (AlltoneCorr >0.1);
scatter(AllspatialCorr(idxSess),AlltoneCorr(idxSess),10,[187/243 86/243 149/243],'filled')

ylim([-0.2 1])
xlim([-0.2 1])
refline(1)
title('ACgN')
xlabel('Space correlation')
ylabel('Tone correlation')

subplot(2,2,3)
b = bar(fractTune','FaceColor','flat');
b(1).CData = [0.5 0.5 0.5];
b(2).CData = [187/243 86/243 149/243];
ylabel('Fraction of place cells')
xticklabels({'Space','Tone'})

subplot(2,2,4)
idxSess = AllsessType==0 & AllcellType == 1 & ((AllspaceField==1) | (AlllinField==1));
data{1} = AlllinCorr(idxSess);
%data{2} = AlllinlinEndCorr(idxSess);
data{2} = AlltoneNoToneCorr(idxSess);

idxSess = AllsessType==1 & AllcellType == 1 & ((AllretFieldlin | AllretFieldCorrect));
data{3} = AllretlinCorrCorr(idxSess);

idxSess = AllsessType==1 & AllcellType == 1 & ((AllspaceField==1) | (AlllinField==1));
data{4} = AlllinCorr(idxSess);
data{5} = AlltoneNoToneCorr(idxSess);
data{6} = AlltonelinEndCorr(idxSess);
data{7} = AlllinlinEndCorr(idxSess);

idxSess = AllsessType==1 & AllcellType == 1 & ((AllretFieldlin | AllretFieldCorrect));
data{8} = AllretlinCorrCorr(idxSess);
data{9} = AllretCorrInCorrCorr(idxSess);

col = [52/243 52/243 52/243;...
    56/243 61/243 150/243;...
    94/243 60/243 108/243;...
    52/243 52/243 52/243;...
    56/243 61/243 150/243;...
    80/243 91/243 166/243;...
    180/243 180/243 180/243;...
    94/243 60/243 108/243;...
    133/243 128/243 177/243];
    
    
stats = groupStats(data,[],'inAxis',true,'color',col);
ylabel('Spatial correlation')

saveas(gcf,strcat(expPath,'Compiled\Fraction_Correlation.png'));
saveas(gcf,strcat(expPath,'Compiled\Fraction_Correlation.eps'),'epsc');
saveas(gcf,strcat(expPath,'Compiled\Fraction_Correlation.fig'));
save(strcat(expPath,'Compiled\Fraction_Correlation.mat'),'stats'); 

figure
set(gcf,'Color','w')
set(gcf,'renderer','painters')
set(gcf,'Position',[2384 318 617 547])

for ss = 1:2
    idxSess = AllsessType==(ss-1) & AllcellType == 1;
    idxActive = AllsessType==(ss-1) & AllcellType == 1 & nanmean(AllretMapLinInit,2)>0.2;
    idxMaps = idxSess &  (AllretFieldlin | AllretFieldCorrect);
    if ss == 1
        iRange = [1 2];
    else
        iRange = [1 2 3];
    end
    for ii = iRange       
        selectedlinMap = AllretMapLinInit(idxMaps,:);
        selectedcorrMap = AllretMapCorr(idxMaps,:);
        selectedincorrMap = AllretMapIncorr(idxMaps,:);
        selectedsessID = AllsessID(idxMaps);
        selectedcellID = AllcellID(idxMaps);
        
        [maxLin,idxLin] = max(selectedlinMap,[],2);
        [maxCorr,idxCorr] = max(selectedcorrMap,[],2);        
        [maxInCorr,idxInCorr] = max(selectedincorrMap,[],2);
        
        %Z score normalize
        normlinMap = (selectedlinMap-nanmean(selectedlinMap,2))./nanstd(selectedlinMap,[],2);
        normcorrMap = (selectedcorrMap-nanmean(selectedcorrMap,2))./nanstd(selectedcorrMap,[],2);
        normincorrMap = (selectedincorrMap-nanmean(selectedincorrMap,2))./nanstd(selectedincorrMap,[],2);
        
        if ii ==1            
            [~,sortidx] = sort(idxLin,'ascend');
        elseif ii == 2            
            [~,sortidx] = sort(idxCorr,'ascend');
        elseif ii == 3            
            [~,sortidx] = sort(idxInCorr,'ascend');
        end
 
        colormap(YlGnBu)
        subplot(3,7,3*(ss-1)+7*(ii-1)+1)
        imagesc(linPos, 1:length(sortidx),normlinMap(sortidx,:))
        ylabel('Cell ID')
        xlabel('Position')
        caxis([-1 4])
        title('Return no tone first half')
        
        subplot(3,7,3*(ss-1)+7*(ii-1)+2)
        imagesc(linPos, 1:length(sortidx),normcorrMap(sortidx,:))
        caxis([-1 4])
        ylabel('Cell ID')
        xlabel('Position')
        title('Return tone correct')
        
        subplot(3,7,3 *(ss-1)+7*(ii-1)+3)
        imagesc(linPos, 1:length(sortidx),normincorrMap(sortidx,:))
        caxis([-1 4])
        ylabel('Cell ID')
        xlabel('Position')
        title('Return tone incorrect')
        
    end    
end

subplot(3,7,14)
idxSess = AllsessType==0 & AllcellType == 1 & ((AllretFieldlin | AllretFieldCorrect));
data = [];

data{1} = AllretlinCorrCorr(idxSess);

idxSess = AllsessType==1 & AllcellType == 1 & ((AllretFieldlin | AllretFieldCorrect));
data{2} = AllretlinCorrCorr(idxSess);
data{3} = AllretCorrInCorrCorr(idxSess);

stats = groupStats(data,[],'inAxis',true,'color',col);
ylabel('Spatial correlation')

saveas(gcf,strcat(expPath,'Compiled\returnpopulationMaps.png'));
saveas(gcf,strcat(expPath,'Compiled\returnpopulationMaps.eps'),'epsc');
saveas(gcf,strcat(expPath,'Compiled\returnpopulationMaps.fig'));

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