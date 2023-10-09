function Summary = compileAvgTonePlaceMaps(varargin)
warning off

sess= {'IZ41\Final\IZ41_220624_sess5','IZ41\Final\IZ41_220701_sess9',...
    'IZ41\Final\IZ41_220704_sess11','IZ41\Final\IZ41_220708_sess13',...
    'IZ41\Final\IZ41_220629_sess7','IZ41\Final\IZ41_220714_sess14',...
    'IZ45\Final\IZ45_230410_sess11','IZ45\Final\IZ45_230420_sess19',...
    'IZ45\Final\IZ45_230414_sess15','IZ45\Final\IZ45_230417_sess16','IZ45\Final\IZ45_230424_sess21',... 
    'IZ46\Final\IZ46_230406_sess9','IZ46\Final\IZ46_230407_sess10',...
    'IZ46\Final\IZ46_230410_sess11','IZ46\Final\IZ46_230413_sess14','IZ46\Final\IZ46_230420_sess19'...    16
    'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
    'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
    'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   23
    'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...
    'IZ40\Final\IZ40_220708_sess17','IZ40\Final\IZ40_220714_sess18',...27
    'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
    'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
    'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
    'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    36
    'IZ44\Final\IZ44_220827_sess4', 'IZ44\Final\IZ44_220828_sess5',...
    'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...
    'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
    'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',... 45
    'IZ47\Final\IZ47_230626_sess15','IZ47\Final\IZ47_230707_sess24',...
    'IZ47\Final\IZ47_230710_sess25','IZ47\Final\IZ47_230712_sess27',...49
    'IZ48\Final\IZ48_230628_sess17','IZ48\Final\IZ48_230703_sess21',...
    'IZ48\Final\IZ48_230705_sess22','IZ48\Final\IZ48_230714_sess28'};  
  

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';

p = inputParser;
addParameter(p,'plotfig',true,@islogical);
addParameter(p,'savefig',false,@islogical);

parse(p,varargin{:});
plotfig = p.Results.plotfig;
savefig = p.Results.savefig;

Summary.AllcellType = [];
Summary.AllsessType = [];
Summary.AlllinCorr = [];
Summary.AllspatialCorr = [];
Summary.AlltoneCorr = [];
Summary.AlllinMapInit = [];
Summary.AlllinMapEnd = [];
Summary.AllspaceMap = [];
Summary.AllspaceMapAvg = [];

Summary.AllretMapLinInit = [];
Summary.AllretMapCorr = [];
Summary.AllretMapIncorr = [];

Summary.AlltoneMap = [];
Summary.AlltoneMapError1 = [];
Summary.AlltoneMapError2 = [];

Summary.AlltoneCorrError1 = [];
Summary.AlltoneCorrError2 = [];

Summary.AlllinInfo = [];
Summary.AllspatialInfo = [];
Summary.AlltoneInfo = [];
Summary.AlllinpeakRate = [];
Summary.AllspacepeakRate = [];
Summary.AlltonepeakRate = [];
Summary.AlllinField = [];
Summary.AlllinField1 = [];
Summary.AlllinField2 = [];
Summary.AlllinEndField = [];
Summary.AllspaceField = [];
Summary.AlltoneField = [];
Summary.AllretFieldlin = [];
Summary.AllretFieldCorrect = [];
Summary.AllretFieldlinEnd = [];

Summary.AllsessID = [];
Summary.AllcellID = [];

Summary.AlltoneNoToneCorr = [];
Summary.AlllinlinEndCorr = [];
Summary.AlltonelinEndCorr = [];
Summary.AllretlinToneCorr = [];
Summary.AllretlinlinCorr = [];
Summary.AllretlinEndToneCorr = [];
Summary.AllretlinlinEndCorr = [];

for ii = 1:length(sess)
    %% Load files
    cd(strcat(expPath,sess{ii}))    
    % This variable has correct and incorrect trials averaged by target
    % port
    file = dir(['*.rateMapsAvg.cellinfo.mat']);
    load(file(1).name);
    
    % Exception for control mice because they don't have the separate error
    % matrix calculated
    if strcmp(sess{ii}(1:4),'IZ41')==1 || strcmp(sess{ii}(1:4),'IZ36')==1 || strcmp(sess{ii}(1:4),'IZ45')==1 || strcmp(sess{ii}(1:4),'IZ46')==1
        errorMaps = firingMaps;
    else
        % This variable has correct and incorrect trials averaged by lick
        % location
        file = dir(['*.rateMapsAvgError.cellinfo.mat']);
        Maps = load(file(1).name);
        errorMaps = Maps.firingMaps;
    end
    
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
    spaceMapAvg(1:length(cell_metrics.UID),1:50) = nan;
    toneMap(1:length(cell_metrics.UID),1:50) = nan;   
    toneMapError1(1:length(cell_metrics.UID),1:50) = nan; 
    toneMapError2(1:length(cell_metrics.UID),1:50) = nan; 
    
    retMapLinInit(1:length(cell_metrics.UID),1:50) = nan;
    retMapCorr(1:length(cell_metrics.UID),1:50) = nan;
    retMapIncorr(1:length(cell_metrics.UID),1:50) = nan;
    retMapLinEnd(1:length(cell_metrics.UID),1:50) = nan;
    
    linInfo(1:length(cell_metrics.UID),1) = nan;
    spatialInfo(1:length(cell_metrics.UID),1) = nan;
    toneInfo(1:length(cell_metrics.UID),1) = nan;    
    
    linpeakRate(1:length(cell_metrics.UID),1) = nan;
    spacepeakRate(1:length(cell_metrics.UID),1) = nan;
    tonepeakRate(1:length(cell_metrics.UID),1) = nan;        

    linField(1:length(cell_metrics.UID),1) = nan;
    linField1(1:length(cell_metrics.UID),1) = nan;
    linField2(1:length(cell_metrics.UID),1) = nan;    
    linEndField(1:length(cell_metrics.UID),1) = nan;
    spaceField(1:length(cell_metrics.UID),1) = nan;
    toneField(1:length(cell_metrics.UID),1) = nan; 
    retFieldlin(1:length(cell_metrics.UID),1) = nan;
    retFieldCorrect(1:length(cell_metrics.UID),1) = nan;
    retFieldlinEnd(1:length(cell_metrics.UID),1) = nan;

    sessID(1:length(cell_metrics.UID),1) = nan;
    cellID(1:length(cell_metrics.UID),1) = nan;  
    
    toneNoToneCorr(1:length(cell_metrics.UID),1) = nan;
    linlinEndCorr(1:length(cell_metrics.UID),1) = nan;
    tonelinEndCorr(1:length(cell_metrics.UID),1) = nan;
    
    retlinToneCorr(1:length(cell_metrics.UID),1) = nan;
    retlinlinCorr(1:length(cell_metrics.UID),1) = nan;
    retlinEndToneCorr(1:length(cell_metrics.UID),1) = nan;
    retlinlinEndCorr(1:length(cell_metrics.UID),1) = nan;    

    toneCorrError1(1:length(cell_metrics.UID),1) = nan; 
    toneCorrError2(1:length(cell_metrics.UID),1) = nan;
    
    for kk=1:length(cell_metrics.UID)
        sessID(kk,1) = ii;
        cellID(kk,1) = kk;
    %% Extract relevant rate maps
        linMapInit(kk,:) = firingMaps.forward.rateMaps{kk}{1};
        retMapLinInit(kk,:) = firingMaps.reverse.rateMaps{kk}{1};
        retMapCorr(kk,:) = firingMaps.reverse.rateMaps{kk}{2};
        retMapIncorr(kk,:) = firingMaps.reverse.rateMaps{kk}{3};
        retMapLinEnd(kk,:) = firingMaps.reverse.rateMaps{kk}{9};
        
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
        retlinToneCorr(kk,1) = corrLin(1,2);
        
        % Average return track correlation (2 halves)
        corrLin = corrcoef(firingMaps.reverse.rateMaps{kk}{7}',firingMaps.reverse.rateMaps{kk}{8}','rows','complete');
        retlinlinCorr(kk,1) = corrLin(1,2);
        
        % Average return track correlation (no tone and no tone end)
        corrLin = corrcoef(firingMaps.reverse.rateMaps{kk}{1}',firingMaps.reverse.rateMaps{kk}{9}','rows','complete');
        retlinlinEndCorr(kk,1) = corrLin(1,2);
        
        % Average return track correlation (tone and no tone end)
        corrLin = corrcoef(firingMaps.reverse.rateMaps{kk}{2}',firingMaps.reverse.rateMaps{kk}{9}','rows','complete');
        retlinEndToneCorr(kk,1) = corrLin(1,2);        
        
        % Average spatial and tone correlation
        dataMat = [];
        dataMatTone = [];
        dataMatToneError1 = [];
        dataMatToneError2 = [];
        info = [];
        infoTone = [];
        for jj = 2:7           
            dataMat = [dataMat;firingMaps.forward.rateMaps{kk}{jj}];
            info = [info;firingMaps.forward.information{kk}{jj}];
            a = fillmissing(firingMaps.tone.rateMaps{kk}{jj},'linear');
            dataMatTone = [dataMatTone;a];
            infoTone = [infoTone;firingMaps.tone.information{kk}{jj}];
            
            dataMatToneError1 = [dataMatToneError1;firingMaps.tone.rateMaps{kk}{jj+6}];           
            dataMatToneError2 = [dataMatToneError2;errorMaps.tone.rateMaps{kk}{jj+6}];
            
        end 
        
        spaceMap(kk,:) = firingMaps.forward.rateMaps{kk}{7};      
        spaceMapAvg(kk,:) = nanmean(dataMat,1);
        toneMap(kk,:) = nanmean(dataMatTone,1);      
        toneMapError1(kk,:) = nanmean(dataMatToneError1,1);
        toneMapError2(kk,:) = nanmean(dataMatToneError2,1);
        
        a = corrcoef(toneMap(kk,:),toneMapError1(kk,:),'rows','complete'); 
        toneCorrError1(kk,1) = a(1,2); 
        
        a = corrcoef(toneMap(kk,:),toneMapError2(kk,:),'rows','complete'); 
        toneCorrError2(kk,1) = a(1,2);  
        
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
        
        Field_Info = detectFields(firingMaps.forward.rateMaps{kk}{26});
        if isempty(Field_Info)
            linField1(kk,1) = 0;
        else 
            linField1(kk,1) = 1;
        end
        
        Field_Info = detectFields(firingMaps.forward.rateMaps{kk}{27});
        if isempty(Field_Info)
            linField2(kk,1) = 0;
        else 
            linField2(kk,1) = 1;
        end
        
        if length(firingMaps.forward.rateMaps{kk})==28
            Field_Info = detectFields(linMapEnd(kk,:));
            if isempty(Field_Info)
                linEndField(kk,1) = 0;
            else 
                linEndField(kk,1) = 1;
            end
            Field_Info = detectFields(retMapLinEnd(kk,:));
            if isempty(Field_Info)
                retFieldlinEnd(kk,1) = 0;
            else 
                retFieldlinEnd(kk,1) = 1;
            end
        else
            linEndField(kk,1) = 0;
            retFieldlinEnd(kk,1) = 0;
        end
        
        Field_Info = detectFields(spaceMapAvg(kk,:));
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
    Summary.AllcellType = [Summary.AllcellType;cellType];
    Summary.AllsessType = [Summary.AllsessType;sessType];
    Summary.AlllinCorr = [Summary.AlllinCorr;linCorr];
    Summary.AllspatialCorr = [Summary.AllspatialCorr;spatialCorr];
    Summary.AlltoneCorr = [Summary.AlltoneCorr;toneCorr];
    Summary.AlllinMapInit = [Summary.AlllinMapInit;linMapInit];
    Summary.AlllinMapEnd = [Summary.AlllinMapEnd;linMapEnd];
    Summary.AllspaceMap = [Summary.AllspaceMap;spaceMap];
    Summary.AllspaceMapAvg = [Summary.AllspaceMapAvg;spaceMapAvg];
    Summary.AlltoneMap = [Summary.AlltoneMap;toneMap];
    Summary.AlltoneMapError1 = [Summary.AlltoneMapError1;toneMapError1];
    Summary.AlltoneMapError2 = [Summary.AlltoneMapError2;toneMapError2];
    Summary.AlllinInfo = [Summary.AlllinInfo;linInfo];
    Summary.AllspatialInfo = [Summary.AllspatialInfo;spatialInfo];
    Summary.AlltoneInfo = [Summary.AlltoneInfo;toneInfo];
    Summary.AlllinpeakRate = [Summary.AlllinpeakRate;linpeakRate];
    Summary.AllspacepeakRate = [Summary.AllspacepeakRate;spacepeakRate];
    Summary.AlltonepeakRate = [Summary.AlltonepeakRate;tonepeakRate];
    Summary.AlllinField = [Summary.AlllinField;linField];
    Summary.AlllinField1 = [Summary.AlllinField1;linField1];
    Summary.AlllinField2 = [Summary.AlllinField2;linField2];
    Summary.AlllinEndField = [Summary.AlllinEndField;linEndField];
    Summary.AllspaceField = [Summary.AllspaceField;spaceField];
    Summary.AlltoneField = [Summary.AlltoneField;toneField];    
    Summary.AllsessID = [Summary.AllsessID;sessID];
    Summary.AllcellID = [Summary.AllcellID;cellID];    
    Summary.AlltoneNoToneCorr = [Summary.AlltoneNoToneCorr;toneNoToneCorr];
    Summary.AlllinlinEndCorr = [Summary.AlllinlinEndCorr;linlinEndCorr];
    Summary.AlltonelinEndCorr = [Summary.AlltonelinEndCorr; tonelinEndCorr];
    
    Summary.AlltoneCorrError1 = [Summary.AlltoneCorrError1; toneCorrError1];
    Summary.AlltoneCorrError2 = [Summary.AlltoneCorrError2; toneCorrError2];
    
    Summary.AllretMapLinInit = [Summary.AllretMapLinInit;retMapLinInit];
    Summary.AllretMapCorr = [Summary.AllretMapCorr;retMapCorr];
    Summary.AllretMapIncorr = [Summary.AllretMapIncorr;retMapIncorr];
    Summary.AllretFieldlin = [Summary.AllretFieldlin;retFieldlin];
    Summary.AllretFieldCorrect = [Summary.AllretFieldCorrect;retFieldCorrect];
    Summary.AllretFieldlinEnd = [Summary.AllretFieldlinEnd;retFieldlinEnd];
      
    Summary.AllretlinToneCorr = [Summary.AllretlinToneCorr;retlinToneCorr];
    Summary.AllretlinlinCorr = [Summary.AllretlinlinCorr;retlinlinCorr];
    Summary.AllretlinlinEndCorr = [Summary.AllretlinlinEndCorr;retlinlinEndCorr];
    Summary.AllretlinEndToneCorr = [Summary.AllretlinEndToneCorr;retlinEndToneCorr];    
       
    clear cellType sessType linCorr spatialCorr toneCorr linMapInit linMapEnd spaceMap spaceMapAvg toneMap toneMapError1 toneMapError2 toneNoToneCorr linlinEndCorr linInfo ...
        spatialInfo toneInfo linpeakRate tonelinEndCorr spacepeakRate tonepeakRate linField linField1 linField2 linEndField spaceField toneField sessID cellID ...
        retMapLinInit retMapCorr retMapIncorr retFieldlin retFieldCorrect retlinToneCorr retlinlinCorr retlinlinEndCorr retlinEndToneCorr ...
        toneCorrError1 toneCorrError2 retFieldlinEnd retMapLinEnd
end

if plotfig
    YlGnBu=cbrewer('seq', 'YlGnBu', 11);
    fractTune = [];
    linPos = linspace(1,122,50);
    linTone = linspace(2000,22000,50);
    tag = {'Control','ACgN'};

    for ss = 1:2
        idxSess = Summary.AllsessType==(ss-1) & Summary.AllcellType == 1;
        idxActive = Summary.AllsessType==(ss-1) & Summary.AllcellType == 1 & nanmean(Summary.AllspaceMap,2)>0.2;
        idxMaps{1} = idxSess & Summary.AlllinField & Summary.AlllinCorr>0.1;
        idxMaps{2} = idxSess & Summary.AllspaceField & Summary.AllspatialCorr>0.1;
        idxMaps{3} = idxSess & Summary.AlltoneField & Summary.AlltoneCorr>0.1;    
        figure
        set(gcf,'Color','w')
        set(gcf,'renderer','painters')    
        set(gcf,'Position',[2384 318 617 547])
        for ii = 1:3
            fractTune(ss,ii) = sum(idxMaps{ii})/sum(idxActive);
            selectedlinMap = Summary.AlllinMapInit(idxMaps{ii},:);
            selectedlinMapEnd = Summary.AlllinMapEnd(idxMaps{ii},:);
            selectedtoneMap = Summary.AlltoneMap(idxMaps{ii},:);
            selectedspaceMap = Summary.AllspaceMapAvg(idxMaps{ii},:);
            selectedsessID = Summary.AllsessID(idxMaps{ii});
            selectedcellID = Summary.AllcellID(idxMaps{ii});

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
        if savefig
            saveas(gcf,strcat(expPath,'Compiled\populationMaps',tag{ss},'.png'));
            saveas(gcf,strcat(expPath,'Compiled\populationMaps',tag{ss},'eps'),'epsc');
            saveas(gcf,strcat(expPath,'Compiled\populationMaps',tag{ss},'fig'));
        end
    end

    figure
    set(gcf,'Color','w')
    set(gcf,'renderer','painters')
    set(gcf,'Position',[2384 318 617 547])
    subplot(2,2,1)
    idxSess = Summary.AllsessType==0 & Summary.AllcellType == 1 & (Summary.AlltoneField==1) & (Summary.AlltoneCorr >0.1);
    scatter(Summary.AllspatialCorr(idxSess),Summary.AlltoneCorr(idxSess),10,[187/243 86/243 149/243],'filled')
    hold on
    idxSess = Summary.AllsessType==0 & Summary.AllcellType == 1 & (Summary.AllspaceField==1) & (Summary.AllspatialCorr >0.1);
    scatter(Summary.AllspatialCorr(idxSess),Summary.AlltoneCorr(idxSess),10,[0.5 0.5 0.5],'filled')
    ylim([-0.2 1])
    xlim([-0.2 1])
    refline(1)
    title('Control')
    xlabel('Space correlation')
    ylabel('Tone correlation')

    subplot(2,2,2)
    idxSess = Summary.AllsessType==1 & Summary.AllcellType == 1 & (Summary.AllspaceField==1);% & (AllspatialCorr >0.1);
    scatter(Summary.AllspatialCorr(idxSess),Summary.AlltoneCorr(idxSess),10,[0.5 0.5 0.5],'filled')
    hold on
    idxSess = Summary.AllsessType==1 & Summary.AllcellType == 1 & (Summary.AlltoneField==1);% & (AlltoneCorr >0.1);
    scatter(Summary.AllspatialCorr(idxSess),Summary.AlltoneCorr(idxSess),10,[187/243 86/243 149/243],'filled')

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
    idxSess = Summary.AllsessType==0 & Summary.AllcellType == 1 & ((Summary.AllspaceField==1) | (Summary.AlllinField==1));
    data{1} = Summary.AlllinCorr(idxSess);
    %data{2} = AlllinlinEndCorr(idxSess);
    data{2} = Summary.AlltoneNoToneCorr(idxSess);

    idxSess = Summary.AllsessType==0 & Summary.AllcellType == 1 & ((Summary.AllretFieldlin | Summary.AllretFieldCorrect));
    data{3} = Summary.AllretlinlinCorr(idxSess);

    idxSess = Summary.AllsessType==1 & Summary.AllcellType == 1 & ((Summary.AllspaceField==1) | (Summary.AlllinField==1));
    data{4} = Summary.AlllinCorr(idxSess);
    data{5} = Summary.AlltoneNoToneCorr(idxSess);
    data{6} = Summary.AlltonelinEndCorr(idxSess);
    data{7} = Summary.AlllinlinEndCorr(idxSess);

    idxSess = Summary.AllsessType==1 & Summary.AllcellType == 1 & ((Summary.AllretFieldlin | Summary.AllretFieldCorrect));
    data{8} = Summary.AllretlinlinCorr(idxSess);
    data{9} = Summary.AllretlinToneCorr(idxSess);
    
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

    if savefig
        saveas(gcf,strcat(expPath,'Compiled\Fraction_Correlation.png'));
        saveas(gcf,strcat(expPath,'Compiled\Fraction_Correlation.eps'),'epsc');
        saveas(gcf,strcat(expPath,'Compiled\Fraction_Correlation.fig'));
        save(strcat(expPath,'Compiled\Fraction_Correlation.mat'),'stats'); 
    end

    figure
    set(gcf,'Color','w')
    set(gcf,'renderer','painters')
    set(gcf,'Position',[2384 318 617 547])

    for ss = 1:2
        idxSess = Summary.AllsessType==(ss-1) & Summary.AllcellType == 1;
        idxActive = Summary.AllsessType==(ss-1) & Summary.AllcellType == 1 & nanmean(Summary.AllretMapLinInit,2)>0.2;
        idxMaps = idxSess &  (Summary.AllretFieldlin | Summary.AllretFieldCorrect);
        if ss == 1
            iRange = [1 2];
        else
            iRange = [1 2 3];
        end
        for ii = iRange       
            selectedlinMap = Summary.AllretMapLinInit(idxMaps,:);
            selectedcorrMap = Summary.AllretMapCorr(idxMaps,:);
            selectedincorrMap = Summary.AllretMapIncorr(idxMaps,:);
            selectedsessID = Summary.AllsessID(idxMaps);
            selectedcellID = Summary.AllcellID(idxMaps);

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
    idxSess = Summary.AllsessType==0 & Summary.AllcellType == 1 & ((Summary.AllretFieldlin | Summary.AllretFieldCorrect));
    data = [];

    data{1} = Summary.AllretlinlinCorr(idxSess);
    data{2} = Summary.AllretlinToneCorr(idxSess);
    data{3} = Summary.AllretlinEndToneCorr(idxSess);
    data{4} = Summary.AllretlinlinEndCorr(idxSess);
    
    idxSess = Summary.AllsessType==1 & Summary.AllcellType == 1 & ((Summary.AllretFieldlin | Summary.AllretFieldCorrect));
    data{5} = Summary.AllretlinlinCorr(idxSess);
    data{6} = Summary.AllretlinToneCorr(idxSess);
    data{7} = Summary.AllretlinEndToneCorr(idxSess);
    data{8} = Summary.AllretlinlinEndCorr(idxSess);

    stats = groupStats(data,[],'inAxis',true,'color',col);
    ylabel('Spatial correlation')
    
    if savefig
        saveas(gcf,strcat(expPath,'Compiled\returnpopulationMaps.png'));
        saveas(gcf,strcat(expPath,'Compiled\returnpopulationMaps.eps'),'epsc');
        saveas(gcf,strcat(expPath,'Compiled\returnpopulationMaps.fig'));
    end
end

end
