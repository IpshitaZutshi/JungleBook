function Summary = compileProbeTrials(varargin)

sesstoAnalyze = {'IZ47\Final\IZ47_230626_sess15','IZ47\Final\IZ47_230707_sess24',...
    'IZ47\Final\IZ47_230710_sess25','IZ47\Final\IZ47_230712_sess27',...
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
Summary.AllspatialCorr = [];
Summary.AlltoneCorr = [];
Summary.AllspatialCorrProbe = [];
Summary.AlltoneCorrProbe = [];
Summary.AllspaceMap = [];
Summary.AlltoneMap = [];
Summary.AllspaceMapProbe = [];
Summary.AlltoneMapProbe = [];
Summary.AllspaceField = [];
Summary.AlltoneField = [];
Summary.AllspaceFieldProbe = [];
Summary.AlltoneFieldProbe = [];
Summary.AllspatialInfo = [];
Summary.AlltoneInfo = [];
Summary.AllspatialInfoProbe = [];
Summary.AlltoneInfoProbe = [];

for ii = 1:length(sesstoAnalyze)
 %% Load files
    cd(strcat(expPath,sesstoAnalyze{ii}))    
    file = dir(['*.rateMapsAvgLickLocProbe.cellinfo.mat']);
    load(file(1).name);
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
    spaceMapProbe(1:length(cell_metrics.UID),1:50) = nan;
    toneMapProbe(1:length(cell_metrics.UID),1:50) = nan;
    spatialCorr(1:length(cell_metrics.UID),1) = nan;
    toneCorr(1:length(cell_metrics.UID),1) = nan;    
    spatialCorrProbe(1:length(cell_metrics.UID),1) = nan;
    toneCorrProbe(1:length(cell_metrics.UID),1) = nan;    
    spaceField(1:length(cell_metrics.UID),1) = nan;
    toneField(1:length(cell_metrics.UID),1) = nan;   
    spaceFieldProbe(1:length(cell_metrics.UID),1) = nan;
    toneFieldProbe(1:length(cell_metrics.UID),1) = nan;       
    spatialInfo(1:length(cell_metrics.UID),1) = nan;
    toneInfo(1:length(cell_metrics.UID),1) = nan; 
    spatialInfoProbe(1:length(cell_metrics.UID),1) = nan;
    toneInfoProbe(1:length(cell_metrics.UID),1) = nan; 
    
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
        dataMatProbe =[];
        dataMatToneProbe = [];   
        
        info = [];
        infoTone = [];
        infoProbe = [];
        infoToneProbe = [];
        
        for jj = 1:6                
            dataMat = [dataMat;firingMaps.forward.rateMaps{kk}{jj}];
            info = [info;firingMaps.forward.information{kk}{jj}];
            
            a = fillmissing(firingMaps.tone.rateMaps{kk}{jj},'linear');
            dataMatTone = [dataMatTone;a];
            infoTone = [infoTone;firingMaps.tone.information{kk}{jj}];
              
            dataMatProbe = [dataMatProbe;firingMaps.forward.rateMaps{kk}{jj+6}];
            infoProbe = [infoProbe;firingMaps.forward.information{kk}{jj+6}];
            
            a = fillmissing(firingMaps.tone.rateMaps{kk}{jj+6},'linear');
            dataMatToneProbe = [dataMatToneProbe;a];
            infoToneProbe = [infoToneProbe;firingMaps.tone.information{kk}{jj+6}];
        end
        
        spaceMap(kk,:) = nanmean(dataMat,1);
        toneMap(kk,:) = nanmean(dataMatTone,1);        
        
        spaceMapProbe(kk,:) = nanmean(dataMatProbe,1);
        toneMapProbe(kk,:) = nanmean(dataMatToneProbe,1);                
        
        corrTone = []; 
        corrSpace = [];
        corrToneProbe = []; 
        corrSpaceProbe = [];
        
        for pp = 1:6
           for jj = (pp+1):6
               a = corrcoef(dataMat(pp,:),dataMat(jj,:),'rows','complete');
               corrSpace = [corrSpace a(1,2)];               
               a = corrcoef(dataMatTone(pp,:),dataMatTone(jj,:),'rows','complete');         
               corrTone = [corrTone a(1,2)];

               a = corrcoef(dataMatProbe(pp,:),dataMatProbe(jj,:),'rows','complete');
               corrSpaceProbe = [corrSpaceProbe a(1,2)];               
               a = corrcoef(dataMatToneProbe(pp,:),dataMatToneProbe(jj,:),'rows','complete');         
               corrToneProbe = [corrToneProbe a(1,2)];               
           end
        end     
        toneCorr(kk,1) = nanmean(corrTone);  
        spatialCorr(kk,1) = nanmean(corrSpace);
        
        toneCorrProbe(kk,1) = nanmean(corrToneProbe);  
        spatialCorrProbe(kk,1) = nanmean(corrSpaceProbe);
        
        spatialInfo(kk,1) = nanmean(info);
        toneInfo(kk,1) = nanmean(infoTone);           
        
        spatialInfoProbe(kk,1) = nanmean(infoProbe);
        toneInfoProbe(kk,1) = nanmean(infoToneProbe);           
        
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
        
        Field_Info = detectFields(spaceMapProbe(kk,:));
        if isempty(Field_Info)
            spaceFieldProbe(kk,1) = 0;
        else 
            spaceFieldProbe(kk,1) = 1;
        end
        
        Field_Info = detectFields(toneMapProbe(kk,:));
        if isempty(Field_Info)
            toneFieldProbe(kk,1) = 0;
        else 
            toneFieldProbe(kk,1) = 1;
        end           
         
        dataMat =[];
        dataMatTone = [];
        dataMatProbe =[];
        dataMatToneProbe = [];   
        
        for jj = 3:6                
            dataMat = [dataMat;firingMaps.forward.rateMaps{kk}{jj}];
            
            a = fillmissing(firingMaps.tone.rateMaps{kk}{jj},'linear');
            dataMatTone = [dataMatTone;a];
              
            dataMatProbe = [dataMatProbe;firingMaps.forward.rateMaps{kk}{jj+6}];
            
            a = fillmissing(firingMaps.tone.rateMaps{kk}{jj+6},'linear');
            dataMatToneProbe = [dataMatToneProbe;a];
        end
        
        spaceMap(kk,:) = nanmean(dataMat,1);
        toneMap(kk,:) = nanmean(dataMatTone,1);        
        
        spaceMapProbe(kk,:) = nanmean(dataMatProbe,1);
        toneMapProbe(kk,:) = nanmean(dataMatToneProbe,1);       
        
    end
    
    Summary.AllcellType = [Summary.AllcellType;cellType];
    Summary.AllspatialCorr = [Summary.AllspatialCorr;spatialCorr];
    Summary.AlltoneCorr = [Summary.AlltoneCorr;toneCorr];
    Summary.AllspatialCorrProbe = [Summary.AllspatialCorrProbe;spatialCorrProbe];
    Summary.AlltoneCorrProbe = [Summary.AlltoneCorrProbe;toneCorrProbe];
    Summary.AllspaceMap = [Summary.AllspaceMap;spaceMap];
    Summary.AlltoneMap = [Summary.AlltoneMap;toneMap];
    Summary.AllspaceMapProbe = [Summary.AllspaceMapProbe;spaceMapProbe];
    Summary.AlltoneMapProbe = [Summary.AlltoneMapProbe;toneMapProbe];
    Summary.AllspaceField = [Summary.AllspaceField;spaceField];
    Summary.AlltoneField = [Summary.AlltoneField;toneField];     
    Summary.AllspaceFieldProbe = [Summary.AllspaceFieldProbe;spaceFieldProbe];
    Summary.AlltoneFieldProbe = [Summary.AlltoneFieldProbe;toneFieldProbe];    
    Summary.AllspatialInfo = [Summary.AllspatialInfo;spatialInfo];
    Summary.AlltoneInfo = [Summary.AlltoneInfo;toneInfo];
    Summary.AllspatialInfoProbe = [Summary.AllspatialInfoProbe;spatialInfoProbe];
    Summary.AlltoneInfoProbe = [Summary.AlltoneInfoProbe;toneInfoProbe];
    
    clear cellType spatialCorr toneCorr spatialCorrProbe toneCorrProbe spaceMap toneMap spaceMapProbe toneMapProbe ...
        spaceField toneField spaceFieldProbe toneFieldProbe spatialInfo toneInfo spatialInfoProbe toneInfoProbe
    
end

if plotfig
    figure
    set(gcf,'Color','w')
    set(gcf,'renderer','painters')    
    set(gcf,'Position',[408 38 1107 946])

    YlGnBu=cbrewer('seq', 'YlGnBu', 11);
    fractTune = [];
    linPos = linspace(1,122,50);
    linTone = linspace(2000,22000,50);

    activeCells = Summary.AllcellType & (nanmean(Summary.AllspaceMap,2)>0.2);
    idxMaps{1} = Summary.AllcellType & Summary.AllspaceField & Summary.AllspatialCorr>0.1; % Idx for spatial cells
    idxMaps{2} = Summary.AllcellType & Summary.AllspaceFieldProbe & Summary.AllspatialCorrProbe>0.1; % Idx for spatial cells
    idxMaps{3} = Summary.AllcellType & Summary.AlltoneField & Summary.AlltoneCorr>0.1;  % Idx for tone cells 
    idxMaps{4} = Summary.AllcellType & Summary.AlltoneFieldProbe & Summary.AlltoneCorrProbe>0.1;  % Idx for tone cells 
    fractTune = round(([sum(idxMaps{1}) sum(idxMaps{2}) sum(idxMaps{3}) sum(idxMaps{4})]/sum(activeCells))*100); 

    subplot(5,4,1)
    bar(fractTune)
    text(1:length(fractTune),fractTune,num2str(fractTune'),'vert','bottom','horiz','center'); 
    box off

    clear idxMaps
    % idxMaps{1} = AllcellType & ((AllspaceField & AllspatialCorr>0.1) | (AllspaceFieldProbe & AllspatialCorrProbe>0.1)); % Idx for spatial cells
    % idxMaps{2} = AllcellType & ((AlltoneField & AlltoneCorr>0.1) | (AlltoneFieldProbe & AlltoneCorrProbe>0.1));  % Idx for tone cells 

    idxMaps{1} = Summary.AllcellType & Summary.AllspaceField & Summary.AllspatialCorr>0.1; % Idx for spatial cells
    idxMaps{2} = Summary.AllcellType & Summary.AllspaceFieldProbe & Summary.AllspatialCorrProbe>0.1; % Idx for spatial cells
    idxMaps{3} = Summary.AllcellType & Summary.AlltoneField & Summary.AlltoneCorr>0.1;  % Idx for tone cells 
    idxMaps{4} = Summary.AllcellType & Summary.AlltoneFieldProbe & Summary.AlltoneCorrProbe>0.1;  % Idx for tone cells 


    for ii = 1:4

       if ii == 1
           infoData{1} = Summary.AllspatialInfo(idxMaps{ii});
           infoData{2} = Summary.AllspatialInfoProbe(idxMaps{ii+1});

           subplot(5,4,2)
           stats.spatialInfo = groupStats(infoData,[],'inAxis',true,'labelSummary',false);%,'repeatedMeasures',true);       
           title('Spatial Information')
       elseif ii == 2
           infoData{1} = Summary.AlltoneInfo(idxMaps{ii+1});
           infoData{2} = Summary.AlltoneInfoProbe(idxMaps{ii+2});

           subplot(5,4,3)
           stats.toneInfo = groupStats(infoData,[],'inAxis',true,'labelSummary',false);%,'repeatedMeasures',true);         
           title('Tone Information')
       end
       selectedtoneMap = Summary.AlltoneMap(idxMaps{ii},:);
       selectedspaceMap = Summary.AllspaceMap(idxMaps{ii},:);
       selectedtoneMapProbe = Summary.AlltoneMapProbe(idxMaps{ii},:);
       selectedspaceMapProbe = Summary.AllspaceMapProbe(idxMaps{ii},:);

       [maxSpace,idxSpace] = max(selectedspaceMap,[],2);        
       [maxTone,idxTone] = max(selectedtoneMap,[],2);

       [maxSpaceProbe,idxSpaceProbe] = max(selectedspaceMapProbe,[],2);        
       [maxToneProbe,idxToneProbe] = max(selectedtoneMapProbe,[],2);   

        %Z score normalize
        normtoneMap = (selectedtoneMap-nanmean(selectedtoneMap,2))./nanstd(selectedtoneMap,[],2);
        normspaceMap = (selectedspaceMap-nanmean(selectedspaceMap,2))./nanstd(selectedspaceMap,[],2);
        normtoneMapProbe = (selectedtoneMapProbe-nanmean(selectedtoneMapProbe,2))./nanstd(selectedtoneMapProbe,[],2);
        normspaceMapProbe = (selectedspaceMapProbe-nanmean(selectedspaceMapProbe,2))./nanstd(selectedspaceMapProbe,[],2);

        sortidx = [];
        sortidxProbe = [];
        if ii == 1            
            [~,sortidx] = sort(idxSpace,'ascend');
        elseif ii == 2
            [~,sortidx] = sort(idxSpaceProbe,'ascend');
        elseif ii == 3            
            [~,sortidx] = sort(idxTone,'ascend'); 
        else
            [~,sortidx] = sort(idxToneProbe,'ascend'); 
        end

        colormap(YlGnBu)
        subplot(5,4,4*(ii-1)+5)
        imagesc(linPos, 1:length(sortidx),normspaceMap(sortidx,:))
        ylabel('Cell ID')
        xlabel('Frequency')  
        title('Task spatial')
        caxis([-1 4])

        subplot(5,4,4*(ii-1)+6)
        imagesc(linPos, 1:length(sortidx),normspaceMapProbe(sortidx,:))
        ylabel('Cell ID')
        xlabel('Position')
        caxis([-1 4])
        title('Task spatial Probe')  

        subplot(5,4,4*(ii-1)+7)
        imagesc(linTone, 1:length(sortidx),normtoneMap(sortidx,:))
        ylabel('Cell ID')
        xlabel('Position')
        caxis([-1 4])
        title('Task auditory')  

        subplot(5,4,4*(ii-1)+8)
        imagesc(linTone, 1:length(sortidx),normtoneMapProbe(sortidx,:))
        ylabel('Cell ID')
        xlabel('Position')
        caxis([-1 4])
        title('Task auditory Probe')  

    %     subplot(5,4,8*(ii-1)+9)
    %     imagesc(linPos, 1:length(sortidxProbe),normspaceMap(sortidxProbe,:))
    %     ylabel('Cell ID')
    %     xlabel('Frequency')  
    %     title('Task spatial')
    %     caxis([-1 4])
    %     
    %     subplot(5,4,8*(ii-1)+10)
    %     imagesc(linPos, 1:length(sortidxProbe),normspaceMapProbe(sortidxProbe,:))
    %     ylabel('Cell ID')
    %     xlabel('Position')
    %     caxis([-1 4])
    %     title('Task spatial Probe')  
    %     
    %     subplot(5,4,8*(ii-1)+11)
    %     imagesc(linTone, 1:length(sortidxProbe),normtoneMap(sortidxProbe,:))
    %     ylabel('Cell ID')
    %     xlabel('Position')
    %     caxis([-1 4])
    %     title('Task auditory')  
    %     
    %     subplot(5,4,8*(ii-1)+12)
    %     imagesc(linTone, 1:length(sortidxProbe),normtoneMapProbe(sortidxProbe,:))
    %     ylabel('Cell ID')
    %     xlabel('Position')
    %     caxis([-1 4])
    %     title('Task auditory Probe')      
    end
    %
    if savefig
        saveas(gcf,strcat(expPath,'Compiled\probeTrialMaps.png'));
        saveas(gcf,strcat(expPath,'Compiled\probeTrialMaps.eps'),'epsc');
        saveas(gcf,strcat(expPath,'Compiled\probeTrialMaps.fig'));
    end
end
end