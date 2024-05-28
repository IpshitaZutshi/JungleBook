function Summary  = compareActionNoToneTrials(varargin)

sess={'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
    'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
    'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   
    'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...
    'IZ40\Final\IZ40_220708_sess17','IZ40\Final\IZ40_220714_sess18',...27
    'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
    'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
    'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
    'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    
    'IZ44\Final\IZ44_220827_sess4', 'IZ44\Final\IZ44_220828_sess5',...
    'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...40
    'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
    'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',...
    'IZ47\Final\IZ47_230626_sess15','IZ47\Final\IZ47_230707_sess24',...
    'IZ47\Final\IZ47_230710_sess25','IZ47\Final\IZ47_230712_sess27',...
    'IZ48\Final\IZ48_230628_sess17','IZ48\Final\IZ48_230703_sess21',...50
    'IZ48\Final\IZ48_230705_sess22','IZ48\Final\IZ48_230714_sess28'};  

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';

p = inputParser;
addParameter(p,'plotfig',true,@islogical);
addParameter(p,'savefig',false,@islogical);

parse(p,varargin{:});
plotfig = p.Results.plotfig;
savefig = p.Results.savefig;

Summary.AllrateMapLin = [];
Summary.AllrateMapLin1 = [];
Summary.AllrateMapLin2 = [];
Summary.AllrateMapLinEnd = [];
Summary.AllrateMapSpace = [];
Summary.AllrateMaptone = [];
Summary.AlllinField = [];
Summary.AlllinFieldEnd = [];
Summary.AllsessID = [];
Summary.AllcellID = [];

for ii = 1:length(sess)
    %% Load files
    cd(strcat(expPath,sess{ii}))    
    file = dir(['*.rateMapsAvgnotLog.cellinfo.mat']);
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
    sessID = [];
    cellID = [];    
    
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
        if cellType == 1 && (toneField == 1) && (toneCorr > 0.1) %&& idxMax>15
            
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
            
            sessID = [sessID; ii];
            cellID = [cellID; kk];
        end
            
    end
    
    Summary.AllrateMapLin = [Summary.AllrateMapLin;linMapInit];
    Summary.AllrateMapLin1 = [Summary.AllrateMapLin1;linMapInit1];
    Summary.AllrateMapLin2 = [Summary.AllrateMapLin2;linMapInit2];
    Summary.AllrateMapLinEnd = [Summary.AllrateMapLinEnd;linMapEnd];
    Summary.AllrateMapSpace = [Summary.AllrateMapSpace;spaceMap1];
    Summary.AllrateMaptone = [Summary.AllrateMaptone;toneMap1];
    Summary.AlllinField = [Summary.AlllinField;linField];
    Summary.AlllinFieldEnd = [Summary.AlllinFieldEnd;linFieldEnd];
    Summary.AllsessID = [Summary.AllsessID;sessID];
    Summary.AllcellID = [Summary.AllcellID;cellID];
    
    clear linMapInit linMapEnd spaceMap1 toneMap1 linField linFieldEnd linMapInit1 linMapInit2 sessID cellID
end

if plotfig
    normlinMap = (Summary.AllrateMapLin-nanmean(Summary.AllrateMapLin,2))./nanstd(Summary.AllrateMapLin,[],2);
    normlinMap1 = (Summary.AllrateMapLin1-nanmean(Summary.AllrateMapLin1,2))./nanstd(Summary.AllrateMapLin1,[],2);
    normlinMap2 = (Summary.AllrateMapLin2-nanmean(Summary.AllrateMapLin2,2))./nanstd(Summary.AllrateMapLin2,[],2);
    normlinMapEnd = (Summary.AllrateMapLinEnd-nanmean(Summary.AllrateMapLinEnd,2))./nanstd(Summary.AllrateMapLinEnd,[],2);
    normspaceMap = (Summary.AllrateMapSpace-nanmean(Summary.AllrateMapSpace,2))./nanstd(Summary.AllrateMapSpace,[],2);
    normtoneMap = (Summary.AllrateMaptone-nanmean(Summary.AllrateMaptone,2))./nanstd(Summary.AllrateMaptone,[],2);

    selectsessID = Summary.AllsessID;
    selectcellID = Summary.AllcellID;
    
    sortsessID = selectsessID(sortidx);
    sortcellID = selectcellID(sortidx);
        
    figure
    set(gcf,'Renderer','painters')
    set(gcf,'Color','w')
    linPos = linspace(1,122,50);
    YlGnBu=cbrewer('seq', 'YlGnBu', 11);
    colormap(YlGnBu)

    subplot(3,3,1)
    [maxLin,idxLin] = max(Summary.AllrateMapLin,[],2);
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
    bar([sum(Summary.AlllinField)./length(Summary.AlllinField) sum(Summary.AlllinFieldEnd)./length(Summary.AlllinFieldEnd)])

    selecttoneMap = normtoneMap(Summary.AlllinField==1,:);
    selectspaceMap = normspaceMap(Summary.AlllinField==1,:);
    selectlinMap = normlinMap(Summary.AlllinField==1,:);
    selectlinMap1 = normlinMap1(Summary.AlllinField==1,:);
    selectlinMap2 = normlinMap2(Summary.AlllinField==1,:);
    
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
    if savefig
        saveas(gcf,strcat(expPath,'Compiled\ToneCellsLinearTrials.png'));
        saveas(gcf,strcat(expPath,'Compiled\ToneCellsLinearTrials.eps'),'epsc');
        saveas(gcf,strcat(expPath,'Compiled\ToneCellsLinearTrials.fig'));
    end
end
end
