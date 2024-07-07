function Summary = linkToneCellFiringtoDeceleration(varargin)

p = inputParser;
addParameter(p,'plotfig',true,@islogical);

parse(p,varargin{:});
plotfig = p.Results.plotfig;

sess = {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
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

Summary.psthReward{1} = [];
Summary.psthReward{2} = [];
Summary.psthReward{3} = [];
Summary.sessID = [];

for ii = 1:length(sess)

    %% Load files
    cd(strcat(expPath,sess{ii}))    
    % Find deceleration points
    Dec = findDecelerationPoints('plotfig',false);

    file = dir(['*.rateMapsAvgnotLog.cellinfo.mat']);
    load(file(1).name); 
    file = dir(['*.Tracking.Behavior.mat']);
    load(file(1).name);
    file = dir(['*TrialBehavior.Behavior.mat']);
    load(file(1).name);    
    file = dir(['*spikes.cellinfo.mat']);
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
        Field_Info = detectFields(toneMap,'maxFieldSize',40);
        if isempty(Field_Info)
            toneField = 0;
        else 
            toneField = 1;
        end        

        [~,idxMax] = max(toneMap);     
        % If its a lick cell, calculate its psth around these two types of
        % deceleration
        if cellType == 1 && (toneField == 1) && (toneCorr > 0.1) && idxMax>40
            for dd = 1:3 % Extract end licks, and "middle" licks
                if dd==1
                    st = Dec.ts(Dec.decType==1);
                elseif dd == 2
                    st = Dec.ts(Dec.decType==2);
                else
                    st = Dec.ts(Dec.decType>=3);
                end

                if ~isempty(st)
                    [stccg, tPSTH] = CCG({spikes.times{kk} st},[],'binSize',0.05,'duration',2,'norm','rate');                
                    Summary.psthReward{dd} = [Summary.psthReward{dd}; stccg(:,2,1)'];               
                    if dd == 1                    
                        Summary.sessID = [Summary.sessID; ii]; 
                    end
                else
                    fillArr(1,1:length(tPSTH)) = nan;
                    Summary.psthReward{dd} = [Summary.psthReward{dd}; fillArr]; 
                    if dd == 1                    
                        Summary.sessID = [Summary.sessID; ii]; 
                    end
                end
            end            
        end
    end
end

save('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\DeliberationPSTHSummary.mat', 'Summary'); 

if plotfig
    fig2 = figure;
    
    col = [83/255 0/255 0/255;...
        184/255 15/255 10/255;...
        241/255 114/255 42/255];
    
    spec=cbrewer('seq', 'Blues', 20);
    spec(spec>1) = 1;
    spec(spec<0) = 0;
    
    [maxRate,maxRateIdx] = max(Summary.psthReward{1},[],2);
    [~,idxmax] = sort(maxRateIdx,'ascend'); 
    
    for ii = 1:3
        ax1 = subplot(1,4,ii);
        norm = zscore(Summary.psthReward{ii},[],2);
        imagesc(tPSTH, 1:size(Summary.psthReward{ii},1),norm(idxmax,:))
        colormap(ax1,spec)
        caxis([-1 3])
        xlim([tPSTH(1) tPSTH(end)])
    
        plotAvgStd(Summary.psthReward{ii},1,4,4,fig2,tPSTH,col(ii,:))
    end
end

end

function plotAvgStd(array,numrows,numcol,subplotlocation,figureHandle,xAxis,col)

    subplot(numrows, numcol, subplotlocation, 'Parent', figureHandle);

    meanpsth = nanmean(array,1);
    stdpsth = nanstd(array,1)./sqrt(size(array,1));
    lArr  = meanpsth-stdpsth;
    uArr = meanpsth+stdpsth;

    fill([xAxis; flipud(xAxis)],[lArr'; flipud(uArr')],col,'linestyle','none','FaceAlpha',0.5);                    
    hold on
    hi = line(xAxis,meanpsth,'LineWidth',1,'Color',col);

end