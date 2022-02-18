function compileMiceResponseCA3(varargin)

p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
parse(p,varargin{:});

parentDir = p.Results.parentDir;

mice = {'IZ34','IZ33','IZ33','IZ27','IZ28','IZ29'};
condition = {'Final','Saline'};

CellClass = {'Narrow Interneuron','Pyramidal Cell','Wide Interneuron'};
BrainRegion = {'CA1','DG'};%,'CA3'};

col = [0.9856 0.7372 0.2537;...
       0.1986 0.7214 0.6310;...
       0.2305 0.2510 0.6173];
% 
colormap default
% YlGnBu=cbrewer('seq', 'YlGnBu', 11);
% YlGnBu = YlGnBu(11:-1:1,:);
% colormap(YlGnBu)

for cond = 1:length(condition)

    compiledData.putativeClass = [];
    compiledData.region = [];
    compiledData.firingrate = [];
    
    figure(cond)
    set(gcf,'Position',[100 100 1200 800])
    set(gcf,'renderer','Painters')
    
    for m = 1:length(mice)

        cd(strcat(parentDir, mice{m},'\',condition{cond}));
        allSess = dir('*_sess*');
        
        %% Start collecting data
        for ii = 1:size(allSess,1)
           cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
           [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);                
           
            % Load cell metrics
            if exist([sessionInfo.FileName '.cell_metrics.cellinfo.mat'],'file') 
                load([sessionInfo.FileName '.cell_metrics.cellinfo.mat']);
            else
                continue
            end
            
            if exist([sessionInfo.FileName '.MergePoints.events.mat'],'file') 
                load([sessionInfo.FileName '.MergePoints.events.mat']);
            else
                continue
            end
                    
            % Assume bin size in cell metrics is 180 seconds  - compile FRs
            % Compile firings rates for 30 minutes before, and 1 hour after
            % injection. Assume injection timestamp is the start of 3rd
            % mergepoint
            injTime = MergePoints.timestamps(3,1); 
            [~,idxstart] = min(abs(cell_metrics.general.responseCurves.firingRateAcrossTime.x_bins - injTime));
            
            % Now take 10 bins before, and 20 bins after
            firingResponse = cell2mat(cell_metrics.responseCurves.firingRateAcrossTime)';
            %firingResponse = firingResponse(:,idxstart-15:idxstart+23);
            range = min(idxstart+30,size(firingResponse,2));
            firingResponse = firingResponse(:,idxstart-16:range);
            
            compiledData.putativeClass = [compiledData.putativeClass cell_metrics.putativeCellType];
            compiledData.region = [compiledData.region cell_metrics.brainRegion];
            compiledData.firingrate = catpad(1,compiledData.firingrate,firingResponse);
            
        end
    end
    
    % Now to make some plots
    for jj = 1:length(compiledData.region)  
        % Assign numerical tag to putative class
        if strcmp(compiledData.putativeClass{jj},'Narrow Interneuron') == 1
            putativeClass(jj) = 1;
        elseif strcmp(compiledData.putativeClass{jj},'Pyramidal Cell') == 1
            putativeClass(jj) = 2;
        elseif strcmp(compiledData.putativeClass{jj},'Wide Interneuron') == 1
            putativeClass(jj) = 3;
        else 
            putativeClass(jj) = 4;                       
        end                    

        % Assign numerical tag to region
        if strcmp(compiledData.region{jj},'CA1') == 1
            region(jj) = 1;
        elseif strcmp(compiledData.region{jj},'DG') == 1
            region(jj) = 2;
        elseif strcmp(compiledData.region{jj},'CA3') == 1
            region(jj) = 2;
        else 
            region(jj) = 4;                       
        end
        
        firingrate(jj,:) = compiledData.firingrate(jj,:)./max(compiledData.firingrate(jj,:));%
        baselineFR(jj) = nanmean(compiledData.firingrate(jj,1:14));
        %rateRatio(jj) = nanmean(compiledData.firingrate(jj,36:38))./nanmean(compiledData.firingrate(jj,1:3));
        rateRatio(jj) = nanmean(compiledData.firingrate(jj,23:38))./nanmean(compiledData.firingrate(jj,1:14));
    end
    
    for c = 1:3
        statData = [];
        statGroup = [];
        for r = 1:2%3
            subplot(3,5,5*(c-1)+r)
            idx = region == r & putativeClass == c & baselineFR>0.2; 
            %[~,idxsort] = sort(rateRatio(idx));
            [~,idxsort] = sort(baselineFR(idx));
            tempFR = firingrate(idx,:);
            %tempFR = zscore(tempFR,[],2);
            imagesc(tempFR(idxsort,:),'AlphaData',~isnan(tempFR(idxsort,:)))
            caxis([0 1])
            line([15 15],[0 size(tempFR,1)+0.5],'Color','white','LineWidth',2)
            title(strcat(BrainRegion{r},num2str(length(idxsort))),'fontweight','bold')
            xlim([0 43])
            if r == 1            
                ylabel(CellClass{c},'fontweight','bold')           
            end
            
            % Plot normalized rates
            subplot(3,5,5*(c-1)+4)
            a = nanmean(tempFR,1);
            meanpsth =a/mean(a(1:14));
            %meanpsth =a;
            stdpsth = nanstd(tempFR,[],1)./sqrt(size(tempFR,1));
            x = 1:1:size(stdpsth,2);
            hold on
            fill([x'; fliplr(x)'],[meanpsth'-stdpsth';flipud(meanpsth'+stdpsth')],col(r,:),'linestyle','none','FaceAlpha',0.5);
            line(x,meanpsth,'Color',col(r,:),'LineWidth',1.5)   
            line([15 15],[0 1.5],'Color','red','LineWidth',2)
            %xlim([0 39])
            ylim([0 1.5]) 
            
            statData = [statData rateRatio(idx)];
            statGroup = [statGroup ones(1,length(rateRatio(idx)))*r];           
        end
        
        subplot(3,5,5*(c-1)+5)
        h = cdfplot(statData(statGroup==1));
        set(h,'Color',[0.5 0.5 0.5],'LineWidth',2)
        hold on
        h = cdfplot(statData(statGroup==2));
        set(h,'Color','r','LineWidth',2)
%         h = cdfplot(statData(statGroup==3));
%         set(h,'Color','c','LineWidth',2)
        xlim([0 2])
        stats{c} = groupStats(statData,statGroup,'inAxis',true,'doPlot',false);         
        for x = 1:2%3
            [stats{c}.signrank.p(x), ~,stats{c}.signrank.stats{x}] = signrank(statData(statGroup==x));
        end
    end
   
    saveas(figure(cond),strcat(parentDir,'Compiled\Unit PSTH\CA3UnitResponse_',condition{cond},'.png'))
    saveas(figure(cond),strcat(parentDir,'Compiled\Unit PSTH\CA3UnitResponse_',condition{cond},'.fig'),'fig')    
    saveas(figure(cond),strcat(parentDir,'Compiled\Unit PSTH\CA3UnitResponse_',condition{cond},'.eps'),'epsc') 
    save(strcat(parentDir,'Compiled\Unit PSTH\UnitPSTHStatsCA3_',condition{cond},'.mat'),'stats')
    clear putativeClass region firingRate baselineFR rateRatio
end





  