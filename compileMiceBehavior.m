% Compile data across all sessions

function pVal = compileMiceBehavior(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
addParameter(p,'avgSessions',false,@islogical);
parse(p,varargin{:});

parentDir = p.Results.parentDir;
avgSessions = p.Results.avgSessions;

tag = 'CA3Saline';

if strcmp(tag,'CA1') == 1
    mice = {'IZ15\Final','IZ18\Final','IZ20\Final','IZ21\Final','IZ30\Final','IZ31\Final'};
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ11\Final','IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final',...
        'IZ21\Final','IZ24\Final', 'IZ25\Final', 'IZ26\Final','IZ27\Final','IZ28\Final',...
        'IZ29\Final','IZ30\Final','IZ31\Final','IZ32\Saline','IZ33\Saline','IZ34\Saline'};  % To add: IZ23, IZ32, IZ33, IZ34
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ27\Final','IZ28\Final','IZ29\Final','IZ32\Final','IZ33\Final','IZ34\Final'}; 
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'CA3Saline') == 1
    mice = {'IZ27\Saline','IZ28\Saline','IZ29\Saline','IZ32\Saline','IZ33\Saline','IZ34\Saline'}; % To add: IZ32, IZ33, IZ34
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'mECBilateral') == 1 
    mice = {'IZ24\Final','IZ25\Final','IZ26\Final'};    %
    reg = {'contramEC','ipsimEC','Both'};
end


for rr = 1:3
    for cc = 1:2
        BehavData.data{rr,cc} = [];
        BehavData.mice{rr,cc} = [];
        BehavDataavg{rr,cc} = [];
    end
end

%% Loop through the mice
for m = 1:length(mice)
    
    pathdir = strcat(parentDir, mice{m});
    [~ ,BehavDataSumm] = SessBehavior('expPath',pathdir,'plotfig',false);
    for rr = 1:3
        for cc = 1:2
            BehavData.data{rr,cc} = [BehavData.data{rr,cc}; BehavDataSumm{rr,cc}];
            mouse(1:size(BehavDataSumm{rr,cc},1),1) = m;
            BehavData.mice{rr,cc} = [BehavData.mice{rr,cc}; mouse];
            clear mouse
        end
    end
end

%% Now plot the data

[colormap] = cbrewer('seq','YlGnBu',length(mice)+4);
reg = {'CA1','mEC','Both'};
zone = {'Stem','Return'};


if strcmp(tag,'CA3')== 1 || strcmp(tag,'CA3Saline') == 1
    figure
    set(gcf,'Renderer','painters')
    set(gcf,'Position',[10 10 500 500])
    if avgSessions
        for ii = 2:3
            for jj = 1:2
                subplot(3,2,2*(ii-1)+jj)             
                for mm = 1:length(mice)
                    idxMouse = BehavData.mice{ii,jj}==mm;
                    BehavDataavg{ii,jj} = [BehavDataavg{ii,jj}; nanmean(BehavData.data{ii,jj}(idxMouse,:),1)];             
                    %err = std(BehavData.data{ii,jj}(idxMouse,:),0,1)/sqrt(sum(idxMouse));
                    %errorbar(nanmean(BehavData.data{ii,jj}(idxMouse,:),1),err,'LineWidth',2,'Color',colormap(mm+4,:)); 
                    plot([1:1:2],nanmean(BehavData.data{ii,jj}(idxMouse,:),1),'LineWidth',2,'Color',colormap(mm+4,:));
                    hold on
                    mousePoint = find(idxMouse);
                    for kk = 1:length(mousePoint)
                        f1 = scatter([1:1:2],BehavData.data{ii,jj}(mousePoint(kk),:),'.','MarkerFaceColor',colormap(mm+4,:),'MarkerEdgeColor',colormap(mm+4,:));
                        f1.MarkerFaceAlpha = 0.5;
                    end
                    xlabel('Trial Blocks')
                    ylabel('Performance')
                    ylim([0.5 1.1])
                    xlim([0 3])
                    title(strcat('Region: ',reg{ii},' Zone:',zone{jj}));
                end
                err = std(BehavDataavg{ii,jj},0,1)/sqrt(length(mice));
                errorbar(nanmean(BehavDataavg{ii,jj},1),err,'LineWidth',2,'Color','k');               
                [pVal(ii,jj)] = signrank(BehavData.data{ii,jj}(:,1),BehavData.data{ii,jj}(:,2));
            end
        end
    else              
        data = [];
        dataid = [];
        col = [85/243 85/243 85/243;8/243 133/243 161/243;224/243 163/243 46/243;56/243 61/243 150/243];
        data = reshape(BehavData.data{2,1},[2*size(BehavData.data{2,1},1) 1]);
        dataid = [ones(size(BehavData.data{2,1},1),1)*1; ones(size(BehavData.data{2,1},1),1)*2];
        data = [data; reshape(BehavData.data{3,1},[2*size(BehavData.data{3,1},1) 1])];
        dataid = [dataid; ones(size(BehavData.data{3,1},1),1)*3; ones(size(BehavData.data{3,1},1),1)*4];       
        behavStats = groupStats(data,dataid,'inAxis',true,'color',col,'repeatedMeasures',true);
        hold on
        plot([1:1:4],[BehavData.data{2,1} BehavData.data{3,1}],'LineWidth',1,'Color',[0.5 0.5 0.5]);
        ylabel('Performance')
        %ylim([0.5 1])
        xlim([0 5])
        for kk = unique(dataid')
            [behavStats.signrank.p(kk),~, behavStats.signrank.stats{kk}] = signrank(data(dataid==kk),0.5);
        end
    end    
else
    figure
    set(gcf,'Renderer','painters')
    set(gcf,'Position',[10 10 500 1200])
    for ii = 1:3
        if ii == 1 
            col = [85/243 85/243 85/243;224/243 163/243 46/243];
        elseif ii == 2
            col = [85/243 85/243 85/243; 8/243 133/243 161/243];
        elseif ii ==3
            col = [85/243 85/243 85/243;56/243 61/243 150/243];
        end
        for jj = 1%:2
            subplot(3,2,2*(ii-1)+jj)  
            if avgSessions
                for mm = 1:length(mice)
                    idxMouse = BehavData.mice{ii,jj}==mm;
                    BehavDataavg{ii,jj} = [BehavDataavg{ii,jj}; nanmean(BehavData.data{ii,jj}(idxMouse,:),1)];             
                    %err = std(BehavData.data{ii,jj}(idxMouse,:),0,1)/sqrt(sum(idxMouse));
                    %errorbar(nanmean(BehavData.data{ii,jj}(idxMouse,:),1),err,'LineWidth',2,'Color',colormap(mm+4,:)); 
                    plot([1:1:2],nanmean(BehavData.data{ii,jj}(idxMouse,:),1),'LineWidth',2,'Color',colormap(mm+4,:));
                    hold on
                    mousePoint = find(idxMouse);
                    for kk = 1:length(mousePoint)
                        f1 = scatter([1:1:2],BehavData.data{ii,jj}(mousePoint(kk),:),'.','MarkerFaceColor',colormap(mm+4,:),'MarkerEdgeColor',colormap(mm+4,:));
                        f1.MarkerFaceAlpha = 0.5;
                    end
                    xlabel('Trial Blocks')
                    ylabel('Performance')
                    ylim([0.5 1.1])
                    xlim([0 3])
                    title(strcat('Region: ',reg{ii},' Zone:',zone{jj}));
                end
                err = std(BehavDataavg{ii,jj},0,1)/sqrt(length(mice));
                errorbar(nanmean(BehavDataavg{ii,jj},1),err,'LineWidth',2,'Color','k');               
                [pVal(ii,jj)] = signrank(BehavData.data{ii,jj}(:,1),BehavData.data{ii,jj}(:,2));
            else
                data = [];
                dataid = [];
                data = reshape(BehavData.data{ii,jj},[2*size(BehavData.data{ii,jj},1) 1]);
                dataid = [ones(size(BehavData.data{ii,jj},1),1)*1; ones(size(BehavData.data{ii,jj},1),1)*2];
                behavStats{ii,jj} = groupStats(data,dataid,'inAxis',true,'color',col,'repeatedMeasures',true,'labelSummary',false,'sigStar',false);
                hold on
                plot([1:1:2],BehavData.data{ii,jj},'LineWidth',1,'Color',[0.5 0.5 0.5]);
                ylabel('Performance')
                ylim([0.5 1])
                xlim([0 3])
                title(strcat('Region: ',reg{ii},' Zone:',zone{jj}));
                [behavStats{ii,jj}.signrank.pVal, ~ ,behavStats{ii,jj}.signrank.stats] = signrank(BehavData.data{ii,jj}(:,1),BehavData.data{ii,jj}(:,2));
                sigstar({[1,2]},[behavStats{ii,jj}.signrank.pVal])
            end
        end
    end
end

 saveas(gcf,strcat(parentDir,'\Compiled\Behavior\Behavior',tag,'.png'));
 saveas(gcf,strcat(parentDir,'\Compiled\Behavior\Behavior',tag,'.fig'));
 saveas(gcf,strcat(parentDir,'\Compiled\Behavior\Behavior',tag,'.eps'),'epsc');
 save(strcat(parentDir,'\Compiled\Behavior\Behavior',tag,'.mat'),'behavStats');
 
end