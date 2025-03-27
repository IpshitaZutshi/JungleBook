function examinePokePatchPSTH

% Look at different types of pokes - 
%   1) Rewarded vs unrewarded trials
%   2) Pokes at each of the 7 different ports

sess= {'N7\N7_241216_sess23',...
    }; 

% sess= {'N11\Final\N11_250314_sess18',...
%     };

expPath = 'Z:\Buzsakilabspace\LabShare\ZutshiI\patchTask\';

mode = 7; % 0 = 'pyr_vs_int' or 7 = 'ports'

if mode == 7
    end_pt = 7;
else
    end_pt = 2;
end

for tt = 1:end_pt
    psthReward{tt} = [];
    psthReward_i{tt} = [];
    pearson_corr_per_cell{tt} = [];
end

for ii = 1:length(sess)
    %% Load files
    cd(strcat(expPath,sess{ii}))    
    file = dir(['*TrialBehavior.mat']);
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
        
        % calculate PSTH in response to each trial type, correct or incorrect

        if cellType == 1 % pyramidal cell     
            for tt = 1:end_pt

                if mode == 7
                    start_time = behavTrials.timestamps(96); %61
                    end_time = behavTrials.timestamps(146);%95
                    if tt == 1 
                        st = behavTrials.timestamps(behavTrials.port == 1 & ...
                            behavTrials.timestamps >= start_time & ...
                            behavTrials.timestamps <= end_time);
                    elseif tt == 2 
                        st = behavTrials.timestamps(behavTrials.port == 2 & ... % behavTrials.timestamps(behavTrials.reward_outcome == 1);
                            behavTrials.timestamps >= start_time & ...
                            behavTrials.timestamps <= end_time);
                    elseif tt == 3 
                        st = behavTrials.timestamps(behavTrials.port == 3 & ...
                            behavTrials.timestamps >= start_time & ...
                            behavTrials.timestamps <= end_time);
                    elseif tt == 4 
                        st = behavTrials.timestamps(behavTrials.port == 4 & ...
                            behavTrials.timestamps >= start_time & ...
                            behavTrials.timestamps <= end_time);
                    elseif tt == 5 
                        st = behavTrials.timestamps(behavTrials.port == 5 & ...
                            behavTrials.timestamps >= start_time & ...
                            behavTrials.timestamps <= end_time);
                    elseif tt == 6 
                        st = behavTrials.timestamps(behavTrials.port == 6 & ...
                            behavTrials.timestamps >= start_time & ...
                            behavTrials.timestamps <= end_time);
                    elseif tt == 7 
                        st = behavTrials.timestamps(behavTrials.port == 7 & ...
                            behavTrials.timestamps >= start_time & ...
                            behavTrials.timestamps <= end_time);
                    end
                elseif mode == 0
                    if tt == 1 
                        st = behavTrials.timestamps(behavTrials.reward_outcome == 0); % NOT SURE ABOUT 0 and 1
                    elseif tt == 2 
                        st = behavTrials.timestamps(behavTrials.reward_outcome == 1);                
                    end
                end
                    
                if ~isempty(st)
                    [stccg, tPSTH] = CCG({spikes.times{kk} st},[],'binSize',0.1,'duration',4,'norm','rate');
                    psthReward{tt} = [psthReward{tt}; stccg(:,2,1)']; 
                    %disp(strjoin(string(tt)))
                else
                    fillArr(1,1:41) = nan;
                    psthReward{tt} = [psthReward{tt}; fillArr];               
                end
            end

        elseif cellType == 0 % interneuron      
            if mode == 0
                for tt = 1:2

                    if tt ==1 
                        st_i = behavTrials.timestamps(behavTrials.reward_outcome == 0);
                    elseif tt == 2 
                        st_i = behavTrials.timestamps(behavTrials.reward_outcome == 1);                
                    end
                    
                    if ~isempty(st_i)
                        [stccg_i, tPSTH_i] = CCG({spikes.times{kk} st_i},[],'binSize',0.1,'duration',4,'norm','rate');
                        psthReward_i{tt} = [psthReward_i{tt}; stccg_i(:,2,1)'];               
                    else
                        fillArr(1,1:41) = nan;
                        psthReward_i{tt} = [psthReward_i{tt}; fillArr];               
                    end
                end
            else
                continue
            end
        end
    end
end


%% PLOT

figure
set(gcf,'Renderer','painters')
set(gcf,'Color','w')
YlGnBu=cbrewer('seq', 'YlGnBu', 11);
colormap(YlGnBu)


if mode == 0
    % pyramidal sorting
    idxT = tPSTH<0 & tPSTH>=-0.2;
    avgRate1 = nanmean(psthReward{1}(:,idxT),2);
    avgRate2 = nanmean(psthReward{2}(:,idxT),2);
    avgRate = nanmean([avgRate1 avgRate2],2);
    for tt = 1:2
        newpsth{tt} = psthReward{tt}(avgRate>0.5,:);
    end
    idxT = tPSTH<0 & tPSTH>=-0.5;
    [~,idxMax2] = max(newpsth{1}(:,idxT),[],2); % the 1 keeps sorting the same for all plots. change to tt for individually sorted
    [~,idxMax] = sort(idxMax2);

 
    % interneuron sorting
    idxT_i = tPSTH_i<0 & tPSTH_i>=-0.2;
    avgRate1_i = nanmean(psthReward_i{1}(:,idxT_i),2);
    avgRate2_i = nanmean(psthReward_i{2}(:,idxT_i),2);
    avgRate_i = nanmean([avgRate1_i avgRate2_i],2);
    for tt = 1:2
        newpsth_i{tt} = psthReward_i{tt}(avgRate_i>1,:);
    end
    idxT_i = tPSTH_i<0 & tPSTH_i>=-0.5;
    [~,idxMax2_i] = max(newpsth_i{1}(:,idxT_i),[],2);
    [~,idxMax_i] = sort(idxMax2_i);

    allData = [];
       % Combine pyramidal cells data
    for tt = 1:2
        allData = [allData; zscore(newpsth{tt}, [], 2)];
        allData = [allData; zscore(newpsth_i{tt}, [], 2)];
    end
    
    % Get global color scale limits
    collim = [nanmin(allData(:)), nanmax(allData(:))];


    % plot heatmaps 
    for tt = 1:2
        subplot(3,2,tt)
        temp = zscore(newpsth{tt},[],2);
        h = imagesc(tPSTH, 1:size(newpsth{tt},1),temp(idxMax,:));
        set(h, 'AlphaData', ~isnan(temp))
        %clim([-0.5 2])
        clim(collim);
        colorbar
        %ylim([125 350])
        title('Pyramidal Cells')
    
        subplot(3,2,tt+2)
        temp = zscore(newpsth_i{tt},[],2);
        h = imagesc(tPSTH_i, 1:size(newpsth_i{tt},1),temp(idxMax_i,:));
        set(h, 'AlphaData', ~isnan(temp))
        %clim([-0.5 2])
        clim(collim);
        colorbar
        %ylim([125 350])
        title('Interneurons')
            
        subplot(3,2,tt+4)
        % pyramidal cells
        col = [0.5 0.5 0.5];
        col_2 = [0.980392156862745 0.647058823529412 0.980392156862745];
        meanpsth = nanmean(newpsth{tt},1);
        stdpsth = nanstd(newpsth{tt},1)./sqrt(size(newpsth{tt},1));         
        hold on
        fill([tPSTH; flipud(tPSTH)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],col,'linestyle','none','FaceAlpha',0.5);                    
        hi = line(tPSTH,meanpsth,'LineWidth',1.5,'Color',col);
        % interneurons
        meanpsth_i = nanmean(newpsth_i{tt},1);
        stdpsth_i = nanstd(newpsth_i{tt},1)./sqrt(size(newpsth_i{tt},1));         
        fill([tPSTH_i; flipud(tPSTH_i)],[meanpsth_i'-stdpsth_i'; flipud(meanpsth_i'+stdpsth_i')],col_2,'linestyle','none','FaceAlpha',0.5);                    
        hi = line(tPSTH_i,meanpsth_i,'LineWidth',1.5,'Color',col_2);
        line([0 0],[1 20],'Color','r')
        %ylim([1 11])
        hold off
    end

elseif mode == 7
    idxT = tPSTH<0 & tPSTH>=-0.2;
    avgRate1 = nanmean(psthReward{1}(:,idxT),2);
    avgRate2 = nanmean(psthReward{2}(:,idxT),2);
    avgRate3 = nanmean(psthReward{3}(:,idxT),2);
    avgRate4 = nanmean(psthReward{4}(:,idxT),2);
    avgRate5 = nanmean(psthReward{5}(:,idxT),2);
    avgRate6 = nanmean(psthReward{6}(:,idxT),2);
    avgRate7 = nanmean(psthReward{7}(:,idxT),2);
    avgRate = nanmean([avgRate1 avgRate2 avgRate3 avgRate4 avgRate5 avgRate6 avgRate7],2);
    
    for tt = 1:7
        newpsth{tt} = psthReward{tt}(avgRate>0.5,:);
    end
    idxT = tPSTH<0 & tPSTH>=-0.5;
    [~,idxMax2] = max(newpsth{1}(:,idxT),[],2); % 1 for consistent sorting. tt for individually sorted
    [~,idxMax] = sort(idxMax2);

    allData = [];
       % Combine pyramidal cells data
    for tt = 1:7
        allData = [allData; zscore(newpsth{tt}, [], 2)];
    end
    
    % Get global color scale limits
    collim = [nanmin(allData(:)), nanmax(allData(:))];

    % plot heatmaps
    for tt = 1:7
        subplot(2,7,tt)
        temp = zscore(newpsth{tt},[],2);
        h = imagesc(tPSTH, 1:size(newpsth{tt},1),temp(idxMax,:));
        set(h, 'AlphaData', ~isnan(temp))
        %clim([-0.5 2])
        clim(collim);
        colorbar
        %ylim([125 350])
        title('Pyramidal Cells')
            
        subplot(2,7,tt+7)
        col = [0.047058823529412   0.419607843137255   0.121568627450980];
        meanpsth = nanmean(newpsth{tt},1);
        stdpsth = nanstd(newpsth{tt},1)./sqrt(size(newpsth{tt},1));         
        hold on
        fill([tPSTH; flipud(tPSTH)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],col,'linestyle','none','FaceAlpha',0.5); 
        hi = line(tPSTH,meanpsth,'LineWidth',1.5,'Color',col);
        line([0 0],[1 9],'Color','r')
        ylim([1 9])    
    end

end


%% heatmap corrolations

%pearson_corr = corr(newpsth{1}(:), newpsth{2}(:), 'Type', 'Pearson');

for kk = 1:7
    num_cells = size(newpsth{kk}, 1);
    pearson_corr_per_cell{kk} = corr(newpsth{1}(:), newpsth{kk}(:), 'Type', 'Pearson');
    %pearson_corr_per_cell{kk} = zeros (num_cells, 1);
    %for i = 1:num_cells
        %pearson_corr_per_cell{kk}(i) = corr(newpsth{1}(i, :)', newpsth{kk}(i, :)', 'Type', 'Pearson');
    %end
    % max(pearson_corr_per_cell);
    % mean(pearson_corr_per_cell(60:79));
    % nanmean(pearson_corr_per_cell);
end

figure
hold on
x = linspace(0, 10, 100);
for ii = 2:7
    m=pearson_corr_per_cell{ii};
    y = m*x;
    plot(x, y, 'LineWidth', 2)
    legend
end
hold off

end

