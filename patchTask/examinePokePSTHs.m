function examinePokePSTHsCorrectIncorrect

% Look at 2 types of pokes - 
%   1) Return run, in any of the middle ports while facing forward
%   2) Return run, in any of the middle ports while facing backwards

sess= {'N7\N7_241216_sess23',...
    }; 


expPath = 'Z:\Buzsakilabspace\LabShare\ZutshiI\patchTask\';

mode = 0; % 0 = 'pyr_vs_int' or 7 = 'ports'

for tt = 1:7
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
    
    if mode == 0
        if size(digitalIn.ints{4},1)==2
            tsLick = [digitalIn.ints{4}';digitalIn.ints{5}';digitalIn.ints{6}';...
                digitalIn.ints{7}';digitalIn.ints{8}'];        
        else
            tsLick = [digitalIn.ints{4};digitalIn.ints{5};digitalIn.ints{6};...
                digitalIn.ints{7};digitalIn.ints{8}];        
        end
    
        intsPeriods =[];
        intsPeriods(1,1) = tsLick(1,1); % find stimulation intervals
        intPeaks =find(diff(tsLick(:,1))>1);
        for jj = 1:length(intPeaks)
            intsPeriods(jj,2) = tsLick(intPeaks(jj),2);
            intsPeriods(jj+1,1) = tsLick(intPeaks(jj)+1,1);
        end
        intsPeriods(end,2) = tsLick(end,2);  
        tsLick2 = intsPeriods;
    end

    for kk=1:length(cell_metrics.UID)
    %% Extract pyramidal cells and calculate rate map correlations and information 
        if strcmp(cell_metrics.putativeCellType{kk},'Pyramidal Cell') == 1 
            cellType = 1;
        else
            cellType = 0;
        end
        
 
        % If its a trajectory cell, calculate its PSTH in response to each
        % trial type, correct, incorrect, return

        if cellType == 1 % pyramidal cell      && (toneField == 1) && (toneCorr > 0.1) && idxMax>40
            for tt = 1:7

                if tt == 1 
                    st = behavTrials.timestamps(behavTrials.port == 1);
                elseif tt == 2 
                    st = behavTrials.timestamps(behavTrials.port == 2); % behavTrials.timestamps(behavTrials.reward_outcome == 1);
                elseif tt == 3 
                    st = behavTrials.timestamps(behavTrials.port == 3);
                elseif tt == 4 
                    st = behavTrials.timestamps(behavTrials.port == 4);
                elseif tt == 5 
                    st = behavTrials.timestamps(behavTrials.port == 5);
                elseif tt == 6 
                    st = behavTrials.timestamps(behavTrials.port == 6);
                elseif tt == 7 
                    st = behavTrials.timestamps(behavTrials.port == 7);
                end
                    
                if ~isempty(st)
                    [stccg, tPSTH] = CCG({spikes.times{kk} st},[],'binSize',0.1,'duration',4,'norm','rate');
                    psthReward{tt} = [psthReward{tt}; stccg(:,2,1)'];               
                else
                    fillArr(1,1:41) = nan;
                    psthReward{tt} = [psthReward{tt}; fillArr];               
                end
            end

        end
%{
        if cellType == 0 % interneuron      && (toneField == 1) && (toneCorr > 0.1) && idxMax>40
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

        end
      %}
    end
end


figure
set(gcf,'Renderer','painters')
set(gcf,'Color','w')
YlGnBu=cbrewer('seq', 'YlGnBu', 11);
colormap(YlGnBu)


idxT = tPSTH<0 & tPSTH>=-0.2;
avgRate1 = nanmean(psthReward{1}(:,idxT),2);
avgRate2 = nanmean(psthReward{2}(:,idxT),2);
avgRate = nanmean([avgRate1 avgRate2],2);
for tt = 1:7
    newpsth{tt} = psthReward{tt}(avgRate>0.5,:);
end
idxT = tPSTH<0 & tPSTH>=-0.5;
[~,idxMax2] = max(newpsth{1}(:,idxT),[],2); % the 1 keeps sorting the same for all plots
[~,idxMax] = sort(idxMax2);


% interneuron
% idxT_i = tPSTH_i<0 & tPSTH_i>=-0.2;
% avgRate1_i = nanmean(psthReward_i{1}(:,idxT_i),2);
% avgRate2_i = nanmean(psthReward_i{2}(:,idxT_i),2);
% avgRate_i = nanmean([avgRate1_i avgRate2_i],2);
% for tt = 1:2
%     newpsth_i{tt} = psthReward_i{tt}(avgRate_i>1,:);
% end
% idxT_i = tPSTH_i<0 & tPSTH_i>=-0.5;
% [~,idxMax2_i] = max(newpsth_i{1}(:,idxT_i),[],2);
% [~,idxMax_i] = sort(idxMax2_i);

%{
for tt = 1:2
    subplot(3,2,tt)
    temp = zscore(newpsth{tt},[],2);
    h = imagesc(tPSTH, 1:size(newpsth{tt},1),temp(idxMax,:));
    set(h, 'AlphaData', ~isnan(temp))
    %clim([-0.5 2])
    colorbar
    %ylim([125 350])
    title('Pyramidal Cells')

    subplot(3,2,tt+2)
    temp = zscore(newpsth_i{tt},[],2);
    h = imagesc(tPSTH_i, 1:size(newpsth_i{tt},1),temp(idxMax_i,:));
    set(h, 'AlphaData', ~isnan(temp))
    %clim([-0.5 2])
    colorbar
    %ylim([125 350])
    title('Interneurons')
        
    subplot(3,2,tt+4)
    col = [0.5 0.5 0.5];
    col_2 = [0.980392156862745 0.647058823529412 0.980392156862745];
    meanpsth = nanmean(newpsth{tt},1);
    stdpsth = nanstd(newpsth{tt},1)./sqrt(size(newpsth{tt},1));         
    hold on
    fill([tPSTH; flipud(tPSTH)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],col,'linestyle','none','FaceAlpha',0.5);                    
    hold on
    hi = line(tPSTH,meanpsth,'LineWidth',1.5,'Color',col);

    % interneurons
    meanpsth_i = nanmean(newpsth_i{tt},1);
    stdpsth_i = nanstd(newpsth_i{tt},1)./sqrt(size(newpsth_i{tt},1));         
    hold on
    fill([tPSTH_i; flipud(tPSTH_i)],[meanpsth_i'-stdpsth_i'; flipud(meanpsth_i'+stdpsth_i')],col_2,'linestyle','none','FaceAlpha',0.5);                    
    hold on
    hi = line(tPSTH_i,meanpsth_i,'LineWidth',1.5,'Color',col_2);

    line([0 0],[1 20],'Color','r')
    %ylim([1 11])    
end
%}
%{
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
%}

for tt = 1:7
%{
    idxT = tPSTH<0 & tPSTH>=-0.5;
    [~,idxMax2] = max(newpsth{tt}(:,idxT),[],2);
    [~,idxMax] = sort(idxMax2);
%}
    subplot(2,7,tt)
    temp = zscore(newpsth{tt},[],2);
    h = imagesc(tPSTH, 1:size(newpsth{tt},1),temp(idxMax,:));
    set(h, 'AlphaData', ~isnan(temp))
    %clim([-0.5 2])
    colorbar
    %ylim([125 350])
    title('Pyramidal Cells')
        
    subplot(2,7,tt+7)
    col = [0.5 0.5 0.5];
    meanpsth = nanmean(newpsth{tt},1);
    stdpsth = nanstd(newpsth{tt},1)./sqrt(size(newpsth{tt},1));         
    hold on
    fill([tPSTH; flipud(tPSTH)],[meanpsth'-stdpsth'; flipud(meanpsth'+stdpsth')],col,'linestyle','none','FaceAlpha',0.5);                    
    hold on
    hi = line(tPSTH,meanpsth,'LineWidth',1.5,'Color',col);

    line([0 0],[1 9],'Color','r')
    ylim([1 9])    
end

end

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
    plot(x, y)
    legend
end
hold off

