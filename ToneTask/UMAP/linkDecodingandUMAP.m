function linkDecodingandUMAP

fig2 = figure;
set(gcf,'Color','w')
set(gcf,'Position',[88 4 1700 927])
set(gcf,'Renderer','painters')
    
%% Extract Decoding and UMAP for a single session

sess  = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230626_sess15';
%sess = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ48\Final\IZ48_230714_sess28';
cd(sess)

plotindividualTrials = 0;

posterior_goal =  ncread('causal_posterior_lickLoc_y.nc','x_position') ;
posterior_pos =  ncread('causal_posterior_lickLoc_y.nc','y_position') ;
post_time =  ncread('causal_posterior_lickLoc_y.nc','time') ;
post_pos =  ncread('causal_posterior_lickLoc_y.nc','y_position_value') ;
post_goal =  ncread('causal_posterior_lickLoc_y.nc','x_position_value') ;
% 
% posterior_goal =  ncread('IZ48_230714_sess28.causal_posterior_lickLoc_y.nc','x_position') ;
% posterior_pos =  ncread('IZ48_230714_sess28.causal_posterior_lickLoc_y.nc','y_position') ;
% post_time =  ncread('IZ48_230714_sess28.causal_posterior_lickLoc_y.nc','time') ;
% post_pos =  ncread('IZ48_230714_sess28.causal_posterior_lickLoc_y.nc','y_position_value') ;
% post_goal =  ncread('IZ48_230714_sess28.causal_posterior_lickLoc_y.nc','x_position_value') ;

umap_name = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ47_230626_sess15\manifold\Umap_behavior_speed_1_smooth_5_bin_0.1.csv';
behav_file = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ47_230626_sess15\manifold\IZ47_230626_sess15.position_behavior_speed_1_smooth_5_bin_0.1.mat';
% umap_name = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ48_230714_sess28\manifold\Umap_behavior_speed_1_smooth_5.csv';
% behav_file = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ48_230714_sess28\manifold\IZ48_230714_sess28.position_behavior_speed_1_smooth_5_bin_0.1.mat';
A = -1.17; E = -31.94;
% A = 23;
% E = 1.55;

%% Only extract goal decoding that is consistent for 250 ms
[~,goal_dec] = max(posterior_goal,[],1);
result = NaN(size(goal_dec));  % Initialize result with NaNs
last_value = nan;
count = 0;

for i = 1:length(goal_dec)
    if goal_dec(i) ~= last_value
        if count >= 10
            result(i-count:i-1) = last_value;  % Mark previous values as last_value
        end
        last_value = goal_dec(i);
        count = 1;
    else
        count = count + 1;
    end
end
if count >= 10
    result(end-count+1:end) = last_value;  % Mark last batch as last_value
end

file = dir('*.TrialBehavior.Behavior.mat');
load(file.name)
file = dir('*.Tracking.Behavior.mat');
load(file.name)

%% Within each trial, find the earliest time that the goal is decoded consistently for 250 ms - do this separately for probe and non-probe trials

for pp = 1:2
    ts_dec = [];
    trial_dec = [];
    
    for tt = 1:length(behavTrials.lickLoc)
        if behavTrials.linTrial(tt)==0 && behavTrials.probe(tt) == pp-1 
            [~,idxstart] = min(abs(post_time-behavTrials.timestamps(tt,1)));
            if post_time(idxstart)<behavTrials.timestamps(tt,1) %Take the next index
                idxstart = idxstart+1;
            end        
            [~,idxend] = min(abs(post_time-behavTrials.timestamps(tt,2)));
            if post_time(idxend)>behavTrials.timestamps(tt,2) %Take the previous index
                idxend = idxend-1;
            end   
            idxGoal = find(result(idxstart:idxend)==(behavTrials.lickLoc(tt)+1),1,'first');
            if ~isempty(idxGoal)
                ts_dec = [ts_dec post_time(idxGoal+idxstart-1)];
                trial_dec = [trial_dec tt];
            end
        end
    end
    
    %% load Umap result
    Umap_results = readtable(umap_name);
    Umap_results = table2array(Umap_results);
    
    %% load position direction and other information
    load(behav_file);
    TRIAL_TYPE = 0:5;
    plot_ind = [];
    
    for tt = 1:length(TRIAL_TYPE)
        plot_ind =  [plot_ind,find(lick_loc_ds==TRIAL_TYPE(tt) & probe_ds==pp-1 & position_y_all>4 & speed_all'>2)];   
    end
    
    gain = [122/9, 122/32 122/55.53 122/79.62 122/102.79 122/122];
    freqExp = log10(22000/2000);
    
    for tt = 1:length(position_y_all)
        if trial_type_ds(tt)<6
            freq = (position_y_all(tt)*gain(trial_type_ds(tt)+1))/122;
            tonepos_all(tt) = 2000*(10.^(freqExp*freq));
        else
            tonepos_all(tt) = 0;
        end
    end
    
    pos_plot = position_y_all;
    pos_plot(isnan(pos_plot))=0; % deal with nan
    pos_plot = pos_plot(plot_ind);
    
    %tonepos_all(tonepos_all>25000) = 0;
    freq_plot = tonepos_all;
    freq_plot(isnan(freq_plot))=0; % deal with nan
    freq_plot = freq_plot(plot_ind);
    
    arm_plot = trial_type_ds;
    arm_plot(isnan(arm_plot))=0; % deal with nan
    arm_plot = arm_plot(plot_ind);
    
    lick_plot = lick_loc_ds; 
    lick_plot(isnan(lick_plot))=0; % deal with nan
    lick_plot = lick_plot(plot_ind);
    
    % ax2 = subplot(2,4,pp);
    % scatter3(Umap_results(plot_ind,1),Umap_results(plot_ind,2),Umap_results(plot_ind,3),5,freq_plot,'filled','MarkerFaceAlpha',1);
    % colormap(ax2,"viridis")
    % view(A,E)
    % grid off;
    % axis off;
    % axis tight
    % clim([2000 23000])
    % colorbar
    % set(gca,'ColorScale','log')
    
    col = [83/255 0/255 0/255;...
        184/255 15/255 10/255;...
        241/255 114/255 42/255;...
        249/255 197/255 81/255;...
        143/255 189/255 107/255;...
        87/255 116/255 144/255];
    
    ax2 = subplot(2,4,pp);
    scatter3(Umap_results(plot_ind,1),Umap_results(plot_ind,2),Umap_results(plot_ind,3),5,lick_plot,'filled','MarkerFaceAlpha',1);
    colormap(ax2,col);
    view(A,E)
    grid off;
    axis off;
    axis tight
    %colorbar
    
    subplot(2,4,pp+2)
    scatter3(Umap_results(plot_ind,1),Umap_results(plot_ind,2),Umap_results(plot_ind,3),5,[0.9 0.9 0.9],'filled','MarkerFaceAlpha',1);
    view(A,E)
    grid off;
    axis off;
    axis tight
    hold on   
    
    for ii = 1:length(ts_dec)
        [~,idxUMAP] = min(abs(timestamp_beh-ts_dec(ii)));
        colid = behavTrials.lickLoc(trial_dec(ii))+1;
        subplot(2,4,pp+2)
        scatter3(Umap_results(idxUMAP,1),Umap_results(idxUMAP,2),Umap_results(idxUMAP,3),35,col(colid,:),'filled');
    end
    %colorbar
end

%% Get decoding stats across mice
Dec = compileGoalDecodingStats;

t = linspace(-2,2,121);

subplot(2,4,5)
nhist(Dec.freqtoGoal,'samebins','binfactor',3,'color','sequential','proportion','smooth')
xscale log
xlim([2000 23000])

% subplot(2,4,4)
% fractDec = Dec.countDec./Dec.totalDec;
% boxplot(fractDec)
% 
% plotAvgStd(Dec.nose,2,4,5,fig2,t',[0 0 1], 0)
% title('Nose position')

plotAvgStd(Dec.distToGoal,2,4,6,fig2,t',[0 0 1], 0)
title('distancetoGoal')

plotAvgStd(Dec.speed,2,4,7,fig2,t',[0 0 1], 0)
title('speed')

plotAvgStd(Dec.hd,2,4,8,fig2,t',[0 0 1], 0)
title('headdirection')


if plotindividualTrials
    for ii = 1:length(ts_dec)
        figure
        set(gcf,'Color','w')
    
        subplot(1,4,1)
        scatter3(Umap_results(plot_ind,1),Umap_results(plot_ind,2),Umap_results(plot_ind,3),5,freq_plot,'filled','MarkerFaceAlpha',1);
        colormap jet
        view(A,E)
        grid off;
        axis off;
        axis tight
        clim([2000 23000])
        
        subplot(1,4,2)
        scatter3(Umap_results(plot_ind,1),Umap_results(plot_ind,2),Umap_results(plot_ind,3),5,[0.9 0.9 0.9],'filled','MarkerFaceAlpha',1);
        view(A,E)
        grid off;
        axis off;
        axis tight
        hold on   
    
        [~,idxUMAP1] = min(abs(timestamp_beh-behavTrials.timestamps(trial_dec(ii),1)));
        [~,idxUMAP2] = min(abs(timestamp_beh-behavTrials.timestamps(trial_dec(ii),2)));
        subplot(1,4,2)
        tsLine = idxUMAP1:idxUMAP2;
        scatter3(Umap_results(idxUMAP1:idxUMAP2,1),Umap_results(idxUMAP1:idxUMAP2,2),Umap_results(idxUMAP1:idxUMAP2,3),10,1:length(tsLine),'filled');
        colormap jet
    
        [~,idxUMAP] = min(abs(timestamp_beh-ts_dec(ii)));
        scatter3(Umap_results(idxUMAP,1),Umap_results(idxUMAP,2),Umap_results(idxUMAP,3),50,'m','filled');
    
        subplot(1,4,3)
        [minTS,idxDec1] = min(abs(post_time-behavTrials.timestamps(trial_dec(ii),1)));
        if post_time(idxDec1)<behavTrials.timestamps(trial_dec(ii),1) %Take the next index %Take the next index
            idxDec1 = idxDec1+1;
        end
        [minTS,idxDec2] = min(abs(post_time-behavTrials.timestamps(trial_dec(ii),2)));
        if post_time(idxDec2)>behavTrials.timestamps(trial_dec(ii),2) %Take the previous index
            idxDec2 = idxDec2-1;
        end
        imagesc(post_time(idxDec1:idxDec2),1:6,posterior_goal(:,idxDec1:idxDec2));
        set(gca,'YDir','normal')
        title(num2str(behavTrials.lickLoc(trial_dec(ii))+1))
        xlim([post_time(idxDec1) post_time(idxDec2)])
    
    
        subplot(1,4,4)
        plot(post_time(idxDec1:idxDec2),goal_dec(idxDec1:idxDec2))
        ylim([0.5 6.5])
        hold on
        plot(post_time(idxDec1:idxDec2),result(idxDec1:idxDec2))
        title(num2str(trial_dec(ii)))
        xlim([post_time(idxDec1) post_time(idxDec2)])
    
    end
end

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\MainFigures\';
saveas(gcf,strcat(expPath,'Figure4B_decodingUMAP.png'));
saveas(gcf,strcat(expPath,'Figure4B_decodingUMAP.eps'),'epsc');
saveas(gcf,strcat(expPath,'Figure4B_decodingUMAP.fig'));

end

function plotAvgStd(array,numrows,numcol,subplotlocation,figureHandle,xAxis,col, useMedian)

    subplot(numrows, numcol, subplotlocation, 'Parent', figureHandle);
    if ~useMedian
        meanpsth = nanmean(array,1);
        stdpsth = nanstd(array,1)./sqrt(size(array,1));
        lArr  = meanpsth-stdpsth;
        uArr = meanpsth+stdpsth;
    else
        meanpsth = nanmedian(array,1);
        a = quantile(array,4,1);
        lArr  = meanpsth-a(2,:);
        uArr = meanpsth+a(3,:);
    end
    fill([xAxis; flipud(xAxis)],[lArr'; flipud(uArr')],col,'linestyle','none','FaceAlpha',0.5);                    
    hold on
    hi = line(xAxis,meanpsth,'LineWidth',1,'Color',col);

end
