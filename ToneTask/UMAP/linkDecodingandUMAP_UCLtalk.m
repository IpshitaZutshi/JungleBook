function linkDecodingandUMAP_UCLtalk

fig2 = figure;
set(gcf,'Color','w')
set(fig2,'Position',[463 2 1204 745]);
set(gcf,'Renderer','painters')
    
numrows = 4;
numcol = 2;

%% Extract Decoding and UMAP for a single session

sess  = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230710_sess25';
%sess = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ48\Final\IZ48_230714_sess28';
cd(sess)

plotindividualTrials = 0;
speedThresh = 1;
decodingWin = 15;

decodingPath = 'Z:\Homes\zz737\ipshita_data\Auditory_Task\IZ47\Final\IZ47_230710_sess25\py_data\theta_decoding_lickLoc_y\up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].nc';
changePointName = 'Z:\Homes\zz737\ipshita_data\Auditory_Task\IZ47\Final\IZ47_230710_sess25\py_data\theta_decoding_lickLoc_y\change_point_posterior_up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].mat';

posterior_goal =  ncread(decodingPath,'x_position') ;
posterior_pos =  ncread(decodingPath,'y_position') ;
post_time =  ncread(decodingPath,'time') ;
post_pos =  ncread(decodingPath,'y_position_value') ;
post_goal =  ncread(decodingPath,'x_position_value') ;

load(changePointName);

umap_name = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ47_230710_sess25\manifold\Umap_behavior_speed_1_smooth_5_bin_0.1.csv';
behav_file = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ47_230710_sess25\manifold\IZ47_230710_sess25.position_behavior_speed_1_smooth_5_bin_0.1.mat';
% umap_name = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ48_230714_sess28\manifold\Umap_behavior_speed_1_smooth_5.csv';
% behav_file = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ48_230714_sess28\manifold\IZ48_230714_sess28.position_behavior_speed_1_smooth_5_bin_0.1.mat';
A = -98; E = 1.69;

file = dir('*.TrialBehavior.Behavior.mat');
load(file.name)
file = dir('*.Tracking.Behavior.mat');
load(file.name)

%% Within each trial, find the earliest time that the goal is decoded consistently for 250 ms - do this separately for probe and non-probe trials

for pp = 1%:2
    ts_dec = [];
    trial_dec = [];
    
    for tt = 1:length(behavTrials.lickLoc)
        if behavTrials.linTrial(tt)==0 && behavTrials.correct(tt) == 1%pp-1 

            [~,idxstart] = min(abs(post_time-behavTrials.timestamps(tt,1)));
            if post_time(idxstart)<behavTrials.timestamps(tt,1) %Take the next index
                idxstart = idxstart+1;
            end        
            [~,idxend] = min(abs(post_time-behavTrials.timestamps(tt,2)));
            if post_time(idxend)>behavTrials.timestamps(tt,2) %Take the previous index
                idxend = idxend-1;
            end   

            [~,decGoal] = max(posterior_goal(:,idxstart:idxend));

            %% Get last change point for that trial
            if sum(trial==(tt-1))>0
                curChanges = change_point{trial==(tt-1)};

                idxGoal = curChanges(end);
                trialDecGoal = mode(decGoal(curChanges(end)+1:end));

                if trialDecGoal==(behavTrials.lickLoc(tt)+1)             
                    ts_dec = [ts_dec post_time(idxGoal+idxstart)];
                    trial_dec = [trial_dec tt];
                end

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
        plot_ind =  [plot_ind,find(lick_loc_ds==TRIAL_TYPE(tt) & correct_ds==1 & position_y_all>2 & speed_all'>speedThresh)];   
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
    
    if pp == 1
        ax2 = subplot(numrows,numcol,numcol+1);
        scatter3(Umap_results(plot_ind,1),Umap_results(plot_ind,2),Umap_results(plot_ind,3),5,freq_plot,'filled','MarkerFaceAlpha',1);
        colormap(ax2,"viridis")
        view(A,E)
        grid off;
        axis off;
        axis tight
        clim([2000 23000])
        colorbar
        set(gca,'ColorScale','log')
    end
    
    col = [83/255 0/255 0/255;...
        184/255 15/255 10/255;...
        241/255 114/255 42/255;...
        249/255 197/255 81/255;...
        143/255 189/255 107/255;...
        87/255 116/255 144/255];
    
    ax2 = subplot(numrows,numcol,numcol*2+1);
    scatter3(Umap_results(plot_ind,1),Umap_results(plot_ind,2),Umap_results(plot_ind,3),5,lick_plot,'filled','MarkerFaceAlpha',1);
    colormap(ax2,col);
    view(A,E)
    grid off;
    axis off;
    axis tight    
    colorbar
    
    ax2 = subplot(numrows,numcol,1);
    scatter3(Umap_results(plot_ind,1),Umap_results(plot_ind,2),Umap_results(plot_ind,3),5,pos_plot,'filled','MarkerFaceAlpha',1);
    colormap(ax2,"jet");
    view(A,E)
    grid off;
    axis off;
    axis tight
    colorbar

    subplot(numrows,numcol,numcol*3+1)
    scatter3(Umap_results(plot_ind,1),Umap_results(plot_ind,2),Umap_results(plot_ind,3),5,[0.9 0.9 0.9],'filled','MarkerFaceAlpha',1);
    view(A,E)
    grid off;
    axis off;
    axis tight
    hold on   
    
    for ii = 1:length(ts_dec)
        [~,idxUMAP] = min(abs(timestamp_beh-ts_dec(ii)));
        colid = behavTrials.lickLoc(trial_dec(ii))+1;
        subplot(numrows,numcol,numcol*3+1)
        scatter3(Umap_results(idxUMAP-1,1),Umap_results(idxUMAP-1,2),Umap_results(idxUMAP-1,3),35,col(colid,:),'filled');
    end
    colorbar
end

%% Plot decoding

decodingPath = 'Z:\Homes\zz737\ipshita_data\Auditory_Task\';
decodingName = 'py_data\theta_decoding_lickLoc_y\up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].nc';
changePointName = 'py_data/theta_decoding_lickLoc_y/change_point_posterior_up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].mat';

sess = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230710_sess25';
cd(sess);
file = dir('*.Tracking.Behavior.mat');
load(file.name)
file = dir('*.TrialBehavior.Behavior.mat');
load(file.name)

% Load decoding data
[posterior_goal, posterior_pos, post_time, post_pos, post_goal,change_point,trial] = loadDecodingData(decodingPath, 'IZ47\Final\IZ47_230710_sess25', decodingName,changePointName);

% Load manifold data
umap_path = strcat('Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ47_230710_sess25\manifold');
cd(umap_path)    
umap_name = 'behavior_speed_1_smooth_5_bin_0.1';
behav_file = 'IZ47_230710_sess25.position_behavior_speed_1_smooth_5_bin_0.1.mat';
A = -98; E = 1.69;

lickport = [2 3];
% Randomly find a trial to that port
for ii = 1%:2
    trialidx = find(behavTrials.linTrial==0 & behavTrials.lickLoc==lickport(ii) & behavTrials.probe==0);
    %curTrial = trialNum(ii);%randsample(trialidx,1);
    curTrial = 60;%randsample(trialidx,1);
    tWin = behavTrials.timestamps(curTrial,:);
    idxVel = find(tracking.position.y(tracking.timestamps>tWin(1) & tracking.timestamps<tWin(2))>4,1,'first');
    startTime = tracking.timestamps(find((tracking.timestamps>tWin(1))==1,1,'first')+idxVel-1);  
    tWin(1) = startTime;        

    success = plotDecoding(tWin, posterior_pos,posterior_goal,post_time,post_pos,post_goal,change_point,trial,numrows,numcol,2,ii+1,fig2,(curTrial-1),behavTrials.timestamps(curTrial,1),tracking);
    if ~success
        ii = ii-1;
    else
        %plot manifold
        manifoldPlot('figHandle',fig2,'umap_path',umap_path,'umap_name',umap_name,'behav_file',behav_file,'dotSize',15,'poscutOff',1,'speedThresh',1,...
    'numrow',numrows,'numcol',numcol,'rowloc',1,'colloc',ii+1,'TRIAL_TYPE', [0:5],'A',A,'E',E,'singleTrial',true,'tsWin',tWin)

        subplot(numrows,numcol, ii)
        title(strcat(num2str(behavTrials.lickLoc(curTrial)+1),',',num2str(behavTrials.correct(curTrial)),',',num2str(curTrial)));
    end
end
end

function success = plotDecoding(tWin, posterior_pos,posterior_goal,post_time,post_pos,post_goal,change_point,trial,numrows,numcol,rowloc,colloc,fighandle,trialNum,trueTStart,tracking)

RdPu=cbrewer('seq', 'RdPu', 11);
Purples=cbrewer('seq', 'PuBu', 11);

[minTS,idxDec1] = min(abs(post_time-tWin(1)));
if minTS>0.5 %Take the next index
    idxDec1 = idxDec1+1;
end
[minTS,idxDec2] = min(abs(post_time-tWin(2)));
if minTS>0.5 %Take the previous index
    idxDec2 = idxDec2-1;
end

if any(diff(post_time(idxDec1:idxDec2))> 0.011)
    success = 0;
    return
else
    success = 1;
end

[~,idxstart] = min(abs(tracking.timestamps-tWin(1)));      
[~,idxend] = min(abs(tracking.timestamps-tWin(2)));

ax1 = subplot(numrows,numcol,numcol*(rowloc-1)+colloc,'Parent',fighandle);
imagesc(post_time(idxDec1:idxDec2),post_pos,posterior_pos(:,idxDec1:idxDec2));
set(gca,'YDir','normal')
colormap(ax1,RdPu)
hold on
plot(tracking.timestamps(idxstart:idxend),tracking.position.y(idxstart:idxend),'Color','b')
box off
clim([0 0.4])
xlim([tWin(1) tWin(2)])
set(gca,'YTickLabel',[]);
%set(gca,'XTickLabel',[]);
axis tight

ax1 = subplot(numrows,numcol,numcol*(rowloc)+colloc,'Parent',fighandle);
imagesc(post_time(idxDec1:idxDec2),1:6,posterior_goal(:,idxDec1:idxDec2));
set(gca,'YDir','normal')
colormap(ax1,Purples)
box off
hold on
set(gca,'YTickLabel',[]);
plotChangePoints(fighandle,numrows,numcol,numcol*(rowloc)+colloc,change_point,trial,tWin(1),tWin(2),post_time,trialNum,posterior_goal,trueTStart)
xlim([tWin(1) tWin(2)])

end

function [posterior_goal, posterior_pos, post_time, post_pos, post_goal,change_point,trial] = loadDecodingData(decodingPath, sess, decodingName,changePointName)
    % Load decoding data
    file_name = strcat(decodingPath, '\', sess, '\', decodingName);
    posterior_goal = ncread(file_name, 'x_position');
    posterior_pos = ncread(file_name, 'y_position');
    post_time = ncread(file_name, 'time');
    post_pos = ncread(file_name, 'y_position_value');
    post_goal = ncread(file_name, 'x_position_value');

    file_name = strcat(decodingPath, '\', sess, '\', changePointName);
    load(file_name)
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

function plotChangePoints(fighandle,numrow,numcol,plotloc,change_points,trial,tsStart,tsEnd, post_time,trialNum,posterior_goal,trueTStart)

[~,idxstart] = min(abs(post_time-trueTStart));
if post_time(idxstart)<tsStart %Take the next index
    idxstart = idxstart+1;
end        
[~,idxend] = min(abs(post_time-tsEnd));
if post_time(idxend)>tsEnd %Take the previous index
    idxend = idxend-1;
end   

[~,decGoal] = max(posterior_goal(:,idxstart:idxend));

[~, goal_dec] = max(posterior_goal, [], 1);

curChanges = change_points{trial==(trialNum)};
subplot(numrow,numcol,plotloc,'Parent',fighandle)
for cr = 1:length(curChanges)
    if cr == length(curChanges)
        trialDecGoal = mode(decGoal(curChanges(cr)+1:end));
        line([post_time(idxstart+curChanges(cr)-1) post_time(idxend)],[trialDecGoal trialDecGoal],'Color','m')
    elseif cr == 1
        trialDecGoal = mode(decGoal(curChanges(cr)+1:curChanges(cr+1)+1));
        line([post_time(idxstart) post_time(idxstart+curChanges(cr+1)-1)],[trialDecGoal trialDecGoal],'Color','m')
    else
        trialDecGoal = mode(decGoal(curChanges(cr)+1:curChanges(cr+1)+1));
        line([post_time(idxstart+curChanges(cr)-1) post_time(idxstart+curChanges(cr+1)-1)],[trialDecGoal trialDecGoal],'Color','m')
    end
    line([post_time(idxstart+curChanges(cr)-1) post_time(idxstart+curChanges(cr)-1)],[0 6],'Color','r','LineWidth',1.5);  
end
xlim([tsStart tsEnd])
end