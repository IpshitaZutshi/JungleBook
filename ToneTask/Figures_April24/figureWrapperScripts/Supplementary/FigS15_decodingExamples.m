function FigS15_decodingExamples

plotAvgs = 1;
plotExamples = 0;

if plotAvgs
    fig2  = figure;
    set(fig2,'Renderer','painters')
    set(fig2,'Color','w')
    set(fig2,'Position',[462 491 1160 237]);
    
    numrows = 1;
    numcol = 5;
    
    %Dec = compileGoalDecodingStats_changePoint('plotfig',false);
    load('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\goalDecodingStats.mat')
    
    t = linspace(-2,2,121);
    col = [83/255 0/255 0/255;...
        184/255 15/255 10/255;...
        241/255 114/255 42/255;...
        249/255 197/255 81/255;...
        143/255 189/255 107/255;...
        87/255 116/255 144/255];
    
    subplot(numrows,numcol,1)
    Stats.propDecoded = groupStats([{Dec.propDecoded},{Dec.propDecodedShuff},{Dec.propDecodedError}, {Dec.propDecodedProbe}],[],'inAxis',true);
    
    subplot(numrows,numcol,2)
    Stats.timetoDecode = groupStats([{Dec.timetoGoal},{Dec.timetoGoalError}],[],'inAxis',true);
    
    plotAvgStd(Dec.distToGoal,numrows,numcol,3,fig2,t',[0 0 1], 0)
    title('distancetoGoal')
    
    plotAvgStd(Dec.acc,numrows,numcol,4,fig2,t',[0 0 1], 0)
    ylabel('acceleration')
    xlabel('Time around goal decoding')
    hold on
    line([-2 2],[0 0])
    line([0 0],[-20 10])
    
    subplot(numrows,numcol,5)
    data = [{Dec.propDecoded(Dec.mouseID==39)},{Dec.propDecoded(Dec.mouseID==40)},{Dec.propDecoded(Dec.mouseID==43)},{Dec.propDecoded(Dec.mouseID==44)},{Dec.propDecoded(Dec.mouseID==47)},{Dec.propDecoded(Dec.mouseID==48)}];
    Stats.mouseID = groupStats(data,[],'inAxis',true);

    expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
    saveas(gcf,strcat(expPath,'SupFigure15_DecodingAverages.png'));
    saveas(gcf,strcat(expPath,'SupFigure15_DecodingAverages.eps'),'epsc');
    saveas(gcf,strcat(expPath,'SupFigure15_DecodingAverages.fig'));
end

%% Plot examples of decoding
if plotExamples

    fig2  = figure;
    set(fig2,'Renderer','painters')
    set(fig2,'Color','w')
    set(fig2,'Position',[463 2 1204 745]);

    numrows = 3;
    numcol = 5;
    % 
    decodingPath = 'Z:\Homes\zz737\ipshita_data\Auditory_Task\';
    decodingName = 'py_data\theta_decoding_lickLoc_y\up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].nc';
    changePointName = 'py_data/theta_decoding_lickLoc_y/change_point_posterior_up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].mat';
    

    % IZ44 
    sess = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ44\Final\IZ44_220830_sess7';
    cd(sess);
    file = dir('*.Tracking.Behavior.mat');
    load(file.name)
    file = dir('*.TrialBehavior.Behavior.mat');
    load(file.name)

    % Load decoding data
    [posterior_goal, posterior_pos, post_time, post_pos, post_goal,change_point,trial] = loadDecodingData(decodingPath, 'IZ44\Final\IZ44_220830_sess7', decodingName,changePointName);

    % Load manifold data
    umap_path = strcat('Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ44_220830_sess7\manifold');
    cd(umap_path)    
    umap_name = 'behavior_speed_1_smooth_5_bin_0.1';
    behav_file = 'IZ44_220830_sess7.position_behavior_speed_1_smooth_5_bin_0.1.mat';
    A = 363.1123;
    E = 6.4;

    trialNum = [35 13 63 48 57];

    % Randomly find a trial to that port
    for ii = 1:5        
        %trialidx = find(behavTrials.linTrial==0 & behavTrials.lickLoc==ii);
        curTrial = trialNum(ii);%randsample(trialidx,1);
        tWin = behavTrials.timestamps(curTrial,:);
        idxVel = find(tracking.position.y(tracking.timestamps>tWin(1) & tracking.timestamps<tWin(2))>5,1,'first');
        startTime = tracking.timestamps(find((tracking.timestamps>tWin(1))==1,1,'first')+idxVel-1);  
        tWin(1) = startTime;

        success = plotDecoding(tWin, posterior_pos,posterior_goal,post_time,post_pos,post_goal,change_point,trial,numrows,numcol,2,ii,fig2,(curTrial-1),behavTrials.timestamps(curTrial,1),tracking);
        if ~success
            ii = ii-1;
        else
             manifoldPlot('figHandle',fig2,'umap_path',umap_path,'umap_name',umap_name,'behav_file',behav_file,'dotSize',15,'poscutOff',1,'speedThresh',1,...
    'numrow',numrows,'numcol',numcol,'rowloc',1,'colloc',ii,'TRIAL_TYPE', [0:5],'A',A,'E',E,'singleTrial',true,'tsWin',tWin)

            subplot(numrows,numcol, ii)
            title(strcat(num2str(behavTrials.lickLoc(curTrial)+1),',',num2str(behavTrials.correct(curTrial)),',',num2str(curTrial)));
        end
    end

    expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
    saveas(gcf,strcat(expPath,'SupFigure15_DecodingIZ44.png'));
    saveas(gcf,strcat(expPath,'SupFigure15_DecodingIZ44.eps'),'epsc');
    saveas(gcf,strcat(expPath,'SupFigure15_DecodingIZ44.fig'));

    % IZ47
    fig2  = figure;
    set(fig2,'Renderer','painters')
    set(fig2,'Color','w')
    set(fig2,'Position',[463 2 1204 745]);

    numrows = 3;
    numcol = 5;

    sess = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230707_sess24';
    cd(sess);
    file = dir('*.Tracking.Behavior.mat');
    load(file.name)
    file = dir('*.TrialBehavior.Behavior.mat');
    load(file.name)

    % Load decoding data
    [posterior_goal, posterior_pos, post_time, post_pos, post_goal,change_point,trial] = loadDecodingData(decodingPath, 'IZ47\Final\IZ47_230707_sess24', decodingName,changePointName);

    % Load manifold data
    umap_path = strcat('Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ47_230707_sess24\manifold');
    cd(umap_path)    
    umap_name = 'behavior_speed_1_smooth_5_bin_0.1';
    behav_file = 'IZ47_230707_sess24.position_behavior_speed_1_smooth_5_bin_0.1.mat';
    A = -74;
    E = 90;

    trialNum = [37 89 33 23 81];

    % Randomly find a trial to that port
    for ii = 1:5
        %trialidx = find(behavTrials.linTrial==0 & behavTrials.lickLoc==ii & behavTrials.probe==0);
        curTrial = trialNum(ii);%randsample(trialidx,1);
        tWin = behavTrials.timestamps(curTrial,:);
        idxVel = find(tracking.position.y(tracking.timestamps>tWin(1) & tracking.timestamps<tWin(2))>5,1,'first');
        startTime = tracking.timestamps(find((tracking.timestamps>tWin(1))==1,1,'first')+idxVel-1);  
        tWin(1) = startTime;        

        success = plotDecoding(tWin, posterior_pos,posterior_goal,post_time,post_pos,post_goal,change_point,trial,numrows,numcol,2,ii,fig2,(curTrial-1),behavTrials.timestamps(curTrial,1),tracking);
        if ~success
            ii = ii-1;
        else
            %plot manifold
            manifoldPlot('figHandle',fig2,'umap_path',umap_path,'umap_name',umap_name,'behav_file',behav_file,'dotSize',15,'poscutOff',1,'speedThresh',1,...
        'numrow',numrows,'numcol',numcol,'rowloc',1,'colloc',ii,'TRIAL_TYPE', [0:5],'A',A,'E',E,'singleTrial',true,'tsWin',tWin)

            subplot(numrows,numcol, ii)
            title(strcat(num2str(behavTrials.lickLoc(curTrial)+1),',',num2str(behavTrials.correct(curTrial)),',',num2str(curTrial)));
        end
    end

    expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
    saveas(gcf,strcat(expPath,'SupFigure15_DecodingIZ47.png'));
    saveas(gcf,strcat(expPath,'SupFigure15_DecodingIZ47.eps'),'epsc');
    saveas(gcf,strcat(expPath,'SupFigure15_DecodingIZ47.fig'));

    % IZ47 probe
    fig2  = figure;
    set(fig2,'Renderer','painters')
    set(fig2,'Color','w')
    set(fig2,'Position',[463 2 1204 745]);

    numrows = 3;
    numcol = 5;

    trialNum = [55 44 26 64 73];

    % Randomly find a trial to that port
    for ii = 1:5
        %trialidx = find(behavTrials.linTrial==0 & behavTrials.lickLoc==ii & behavTrials.probe==1);
        curTrial = trialNum(ii);%randsample(trialidx,1);
        tWin = behavTrials.timestamps(curTrial,:);
        idxVel = find(tracking.position.y(tracking.timestamps>tWin(1) & tracking.timestamps<tWin(2))>5,1,'first');
        startTime = tracking.timestamps(find((tracking.timestamps>tWin(1))==1,1,'first')+idxVel-1);  
        tWin(1) = startTime;

        success = plotDecoding(tWin, posterior_pos,posterior_goal,post_time,post_pos,post_goal,change_point,trial,numrows,numcol,2,ii,fig2,(curTrial-1),behavTrials.timestamps(curTrial,1),tracking);
        if ~success
            ii = ii-1;

        else
            %plot manifold
            manifoldPlot('figHandle',fig2,'umap_path',umap_path,'umap_name',umap_name,'behav_file',behav_file,'dotSize',15,'poscutOff',2,'speedThresh',2,...
        'numrow',numrows,'numcol',numcol,'rowloc',1,'colloc',ii,'TRIAL_TYPE', [0:5],'A',A,'E',E,'singleTrial',true,'tsWin',tWin)

            subplot(numrows,numcol, ii)
            title(strcat(num2str(behavTrials.lickLoc(curTrial)+1),',',num2str(behavTrials.correct(curTrial)),',',num2str(curTrial)));
        end
    end

    expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
    saveas(gcf,strcat(expPath,'SupFigure15_DecodingIZ47_probe.png'));
    saveas(gcf,strcat(expPath,'SupFigure15_DecodingIZ47_probe.eps'),'epsc');
    saveas(gcf,strcat(expPath,'SupFigure15_DecodingIZ47_probe.fig'));

    %% IZ48

    fig2  = figure;
    set(fig2,'Renderer','painters')
    set(fig2,'Color','w')
    set(fig2,'Position',[463 2 1204 745]);

    numrows = 3;
    numcol = 5;

    sess = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ48\Final\IZ48_230628_sess17';
    cd(sess);
    file = dir('*.Tracking.Behavior.mat');
    load(file.name)
    file = dir('*.TrialBehavior.Behavior.mat');
    load(file.name)

    % Load decoding data
    [posterior_goal, posterior_pos, post_time, post_pos, post_goal,change_point,trial] = loadDecodingData(decodingPath, 'IZ48\Final\IZ48_230628_sess17', decodingName,changePointName);

    % Load manifold data
    umap_path = strcat('Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\IZ48_230628_sess17\manifold');
    cd(umap_path)    
    umap_name = 'behavior_speed_1_smooth_5_bin_0.1';
    behav_file = 'IZ48_230628_sess17.position_behavior_speed_1_smooth_5_bin_0.1.mat';
    A = 179;
    E = 90;

    trialNum = [28 71 65 36 53];

    % Randomly find a trial to that port
    for ii = 1:5
        %trialidx = find(behavTrials.linTrial==0 & behavTrials.lickLoc==ii & behavTrials.probe==0);
        curTrial = trialNum(ii);%randsample(trialidx,1);
        tWin = behavTrials.timestamps(curTrial,:);
        idxVel = find(tracking.position.y(tracking.timestamps>tWin(1) & tracking.timestamps<tWin(2))>5,1,'first');
        startTime = tracking.timestamps(find((tracking.timestamps>tWin(1))==1,1,'first')+idxVel-1);  
        tWin(1) = startTime;        

        success = plotDecoding(tWin, posterior_pos,posterior_goal,post_time,post_pos,post_goal,change_point,trial,numrows,numcol,2,ii,fig2,(curTrial-1),behavTrials.timestamps(curTrial,1),tracking);
        if ~success
            ii = ii-1;
        else
            %plot manifold
            manifoldPlot('figHandle',fig2,'umap_path',umap_path,'umap_name',umap_name,'behav_file',behav_file,'dotSize',15,'poscutOff',1,'speedThresh',1,...
        'numrow',numrows,'numcol',numcol,'rowloc',1,'colloc',ii,'TRIAL_TYPE', [0:5],'A',A,'E',E,'singleTrial',true,'tsWin',tWin)

            subplot(numrows,numcol, ii)
            title(strcat(num2str(behavTrials.lickLoc(curTrial)+1),',',num2str(behavTrials.correct(curTrial)),',',num2str(curTrial)));
        end
    end

    expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
    saveas(gcf,strcat(expPath,'SupFigure15_DecodingIZ48.png'));
    saveas(gcf,strcat(expPath,'SupFigure15_DecodingIZ48.eps'),'epsc');
    saveas(gcf,strcat(expPath,'SupFigure15_DecodingIZ48.fig'));

    % IZ48 probe
    fig2  = figure;
    set(fig2,'Renderer','painters')
    set(fig2,'Color','w')
    set(fig2,'Position',[463 2 1204 745]);

    numrows = 3;
    numcol = 5;

    trialNum = [14 76 59 18 54];

    % Randomly find a trial to that port
    for ii = 1:5
        %trialidx = find(behavTrials.linTrial==0 & behavTrials.lickLoc==ii & behavTrials.probe==1);
        curTrial = trialNum(ii);%randsample(trialidx,1);
        tWin = behavTrials.timestamps(curTrial,:);
        idxVel = find(tracking.position.y(tracking.timestamps>tWin(1) & tracking.timestamps<tWin(2))>5,1,'first');
        startTime = tracking.timestamps(find((tracking.timestamps>tWin(1))==1,1,'first')+idxVel-1);  
        tWin(1) = startTime;

        success = plotDecoding(tWin, posterior_pos,posterior_goal,post_time,post_pos,post_goal,change_point,trial,numrows,numcol,2,ii,fig2,(curTrial-1),behavTrials.timestamps(curTrial,1),tracking);
        if ~success
            ii = ii-1;

        else
            %plot manifold
            manifoldPlot('figHandle',fig2,'umap_path',umap_path,'umap_name',umap_name,'behav_file',behav_file,'dotSize',15,'poscutOff',2,'speedThresh',2,...
        'numrow',numrows,'numcol',numcol,'rowloc',1,'colloc',ii,'TRIAL_TYPE', [0:5],'A',A,'E',E,'singleTrial',true,'tsWin',tWin)

            subplot(numrows,numcol, ii)
            title(strcat(num2str(behavTrials.lickLoc(curTrial)+1),',',num2str(behavTrials.correct(curTrial)),',',num2str(curTrial)));
        end
    end

    expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
    saveas(gcf,strcat(expPath,'SupFigure15_DecodingIZ48_probe.png'));
    saveas(gcf,strcat(expPath,'SupFigure15_DecodingIZ48_probe.eps'),'epsc');
    saveas(gcf,strcat(expPath,'SupFigure15_DecodingIZ48_probe.fig'));

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