function pass_noPassExamples_mPFC

%% Plot examples of decoding
numrows = 2;
numcol = 5;
% 
mouseID = 'IZ40';
sessName = 'IZ40_220707_sess16';
decodingPath = 'Z:\Homes\zz737\ipshita_data\Auditory_Task\';
decodingName = 'py_data\theta_decoding_lickLoc_y\up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].nc';
changePointName = 'py_data/theta_decoding_lickLoc_y/change_point_posterior_up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].mat';


% IZ44 
sess = strcat('Z:\Homes\zutshi01\Recordings\Auditory_Task\',mouseID,'\Final\',sessName);
expPath = strcat('Z:\Homes\zutshi01\Recordings\Auditory_Task\mPFC\',mouseID,'\',sessName,'\decoding\');
cd(sess);
file = dir('*.Tracking.Behavior.mat');
load(file.name)
file = dir('*.TrialBehavior.Behavior.mat');
load(file.name)
file = dir('*.thetaLFP.mat');
load(file.name)
Dec = findDecelerationPoints('plotfig',false);
delibTS = Dec.ts(Dec.decType==2);

% Load decoding data
[posterior_goal, posterior_pos, post_time, post_pos, post_goal,change_point,trial] = loadDecodingData(decodingPath, strcat(mouseID,'\Final\',sessName), decodingName,changePointName);

trialNum = [10:5:85];

% Randomly find a trial to that port
for ii = 1:length(trialNum)     
    fig2  = figure;
    set(fig2,'Renderer','painters')
    set(fig2,'Color','w')
    set(fig2,'Position',[463 2 1204 250]);
    
    for jj = 1:5
        curTrial = trialNum(ii)+jj;%randsample(trialidx,1);
        tWin = behavTrials.timestamps(curTrial,:);
        idxVel = find(tracking.position.y(tracking.timestamps>tWin(1) & tracking.timestamps<tWin(2))>10,1,'first');
        startTime = tracking.timestamps(find((tracking.timestamps>tWin(1))==1,1,'first')+idxVel-1);  
        tWin(1) = startTime;
    
        plotDecoding(tWin, posterior_pos,posterior_goal,post_time,post_pos,post_goal,change_point,trial,numrows,numcol,1,jj,fig2,(curTrial-1),behavTrials.timestamps(curTrial,1),tracking);
        plotTheta(tWin,lfp,posterior_pos,post_time,post_pos,numrows,numcol,2,jj,fig2,tracking);
        % If there's a deliberation deceleration, show that
        idxDelib = InIntervals(delibTS,tWin);
        if sum(idxDelib)>0
            line([delibTS(idxDelib) delibTS(idxDelib)],[min(post_pos) max(post_pos)])
        end    
    end    
    if ~exist(expPath)
        mkdir(expPath)
    end
    saveas(gcf,strcat(expPath,'passNoPass_',num2str(trialNum(ii)),'.png'));
    saveas(gcf,strcat(expPath,'passNoPass_',num2str(trialNum(ii)),'.eps'),'epsc');
    saveas(gcf,strcat(expPath,'passNoPass_',num2str(trialNum(ii)),'.fig'));
    close all
end
% 


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

% Extract maxLocation
[~,idxPos] = max(posterior_pos(:,idxDec1:idxDec2),[],1);
decPos = post_pos(idxPos);

ax1 = subplot(numrows,numcol,numcol*(rowloc-1)+colloc,'Parent',fighandle);
imagesc(post_time(idxDec1:idxDec2),post_pos,posterior_pos(:,idxDec1:idxDec2));
set(gca,'YDir','normal')
colormap(ax1,RdPu)
hold on
plot(tracking.timestamps(idxstart:idxend),tracking.position.y(idxstart:idxend),'Color','b')
plot(post_time(idxDec1:idxDec2),decPos,'Color','k')
box off
clim([0 0.3])
xlim([tWin(1) tWin(2)])
set(gca,'YTickLabel',[]);
%set(gca,'XTickLabel',[]);
axis tight

line([tWin(1) tWin(2)],[18 18],'Color',[0.2 0.2 0.2],'LineWidth',0.75,'LineStyle','--')   
hold on
line([tWin(1) tWin(2)],[35 35],'Color',[0.2 0.2 0.2],'LineWidth',0.75,'LineStyle','--')        
line([tWin(1) tWin(2)],[61 61],'Color',[0.2 0.2 0.2],'LineWidth',0.75,'LineStyle','--')
line([tWin(1) tWin(2)],[82 82],'Color',[0.2 0.2 0.2],'LineWidth',0.75,'LineStyle','--')
line([tWin(1) tWin(2)],[108 108],'Color',[0.2 0.2 0.2],'LineWidth',0.75,'LineStyle','--')
line([tWin(1) tWin(2)],[122 122],'Color',[0.2 0.2 0.2],'LineWidth',0.75,'LineStyle','--')

end

function success = plotTheta(tWin,lfp,posterior_pos,post_time,post_pos,numrows,numcol,rowloc,colloc,fighandle,tracking)

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

[~,idxstartLFP] = min(abs(lfp.timestamps-tWin(1)));      
[~,idxendLFP] = min(abs(lfp.timestamps-tWin(2)));

% Extract maxLocation
[~,idxPos] = max(posterior_pos(:,idxDec1:idxDec2),[],1);
decPos = post_pos(idxPos);

ax1 = subplot(numrows,numcol,numcol*(rowloc-1)+colloc,'Parent',fighandle);
hold on
plot(tracking.timestamps(idxstart:idxend),tracking.position.y(idxstart:idxend),'Color','b')
plot(post_time(idxDec1:idxDec2),decPos,'Color','k')
plot(lfp.timestamps(idxstartLFP:idxendLFP),(lfp.filtered(idxstartLFP:idxendLFP)./max(lfp.filtered))*120,'Color','r')
box off
xlim([tWin(1) tWin(2)])
set(gca,'YTickLabel',[]);
%set(gca,'XTickLabel',[]);
axis tight
ylim([-20 125])
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



