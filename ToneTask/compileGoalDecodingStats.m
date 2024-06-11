function Dec = compileGoalDecodingStats;
% 
% sess = {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
%     'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
%     'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   23    
%     'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
%     'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
%     'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
%     'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    36
%     'IZ44\Final\IZ44_220827_sess4', 'IZ44\Final\IZ44_220828_sess5',...
%     'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...
%     'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
%     'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',... 45
%     'IZ47\Final\IZ47_230626_sess15','IZ47\Final\IZ47_230707_sess24',...
%     'IZ47\Final\IZ47_230710_sess25','IZ47\Final\IZ47_230712_sess27',...49
%     'IZ48\Final\IZ48_230628_sess17','IZ48\Final\IZ48_230703_sess21',...
%     'IZ48\Final\IZ48_230705_sess22','IZ48\Final\IZ48_230714_sess28'};    

sess = {'IZ47\Final\IZ47_230626_sess15','IZ48\Final\IZ48_230714_sess28'};
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';

decodingWin = 10; %500 ms
tWin = 60;
gain = [122/9, 122/32 122/55.53 122/79.62 122/102.79 122/122];
freqExp = log10(22000/2000);
plotfig = 0;

Dec.speed = [];
Dec.hd = [];
Dec.nose = [];
Dec.distToGoal = [];

portlabel = {'port1','port2','port3','port4','port5','port6'};

for ff = 1:6
    Dec.freqtoGoal.(portlabel{ff}) = [];
    Dec.portdistToGoal.(portlabel{ff}) = [];
end

for ss = 1:length(sess)

    cd(strcat(expPath,sess{ss}))   

    file = dir(['*.Tracking.Behavior.mat']);
    load(file(1).name);
    file = dir(['*TrialBehavior.Behavior.mat']);
    load(file(1).name);
    [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts', true);

    %% Load decoding
    file = dir(['*causal_posterior_lickLoc_y.nc']);
    posterior_goal =  ncread(file.name,'x_position') ;
    posterior_pos =  ncread(file.name,'y_position') ;
    post_time =  ncread(file.name,'time') ;
    post_pos =  ncread(file.name,'y_position_value') ;
    post_goal =  ncread(file.name,'x_position_value') ;

    %% Load HD
    load(strcat('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\rotationData\',sessionInfo.FileName,'.rotation.mat'))

    %% Only extract goal decoding that is consistent for 500 ms
    [~,goal_dec] = max(posterior_goal,[],1);
    result = NaN(size(goal_dec));  % Initialize result with NaNs
    last_value = nan;
    count = 0;
    
    for i = 1:length(goal_dec)
        if goal_dec(i) ~= last_value
            if count >= decodingWin
                result(i-count:i-1) = last_value;  % Mark previous values as last_value
            end
            last_value = goal_dec(i);
            count = 1;
        else
            count = count + 1;
        end
    end
    if count >= decodingWin
        result(end-count+1:end) = last_value;  % Mark last batch as last_value
    end

    %% Within each trial, find the earliest time that the goal is decoded consistently    
    ts_dec = [];
    trial_dec = [];
    Dec.countDec(ss,1:6) = 0 ;
    for tt = 1:6
        Dec.totalDec(ss,tt) = sum(behavTrials.linTrial ==0 & behavTrials.probe == 0 & behavTrials.lickLoc==tt-1);
    end
    for tt = 1:length(behavTrials.lickLoc)
        if behavTrials.linTrial(tt)==0 && behavTrials.probe(tt) == 0 
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
                Dec.countDec(ss,behavTrials.lickLoc(tt)+1) = Dec.countDec(ss,behavTrials.lickLoc(tt)+1)+1;
            end
        end
    end

    %% Relate this decoding to frequency, speed etc    

    speed_goal = [];
    hd_goal = [];
    nose_x = [];
    distToGoal = [];

    portlabel = {'port1','port2','port3','port4','port5','port6'};
    for ff = 1:6
        freqtoGoal.(portlabel{ff}) = [];
        portdistToGoal.(portlabel{ff}) = [];
    end

    for ii = 1:length(ts_dec)
        [~,idxCurGoal] = min(abs(tracking.timestamps-ts_dec(ii)));
        speed_goal(ii,:) = tracking.position.vy(idxCurGoal-tWin:idxCurGoal+tWin);
        pos = tracking.position.y(idxCurGoal-tWin:idxCurGoal+tWin);
        gainid = behavTrials.lickLoc(trial_dec(ii))+1;
        distToGoal(ii,:) = (1/gain(gainid))*122-pos;
        if behavTrials.probe(trial_dec(ii))==1 || behavTrials.correct(trial_dec(ii))==0
            freq = nan;
            posG = nan;
        else
            freq = tracking.position.y(idxCurGoal)*gain(gainid)/122;
            posG = tracking.position.y(idxCurGoal);
        end
        freqtoGoal.(portlabel{gainid}) = [freqtoGoal.(portlabel{gainid}) 2000*(10.^(freqExp*freq))];
        portdistToGoal.(portlabel{gainid}) = [portdistToGoal.(portlabel{gainid}) (1/gain(gainid))*122-posG];

        % Only do this for the newer mice
        hd_goal(ii,:) = trackRotation.angle(idxCurGoal-tWin:idxCurGoal+tWin);
        nose_x(ii,:) = trackRotation.nosePos(idxCurGoal-tWin:idxCurGoal+tWin,1);
    end

    Dec.speed = [Dec.speed; speed_goal];
    Dec.distToGoal = [Dec.distToGoal; distToGoal];
    for ff = 1:6
        Dec.freqtoGoal.(portlabel{ff}) = [Dec.freqtoGoal.(portlabel{ff}) freqtoGoal.(portlabel{ff})];
        Dec.portdistToGoal.(portlabel{ff}) = [Dec.portdistToGoal.(portlabel{ff}) portdistToGoal.(portlabel{ff})];
    end

    % Only do this for the newer mice
    Dec.hd = [Dec.hd; hd_goal];
    Dec.nose = [Dec.nose; nose_x];
        
end

if plotfig

    fig2 = figure;
    t = linspace(-2,2,121);

    subplot(2,3,1)
    nhist(Dec.freqtoGoal,'samebins','binfactor',3,'color','sequential','proportion','smooth')
    xscale log
    xlim([2000 23000])
    
    subplot(2,3,2)
    nhist(Dec.portdistToGoal,'samebins','binfactor',3,'color','sequential','proportion','smooth')
    xlim([0 50])
    
    plotAvgStd(Dec.nose,2,3,3,fig2,t',[0 0 1], 0)
    title('Nose position')
    
    plotAvgStd(Dec.distToGoal,2,3,4,fig2,t',[0 0 1], 0)
    title('distancetoGoal')
    
    plotAvgStd(Dec.speed,2,3,5,fig2,t',[0 0 1], 0)
    title('speed')
    
    plotAvgStd(Dec.hd,2,3,6,fig2,t',[0 0 1], 0)
    title('headdirection')

end

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
