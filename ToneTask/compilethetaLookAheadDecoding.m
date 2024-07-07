function Dec = compilethetaLookAheadDecoding(varargin)

    p = inputParser;
    addParameter(p,'plotfig',true);
    addParameter(p,'decWin',150); % Window for speed/ distance profile
    
    parse(p,varargin{:});
    plotfig = p.Results.plotfig;
    decWin = p.Results.decWin;

    % Define session data
    sess = {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
        'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
        'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   
        'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
        'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
        'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
        'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    
        'IZ44\Final\IZ44_220827_sess4', 'IZ44\Final\IZ44_220828_sess5',...
        'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...
        'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
        'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',...  
        'IZ47\Final\IZ47_230626_sess15','IZ47\Final\IZ47_230707_sess24',...
        'IZ47\Final\IZ47_230710_sess25','IZ47\Final\IZ47_230712_sess27',...  
        'IZ48\Final\IZ48_230628_sess17','IZ48\Final\IZ48_230703_sess21',...
        'IZ48\Final\IZ48_230705_sess22','IZ48\Final\IZ48_230714_sess28'};  
    %sess = {'IZ47\Final\IZ47_230626_sess15','IZ48\Final\IZ48_230714_sess28'};

     % Define paths and parameters
    expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';
    decodingPath = 'Z:\Homes\zz737\ipshita_data\Auditory_Task\';
    decodingName = 'py_data\theta_decoding_lickLoc_y\up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].nc';

    portPos = [57.53 81.62 104.79 122];
    freqExp = log10(22000/2000);

     % Initialize Dec structure
    Dec = initializeDecStruct();

     % Loop through sessions
    for ss = 1:length(sess)

        cd(strcat(expPath,sess{ss}))   
        [behavTrials, spikes, sessionInfo, tracking] = loadSessionData();

        if spikes.numcells<150
            continue
        end

        % Load decoding data
        [posterior_goal, posterior_pos, post_time, post_pos, post_goal] = loadDecodingData(decodingPath, sess{ss}, decodingName);

        % Load HD
        if str2double(sess{ss}(3:4))==47 || str2double(sess{ss}(3:4))==48
            load(strcat('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\rotationData\',sessionInfo.FileName,'.rotation.mat'))
        end

        %% Cycle through trials to find the look ahead for the target port compared to the previous port

        portName = {'Port3', 'Port4','Port5','Port6'};
        for pp = 1:4
            decGoalDiffTarget.(portName{pp}) = [];
            decPosDiffTarget.(portName{pp}) = [];                
            decGoalDiffPreTarget.(portName{pp}) = [];
            decPosDiffPreTarget.(portName{pp}) = [];
        end

        
        for tt = 1:(length(behavTrials.lickLoc)-1)
            if behavTrials.linTrial(tt)==0 && behavTrials.probe(tt) == 0 && behavTrials.correct(tt) == 1 && behavTrials.lickLoc(tt)>=2           
                % Only look at times when the mouse has crossed the first
                % port
                [~,idxstart] = min(abs(post_time-behavTrials.timestamps(tt,1)));
                if post_time(idxstart)<behavTrials.timestamps(tt,1) %Take the next index
                    idxstart = idxstart+1;
                end        
                [~,idxend] = min(abs(post_time-behavTrials.timestamps(tt,2)));
                if post_time(idxend)>behavTrials.timestamps(tt,2) %Take the previous index
                    idxend = idxend-1;
                end  

                %For each timepoint in the posterior, calculate y position
                %and goal
                truePos = [];
                trueGoal = [];
                for kk = idxstart:idxend
                    [~,idxTrack] = min(abs(tracking.timestamps-post_time(kk)));
                    truePos = [truePos tracking.position.y(idxTrack)];
                    trueGoal = [trueGoal behavTrials.lickLoc(tt)+1];
                end
                
                if isempty(truePos)
                    continue
                end

                % Calculate the difference between current position and
                % decoded position    
                [~,decPosidx] = max(posterior_pos(:,idxstart:idxend));
                decPos = post_pos(decPosidx);
                decPosDiff = decPos-truePos';

                % Calculate the difference between current goal and
                % decoded goal
                [~,decGoalidx] = max(posterior_goal(:,idxstart:idxend));
                decGoalDiff = decGoalidx-trueGoal;                     

                %Find points when the mouse is crossing ports
                curTarget = behavTrials.lickLoc(tt)+1;
                prevTargets = 3:curTarget;
                for pt = 1:length(prevTargets)
                    if prevTargets(pt) ==curTarget
                        if length(decGoalDiff) <=decWin
                            continue
                        end
                        decGoalDiffTarget.(portName{prevTargets(pt)-2}) = [decGoalDiffTarget.(portName{prevTargets(pt)-2});decGoalDiff(end-decWin:end)];
                        decPosDiffTarget.(portName{prevTargets(pt)-2}) = [decPosDiffTarget.(portName{prevTargets(pt)-2});decPosDiff(end-decWin:end)'];
                    else
                        [~,idxPort] = find(truePos>portPos(prevTargets(pt)-2),1,'first');
                        if idxPort <=decWin
                            continue
                        end
                        decGoalDiffPreTarget.(portName{prevTargets(pt)-2}) = [decGoalDiffPreTarget.(portName{prevTargets(pt)-2});decGoalDiff(idxPort-decWin:idxPort)];
                        decPosDiffPreTarget.(portName{prevTargets(pt)-2}) = [decPosDiffPreTarget.(portName{prevTargets(pt)-2});decPosDiff(idxPort-decWin:idxPort)'];
                    end
                end
            end
        end

        for pp = 1:4
            Dec.decGoalDiffTarget.(portName{pp}) = [Dec.decGoalDiffTarget.(portName{pp});nanmean(decGoalDiffTarget.(portName{pp}),1)];
            Dec.decPosDiffTarget.(portName{pp}) = [Dec.decPosDiffTarget.(portName{pp});nanmean(decPosDiffTarget.(portName{pp}),1)];                
            Dec.decGoalDiffPreTarget.(portName{pp}) = [Dec.decGoalDiffPreTarget.(portName{pp});nanmean(decGoalDiffPreTarget.(portName{pp}),1)];
            Dec.decPosDiffPreTarget.(portName{pp}) = [Dec.decPosDiffPreTarget.(portName{pp});nanmean(decPosDiffPreTarget.(portName{pp}),1)];
        end

    end

    if plotfig
    
        fig2 = figure;
        t = linspace(-1.5,0,151);
        
        for kk = 1:3
            plotAvgStd(Dec.decGoalDiffTarget.(portName{kk}),4,2,2*(kk-1)+1,fig2,t',[0 0 1], 0)
            %hold on
            %plotAvgStd(Dec.decGoalDiffPreTarget.(portName{kk}),4,2,2*(kk-1)+1,fig2,t',[0 1 0], 0)
            title('Goal')
            
            plotAvgStd(Dec.decPosDiffTarget.(portName{kk}),4,2,2*(kk-1)+2,fig2,t',[0 0 1], 0)
            hold on
            plotAvgStd(Dec.decPosDiffPreTarget.(portName{kk}),4,2,2*(kk-1)+2,fig2,t',[0 1 0], 0)
            title('Target position')
        end
    
    end

end

function Dec = initializeDecStruct()
    % Initialize Dec structure with empty fields
    portName = {'Port3', 'Port4','Port5','Port6'};
    for pp = 1:4
        Dec.decGoalDiffTarget.(portName{pp}) = [];
        Dec.decPosDiffTarget.(portName{pp}) = [];                
        Dec.decGoalDiffPreTarget.(portName{pp}) = [];
        Dec.decPosDiffPreTarget.(portName{pp}) = [];
    end
end

function [behavTrials, spikes, sessionInfo, tracking] = loadSessionData()
    % Load session data
    file = dir(['*.Tracking.Behavior.mat']);
    load(file(1).name);
    file = dir(['*TrialBehavior.Behavior.mat']);
    load(file(1).name);
    file = dir(['*spikes.cellinfo.mat']);
    load(file(1).name);
    [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts', true);
    if ~isfield(behavTrials,'probe')
        behavTrials.probe(1:length(behavTrials.linTrial),1) = 0;
    end
end

function [posterior_goal, posterior_pos, post_time, post_pos, post_goal] = loadDecodingData(decodingPath, sess, decodingName)
    % Load decoding data
    file_name = strcat(decodingPath, '\', sess, '\', decodingName);
    posterior_goal = ncread(file_name, 'x_position');
    posterior_pos = ncread(file_name, 'y_position');
    post_time = ncread(file_name, 'time');
    post_pos = ncread(file_name, 'y_position_value');
    post_goal = ncread(file_name, 'x_position_value');
end

function result = processGoalDecoding(posterior_goal, decodingWin)

    % Process goal decoding to ensure consistency over a certain bin size
    [~, goal_dec] = max(posterior_goal, [], 1);
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
