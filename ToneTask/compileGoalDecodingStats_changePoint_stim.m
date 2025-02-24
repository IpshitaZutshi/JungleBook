function Dec = compileGoalDecodingStats_changePoint_stim(varargin)

    p = inputParser;
    addParameter(p,'plotfig',true);
    addParameter(p,'tWin',60);
    addParameter(p,'plotIndTrials',false);
    addParameter(p,'selectFirst',false);
    
    parse(p,varargin{:});
    plotfig = p.Results.plotfig;
    tWin = p.Results.tWin;
    plotIndTrials = p.Results.plotIndTrials;

    sess = {'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   
        'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...
        'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    
        'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15'};   

     % Define paths and parameters
    expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';
    decodingPath = 'Z:\Homes\zz737\ipshita_data\Auditory_Task\';
    decodingName = 'py_data\theta_decoding_lickLoc_y\up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].nc';
    changePointName = 'py_data/theta_decoding_lickLoc_y/change_point_posterior_up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].mat';

     % Initialize Dec structure
    Dec = initializeDecStruct();

     % Loop through sessions
    for ss = 1:length(sess)

        cd(strcat(expPath,sess{ss}))   
        [behavTrials, spikes, sessionInfo, tracking] = loadSessionData();

        % Load decoding data
        [posterior_goal, posterior_pos, post_time, post_pos, post_goal,change_point,trial] = loadDecodingData(decodingPath, sess{ss}, decodingName,changePointName);

       
        %% Within each trial, find the earliest time that the goal is decoded consistently    
        ts_dec_nostim = [];
        ts_dec_shuff = [];
        ts_dec_stim = [];

        trial_dec_nostim = [];
        trial_dec_shuff = [];
        trial_dec_stim = []; 

        Dec.countDecNoStim(ss,1:6) = 0 ;
        Dec.countDecStim(ss,1:6) = 0 ;
        Dec.countDecShuffle(ss,1:6) = 0 ;

        for tt = 1:6
            Dec.totalDecNoStim(ss,tt) = sum(behavTrials.linTrial ==0 & behavTrials.stim == 0 & behavTrials.lickLoc==tt-1);
        end

        for tt = 1:6           
            Dec.totalDecStim(ss,tt) = sum(behavTrials.linTrial ==0 & behavTrials.stim == 1 & behavTrials.lickLoc==tt-1);
        end

        % Cycle through trials to find the first detected time 
        for tt = 1:(length(behavTrials.lickLoc)-1)
            if behavTrials.linTrial(tt)==0 && behavTrials.stim(tt) == 0

                [~,idxstart] = min(abs(post_time-behavTrials.timestamps(tt,1)));
                if post_time(idxstart)<behavTrials.timestamps(tt,1) %Take the next index
                    idxstart = idxstart+1;
                end        
                [~,idxend] = min(abs(post_time-behavTrials.timestamps(tt,2)));
                if post_time(idxend)>behavTrials.timestamps(tt,2) %Take the previous index
                    idxend = idxend-1;
                end   

                [~,decGoal] = max(posterior_goal(:,idxstart:idxend));

                %% Get change points for that trial
                 if sum(trial==(tt-1))>0
                    curChanges = change_point{trial==(tt-1)};
    
                    idxGoal = curChanges(end);
                    trialDecGoal = mode(decGoal(curChanges(end)+1:end));
    
                    if trialDecGoal==(behavTrials.lickLoc(tt)+1)             
                        ts_dec_nostim = [ts_dec_nostim post_time(idxGoal+idxstart)];
                        trial_dec_nostim = [trial_dec_nostim tt];
                        Dec.countDecNoStim(ss,behavTrials.lickLoc(tt)+1) = Dec.countDecNoStim(ss,behavTrials.lickLoc(tt)+1)+1;
                    end

                    if trialDecGoal==randsample([behavTrials.lickLoc(tt)+1 6],1)
                        ts_dec_shuff = [ts_dec_shuff post_time(idxGoal+idxstart)];
                        trial_dec_shuff = [trial_dec_shuff tt];
                    end

                 end

                if plotIndTrials
                    figure
                    imagesc(idxstart:idxend,post_goal,posterior_goal(:,idxstart:idxend));
                    title(num2str(behavTrials.lickLoc(tt)))
                    set(gca,'YDir','normal')
                    hold on
                    for cc = 1:length(curChanges)
                        line([idxstart+curChanges(cc)-1 idxstart+curChanges(cc)-1],[0 6],'Color','w','LineWidth',1.5);  
                    end
                    plot(idxstart:idxend,decGoal-1,'Color','w','LineWidth',1.5); 
    
                end

           elseif behavTrials.linTrial(tt)==0 && behavTrials.stim(tt) == 1 

                [~,idxstart] = min(abs(post_time-behavTrials.timestamps(tt,1)));
                if post_time(idxstart)<behavTrials.timestamps(tt,1) %Take the next index
                    idxstart = idxstart+1;
                end        
                [~,idxend] = min(abs(post_time-behavTrials.timestamps(tt,2)));
                if post_time(idxend)>behavTrials.timestamps(tt,2) %Take the previous index
                    idxend = idxend-1;
                end   

                [~,decGoal] = max(posterior_goal(:,idxstart:idxend));

                %% Get change points for that trial
                 if sum(trial==(tt-1))>0
                    curChanges = change_point{trial==(tt-1)};
    
                    idxGoal = curChanges(end);
                    trialDecGoal = mode(decGoal(curChanges(end)+1:end));
    
                    if trialDecGoal==(behavTrials.lickLoc(tt)+1)             
                        ts_dec_stim = [ts_dec_stim post_time(idxGoal+idxstart)];
                        trial_dec_stim = [trial_dec_stim tt];
                        Dec.countDecStim(ss,behavTrials.lickLoc(tt)+1) = Dec.countDecStim(ss,behavTrials.lickLoc(tt)+1)+1;
                    end
                 
                 end

            end
        end


        totalTrials = sum(behavTrials.linTrial==0 & behavTrials.stim == 0);
        Dec.propDecoded = [Dec.propDecoded; length(ts_dec_nostim)./totalTrials];
        Dec.propDecodedShuff = [Dec.propDecodedShuff; length(ts_dec_shuff)./totalTrials];
        Dec.mouseID = [Dec.mouseID;str2double(sess{ss}(3:4))];

        totalTrials = sum(behavTrials.linTrial==0 & behavTrials.stim == 1);
        Dec.propDecodedStim = [Dec.propDecodedStim; length(ts_dec_stim)./totalTrials];
        

        %% Relate this decoding to frequency, speed etc    
    
        timetoGoalShuff = [];
        timetoGoalNoStim = [];
        timetoGoalStim = [];    

        for ii = 1:length(ts_dec_nostim)
            timetoGoalNoStim(ii,1) = behavTrials.timestamps(trial_dec_nostim(ii),2) - ts_dec_nostim(ii);
        end
    
        for ii = 1:length(ts_dec_shuff)
            timetoGoalShuff(ii,1) = behavTrials.timestamps(trial_dec_shuff(ii),2) - ts_dec_shuff(ii);
        end
    
        for ii = 1:length(ts_dec_stim)
            timetoGoalStim(ii,1) = behavTrials.timestamps(trial_dec_stim(ii),2) - ts_dec_stim(ii);
        end

        Dec.timetoGoalNoStim = [Dec.timetoGoalNoStim; timetoGoalNoStim];
        Dec.timetoGoalShuff = [Dec.timetoGoalShuff; timetoGoalShuff];
        Dec.timetoGoalStim = [Dec.timetoGoalStim; timetoGoalStim];        
            
    end


    if plotfig
    
        fig2 = figure;

        subplot(1,2,1)
        Stats.propDecoded = groupStats([{Dec.propDecoded},{Dec.propDecodedShuff},{Dec.propDecodedStim}],[],'inAxis',true);
        title('Fraction of trials decoded')

        subplot(1,2,2)
        Stats.timetoGoal = groupStats([{Dec.timetoGoalNoStim},{Dec.timetoGoalShuff},{Dec.timetoGoalStim}],[],'inAxis',true);        title('time to goal decoding')
    
    end

end

function Dec = initializeDecStruct()
    % Initialize Dec structure with empty fields

    Dec.timetoGoalStim = [];
    Dec.timetoGoalShuff = [];
    Dec.timetoGoalNoStim = [];

    Dec.propDecoded = [];
    Dec.propDecodedShuff = [];
    Dec.propDecodedStim = [];
    Dec.mouseID = [];
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
