function Dec = compileGoalDecodingStats_changePoint(varargin)

    p = inputParser;
    addParameter(p,'plotfig',true);
    addParameter(p,'tWin',60);
    addParameter(p,'plotIndTrials',false);
    addParameter(p,'selectFirst',false);
    
    parse(p,varargin{:});
    plotfig = p.Results.plotfig;
    tWin = p.Results.tWin;
    plotIndTrials = p.Results.plotIndTrials;
    selectFirst = p.Results.selectFirst;

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

     % Define paths and parameters
    expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';
    decodingPath = 'Z:\Homes\zz737\ipshita_data\Auditory_Task\';
    decodingName = 'py_data\theta_decoding_lickLoc_y\up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].nc';
    changePointName = 'py_data/theta_decoding_lickLoc_y/change_point_posterior_up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].mat';

    gain = [122/9, 122/32 122/55.53 122/79.62 122/102.79 122/122];
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
        [posterior_goal, posterior_pos, post_time, post_pos, post_goal,change_point,trial] = loadDecodingData(decodingPath, sess{ss}, decodingName,changePointName);

        % Load HD
        if str2double(sess{ss}(3:4))==47 || str2double(sess{ss}(3:4))==48
            load(strcat('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\rotationData\',sessionInfo.FileName,'.rotation.mat'))
        end
       
        %% Within each trial, find the earliest time that the goal is decoded consistently    
        ts_dec = [];
        ts_dec_shuff = [];
        trial_dec = [];
        trial_dec_shuff = [];
        ts_dec_error = [];
        trial_dec_error = [];

        Dec.countDec(ss,1:6) = 0 ;
        Dec.countDecError(ss,1:6) = 0 ;
        Dec.countDecShuffle(ss,1:6) = 0 ;

        for tt = 1:6
            Dec.totalDec(ss,tt) = sum(behavTrials.linTrial ==0 & behavTrials.probe == 0 & behavTrials.correct == 1 & behavTrials.lickLoc==tt-1);
        end
        for tt = 1:6
            Dec.totalDecError(ss,tt) = sum(behavTrials.linTrial ==0 & behavTrials.probe == 0 & behavTrials.correct == 0 & behavTrials.lickLoc==tt-1);
        end

        % Cycle through trials to find the first detected time 
        for tt = 1:(length(behavTrials.lickLoc)-1)
            if behavTrials.linTrial(tt)==0 && behavTrials.probe(tt) == 0 && behavTrials.correct(tt) == 1

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
                        ts_dec = [ts_dec post_time(idxGoal+idxstart)];
                        trial_dec = [trial_dec tt];
                        Dec.countDec(ss,behavTrials.lickLoc(tt)+1) = Dec.countDec(ss,behavTrials.lickLoc(tt)+1)+1;
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

            elseif behavTrials.linTrial(tt)==0 && behavTrials.probe(tt) == 0 && behavTrials.correct(tt) == 0

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
                        ts_dec_error = [ts_dec_error post_time(idxGoal+idxstart)];
                        trial_dec_error = [trial_dec_error tt];
                        Dec.countDecError(ss,behavTrials.lickLoc(tt)+1) = Dec.countDecError(ss,behavTrials.lickLoc(tt)+1)+1;
                    end
                 
                 end

            end
        end

        totalTrials = sum(behavTrials.linTrial==0 & behavTrials.probe == 0 & behavTrials.correct == 1);
        Dec.propDecoded = [Dec.propDecoded; length(ts_dec)./totalTrials];
        Dec.propDecodedShuff = [Dec.propDecodedShuff; length(ts_dec_shuff)./totalTrials];
        Dec.mouseID = [Dec.mouseID;str2double(sess{ss}(3:4))];
        totalTrials = sum(behavTrials.linTrial==0 & behavTrials.probe == 0 & behavTrials.correct == 0);
        Dec.propDecodedError = [Dec.propDecodedError; length(ts_dec_error)./totalTrials];
        

        %% Relate this decoding to frequency, speed etc    
    
        speed_goal = [];
        hd_goal = [];
        nose_x = [];
        distToGoal = [];
        timetoGoal = [];
        timetoGoalShuff = [];
        timetoGoalError = [];
        decFreq =[];
    
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
            timetoGoal(ii,1) = behavTrials.timestamps(trial_dec(ii),2) - ts_dec(ii);
            if behavTrials.probe(trial_dec(ii))==1 || behavTrials.correct(trial_dec(ii))==0
                freq = nan;
                posG = nan;
            else
                freq = tracking.position.y(idxCurGoal)*gain(gainid)/122;
                posG = tracking.position.y(idxCurGoal);
            end
            freqtoGoal.(portlabel{gainid}) = [freqtoGoal.(portlabel{gainid}) 2000*(10.^(freqExp*freq))];
            portdistToGoal.(portlabel{gainid}) = [portdistToGoal.(portlabel{gainid}) (1/gain(gainid))*122-posG];
            decFreq = [decFreq;2000*(10.^(freqExp*freq))];

            % Only do this for the newer mice
            if str2double(sess{ss}(3:4))==47 || str2double(sess{ss}(3:4))==48
                hd_goal(ii,:) = trackRotation.angle(idxCurGoal-tWin:idxCurGoal+tWin);
                nose_x(ii,:) = trackRotation.nosePos(idxCurGoal-tWin:idxCurGoal+tWin,1);
            end
        end
    
        for ii = 1:length(ts_dec_shuff)
            timetoGoalShuff(ii,1) = behavTrials.timestamps(trial_dec_shuff(ii),2) - ts_dec_shuff(ii);
        end
    
        Dec.speed = [Dec.speed; nanmean(speed_goal,1)];
        Dec.distToGoal = [Dec.distToGoal; distToGoal];
        Dec.timetoGoal = [Dec.timetoGoal; timetoGoal];
        Dec.timetoGoalShuff = [Dec.timetoGoalShuff; timetoGoalShuff];
        
        Dec.decFreq = [Dec.decFreq; decFreq];
        for ff = 1:6
            Dec.freqtoGoal.(portlabel{ff}) = [Dec.freqtoGoal.(portlabel{ff}) freqtoGoal.(portlabel{ff})];
            Dec.portdistToGoal.(portlabel{ff}) = [Dec.portdistToGoal.(portlabel{ff}) portdistToGoal.(portlabel{ff})];
        end
    
        for ii = 1:length(ts_dec_error)
            timetoGoalError(ii,1) = behavTrials.timestamps(trial_dec_error(ii),2) - ts_dec_error(ii);
        end
        Dec.timetoGoalError = [Dec.timetoGoalError; timetoGoalError];
        
        % Only do this for the newer mice
        if str2double(sess{ss}(3:4))==47 || str2double(sess{ss}(3:4))==48
            Dec.hd = [Dec.hd; hd_goal];
            Dec.nose = [Dec.nose; nose_x];
        end
            
    end

    if plotfig
    
        fig2 = figure;
        t = linspace(-2,2,121);
    
        col = [83/255 0/255 0/255;...
        184/255 15/255 10/255;...
        241/255 114/255 42/255;...
        249/255 197/255 81/255;...
        143/255 189/255 107/255;...
        87/255 116/255 144/255];

        % subplot(3,3,1)
        % Stats.propDecoded = groupStats([{Dec.propDecoded},{Dec.propDecodedShuff},{Dec.propDecodedError}],[],'inAxis',true);
        % 
        subplot(3,3,2)
        Stats.propDecoded = groupStats([{Dec.timetoGoal},{Dec.timetoGoalShuff},{Dec.timetoGoalError}],[],'inAxis',true);
    
        subplot(3,3,3)
        data = [{Dec.freqtoGoal.port2'},{Dec.freqtoGoal.port3'},{Dec.freqtoGoal.port4'},{Dec.freqtoGoal.port5'},{Dec.freqtoGoal.port6'}];
        Stats = groupStats(data,[],'inAxis',true,'color',col(2:end,:));
        yscale log
        
        % subplot(3,3,4)
        % histogram(Dec.decFreq,[2000:1000:23000])
        % 
        plotAvgStd(Dec.distToGoal,3,3,5,fig2,t',[0 0 1], 0)
        title('distancetoGoal')
        
        plotAvgStd(Dec.speed,3,3,6,fig2,t',[0 0 1], 0)
        title('speed')

        plotAvgStd(Dec.hd,3,3,7,fig2,t',[0 0 1], 0)
        title('headdirection')

        plotAvgStd(Dec.nose,3,3,8,fig2,t',[0 0 1], 0)
        title('Nose position')
    
        % subplot(3,3,9)
        % data = [{Dec.propDecoded(Dec.mouseID==39)},{Dec.propDecoded(Dec.mouseID==43)},{Dec.propDecoded(Dec.mouseID==44)},{Dec.propDecoded(Dec.mouseID==47)},{Dec.propDecoded(Dec.mouseID==48)}];
        % Stats.mouseID = groupStats(data,[],'inAxis',true);
    end

end

function Dec = initializeDecStruct()
    % Initialize Dec structure with empty fields
    Dec.speed = [];
    Dec.hd = [];
    Dec.nose = [];    
    Dec.distToGoal = [];

    Dec.timetoGoal = [];
    Dec.timetoGoalShuff = [];
    Dec.timetoGoalError = [];

    Dec.decFreq = [];
    Dec.propDecoded = [];
    Dec.propDecodedShuff = [];
    Dec.propDecodedError = [];
    Dec.mouseID = [];

    portlabel = {'port1','port2','port3','port4','port5','port6'};
    for ff = 1:6
        Dec.freqtoGoal.(portlabel{ff}) = [];
        Dec.portdistToGoal.(portlabel{ff}) = [];
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
