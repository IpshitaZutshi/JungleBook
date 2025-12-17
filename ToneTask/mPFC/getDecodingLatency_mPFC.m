function getDecodingLatency_mPFC(varargin)

p = inputParser;
addParameter(p,'plotfig',true);
addParameter(p,'tWin',60);
addParameter(p,'plotIndTrialsBase',false);
addParameter(p,'plotIndTrialsStim',false);

parse(p,varargin{:});
plotfig = p.Results.plotfig;
tWin = p.Results.tWin;
plotIndTrialsBase =p.Results.plotIndTrialsBase;
plotIndTrialsStim = p.Results.plotIndTrialsStim;

sess = {'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   23
        'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...        
        'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    36        
        'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15'};  

 % Define paths and parameters
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';
decodingPath = 'Z:\Homes\zz737\ipshita_data\Auditory_Task\';
decodingName = 'py_data\theta_decoding_lickLoc_y\up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].nc';
changePointName = 'py_data/theta_decoding_lickLoc_y/change_point_posterior_up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].mat';

gain = [122/9, 122/32 122/55.53 122/79.62 122/102.79 122/122];
freqExp = log10(22000/2000);

 % Initialize Dec structure
Dec = initializeDecStruct();

figure

for ss = 1:length(sess)

    cd(strcat(expPath,sess{ss}))   
    [behavTrials, spikes, sessionInfo, tracking] = loadSessionData();
    % 
    % if spikes.numcells<150
    %     continue
    % end

    % Load decoding data
    [posterior_goal, posterior_pos, post_time, post_pos, post_goal,change_point,trial] = loadDecodingData(decodingPath, sess{ss}, decodingName,changePointName);

    %% Within each trial, find the earliest time that the goal is decoded consistently    
    ts_dec = [];
    ts_dec_shuff = [];
    ts_dec_stim = [];
    ts_dec_error = [];

    trial_dec = [];
    trial_dec_shuff = [];
    trial_dec_stim = [];                
    trial_dec_error = [];

    Dec.countDec(ss,1:6) = 0 ;
    Dec.countDecError(ss,1:6) = 0 ;
    Dec.countDecStim(ss,1:6) = 0 ;
    Dec.countDecShuffle(ss,1:6) = 0 ;

    for tt = 1:6
        Dec.totalDec(ss,tt) = sum(behavTrials.linTrial ==0 & behavTrials.stim == 0 & behavTrials.correct == 1 & behavTrials.lickLoc==tt-1);
    end
    for tt = 1:6
        Dec.totalDecError(ss,tt) = sum(behavTrials.linTrial ==0 & behavTrials.stim == 0 & behavTrials.correct == 0 & behavTrials.lickLoc==tt-1);
    end
    for tt = 1:6
        Dec.totalDecStim(ss,tt) = sum(behavTrials.linTrial ==0 & behavTrials.stim == 1 & behavTrials.lickLoc==tt-1);
    end

    % Cycle through trials to find the first detected time 
    for tt = 1:(length(behavTrials.lickLoc)-1)
        if behavTrials.linTrial(tt)==0 && behavTrials.stim(tt) == 0 && behavTrials.correct(tt) == 1

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

             if plotIndTrialsBase
                figure
                imagesc(idxstart:idxend,post_goal,posterior_goal(:,idxstart:idxend));
                title(strcat(num2str(behavTrials.lickLoc(tt)),num2str(tt)))
                set(gca,'YDir','normal')
                hold on
                for cc = 1:length(curChanges)
                    line([idxstart+curChanges(cc)-1 idxstart+curChanges(cc)-1],[0 6],'Color','w','LineWidth',1.5);  
                end
                plot(idxstart:idxend,decGoal-1,'Color','w','LineWidth',1.5); 

            end

        elseif behavTrials.linTrial(tt)==0 && behavTrials.stim(tt) == 0 && behavTrials.correct(tt) == 0

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

             if plotIndTrialsStim
                figure
                imagesc(idxstart:idxend,post_goal,posterior_goal(:,idxstart:idxend));
                title(strcat(num2str(behavTrials.lickLoc(tt)),num2str(tt)))
                set(gca,'YDir','normal')
                hold on
                for cc = 1:length(curChanges)
                    line([idxstart+curChanges(cc)-1 idxstart+curChanges(cc)-1],[0 6],'Color','w','LineWidth',1.5);  
                end
                plot(idxstart:idxend,decGoal-1,'Color','w','LineWidth',1.5); 

            end

        end
    end

    totalTrials = sum(behavTrials.linTrial==0 & behavTrials.stim == 0 & behavTrials.correct == 1);
    Dec.propDecoded = [Dec.propDecoded; length(ts_dec)./totalTrials];
    Dec.propDecodedShuff = [Dec.propDecodedShuff; length(ts_dec_shuff)./totalTrials];
    Dec.mouseID = [Dec.mouseID;str2double(sess{ss}(3:4))];

    totalTrials = sum(behavTrials.linTrial==0 & behavTrials.stim == 0 & behavTrials.correct == 0);
    Dec.propDecodedError = [Dec.propDecodedError; length(ts_dec_error)./totalTrials];

    totalTrials = sum(behavTrials.linTrial==0 & behavTrials.stim == 1);
    Dec.propDecodedStim = [Dec.propDecodedStim; length(ts_dec_stim)./totalTrials];
    
    %% Relate this decoding to frequency, speed etc    

    distToGoal = [];
    timetoGoal = [];
    timetoGoalShuff = [];
    timetoGoalError = [];
    timetoGoalProbe = [];

    portlabel = {'port1','port2','port3','port4','port5','port6'};
    for ff = 1:6
        freqtoGoal.(portlabel{ff}) = [];
        portdistToGoal.(portlabel{ff}) = [];
    end

    for ii = 1:length(ts_dec)
        [~,idxCurGoal] = min(abs(tracking.timestamps-ts_dec(ii)));
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
        portdistToGoal.(portlabel{gainid}) = [portdistToGoal.(portlabel{gainid}) (1/gain(gainid))*122-posG];
    end

    for ii = 1:length(ts_dec_shuff)
        timetoGoalShuff(ii,1) = behavTrials.timestamps(trial_dec_shuff(ii),2) - ts_dec_shuff(ii);
    end

    Dec.distToGoal = [Dec.distToGoal; distToGoal];
    Dec.timetoGoal = [Dec.timetoGoal; timetoGoal];
    Dec.timetoGoalShuff = [Dec.timetoGoalShuff; timetoGoalShuff];

    for ff = 1:6
        Dec.portdistToGoal.(portlabel{ff}) = [Dec.portdistToGoal.(portlabel{ff}) portdistToGoal.(portlabel{ff})];
    end

    for ii = 1:length(ts_dec_error)
        timetoGoalError(ii,1) = behavTrials.timestamps(trial_dec_error(ii),2) - ts_dec_error(ii);
    end
    Dec.timetoGoalError = [Dec.timetoGoalError; timetoGoalError];
    
    for ii = 1:length(ts_dec_stim)
        timetoGoalStim(ii,1) = behavTrials.timestamps(trial_dec_stim(ii),2) - ts_dec_stim(ii);
    end
    Dec.timetoGoalStim = [Dec.timetoGoalStim; timetoGoalStim];
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

    subplot(2,2,1)
    Stats.propDecoded = groupStats([{Dec.propDecoded},{Dec.propDecodedShuff},{Dec.propDecodedError}, {Dec.propDecodedStim}],[],'inAxis',true,'labelSummary',false);
    ylabel('Fraction decoded correctly')
    legend({'Stim','Error','Shuffle','Correct'})

    subplot(2,2,2)
    Stats.timetoGoal = groupStats([{Dec.timetoGoal},{Dec.timetoGoalError},{Dec.timetoGoalStim}],[],'inAxis',true,'labelSummary',false);
    ylabel('Time from lick')
    legend({'Stim','Error','Correct'})

end

end

function Dec = initializeDecStruct()
    % Initialize Dec structure with empty fields
    Dec.distToGoal = [];

    Dec.timetoGoal = [];
    Dec.timetoGoalShuff = [];
    Dec.timetoGoalError = [];
    Dec.timetoGoalStim = [];

    Dec.decFreq = [];
    Dec.propDecoded = [];
    Dec.propDecodedShuff = [];
    Dec.propDecodedError = [];
    Dec.propDecodedStim = [];
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