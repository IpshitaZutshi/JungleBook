function AllDec = linkDecodingtoDeceleration(varargin)

p = inputParser;
addParameter(p,'plotfig',true,@islogical);

parse(p,varargin{:});
plotfig = p.Results.plotfig;

sess = {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
    'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
    'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   23
    'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...
    'IZ40\Final\IZ40_220708_sess17','IZ40\Final\IZ40_220714_sess18',...27
    'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
    'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
    'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
    'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    36
    'IZ44\Final\IZ44_220827_sess4', 'IZ44\Final\IZ44_220828_sess5',...
    'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...
    'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
    'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',... 45
    'IZ47\Final\IZ47_230626_sess15','IZ47\Final\IZ47_230707_sess24',...
    'IZ47\Final\IZ47_230710_sess25','IZ47\Final\IZ47_230712_sess27',...49
    'IZ48\Final\IZ48_230628_sess17','IZ48\Final\IZ48_230703_sess21',...
    'IZ48\Final\IZ48_230705_sess22','IZ48\Final\IZ48_230714_sess28'};  
  
expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';

%sess = {'IZ47\Final\IZ47_230626_sess15'};

for dt = 1:4
    AllDec.pos{dt} = [];
    AllDec.goal{dt} = [];
    AllDec.truegoal{dt} = [];
end

for ii = 1:length(sess)

    %% Load files
    cd(strcat(expPath,sess{ii}))    
    file = dir('*.Tracking.Behavior.mat');
    load(file.name)
    file = dir('*.TrialBehavior.Behavior.mat');
    load(file.name)
    file = dir(['*spikes.cellinfo.mat']);
    load(file(1).name);

    if spikes.numcells<150
        continue
    end


    % Find deceleration points
    Dec = findDecelerationPoints('plotfig',false);

    % Load decoding
    decodingPath1 = 'Z:\Homes\zz737\ipshita_data\Auditory_Task\';
    decodingPath2 = '\py_data\theta_decoding_lickLoc_y\up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].nc';
    
    posterior_goal =  ncread(strcat(decodingPath1,sess{ii},decodingPath2),'x_position') ;
    posterior_pos =  ncread(strcat(decodingPath1,sess{ii},decodingPath2),'y_position') ;
    post_time =  ncread(strcat(decodingPath1,sess{ii},decodingPath2),'time') ;
    post_pos =  ncread(strcat(decodingPath1,sess{ii},decodingPath2),'y_position_value') ;
    post_goal =  ncread(strcat(decodingPath1,sess{ii},decodingPath2),'x_position_value') ;
    
    % For the 3 types of decelerations, find the difference between current
    % position and decoded position around a 1 second window
    for dt = 1:4
                
        decTS = Dec.ts(Dec.decType==dt);            
        
        decAvg.pos{dt} = nan(length(decTS),201);
        decAvg.goal{dt} = nan(length(decTS),201);
        decAvg.truegoal{dt} = nan(length(decTS),201);

        for dd = 1:length(decTS)
            
            %Find the trial corresponding to this deceleration
            [~,idxTrial] = min(abs(behavTrials.timestamps(:,2)-decTS(dd)));
             
            % Find the indexes in tracking, as well as decoding corresponding to a +/-1 window around the
            % timestamps

            % tracking
            [~,frameidx] = min(abs(tracking.timestamps-decTS(dd)));
            startTime = tracking.timestamps(frameidx-30);
            endTime = tracking.timestamps(frameidx+30);

            %decoding
            [diffTimestart,idxstart] = min(abs(post_time-startTime));
            if post_time(idxstart)<startTime %Take the next index
                idxstart = idxstart+1;
            end        


            [diffTimeend,idxend] = min(abs(post_time-endTime));
            if post_time(idxend)>endTime %Take the previous index
                idxend = idxend-1;
            end  

            if diffTimestart>1 ||diffTimeend>1
                continue
            end

            [~,decGoalidx] = max(posterior_goal(:,idxstart:idxend));

            truePos = [];
            trueGoal = [];
            curGoal = [];
            tsOffset = [];
            for kk = idxstart:idxend
                [~,idxTrack] = min(abs(tracking.timestamps-post_time(kk)));
                tsOffset = [tsOffset post_time(kk)-decTS(dd)];
                truePos = [truePos tracking.position.y(idxTrack)];
                curGoal = [curGoal decGoalidx(1)];
                trueGoal = [trueGoal behavTrials.lickLoc(idxTrial)+1];
            end
            tsOffset = round(tsOffset*100)+101;

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
            decGoalDiff = decGoalidx-trueGoal;  
            decGoalDiffCur = decGoalidx-curGoal; 

            decAvg.truegoal{dt}(dd,tsOffset) = decGoalDiff;
            decAvg.pos{dt}(dd,tsOffset) = decPosDiff;
            decAvg.goal{dt}(dd,tsOffset) = decGoalDiffCur;
        end
    end

    for dt = 1:4
        AllDec.pos{dt} = [AllDec.pos{dt}; decAvg.pos{dt}];
        AllDec.goal{dt} = [AllDec.goal{dt}; decAvg.goal{dt}];
        AllDec.truegoal{dt} = [AllDec.truegoal{dt}; decAvg.truegoal{dt}];
    end
       
end

save('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\DeliberationDecodingSummary.mat', 'AllDec'); 

if plotfig

    fig2 = figure;
    col = {'b','k','r','m'};

    t = linspace(-1,1,201);
    for dt = 1:4
        plotAvgStd(AllDec.pos{dt},3,4,dt,fig2,t',col{dt})
        title('Position')  
        ylim([-7 3])
        if dt == 1
            xlim([-1 0.3])
        end

        plotAvgStd(AllDec.goal{dt},3,4,dt+4,fig2,t',col{dt})
        title('Goal')  
        ylim([-1 1])
        if dt == 1
            xlim([-1 0.3])
        end

        plotAvgStd(AllDec.truegoal{dt},3,4,dt+8,fig2,t',col{dt})
        title('Relative to chosen goal')
        ylim([-1 1])
        if dt == 1
            xlim([-1 0.3])
        end

    end

end
end

function plotAvgStd(array,numrows,numcol,subplotlocation,figureHandle,xAxis,col)

    subplot(numrows, numcol, subplotlocation, 'Parent', figureHandle);

    meanpsth = nanmean(array,1);
    rowsWithNaN = all(isnan(array), 2);
    numRowsWithoutNaN = sum(~rowsWithNaN);

    stdpsth = nanstd(array,1)./sqrt(numRowsWithoutNaN);
    lArr  = meanpsth-stdpsth;
    uArr = meanpsth+stdpsth;

    a = ~isnan(lArr);

    fill([xAxis(a); flipud(xAxis(a))],[lArr(a)'; flipud(uArr(a)')],col,'linestyle','none','FaceAlpha',0.5);                    
    hold on
    hi = line(xAxis(a),meanpsth(a),'LineWidth',1,'Color',col);

end