function Summary = calculateTrajectoryVariance(cond);

if cond == 1
    sess= {'IZ41\Final\IZ41_220624_sess5','IZ41\Final\IZ41_220701_sess9',...
        'IZ41\Final\IZ41_220704_sess11','IZ41\Final\IZ41_220708_sess13',...
        'IZ41\Final\IZ41_220629_sess7','IZ41\Final\IZ41_220714_sess14',...
        'IZ45\Final\IZ45_230410_sess11','IZ45\Final\IZ45_230420_sess19',...
        'IZ45\Final\IZ45_230414_sess15','IZ45\Final\IZ45_230417_sess16','IZ45\Final\IZ45_230424_sess21',... 
        'IZ46\Final\IZ46_230406_sess9','IZ46\Final\IZ46_230407_sess10',...
        'IZ46\Final\IZ46_230410_sess11','IZ46\Final\IZ46_230413_sess14','IZ46\Final\IZ46_230420_sess19'
        }; 
elseif cond==2
    sess= {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
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
    'IZ47\Final\IZ47_230707_sess24',...
    'IZ47\Final\IZ47_230710_sess25','IZ47\Final\IZ47_230712_sess27',...49
    'IZ48\Final\IZ48_230628_sess17','IZ48\Final\IZ48_230703_sess21',...
    'IZ48\Final\IZ48_230705_sess22','IZ48\Final\IZ48_230714_sess28'};   

end

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';

% bin y position by 2 cm bins
yBins = 1:2:122;

for ii = 1:length(sess)
    cd(strcat(expPath,sess{ii}))    
    file = dir(['*.Tracking.Behavior.mat']);
    load(file(1).name);
    file = dir(['*TrialBehavior.Behavior.mat']);
    load(file(1).name);
    [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts', true);
    
    % Extract x & y position of each trial
    for pf = 1:(size(behavTrials.timestamps,1))    
        [idx] = InIntervals(tracking.timestamps,behavTrials.timestamps(pf,:));
        yCur = tracking.position.y(idx);
        xCur = (tracking.position.x(idx));
            
        for yy = 1:(length(yBins)-1)
            idxY = yCur>=yBins(yy) & yCur<yBins(yy+1);
            position(pf,yy) = nanmean(xCur(idxY));
        end 
    end
    
    linIdx = find(behavTrials.linTrial==1);
    jumpLin = find(diff(linIdx)>1);  
    if isempty(jumpLin)
        lin_begin = linIdx;
    else
        lin_begin = linIdx(1:jumpLin);
    end

    if isempty(jumpLin)
        lin_end = [];
    else
        lin_end = linIdx(jumpLin+1:end);
    end

%     tone6 = find(behavTrials.linTrial(1:(end-1))==0 & behavTrials.correct(1:(end-1)) ==1 & ...
%                 behavTrials.toneGain(1:(end-1)) == 5);
    tone6 = find(behavTrials.linTrial(1:(end-1))==0 & behavTrials.lickLoc(1:(end-1)) == 5);            

    Summary.stdLin(ii,:) = std(position(lin_begin,:),[],1);
    Summary.posLin(ii,:) = mean(position(lin_begin,:),1);

    Summary.stdTone6(ii,:) = std(position(tone6,:),[],1);
    Summary.posTone6(ii,:) = mean(position(tone6,:),1);
    if ~isempty(lin_end)
        Summary.stdLinEnd(ii,:) = std(position(lin_end,:),[],1);
        Summary.posLinEnd(ii,:) = mean(position(lin_end,:),1);
    else
        Summary.stdLinEnd(ii,1:(length(yBins)-1)) = nan;
        Summary.posLinEnd(ii,1:(length(yBins)-1)) = nan;
    end
end