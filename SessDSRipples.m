function SessDSRipples(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'makePlot',true,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',true,@islogical);
parse(p,varargin{:});

expPath = p.Results.expPath;
makePlot = p.Results.makePlot;
saveMat = p.Results.saveMat;
force = p.Results.force;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');

if exist(strcat('Summ\DSRipples.mat'),'file') && ~force 
    disp('DS Ripples already computed! Loading file.');
    load(strcat('Summ\DSRipples.mat'));
else
    for rr = 1:3
        DSRipples{rr,1} = [];
        DSRipples{rr,2} = [];
    end
    
    for ii = 1:size(allSess,1)
        fprintf(' ** Examining session %3.i of %3.i... \n',ii, size(allSess,1));
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
        load([sessionInfo.FileName '.cell_metrics.cellinfo.mat']);
        file = dir(('*.session.mat'));
        load(file.name);        
        
        %Load ripples
        load([sessionInfo.FileName '.ripples.events.mat']);
        load([sessionInfo.FileName '.pulses.events.mat']);
        
        %Load sleep states
        load([sessionInfo.FileName '.SleepState.states.mat']);
        load([sessionInfo.FileName '.spikes.cellinfo.mat']);        
        
        load([sessionInfo.FileName '.SessionPulses.Events.mat']); 
        load([sessionInfo.FileName '.MergePoints.events.mat']); 
        
        %Find home cage intervals
        efields = fields(sessionPulses);
        ts_hc = [];
        for fn = 1:length(MergePoints.foldernames)
            flag = 1;
            for fx = 1:length(efields)
                if strcmp(efields{fx},MergePoints.foldernames{fn}) == 1 %If folder if a maze folder
                    flag = 0;
                end
            end
            if flag
                ts_hc = [ts_hc;MergePoints.timestamps(fn,:)];
            end
        end
        
        %Build spike matrix for sleep periods
        spikeMat = [];
        timeMat = [];   
        events = SleepState.ints.NREMstate;
        % Build matrix, concatenate sleep epochs
        for trials = 1:size(events,1)         
            spkMat = [];
            spkMat_ts = [];
            %check if interval lies within homecage
            status = InIntervals(events(trials,1),ts_hc);
            if status == 0
                continue
            end
            for unit = 1:length(spikes.times)
                if strcmp(cell_metrics.brainRegion(unit),'CA1')~=1
                    continue
                else 
                   spkData = bz_SpktToSpkmat(spikes.times(unit),'dt', .020, 'win',events(trials,:));
                   spkMat = [spkMat spkData.data];                                        
                   spkMat_ts = [spkData.timestamps];
                end                               
            end
           spikeMat = [spikeMat;spkMat];
           timeMat = [timeMat;spkMat_ts];            
        end
        
        %Find downstates
        A = nanmean(spikeMat,2);
        TF = islocalmin(A,'flatselection','center');
        
        % Assign logicals for the timeMat         
        timeLog = [];               
        timeLog(1,1:length(timeMat)) = 0;
        timeLog(2,1:length(timeMat)) = 0;
        timeLog(3,1:length(timeMat)) = 0;
        
        for rr = 1:3              
            if exist('pulses')
                if rr <= 2
                    pulTr = (pulses.stimComb==rr);
                else
                    pulTr = (pulses.stimPerID'==1 & pulses.stimComb==rr);
                end
            end

            %Only select the pulses that happened in the home cage,
            %i.e., which were 5 seconds long
            homeCagePulse = pulses.intsPeriods(2,:) - pulses.intsPeriods(1,:);
            homeCagePulseidx = homeCagePulse < 5.05 & homeCagePulse > 4.95;
            pulTr = pulTr & homeCagePulseidx;
            events = pulses.intsPeriods(1,pulTr)';
            
            for trials = 1:length(events)
                idx = find(timeMat>= (events(trials)) & timeMat< (events(trials)+5));
                if ~isempty(idx)
                    timeLog(rr,idx) = 2;
                end
                idx = find(timeMat>= (events(trials)-5) & timeMat< (events(trials)));
                if ~isempty(idx)
                    timeLog(rr,idx) = 1;
                end              
            end
            
            %Find "down State" timepoints within timelog 1 or 2
            preIdx = find(timeLog(rr,:)==1 & TF');
            postIdx = find(timeLog(rr,:)==2 & TF');
            
            if isempty(preIdx) || isempty(postIdx)
                continue
            end
            [stccg_pre] = CCG({ripples.peaks timeMat(preIdx)},[],'binSize',0.05,'duration',1);            
            psth_pre = stccg_pre(:,2,1);

            [stccg_post] = CCG({ripples.peaks timeMat(postIdx)},[],'binSize',0.05,'duration',1);            
            psth_post = stccg_post(:,2,1);            
            
            if sum(psth_pre) >100 % Kind of arbitrary, but ensure theres enough sleep
                DSRipples{rr,1}  = [DSRipples{rr,1}; psth_pre'./sum(psth_pre)];
                DSRipples{rr,2}  = [DSRipples{rr,2}; psth_post'./sum(psth_post)];
            end
        end
    end
    
    if saveMat
        save([expPath '\Summ\' 'DSRipples.mat'], 'DSRipples','-v7.3');
    end
    
    if makePlot
        figure
        avgerr = nanmean(DSRipples{2,1});
        stderr = nanstd(DSRipples{2,1})./sqrt(size(DSRipples{2,1},2));
        fill([(1:1:length(avgerr))'; (length(avgerr):-1:1)'],[avgerr'-stderr';flipud(avgerr'+stderr')],'k','linestyle','none','FaceAlpha',0.2);
        hold on
        plot(avgerr,'Color','k','LineWidth',1.5)
        
        avgerr = nanmean(DSRipples{2,2});
        stderr = nanstd(DSRipples{2,2})./sqrt(size(DSRipples{2,2},2));
        fill([(1:1:length(avgerr))'; (length(avgerr):-1:1)'],[avgerr'-stderr';flipud(avgerr'+stderr')],'b','linestyle','none','FaceAlpha',0.2);
        hold on
        plot(avgerr,'Color','b','LineWidth',1.5)        
    end
end
    
end

