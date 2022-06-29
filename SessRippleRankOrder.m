function SessRippleRankOrder(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',true,@islogical);
addParameter(p,'dSample',true,@islogical);
parse(p,varargin{:});

expPath = p.Results.expPath;
saveMat = p.Results.saveMat;
force = p.Results.force;
dSample = p.Results.dSample;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');

if dSample
    strname = 'Summ\RippleRankOrder_DS.mat';
else
    strname = 'Summ\RippleRankOrder.mat';
end

if exist(strname,'file') && ~force
    disp('Ripple rank order already computed! Loading file.');
    load(strname);
else
    rippleRankOrder = zeros(size(allSess,1),3,3);
    numRip = zeros(size(allSess,1),3,3);
    
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
        load([sessionInfo.FileName '.spikes.cellinfo.mat']);        
        
        if ~exist([sessionInfo.FileName '.spkEventTimes.mat'])
            spkEventTimes = bz_getRipSpikes;
        else
            load([sessionInfo.FileName '.spkEventTimes.mat'])
        end
        
        %% calculate ripple rank order
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

            %Generate logicals for ripples in pre versus post
            ripple_pre = [];                
            ripple_post = [];     
            ripple_poststim = [];

            ripple_pre(1:length(ripples.peaks)) = 0;                
            ripple_post(1:length(ripples.peaks)) = 0;                    
            ripple_poststim(1:length(ripples.peaks)) = 0;  

            for pp = 1:length(ripples.peaks)
                tempDiff = ripples.peaks(pp) - events;

                if min(abs(tempDiff)) <=5 % If a ripple occurs within 5 seconds of a stimulus
                   [~,idxmin] =  min(abs(tempDiff));
                   if tempDiff(idxmin) > 0
                       ripple_post(pp) = 1;
                   elseif tempDiff(idxmin) < 0
                       ripple_pre(pp) = 1;
                   end  
                 elseif min(abs(tempDiff)) >5 && min(abs(tempDiff)) <=10% If a ripple occurs within 5 seconds of a stimulus
                   [~,idxmin] =  min(abs(tempDiff));
                   if tempDiff(idxmin) > 0
                       ripple_poststim(pp) = 1;
                   end                     
                else
                    continue
                end
            end   

            if sum(ripple_post) < 1 || sum(ripple_pre)< 1|| sum(ripple_poststim)< 1
                rippleRankOrder(ii,1,rr) = nan;
                rippleRankOrder(ii,2,rr) = nan;
                rippleRankOrder(ii,3,rr) = nan;
            else
                if dSample
                    numPost = sum(ripple_post);
                    numPre = sum(ripple_pre);
                    numPostStim = sum(ripple_poststim);
                    if numPost<numPre %randomly subselect ripples from Pre
                        ripIds = find(ripple_pre);
                        pickIds = randsample(numPre,(numPre-numPost));
                        ripple_pre(ripIds(pickIds)) = 0;

                        ripIds = find(ripple_poststim);
                        pickIds = randsample(numPostStim,(numPostStim-numPost));
                        ripple_poststim(ripIds(pickIds)) = 0;                        
                    else                        
                        ripIds = find(ripple_post);
                        pickIds = randsample(numPost,(numPost-numPre));
                        ripple_post(ripIds(pickIds)) = 0;
                        
                        ripIds = find(ripple_poststim);
                        pickIds = randsample(numPostStim,(numPostStim-numPre));
                        ripple_poststim(ripIds(pickIds)) = 0;                            
                    end
                end
                rankStats_pre = bz_RankOrder('spkEventTimes',spkEventTimes,'eventIDs',ripple_pre,'saveMat',false,'doPlot',false);
                rippleRankOrder(ii,1,rr) = rankStats_pre.corrMean;            
                rankStats_post = bz_RankOrder('spkEventTimes',spkEventTimes,'eventIDs',ripple_post,'saveMat',false,'doPlot',false);
                rippleRankOrder(ii,2,rr) = rankStats_post.corrMean;                       
                rankStats_poststim = bz_RankOrder('spkEventTimes',spkEventTimes,'eventIDs',ripple_poststim,'saveMat',false,'doPlot',false);
                rippleRankOrder(ii,3,rr) = rankStats_poststim.corrMean;     
                numRip(ii,1,rr) = sum(ripple_pre);
                numRip(ii,2,rr) = sum(ripple_post);
                numRip(ii,3,rr) = sum(ripple_poststim);
            end
        end
    end
    
    if saveMat
        if dSample
            save([expPath '\Summ\' 'RippleRankOrder_DS.mat'], 'rippleRankOrder','numRip','-v7.3');
        else
            save([expPath '\Summ\' 'RippleRankOrder.mat'], 'rippleRankOrder','numRip','-v7.3');
        end
    end    
end
    
end

