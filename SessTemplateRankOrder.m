function SessTemplateRankOrder(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',true,@islogical);
parse(p,varargin{:});

expPath = p.Results.expPath;
saveMat = p.Results.saveMat;
force = p.Results.force;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');

if exist('Summ\templateRankOrder','file') && ~force
    disp('rank order already calculated');
    load('Summ\templateRankOrder.mat');
else  
    for rr = 1:3
        templateRankOrder.analysis1{rr} = [];
        templateRankOrder.analysis2{rr,1} = [];% No stim template
        templateRankOrder.analysis2{rr,2} = [];% Stim template
    end
    
    
    for ii = 1:size(allSess,1)
        fprintf(' ** Examining session %3.i of %3.i... \n',ii, size(allSess,1));
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));       
        basepath = pwd;
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
        load([sessionInfo.FileName '.cell_metrics.cellinfo.mat']);
        file = dir(('*.session.mat'));
        load(file.name);        
        
        %Load ripples
        load([sessionInfo.FileName '.ripples.events.mat']);
        load([sessionInfo.FileName '.pulses.events.mat']);
        
        if ~exist([sessionInfo.FileName '.spkEventTimes.mat'])
            spkEventTimes = bz_getRipSpikes;
        else
            load([sessionInfo.FileName '.spkEventTimes.mat'])
        end
        
        file = dir(('*.SessionPulses.Events.mat'));
        load(file.name);
        
        load([sessionInfo.FileName '.MergePoints.events.mat']); 
        
        if ~isempty(dir('*Kilosort*')) &&  ~isempty(dir('summ'))
             spikes = loadSpikes;
        else
            continue;
        end
        
        efields = fieldnames(sessionPulses);    

        for jj = 1:length(efields)
            
            if exist([efields{jj} '\' efields{jj} '.templateRankOrder.mat']) && ~force
                load([efields{jj} '\' efields{jj} '.templateRankOrder.mat']);
            else
                region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
                target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return                        

                if target==2 % If stim was targeted to the return arm, don't analyze
                    continue
                end

                load([efields{jj} '\' efields{jj} '.placeFieldTemplate.mat'])

                %% Make two templates, for the center arm, no stim versus stim
                templateExt = placeFieldTemplate.Peak;
                %exclude return arm section
                for kk = 1:length(templateExt)
                    templateExt{kk} = templateExt{kk}(templateExt{kk}(:,1)>60,:);
                    templateExt{kk}(:,1) = templateExt{kk}(:,1)-60;
                end

                %% Now, find home cage ripple events             
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
                intToUse  = ts_hc(1,:);

                riptoUse = InIntervals(ripples.timestamps,intToUse);

                spkEventTimes_temp = spkEventTimes;
                spkEventTimes_temp.EventRel = spkEventTimes.EventRel(riptoUse);
                spkEventTimes_temp.UnitEventRel = spkEventTimes.UnitEventRel(riptoUse);

                rankStats = bz_RankOrder_IZ('spkEventTimes',spkEventTimes_temp,'templateType','Peak','templateExt',templateExt,'saveMat',false,'doPlot',false);
                corrbase = mean([rankStats.corrMean(1) rankStats.corrMean(3)]);
                corrstim = mean([rankStats.corrMean(2) rankStats.corrMean(4)]);

                %% Now, have the same two templates, for the center arm, no stim versus stim, but look at manipulated ripples
                rankOrderStim_pre = [];
                rankOrderStim_stim = [];
                rankOrderBase_pre = [];
                rankOrderBase_stim = [];

                for rr  = 1:3        
                    if rr <= 2
                        pulTr = (pulses.stimComb==rr);
                    else
                        pulTr = (pulses.stimPerID'==1 & pulses.stimComb==rr);
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

                    ripple_pre(1:length(ripples.peaks)) = 0;                
                    ripple_post(1:length(ripples.peaks)) = 0;                    

                    for pp = 1:length(ripples.peaks)
                        tempDiff = ripples.peaks(pp) - events;

                        if min(abs(tempDiff)) <=5 % If a ripple occurs within 5 seconds of a stimulus
                           [~,idxmin] =  min(abs(tempDiff));
                           if tempDiff(idxmin) > 0
                               ripple_post(pp) = 1;
                           elseif tempDiff(idxmin) < 0
                               ripple_pre(pp) = 1;
                           end                    
                        else
                            continue
                        end
                    end   


                    if sum(ripple_post)>5
                        spkEventTimes_temp = spkEventTimes;
                        spkEventTimes_temp.EventRel = spkEventTimes.EventRel(logical(ripple_post));
                        spkEventTimes_temp.UnitEventRel = spkEventTimes.UnitEventRel(logical(ripple_post));

                        rankStats = bz_RankOrder_IZ('spkEventTimes',spkEventTimes_temp,'templateType','Peak','templateExt',templateExt,'saveMat',false,'doPlot',false);

                        corrbase = mean([rankStats.corrMean(1),rankStats.corrMean(3)]);
                        corrstim = mean([rankStats.corrMean(2),rankStats.corrMean(4)]);

                        rankOrderStim_stim = [rankOrderStim_stim corrstim];
                        rankOrderBase_stim = [rankOrderBase_stim corrbase];                
                    else
                        rankOrderStim_stim = [rankOrderStim_stim nan];
                        rankOrderBase_stim = [rankOrderBase_stim nan];
                    end

                    %Pre rippl

                    if sum(ripple_pre)>5
                        spkEventTimes_temp = spkEventTimes;
                        spkEventTimes_temp.EventRel = spkEventTimes.EventRel(logical(ripple_pre));
                        spkEventTimes_temp.UnitEventRel = spkEventTimes.UnitEventRel(logical(ripple_pre));

                        rankStats = bz_RankOrder_IZ('spkEventTimes',spkEventTimes_temp,'templateType','Peak','templateExt',templateExt,'saveMat',false,'doPlot',false);
                        corrbase = mean([rankStats.corrMean(1) rankStats.corrMean(3)]);
                        corrstim = mean([rankStats.corrMean(2) rankStats.corrMean(4)]);

                        rankOrderStim_pre = [rankOrderStim_pre corrstim];
                        rankOrderBase_pre = [rankOrderBase_pre corrbase];              
                    else
                        rankOrderStim_pre = [rankOrderStim_pre nan];
                        rankOrderBase_pre = [rankOrderBase_pre nan];
                    end
                end

                if saveMat
                    save([efields{jj} '\'  efields{jj} '.templateRankOrder.mat'], 'corrbase','corrstim','rankOrderBase_pre',...
                        'rankOrderBase_stim', 'rankOrderStim_pre', 'rankOrderStim_stim', '-v7.3');
                end 
            end
            
            templateRankOrder.analysis1{region} = [templateRankOrder.analysis1{region}; [corrbase corrstim]];
            templateRankOrder.analysis2{region,1} = [templateRankOrder.analysis2{region,1}; rankOrderBase_pre rankOrderBase_stim];% No stim template
            templateRankOrder.analysis2{region,2} = [templateRankOrder.analysis2{region,2}; rankOrderStim_pre rankOrderStim_stim];% stim template

        end
    end  
    
    if saveMat
        save([expPath '\Summ\' 'templateRankOrder.mat'], 'templateRankOrder','-v7.3');
    end    
end     
    
end

