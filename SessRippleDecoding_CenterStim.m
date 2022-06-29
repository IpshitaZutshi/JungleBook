function SessRippleDecoding_CenterStim(varargin)

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
sleepTag = {'all','analog1','analog2','analog3'};

if exist('Summ\rippleDecoding_Center.mat','file') && ~force
    disp('rippleDecoding for center arm already computed! Loading file.');
    load('Summ\rippleDecoding_Center.mat');
else  
    
    for ii = 1:size(allSess,1)
        fprintf(' ** Examining session %3.i of %3.i... \n',ii, size(allSess,1));
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));       
        basepath = pwd;
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
        load([sessionInfo.FileName '.cell_metrics.cellinfo.mat']);
        load([sessionInfo.FileName '.pulses.events.mat']);
        file = dir(('*.session.mat'));
        load(file.name);        
        file = dir(('*.Behavior.mat'));
        load(file(1).name);
        behavior = correctBehav(behavior);
        
        file = dir(('*.SessionPulses.Events.mat'));
        load(file.name);
        file = dir(('*.SessionArmChoice.Events.mat'));
        load(file.name);    
        file = dir(('*.ripples.events.mat'));
        load(file.name);  
        load([sessionInfo.FileName '.MergePoints.events.mat']); 
        
        if ~isempty(dir('*Kilosort*')) &&  ~isempty(dir('summ'))
             spikes = loadSpikes;
        else
            continue;
        end
        
        efields = fieldnames(sessionPulses);    

        for jj = 1:length(efields)
            if exist([basepath '\' efields{jj} '\' efields{jj} '.bayes_decoding.popInfo.mat'],'file') && ~force
                continue
            end
            behavpath = [basepath '\' efields{jj} '\']; 
            saving_folder = 'BayesDecoding\'; mkdir([behavpath saving_folder]);
            
            region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
            target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return                        
            
            if target==2 % If stim was targeted to the return arm, don't analyze
                continue
            end
            
            rewardTS = sessionArmChoice.(efields{jj}).timestamps; 
            endDelay = sessionArmChoice.(efields{jj}).delay.timestamps(2,:)';    
            intTS = behavior.events.intersection(behavior.trials.recordings==jj);
            if length(intTS)>length(endDelay)
                intTS = intTS(2:end);
            end
            for zz = 1:2           
                switch zz
                    case 1
                        startTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0);        
                        endTS = intTS(sessionPulses.(efields{jj}).stim(1:(end-1))==0);
                        events = [startTS';endTS'];
                        tag = 'noStim_center';
                    case 2
                        startTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1);        
                        endTS = intTS(sessionPulses.(efields{jj}).stim(1:(end-1))==1);
                        events = [startTS';endTS'];    
                        tag = 'Stim_center';
                end
                %% decoding for behavior data
                behavior.samplingRate = 1/(median(diff(behavior.timestamps)));
                [beh_dec{zz}, template_beh{zz}]= bz_decodeBehavior_iz(behavpath,spikes,behavior,'dim_names',{'lin'},'intervals',events','dim_bins',15);
                beh_dec{zz}.mse = nanmean((beh_dec{zz}.decodedPosition.lin - beh_dec{zz}.position.lin').^2);
            
                %% Now, find maze ripple events
                %%find all events on the maze
                for fn = 1:length(MergePoints.foldernames)
                    if strcmp(efields{jj},MergePoints.foldernames{fn}) == 1 %If the folder ifs the maze folder
                        intToUse = MergePoints.timestamps(fn,:);
                    end
                end
                riptoUse = InIntervals(ripples.timestamps,intToUse);
                rips.timestamps = ripples.timestamps(riptoUse,:);
                rips.peaks = ripples.peaks(riptoUse,:);            

                padding = .01; % adding padding
                rips.timestamps(:,1) = rips.timestamps(:,1)-padding;
                rips.timestamps(:,2) = rips.timestamps(:,2)+padding;
                
                if length(rips.peaks)>5 %if < 5 ripple events detected

                    ripple_HSE{zz} = find_rippleHSE_iz(behavpath, rips,'save_evt',false,'useVel',false,'intervals',events');
                    
                    if length(ripple_HSE{zz}.idx)>5 % if < 5 high synchrony events detected
                        ripple_dec{zz} = bz_decodeEvents(spikes,beh_dec{zz}.templates{1},ripple_HSE{zz});
                        bayes_allSpikes.awake{zz} = calculateReplayScore(behavpath,saving_folder,ripple_dec{zz},template_beh{zz},behavior,events,tag);

                        bayes_allSpikes.awake{zz}.region = region;
                        bayes_allSpikes.awake{zz}.target = target;                
                        bayes_allSpikes.awake{zz}.behavior = behavior;
                        bayes_allSpikes.awake{zz}.events = events;
                        bayes_allSpikes.awake{zz}.tag =tag;
                        bayes_allSpikes.awake{zz}.beh_dec = beh_dec{zz};
                        bayes_allSpikes.awake{zz}.template_beh = template_beh{zz};
                        bayes_allSpikes.awake{zz}.ripple_HSE = ripple_HSE{zz};
                        bayes_allSpikes.awake{zz}.ripple_dec = ripple_dec{zz};
                        bayes_allSpikes.awake{zz}.fract_sig = length(bayes_allSpikes.awake{zz}.significant_Imax)/size(bayes_allSpikes.awake{zz}.ripple_score_super,2);
                    else
                        bayes_allSpikes.awake{zz}.region = region;
                        bayes_allSpikes.awake{zz}.target = target;                
                        bayes_allSpikes.awake{zz}.behavior = behavior;
                        bayes_allSpikes.awake{zz}.events = events;
                        bayes_allSpikes.awake{zz}.tag =tag;
                        bayes_allSpikes.awake{zz}.beh_dec = beh_dec{zz};
                        bayes_allSpikes.awake{zz}.template_beh = template_beh{zz};
                    end
                else
                    bayes_allSpikes.awake{zz}.region = region;
                    bayes_allSpikes.awake{zz}.target = target;                
                    bayes_allSpikes.awake{zz}.behavior = behavior;
                    bayes_allSpikes.awake{zz}.events = events;
                    bayes_allSpikes.awake{zz}.tag =tag;
                    bayes_allSpikes.awake{zz}.beh_dec = beh_dec{zz};
                    bayes_allSpikes.awake{zz}.template_beh = template_beh{zz};
                end
                clear rips
                %% Now, find sleep ripple events (but only if the template was defined on the morning session)
                if jj == 1
                    %find all intervals in the homecage
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
                    intToUse  = ts_hc;
                    
                    for ss = 1:4 % all sleep, pre,post for each manipulation
                        if ss == 1
                            riptoUse = InIntervals(ripples.timestamps,intToUse);
                            rips{1}.timestamps = ripples.timestamps(riptoUse,:);
                            rips{1}.peaks = ripples.peaks(riptoUse,:);     
                            padding = .01; % adding padding
                            rips{1}.timestamps(:,1) = rips{1}.timestamps(:,1)-padding;
                            rips{1}.timestamps(:,2) = rips{1}.timestamps(:,2)+padding;    
                            range = 1;
                        else
                            range = [1 2];
                            if (ss-1) <2
                                pulTr = (pulses.stimComb==(ss-1));
                            else                                   
                                pulTr = (pulses.stimPerID'==1 & pulses.stimComb==(ss-1));
                            end

                            %Only select the pulses that happened in the home cage,
                            %i.e., which were 5 seconds long
                            homeCagePulse = pulses.intsPeriods(2,:) - pulses.intsPeriods(1,:);
                            homeCagePulseidx = homeCagePulse < 5.05 & homeCagePulse > 4.95;
                            pulTr = pulTr & homeCagePulseidx;
                            eventPulses = pulses.intsPeriods(1,pulTr)';

                            %Generate logicals for ripples in pre, prestim, post, poststim
                            ripple_pre = [];                
                            ripple_post = [];              

                            ripple_pre(1:length(ripples.peaks)) = 0;                
                            ripple_post(1:length(ripples.peaks)) = 0;                

                            for pp = 1:length(ripples.peaks)
                                tempDiff = ripples.peaks(pp) - eventPulses;

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
                              
                            rips{1}.timestamps = ripples.timestamps(logical(ripple_pre),:);
                            rips{1}.peaks = ripples.peaks(logical(ripple_pre),:);   
                            rips{2}.timestamps = ripples.timestamps(logical(ripple_post),:);
                            rips{2}.peaks = ripples.peaks(logical(ripple_post),:);    
                            
                            padding = .01; % adding padding
                            rips{1}.timestamps(:,1) = rips{1}.timestamps(:,1)-padding;
                            rips{1}.timestamps(:,2) = rips{1}.timestamps(:,2)+padding;  
                            rips{2}.timestamps(:,1) = rips{2}.timestamps(:,1)-padding;
                            rips{2}.timestamps(:,2) = rips{2}.timestamps(:,2)+padding;                               
                        end
    
                        for ripRange = 1:length(range)
                            if length(rips{range(ripRange)}.peaks)<5 %if < 5 ripple events detected
                                continue
                            end 
                            
                            ripple_HSE{zz,range(ripRange)} = find_rippleHSE_iz(behavpath, rips{range(ripRange)},'save_evt',false);
                            if length(ripple_HSE{zz,range(ripRange)}.idx)<5 % if < 5 high synchrony events detected
                                continue
                            else
                                ripple_dec{zz,range(ripRange)} = bz_decodeEvents(spikes,beh_dec{zz}.templates{1},ripple_HSE{zz,range(ripRange)});
                                bayes_allSpikes.(sleepTag{ss}){zz,range(ripRange)} = calculateReplayScore(behavpath,saving_folder,ripple_dec{zz,range(ripRange)},template_beh{zz},behavior,events,tag);

                                bayes_allSpikes.(sleepTag{ss}){zz,range(ripRange)}.region = region;
                                bayes_allSpikes.(sleepTag{ss}){zz,range(ripRange)}.target = target;                
                                bayes_allSpikes.(sleepTag{ss}){zz,range(ripRange)}.behavior = behavior;
                                bayes_allSpikes.(sleepTag{ss}){zz,range(ripRange)}.events = events;
                                bayes_allSpikes.(sleepTag{ss}){zz,range(ripRange)}.tag =tag;
                                bayes_allSpikes.(sleepTag{ss}){zz,range(ripRange)}.beh_dec = beh_dec{zz};
                                bayes_allSpikes.(sleepTag{ss}){zz,range(ripRange)}.template_beh = template_beh{zz};
                                bayes_allSpikes.(sleepTag{ss}){zz,range(ripRange)}.ripple_HSE = ripple_HSE{zz,range(ripRange)};
                                bayes_allSpikes.(sleepTag{ss}){zz,range(ripRange)}.ripple_dec = ripple_dec{zz,range(ripRange)};
                                bayes_allSpikes.(sleepTag{ss}){zz,range(ripRange)}.fract_sig = length(bayes_allSpikes.(sleepTag{ss}){zz,range(ripRange)}.significant_Imax)/size(bayes_allSpikes.(sleepTag{ss}){zz,range(ripRange)}.ripple_score_super,2);  
                            end
                        end

                        clear rips ripple_HSE ripple_dec
                    end     
                end
            end
                            
            save([basepath '/' efields{jj} '/' efields{jj} '.bayes_decoding.popInfo.mat'],'bayes_allSpikes');
            close all
            clear bayes_allSpikes
        end                       
    end        
end     
    
end

function bayes_allSpikes = calculateReplayScore(basepath,saving_folder,ripple_dec,template_beh,behavior, intervals, tag)
    %% set options of caculating distance matrix
    options.t_bin = 0.02; %for sleep time bin is 20ms
    options.dmax = 100; %max distance, my maze is 1.1
    options.alpha = 8;  %exponential's decay speed, won't be used w hen option_y = 'linear_bonded'
    options.delta_t = 3; %bond of time.if |t1-t2| > delta_t, d_time = dmax, won't be used w hen option_y = 'linear_bonded'
    options.Nw = 1; % points with highest Nw points are used in weighted distance correlation
    % option of x/y distance metric
    options.option_x = 'linear';
    options.option_y = 'linear';
    track_length = 100; % unit: m

    % a neuron is considered active if it has at least 2 spike in a candidate event
    options.numActCell_criterion = 1;% is was 2, but most of the cells fire once per event!
    % To be analyzied, an event needs contain at least 7 active neurons and it's
    % number of time bins is larger than 4
    options.leastNumCell = 3; % tried 4
    options.leastNumTimeBin = 2; % 4, 3 is 60ms bins

    % correlation metrics to be used
       % (1)'linear_weighted'
       % (2)'weighted_distance'
       % (3)'both'
    options.correlation_metrics = 'both';

    % type of shuffling to be used
       % (1)'neuron'
       % (2)'time'
       % (3)'both'
       % (4)'all' ==> shurrle neuron, shuffle time, shuffle both
    options.shuffle = 'both';

    %% 4.2 Get the decoding result with p value for replay score
    ripple_score_super = decode_replay_w_score_wrapper_iz(options,basepath,saving_folder,ripple_dec,template_beh,behavior,intervals',tag);
    significant_Imax = find_significant_repaly(basepath,saving_folder,ripple_dec,ripple_score_super,tag,'null_type','both');

    bayes_allSpikes.options = options;
    bayes_allSpikes.ripple_score_super = ripple_score_super;
    bayes_allSpikes.significant_Imax = significant_Imax;


end

function behavior = correctBehav(behavior)

    linMap = behavior.position.lin;
    linMapDouble = linMap+170;
    idxMap = (250-linMapDouble)>0; %min([abs(linMap-120),abs(linMapDouble-120)],[],2);
    linMapNew = linMap;
    linMapNew(idxMap) =  linMapDouble(idxMap);           
    linMapNew = linMapNew-75;
    %remove sections when the mouse went backwards
    a = diff(linMapNew);
    idxA= find(a>150);
    idxB= find(a<-150);
    idxtoKeep = ones(length(linMapNew),1);
    %Remove the sections between idxA and idxB
    for pp = 1:length(idxA)
        temp = idxB((idxB-idxA(pp))>0);
        idxEnd = min(temp);
        idxtoKeep(idxA(pp):idxEnd) = 0;
    end
    idxtoKeep = logical(idxtoKeep);
    behavior.position.lin = linMapNew;
    behavior.position.lin(~idxtoKeep) = nan;
end