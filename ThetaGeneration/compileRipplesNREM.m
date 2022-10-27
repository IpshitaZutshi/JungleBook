% Compile data across all sessions
% Final ripple analysis 
% By IZ - 2/18/2022

% RippleMasterDetector MUST have been run already!

function compiledRipples = compileRipplesNREM(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'analogEv',64,@isnumeric);
addParameter(p,'numAnalog',2,@isnumeric);
addParameter(p,'minEvents',5,@isnumeric); % each session should have atleast 50 ripples to be included in the stats
addParameter(p,'plotfig',false,@islogical);
parse(p,varargin{:});

expPath = p.Results.expPath;
analogEv = p.Results.analogEv;
numAnalog = p.Results.numAnalog;
minEvents = p.Results.minEvents;
plotfig = p.Results.plotfig;

%%
if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end
allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');

if ~isempty(analogEv)
    for ii = 1:numAnalog
        analogCh(ii) = (analogEv-1)+ii;
    end
end

for ww = 1:2
    compiledRipples{ww}.rate_pre = [];
    compiledRipples{ww}.rate_post = [];
    compiledRipples{ww}.rateAnalog = [];
end

%% Start collecting data
for ii = 1:size(allSess,1)

    cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
    [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);    

    if exist([sessionInfo.FileName '.ripples.events.mat'],'file') 
       load([sessionInfo.FileName '.ripples.events.mat']);
    else
       disp('First calculate .ripple file! Skipping');
       return
    end
    
    if exist([sessionInfo.FileName '.sharpwaves.events.mat'],'file') 
       load([sessionInfo.FileName '.sharpwaves.events.mat']);
    else
       disp('No sharp waves associated with this session');
    end

    
    if exist([sessionInfo.FileName '.SleepState.states.mat'],'file') 
       load([sessionInfo.FileName '.SleepState.states.mat']);
    else
       disp('No sleep states associated with this session');
    end
    
    % Load pulses
    disp('Getting analog-in inputs...');
    [pulses] = bz_getAnalogPulsesSine('analogCh',analogCh);

    %Skip session if too few ripples
    if length(ripples.peaks)<minEvents
        continue
    end
    
    for i = 1:(numAnalog+1)
        if exist('pulses')
            if i<=numAnalog
                pulTr = (pulses.stimComb==i);
            else
                pulTr = (pulses.stimPerID'==1 & pulses.stimComb==i);
            end
        end
        
        %Only select the pulses that happened in the home cage,
        %i.e., which were 5 seconds long
        homeCagePulse = pulses.intsPeriods(2,:) - pulses.intsPeriods(1,:);
        homeCagePulseidx = homeCagePulse < 5.05 & homeCagePulse > 4.95;
        pulTr = pulTr & homeCagePulseidx;
        events = pulses.intsPeriods(1,pulTr)';
        
        % Only select pulses during NREM OR only select pulses during Wake
        eventsState{1} = InIntervals(events,SleepState.ints.WAKEstate);
        eventsState{2} = InIntervals(events,SleepState.ints.NREMstate);
        
        for ww = 1:2
            %Generate logicals for ripples in pre versus post
            ripple_pre = [];
            ripple_post = [];

            psthAnalogEv= [];

            ripple_pre(1:length(ripples.peaks)) = 0;
            ripple_post(1:length(ripples.peaks)) = 0;

            if sum(eventsState{ww})>0
                for pp = 1:length(ripples.peaks)
                    tempDiff = ripples.peaks(pp) - events(eventsState{ww});

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
            end

            ripple_pre = logical(ripple_pre);
            ripple_post = logical(ripple_post);

            compiledRipples{ww}.rate_pre = [compiledRipples{ww}.rate_pre sum(ripple_pre)./(sum(eventsState{ww})*5)];
            compiledRipples{ww}.rate_post = [compiledRipples{ww}.rate_post sum(ripple_post)./(sum(eventsState{ww})*5)];
            compiledRipples{ww}.rateAnalog = [compiledRipples{ww}.rateAnalog i];
        end
    end
end

save([expPath '\Summ\' 'compiledRipplesNREM.mat'], 'compiledRipples');
end
