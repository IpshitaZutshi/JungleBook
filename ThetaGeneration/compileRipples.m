% Compile data across all sessions
% Final ripple analysis 
% By IZ - 2/18/2022

% RippleMasterDetector MUST have been run already!

function compiledRipples = compileRipples(varargin)

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

compiledRipples.ripple_pre = [];
compiledRipples.ripple_prestim = [];
compiledRipples.ripple_post = [];
compiledRipples.ripple_poststim = [];
compiledRipples.peakAmplitude = [];
compiledRipples.duration = [];
compiledRipples.peakFrequency = [];
compiledRipples.psthAnalogEv = [];

compiledRipples.sWmag = [];
compiledRipples.sWprop= [];

compiledRipples.rate_pre = [];
compiledRipples.rate_prestim = [];
compiledRipples.rate_post = [];
compiledRipples.rate_poststim = [];
compiledRipples.rate_postpoststim = [];
compiledRipples.ISI = [];

compiledRipples.rate_poststim1 = [];
compiledRipples.rate_poststim2 = [];
compiledRipples.rate_poststim3 = [];
compiledRipples.rate_poststim4 = [];
compiledRipples.rate_poststim5 = [];
compiledRipples.rateAnalog = [];

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

        %Generate logicals for ripples in pre versus post
        ripple_pre = [];
        ripple_prestim = [];
        ripple_post = [];
        ripple_poststim = [];  
        ripple_postpoststim = [];  
        
        ripple_poststim1 = [];  
        ripple_poststim2 = [];  
        ripple_poststim3 = [];  
        ripple_poststim4 = [];  
        ripple_poststim5 = [];  
        psthAnalogEv= [];
        
        ripple_pre(1:length(ripples.peaks)) = 0;
        ripple_prestim(1:length(ripples.peaks)) = 0;
        ripple_post(1:length(ripples.peaks)) = 0;
        ripple_poststim(1:length(ripples.peaks)) = 0;
        ripple_postpoststim(1:length(ripples.peaks)) = 0;
        
        ripple_poststim1(1:length(ripples.peaks)) = 0;            
        ripple_poststim2(1:length(ripples.peaks)) = 0;      
        ripple_poststim3(1:length(ripples.peaks)) = 0;      
        ripple_poststim4(1:length(ripples.peaks)) = 0;      
        ripple_poststim5(1:length(ripples.peaks)) = 0;      

        for pp = 1:length(ripples.peaks)
            tempDiff = ripples.peaks(pp) - events;

            if min(abs(tempDiff)) <=5 % If a ripple occurs within 5 seconds of a stimulus
               [~,idxmin] =  min(abs(tempDiff));
               if tempDiff(idxmin) > 0
                   ripple_post(pp) = 1;
               elseif tempDiff(idxmin) < 0
                   ripple_pre(pp) = 1;
               end
            elseif min(abs(tempDiff)) >=5 & min(abs(tempDiff)) < 10
                [~,idxmin] =  min(abs(tempDiff));
                if tempDiff(idxmin) > 0
                    ripple_poststim(pp) = 1;
                    if min(abs(tempDiff)) >=5 & min(abs(tempDiff)) < 6
                        ripple_poststim1(pp) = 1;
                    elseif min(abs(tempDiff)) >=6 & min(abs(tempDiff)) < 7
                        ripple_poststim2(pp) = 1;
                    elseif min(abs(tempDiff)) >=7 & min(abs(tempDiff)) < 8
                        ripple_poststim3(pp) = 1;
                    elseif min(abs(tempDiff)) >=8 & min(abs(tempDiff)) < 9
                        ripple_poststim4(pp) = 1;
                    elseif min(abs(tempDiff)) >=9 & min(abs(tempDiff)) < 10
                        ripple_poststim5(pp) = 1;
                    end
               elseif tempDiff(idxmin) < 0
                   ripple_prestim(pp) = 1;
                end
            elseif min(abs(tempDiff)) >=10 & min(abs(tempDiff)) < 15
                [~,idxmin] =  min(abs(tempDiff));
                if tempDiff(idxmin) > 0
                    ripple_postpoststim(pp) = 1;
                end
            else
                continue
            end
        end

        ripple_pre = logical(ripple_pre);
        ripple_post = logical(ripple_post);
        ripple_poststim = logical(ripple_poststim);
        ripple_postpoststim = logical(ripple_postpoststim);
        ripple_prestim = logical(ripple_prestim);
        
        ripple_poststim1 = logical(ripple_poststim1);
        ripple_poststim2 = logical(ripple_poststim2);
        ripple_poststim3 = logical(ripple_poststim3);
        ripple_poststim4 = logical(ripple_poststim4);
        ripple_poststim5 = logical(ripple_poststim5);
        
        psthAnalogEv(1:length(ripples.peaks)) = i;
        
        compiledRipples.ripple_pre = [compiledRipples.ripple_pre ripple_pre];
        compiledRipples.ripple_prestim = [compiledRipples.ripple_prestim ripple_prestim];
        compiledRipples.ripple_post = [compiledRipples.ripple_post ripple_post];
        compiledRipples.ripple_poststim = [compiledRipples.ripple_poststim ripple_poststim];
        compiledRipples.peakAmplitude = [compiledRipples.peakAmplitude ripples.data.peakAmplitude'];
        compiledRipples.duration = [compiledRipples.duration ripples.data.duration'];
        compiledRipples.peakFrequency = [compiledRipples.peakFrequency ripples.data.peakFrequency'];
        compiledRipples.sWmag = [compiledRipples.sWmag SW.peakZScore'];
        compiledRipples.ISI = [compiledRipples.ISI ripples.peaks'];
        
        %Find proportion of ripples with SWs
        swmag_pre = SW.peakZScore(ripple_pre);
        swmag_post = SW.peakZScore(ripple_post);
        swmag_prestim = SW.peakZScore(ripple_prestim);
        swmag_poststim = SW.peakZScore(ripple_poststim);
        
        compiledRipples.sWprop = [compiledRipples.sWprop; sum(~isnan(swmag_prestim))./sum(ripple_prestim) ...
            sum(~isnan(swmag_pre))./sum(ripple_pre) sum(~isnan(swmag_post))./sum(ripple_post) sum(~isnan(swmag_poststim))./sum(ripple_poststim)];
        compiledRipples.psthAnalogEv = [compiledRipples.psthAnalogEv psthAnalogEv];
        compiledRipples.rate_pre = [compiledRipples.rate_pre sum(ripple_pre)./(numel(events)*5)];
        compiledRipples.rate_prestim = [compiledRipples.rate_prestim sum(ripple_prestim)./(numel(events)*5)];
        compiledRipples.rate_post = [compiledRipples.rate_post sum(ripple_post)./(numel(events)*5)];
        compiledRipples.rate_poststim = [compiledRipples.rate_poststim sum(ripple_poststim)./(numel(events)*5)];
        compiledRipples.rate_postpoststim = [compiledRipples.rate_postpoststim sum(ripple_postpoststim)./(numel(events)*5)];
        
        compiledRipples.rate_poststim1 = [compiledRipples.rate_poststim1 sum(ripple_poststim1)./(numel(events))];
        compiledRipples.rate_poststim2 = [compiledRipples.rate_poststim2 sum(ripple_poststim2)./(numel(events))];
        compiledRipples.rate_poststim3 = [compiledRipples.rate_poststim3 sum(ripple_poststim3)./(numel(events))];
        compiledRipples.rate_poststim4 = [compiledRipples.rate_poststim4 sum(ripple_poststim4)./(numel(events))];
        compiledRipples.rate_poststim5 = [compiledRipples.rate_poststim5 sum(ripple_poststim5)./(numel(events))];
        compiledRipples.rateAnalog = [compiledRipples.rateAnalog i];
    end
end

save([expPath '\Summ\' 'compiledRipples.mat'], 'compiledRipples');

if plotfig
    figure
    for i = 1:(numAnalog+1)
        idxAnalog = (compiledRipples.psthAnalogEv ==i);
        subplot((numAnalog+1),4,4*(i-1)+1) 
        bar([sum(compiledRipples.ripple_prestim(idxAnalog)) sum(compiledRipples.ripple_pre(idxAnalog)) sum(compiledRipples.ripple_post(idxAnalog)) sum(compiledRipples.ripple_poststim(idxAnalog))]);
        ylabel('Number of ripples')

        dataAll = [];
        subplot((numAnalog+1),4,4*(i-1)+2)
        title('Ripple power')
        dataAll.prestim = compiledRipples.peakAmplitude(compiledRipples.ripple_prestim & idxAnalog);
        dataAll.pre = compiledRipples.peakAmplitude(compiledRipples.ripple_pre & idxAnalog);
        dataAll.post = compiledRipples.peakAmplitude(compiledRipples.ripple_post & idxAnalog);
        dataAll.poststim = compiledRipples.peakAmplitude(compiledRipples.ripple_poststim & idxAnalog);
        nhist(dataAll,'samebins','smooth','pdf','median','median','sem','linewidth',1.5)
        ylabel('Proportion')
        xlabel('Peak ripple power')

        dataAll = [];
        subplot((numAnalog+1),4,4*(i-1)+3)
        title('Ripple duration')
        dataAll.prestim = compiledRipples.duration(compiledRipples.ripple_prestim & idxAnalog);
        dataAll.pre = compiledRipples.duration(compiledRipples.ripple_pre & idxAnalog);
        dataAll.post = compiledRipples.duration(compiledRipples.ripple_post & idxAnalog);
        dataAll.poststim = compiledRipples.duration(compiledRipples.ripple_poststim & idxAnalog);
        nhist(dataAll,'samebins','smooth','pdf','median','median','sem','linewidth',1.5)
        ylabel('Proportion')
        xlabel('Ripple duration')

        dataAll = [];
        subplot((numAnalog+1),4,4*(i-1)+4)
        title('Ripple frequency')
        dataAll.prestim = compiledRipples.peakFrequency(compiledRipples.ripple_prestim & idxAnalog);
        dataAll.pre = compiledRipples.peakFrequency(compiledRipples.ripple_pre & idxAnalog);
        dataAll.post = compiledRipples.peakFrequency(compiledRipples.ripple_post & idxAnalog);
        dataAll.poststim = compiledRipples.peakFrequency(compiledRipples.ripple_poststim & idxAnalog);
        nhist(dataAll,'samebins','smooth','pdf','median','median','sem','linewidth',1.5)
        ylabel('Proportion')
        xlabel('Peak ripple frequency')

    end
end
