
function [CCGAll,CCGAllXtile,synchInfo,rippleAllDOWNs,rippleAllUPs,rippleAll,DOWNInfo]=ContinuousDataByDOWNXtile(basePaths,MUAIn,SlowWavesIn,trigger,timelag,varargin)
% Description: Collects specified continuous data (eg MUA) across sessions
% in basePaths (and if provided d/l's ripples from corresponding
% alternative filePath to session w ripples file), and plots wrt DOWN states (or intervals specified)
% Note: if MUA and SlowWaves = strings, MUA and SlowWaves must be in same basePath.

% Input:
% basePaths - can be one or mulitiple basePaths (cell array of strings).
%       If multiple, will combine needed data across basePaths. Assumes
%       MUA, SlowWaves, and optional ripples structs for each basePath are
%       in the location specified by basePath, not subfolders within
%       basePath.
% MUAIn - can be a string that refers to file in basePath eg
%       'MUA.lfp.dat', or a struct named MUA with fields .data,.timestamps, and
%       .samplingRate (eg in this case MUA.data could be pop. rate)
% SlowWaves - can be a string that refers to file in basePath eg
%       .SlowWavesMUA.events.mat or a SlowWaves struct with following format:
%       SlowWaves.ints.DOWN, SlowWaves.ints.UP (must have both)
%       if string, assumes name leads to struct 'SlowWaves' with field
%       SlowWaves.ints.DOWN
% trigger - what to center the data around wrt DOWN state; 'DOWNtoUP', 'UPtoDOWN','DOWN' (temporal midpoint)
% timelag - timelag interval eg [1 1] for -1 to 1 for CCG
%
% Optional:
% 'restrictIn' - can be strings: 'NREM', 'all', 'WAKE' or be list of
%               intervals to restrict DOWNs to. If strings, assume has
%               variable *.SleepStates.state.mat in basePath
% 'xtile' - also breaks up data into xtiles; for now - can only enter 5
%           (come back to this).
% 'lowerBoundCutoff' - cutoff DOWN states of this duration and below (in
%           sec)
% 'upperBoundCutoff' - cutoff DOWN states of this duration and above (in
%           sec)
% 'ripplesIn' - string for name (or partial name) of ripples file
% 'sharpWavesIn' - string for SharpWaves.events.mat file
% 'altPathSessionGeneral' - alternative filepath to ripples/general session
%       variables (ripples,state,PSS,etc). This is useful if filepath provided above is
%       region specific (particular to Stringer dataset; phase out later).
% 'SWRFeatureSelection' - Which feature to look at as fx of DOWN duration;
%       leave empty if just want ripple raster surrounding DOWN trigger of interest;
%       options: SWamp,SWRamp,SWRdur; if want to breakdown by SWR feature
% 'savePathMouse' - path to save figures in; for now not saving figures
%   (add functionality)
% 'includeState' - include state values for each ripple, specify which (for
%                   now only valid input is 'PSS'
%
% Output:
% CCGAll - struct that contains output specified for each DOWN sorted by duration DOWN
% CCGAllXtile - same as CCGAll, but split up by DOWN duration Xtile
% synchInfo - for each DOWN, timing on/off, duration, synchrony post, dur
%           post UP, duration pre UP, PSS
% rippleAllDOWNs - dur, amp, sw mag, timing wrt surrounding UP/DOWN, PSS
% (ripples during DOWN)
% rippleAllUP
% DOWNInfo - all DOWNs sorted by duration. RippleIndex refers to the
% ripples in that [tlag] window surround the specified DOWN; have to find
% the index in rippleAllUP to find which ripple category it belongs to
% MUAHPCin - can be a string that refers to file in basePath eg
%       'MUAHPC.lfp.dat', or a struct named MUA with fields .data,.timestamps, and
%       .samplingRate (eg in this case MUA.data could be pop. rate)


% To do
% - change so MUA can be any continuous variable
% - add in durUP for synchrony calculations
% - remove the close all
% - if don't include xtiles will fail in figure saving section, fix later
% - update state variable to be combine-able across basePaths (for now
% isn't)
% - confirm UP and DOWN ripple state stuff is combining well across
% basePaths
% - add in option to use SWR peak/onset; for now just using peak
% - add back in infraslow in a bit

%% Optional inputs
p = inputParser;
addParameter(p,'restrictIn',[])
addParameter(p,'xtile',[])
addParameter(p,'lowerBoundCutoff',[])
addParameter(p,'upperBoundCutoff',[])
addParameter(p,'ripplesIn',[])
addParameter(p,'altPathSessionGeneral',[])
addParameter(p,'SWRFeatureSelection',[])
addParameter(p,'sharpWavesIn',[])
addParameter(p,'savePathMouse',[])
addParameter(p,'includeState',[])
addParameter(p,'MUAHPCin',[])

parse(p,varargin{:})
restrictIn = p.Results.restrictIn;
xtile = p.Results.xtile;
lowerBoundCutoff = p.Results.lowerBoundCutoff;
upperBoundCutoff = p.Results.upperBoundCutoff;
ripplesIn = p.Results.ripplesIn;
altPathSessionGeneral = p.Results.altPathSessionGeneral;
SWRFeatureSelection = p.Results.SWRFeatureSelection;
sharpWavesIn = p.Results.sharpWavesIn;
savePathMouse = p.Results.savePathMouse;
includeState = p.Results.includeState;
MUAHPCin = p.Results.MUAHPCin;

% Organize inputs
if ischar(basePaths)
    basePathString=basePaths;
    basePaths={};
    basePaths{1} = basePathString; %if single string, put in cell array
end

close all
%% Loop through each basePath
for i = 1:length(basePaths)
    
    % Vars
    basePath = basePaths{i};
    cd(basePath)
    baseName = bz_BasenameFromBasepath(basePath);
    
    % Deal with strings (writing more flexibly bc may not have buzcode
    % formatted data naming, ie baseName.*)
    % MUAHPC
    if ischar(MUAHPCin)
        temp=dir(fullfile(basePath,['*' MUAHPCin]));
        load(temp.name);
        MUAHPC = MUA;
    end
    % MUA
    if ischar(MUAIn)
        temp=dir(fullfile(basePath,['*' MUAIn]));
        load(temp.name);
    end
    % SlowWaves
    if ischar(SlowWavesIn)
        temp=dir(fullfile(basePath,['*' SlowWavesIn]));
        load(temp.name);
        SlowWaves_noneCut = SlowWaves;
    else
        SlowWaves =  SlowWavesIn;
        SlowWaves_noneCut = SlowWaves;
    end
    % Ripples
    if ischar(ripplesIn)
        if ~isempty(altPathSessionGeneral)
            temp=dir(fullfile(altPathSessionGeneral,['*' ripplesIn]));
        else
            temp=dir(fullfile(basePath,['*' ripplesIn]));
        end
        load(fullfile(temp.folder,temp.name));
    end
    % SharpWaves
    if ischar(sharpWavesIn)
        if ~isempty(altPathSessionGeneral)
            temp=dir(fullfile(altPathSessionGeneral,['*' sharpWavesIn]));
        else
            temp=dir(fullfile(basePath,['*' sharpWavesIn]));
        end
        load(fullfile(temp.folder,temp.name));
        % Sharp waves
        SW.magnitude = SW.magnitude.*-1;
        SW.timestamps = SW.timePeak;
    end
    % Sleep state
    if ischar(restrictIn)
        if ~isempty(altPathSessionGeneral)
            temp=dir(fullfile(altPathSessionGeneral,'*.SleepState.states.mat'));
        else
            temp=dir(fullfile(basePath,'*.SleepState.states.mat'));
        end
        if isempty(temp)
            disp('Selected to restrict but have no SleepState.states file. Use buzcode state editor to state score.')
        else
            if ~isempty(altPathSessionGeneral)
                load(fullfile(altPathSessionGeneral,temp.name))
            else
                load(fullfile(basePath,temp.name))
            end
            if strcmp(restrictIn,'NREM')
                restrict=SleepState.ints.NREMstate;
            elseif strcmp(restrictIn,'WAKE')
                restrict=SleepState.ints.WAKEstate;
            elseif strcmp(restrictIn,'all')
                restrict = [0 inf];
            end
        end
    else
        restrictIn = 'all';
    end
    % PSS
    if ischar(includeState)
        if ~isempty(altPathSessionGeneral)
            temp=dir(fullfile(altPathSessionGeneral,'*.SleepState.states.mat'));
            %temp2 = dir(fullfile(altPathSessionGeneral,'BloodModAnalysis','infraslowROI.mat'));
        else
            temp=dir(fullfile(basePath,'*.SleepState.states.mat'));
            %temp2 = dir(fullfile(basePath,'BloodModAnalysis','infraslowROI.mat'));
        end
        if isempty(temp)
            disp('Selected to restrict but have no SleepState.states file. Use buzcode state editor to state score.')
        else
            if ~isempty(altPathSessionGeneral)
                load(fullfile(altPathSessionGeneral,temp.name))
                %load(fullfile(altPathSessionGeneral,'BloodModAnalysis',temp2.name)) %infraslowROI
            else
                load(fullfile(basePath,temp.name))
                %load(fullfile(basePath,'BloodModAnalysis',temp2.name))
            end
            % Pull out PSS
            timePSS = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.t_clus;
            PSS = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.broadbandSlowWave;
            %timeInfra = infraslowROI.timestamps;
            %phaseInfra = infraslowROI.phase;
        end
    end
    
    %% Deal with SlowWaves/Up states
    
    % 1. Align UPs and DOWNs so one UP following every DOWN
    if length(SlowWaves.ints.UP) > length(SlowWaves.ints.DOWN)
        if SlowWaves.ints.DOWN(1,1) == SlowWaves.ints.UP(1,end) %Down follows UP and are more UPs, drop remaining UPs
            amtRemove=length(SlowWaves.ints.UP)-length(SlowWaves.ints.DOWN)-1; % -1 bc dropping first UP already
            SlowWaves.ints.UP = SlowWaves.ints.UP(2:end-amtRemove,:);
        elseif SlowWaves.ints.DOWN(1,end) == SlowWaves.ints.UP(1,1)
            amtRemove=length(SlowWaves.ints.UP)-length(SlowWaves.ints.DOWN);
            SlowWaves.ints.UP = SlowWaves.ints.UP(1:end-amtRemove,:);
        end
    elseif length(SlowWaves.ints.UP) < length(SlowWaves.ints.DOWN)
        if SlowWaves.ints.DOWN(1,1) == SlowWaves.ints.UP(1,end) % DOWN follows UP
            amtRemove=length(SlowWaves.ints.DOWN)-length(SlowWaves.ints.UP)+1;
            SlowWaves.ints.UP = SlowWaves.ints.UP(2:end,:);
            SlowWaves.ints.DOWN = SlowWaves.ints.DOWN(1:end-amtRemove,:);
        elseif SlowWaves.ints.DOWN(1,end) == SlowWaves.ints.UP(1,1) %UP follows DOWN
            amtRemove=length(SlowWaves.ints.DOWN)-length(SlowWaves.ints.UP);
            SlowWaves.ints.DOWN = SlowWaves.ints.DOWN(1:end-amtRemove,:);
        end
    elseif length(SlowWaves.ints.UP) == length(SlowWaves.ints.DOWN) %same number but UP first (need DOWN first)
        if SlowWaves.ints.DOWN(1,1) == SlowWaves.ints.UP(1,end) % DOWN follows UP
            % Drop first UP, drop a DOWN at the end
            SlowWaves.ints.UP = SlowWaves.ints.UP(2:end,:);
            SlowWaves.ints.DOWN = SlowWaves.ints.DOWN(1:end-1,:);
        elseif SlowWaves.ints.DOWN(1,end) == SlowWaves.ints.UP(1,1) % UP follows DOWN
        end
    end
    
    allUP = SlowWaves.ints.UP; % Use all UP for boolean so boundary clean
    
    % Select subset slowwaves
    % i. Restrict to interval specified
    if strcmp(restrictIn,'all')==0
        [status,~,~] = InIntervals(SlowWaves.ints.DOWN(:,1),restrict);
        SlowWaves.ints.DOWN = SlowWaves.ints.DOWN(status,:);
 %       SlowWaves.ints.UP = SlowWaves.ints.UP(status,:);
    end
    % ii. Cut lower bound
    if lowerBoundCutoff
        temp =[];
        temp = SlowWaves.ints.DOWN(:,2)- SlowWaves.ints.DOWN(:,1);
        SlowWaves.ints.DOWN = SlowWaves.ints.DOWN(temp>=lowerBoundCutoff,:);
        %SlowWaves.ints.UP = SlowWaves.ints.UP(temp>=lowerBoundCutoff,:);
    end
    
    if upperBoundCutoff
        temp =[];
        temp = SlowWaves.ints.DOWN(:,2)- SlowWaves.ints.DOWN(:,1);
        SlowWaves.ints.DOWN = SlowWaves.ints.DOWN(temp<=upperBoundCutoff,:);
        %SlowWaves.ints.UP = SlowWaves.ints.UP(temp<=upperBoundCutoff,:);
    end
    
    % Ripple and UP boolean stuff
    [status,~,~] = InIntervals(MUA.timestamps,allUP);
    %[status,~,~] = InIntervals(MUAIn.timestamps,allUP);
    %[status,~,~] = InIntervals(MUA.timestamps,SlowWaves.ints.UP);
    SlowWaveBool = double(status);
    if ischar(ripplesIn)
        [status,intervalSWR,~] = InIntervals(MUA.timestamps,ripples.timestamps);
        rippleBool = double(status);
        ampR = ripples.peakNormedPower;
        intervalSWR(intervalSWR==0)=nan; % make nan so can add prev indices
        % To remove nans  x(isnan(x)==1)=[];
    end
    
    %% JOIN across basePaths
    if i==1
        MUAAll.timestamps{1} = MUA.timestamps;
        MUAAll.data{1} = MUA.data;
        if ischar(MUAHPCin)
            MUAHPCAll.timestamps{1} = MUAHPC.timestamps;
            MUAHPCAll.data{1} = MUAHPC.data;
        end
        SlowWavesAll{1} = SlowWaves.ints.DOWN;
        SlowWaves_noneCutAll.ints.UP{1} = SlowWaves_noneCut.ints.UP;
        SlowWaves_noneCutAll.ints.DOWN{1} = SlowWaves_noneCut.ints.DOWN;
        UPAll{1} = SlowWaves.ints.UP;
        SlowWaveBoolAll{1} = SlowWaveBool;
        if ischar(ripplesIn)
            rippleBoolAll{1} = rippleBool;
            intervalSWRAll{1} = intervalSWR;
            ampRAll{1} = ampR;
            ripplesAll.timestamps{1} = ripples.timestamps;
            ripplesAll.peaks{1} = ripples.peaks;
            numRipplesPre{1}=0;
        end
        if ischar(sharpWavesIn) % ripple features
            SWAll.magnitude{1} = SW.magnitude;
            SWAll.timestamps{1} = SW.timePeak;
        end
    else
        timeZero = (MUAAll.timestamps{i-1}(end))+1/MUA.samplingRate;
        MUAAll.timestamps{i} = MUA.timestamps+timeZero;
        % Add other stuff
        MUAAll.data{i} = MUA.data;
        if ischar(MUAHPCin)
            MUAHPCAll.timestamps{i} = MUAHPC.timestamps+timeZero;
            MUAHPCAll.data{i} = MUAHPC.data;
        end
        SlowWavesAll{i} = SlowWaves.ints.DOWN+timeZero;
        UPAll{i} = SlowWaves.ints.UP+timeZero;
        SlowWaveBoolAll{i} = SlowWaveBool;
        SlowWaves_noneCutAll.ints.UP{i} = SlowWaves_noneCut.ints.UP+timeZero;
        SlowWaves_noneCutAll.ints.DOWN{i} = SlowWaves_noneCut.ints.DOWN+timeZero;
        if ischar(ripplesIn)
            rippleBoolAll{i} = rippleBool;
            numRipplesPre{i} = length(ripplesAll.timestamps{i-1});
            intervalSWRAll{i} = intervalSWR+numRipplesPre{i};
            ampRAll{i} = ampR;
            ripplesAll.timestamps{i} = ripples.timestamps+timeZero;
            ripplesAll.peaks{i} = ripples.peaks+timeZero;
        end
        if ischar(sharpWavesIn)
            SWAll.magnitude{i} = SW.magnitude;
            SWAll.timestamps{i} = SW.timePeak+timeZero;
        end
    end
    
end

% Concatenate vars across basePaths
MUAAll.timestamps=vertcat(MUAAll.timestamps{:});
MUAAll.data=vertcat(MUAAll.data{:});
if ischar(MUAHPCin)
    MUAHPCAll.timestamps=vertcat(MUAHPCAll.timestamps{:});
    MUAHPCAll.data=vertcat(MUAHPCAll.data{:});
end
SlowWavesAll=vertcat(SlowWavesAll{:});
UPAll=vertcat(UPAll{:});
SlowWaveBoolAll=vertcat(SlowWaveBoolAll{:});
intervalSWR_separate = intervalSWRAll; % hold onto these as separate and subtract numRipplesPre to get indices into og ripples struct for each basePath
intervalSWRAll=vertcat(intervalSWRAll{:});
SlowWaves_noneCutAll.ints.UP=vertcat(SlowWaves_noneCutAll.ints.UP{:});
SlowWaves_noneCutAll.ints.DOWN=vertcat(SlowWaves_noneCutAll.ints.DOWN{:});
numRipplesPre = [numRipplesPre{:}];

if ischar(ripplesIn)
    rippleBoolAll=vertcat(rippleBoolAll{:});
    ampRAll = vertcat(ampRAll{:});
    ripplesAll.peaks = vertcat(ripplesAll.peaks{:});
    ripplesAll.timestamps = vertcat(ripplesAll.timestamps{:});
    ampR = ampRAll; % Is this amp estimate i wanna use? Inspect how calculated
    %     durR = data.duration; % Inspect how calculated
    [~,ampRsort]=sort(ampR,'descend');
    %     [~,durRsort]=sort(durR,'descend');
    [xtileIndsA,~,~,~,~] = Xtiles(ampR,xtile,'log');
end
if ischar(sharpWavesIn) % ripple features
    SWAll.magnitude = vertcat(SWAll.magnitude{:});
    SWAll.timestamps = vertcat(SWAll.timestamps{:});
    magSW = SWAll.magnitude; % decide how want to deal with SW magnitude
    [xtileInds,~,~,~,~] = Xtiles(magSW,xtile,'lin');
    [~,magsort]=sort(magSW,'descend');
    [xtileIndsSW,~,~,~,~] = Xtiles(magSW,xtile,'lin');
end

% Choose trigger
if strcmp(trigger,'DOWN') %align to midpoint of DOWN state
    triggers = SlowWavesAll(:,1)+((SlowWavesAll(:,2)- SlowWavesAll(:,1))./2); %midpoints
elseif strcmp(trigger,'DOWNtoUP')
    triggers = SlowWavesAll(:,2); %end of DOWN state
elseif strcmp(trigger,'UPtoDOWN')
    triggers = SlowWavesAll(:,1);
end

% DOWN sort
durations = SlowWavesAll(:,2)- SlowWavesAll(:,1);
durUP = UPAll(:,2)- UPAll(:,1);
[~,sortUse]=sort(durations,'descend');
downSort = durations(sortUse);
if ~isempty(xtile)
    [xtileInds,~,~,~,~] = Xtiles(durations,xtile,'log'); % refer to OG ordering
end

% Additional ripple features stuff - comment out for now
%     maps=[]; data=[]; stats = [];
%     [maps,data,stats] = bz_RippleStats(double(rippleLFP.data),hpc.timestamps,ripples,'durations',[-.1 .1]); %-.3 to .
%     [xtileIndsD,~,~,~,~] = Xtiles(durR,xtile,'log');

%% Calculate CCGs between trigger and continuous variable
% MUA
ccgDataCont = MUAAll.data;
ccgTimeCont = MUAAll.timestamps;
[~,t_lag,CCGstd,CCGall,~ ] = EventVsContinousCCG2( ccgDataCont,ccgTimeCont,triggers,timelag); %CCG all dims: time x slow waves (trials); each cell = basePath

% MUA HPC
if ischar(MUAHPCin)
    ccgDataCont = MUAHPCAll.data;
    ccgTimeCont = MUAHPCAll.timestamps;
    [~,t_lag,CCGstd_HPC,CCGall_HPC,~ ] = EventVsContinousCCG2( ccgDataCont,ccgTimeCont,triggers,timelag); %CCG all dims: time x slow waves (trials); each cell = basePath
end

% UP bool
ccgDataCont = SlowWaveBoolAll;
ccgTimeCont = MUAAll.timestamps;
[~,~,~,CCGall_SlowWaveBool,~ ] = EventVsContinousCCG2( ccgDataCont,ccgTimeCont,triggers,timelag); %CCG all dims: time x slow waves (trials); each cell = basePath

% SWR bool
if ischar(ripplesIn)
    ccgDataCont = rippleBoolAll;
    ccgTimeCont = MUAAll.timestamps;
    [~,~,~,CCGall_rippleBool,~ ] = EventVsContinousCCG2( ccgDataCont,ccgTimeCont,triggers,timelag); %CCG all dims: time x slow waves (trials); each cell = basePath
    
    ccgDataCont = intervalSWRAll;
    ccgDataCont(isnan(ccgDataCont)==1)=0; % turn nans to 0s before calc mean
    ccgTimeCont = MUAAll.timestamps;
    [~,~,~,CCGall_rippleIndex,~ ] = EventVsContinousCCG2( ccgDataCont,ccgTimeCont,triggers,timelag); %CCG all dims: time x slow waves (trials); each cell = basePath
    
end
%% Operations for plots - add/refine this later
% i. Identify features synchrony and UP state post DOWN
if strcmp(trigger,'DOWNtoUP')
    % Choosing window for summing
    ccgSort = squeeze(CCGall(:,sortUse))'; %trials x time; long to short DOWN
    figure(20);
    clf
    subplot 141
    imagesc(t_lag,1:size(ccgSort,2),ccgSort);
    %caxis([.05e2 3e2])
    
    subplot 142
    imagesc(t_lag,1:size(ccgSort,2),ccgSort);
    %caxis([.05e2 3e2])
    hold on
    plot([0 0],get(gca,'ylim'),'y')
    plot([.05 .05],get(gca,'ylim'),'y')
    plot([.008 .008],get(gca,'ylim'),'r')
    plot([.02 .02],get(gca,'ylim'),'r')
    
    subplot 143
    windowDOWNtoUP = [0 .05]; % first 50 ms of up onset
    downSort = durations(sortUse);
    ccgSort = squeeze(CCGall(:,sortUse))'; %trials x time; long to short DOWN
    windowCCG_50ms = ccgSort(:,t_lag>=windowDOWNtoUP(1) & t_lag<windowDOWNtoUP(2));
    windowTime_50ms = t_lag((t_lag>=windowDOWNtoUP(1) & t_lag<windowDOWNtoUP(2)));
    imagesc(windowTime_50ms,1:size(windowCCG_50ms,1),windowCCG_50ms);
    %caxis([.05e2 3e2])
    title('Window: 0 to 50 ms')
    hold on
    plot([.008 .008],get(gca,'ylim'),'r')
    plot([.02 .02],get(gca,'ylim'),'r')
    
    subplot 144
    windowDOWNtoUP = [0.008 .02]; %8 ms from UP onset to 20 ms
    downSort = durations(sortUse);
    ccgSort = squeeze(CCGall(:,sortUse))'; %trials x time; long to short DOWN
    windowCCG_20ms = ccgSort(:,t_lag>=windowDOWNtoUP(1) & t_lag<windowDOWNtoUP(2));
    windowTime_20ms = t_lag(t_lag>=windowDOWNtoUP(1) & t_lag<windowDOWNtoUP(2));
    imagesc(windowTime_20ms,1:size(windowCCG_20ms,1),windowCCG_20ms);
    %caxis([.05e2 3e2])
    title('Window: 8 to 20 ms')
    subtitle('Selection of window for synchrony var')
    
    % Save for later plotting **RETURN TO THIS LATER**
    synchInfo(:).timeDOWNonset = SlowWavesAll(sortUse,1);
    synchInfo(:).timeDOWNoffset = SlowWavesAll(sortUse,2);
    synchInfo(:).durationDOWN = downSort; %duration down; sorted
    synchInfo(:).synchronyDOWNtoUP_50ms = sum(windowCCG_50ms,2); %50 ms post UP
    synchInfo(:).synchronyDOWNtoUP_20ms  = sum(windowCCG_20ms,2); %20 ms post UP
    synchInfo(:).durPostUP = durUP(sortUse); %duration of following UP; sort by same thing as DOWN, works bc = UP following the DOWN
    temp=sortUse-1;
    temp=temp(temp>=1);
    synchInfo(:).durPreUP = [nan; durUP(temp)];%duration preceding UP;
    %PSS value at DOWN->UP transition
    if ischar(includeState)
        triggersPSS=interp1(timePSS,PSS,triggers);
        synchInfo(:).PSS = triggersPSS(sortUse); %PSS value at each DOWN state (too short timescale to be relevant)
        %triggersInfra = interp1(timeInfra,phaseInfra,triggers);
        %synchInfo(:).infra = triggersInfra(sortUse);
    end
    %windowTime_50ms; windowTime_20ms;
else
    synchInfo = [];
end

%%%%%%%%%%%%%%%%%%%
%% Timing of SWRs wrt UP/DOWN

if ischar(ripplesIn)
    
    % Steps:
    % - SWR indices for ripples in subset DOWNs actually used
    % - For each of these ripples, find DOWN/UP info from list of ALL
    % slowwaves (no restriction so can get accurate timing) [is this
    % correct decision?]
    % - Get info re proximity to transitions
    % - leave out SWRs w excluded DOWN state times
    
    % Ripple indices for every DOWN state of subset plotted: SlowWaves var
    % CCG  sorted by DOWN dur
    ripplesDOWN=[];
    sort_ccgRippleIndex = CCGall_rippleIndex(:,sortUse)'; % index into rippleAll
    for dd=1:size(sort_ccgRippleIndex,1)
        rippleInds=unique(sort_ccgRippleIndex(dd,:));
        ripplesDOWN(dd).rippleIndex=rippleInds(rippleInds~=0);
        ripplesDOWN(dd).numberRipples = length(rippleInds(rippleInds~=0));
        ripplesDOWN(dd).duration = durations(dd);
        ripplesDOWN(dd).xtileInds = xtileInds(sortUse(dd));
        % use below if want unsorted
        % ripplesDOWN(dd).sortByDurIndex = sortUse(dd);
        %ripplesDOWN(dd).xtileInds = xtileInds(dd);
    end
    DOWNInfo = ripplesDOWN; %just renaming output bc too lazy to change var names, fix later
    
    % Find prox SWRs to UP/DOWN
    rippleIndices=unique([ripplesDOWN(:).rippleIndex]); % ripple indices used
    rippleIndices(isnan(rippleIndices)) = [];
    numRipples = length(rippleIndices);
    % In DOWN
    [statusDOWN,intervalDOWN,~] = InIntervals(ripplesAll.peaks(rippleIndices),SlowWaves_noneCutAll.ints.DOWN);
    rippleAllDOWN(:,1) = find(statusDOWN==1); % index into rippleIndices (NOT ripplesAll)
    rippleAllDOWN(:,1) = rippleIndices(rippleAllDOWN(:,1)); % index into ripplesAll
    rippleAllDOWN(:,2) = intervalDOWN(statusDOWN); % index into DOWN ints in SlowWaves_noneCutAll
    % In UP
    [statusUP,intervalUP,~] = InIntervals(ripplesAll.peaks(rippleIndices),SlowWaves_noneCutAll.ints.UP);
    rippleAllUP(:,1) = find(statusUP==1); % index into rippleIndices (NOT ripplesAll)
    rippleAllUP(:,1) = rippleIndices(rippleAllUP(:,1)); % index into ripplesAll
    rippleAllUP(:,2) = intervalUP(statusUP); % index into UP ints
    % Save some info
    % SWRs during DOWN
    for rr=1:length(rippleAllDOWN)
        rippleAllDOWNs(rr).indexRipplesAll = rippleAllDOWN(rr,1);
        rippleAllDOWNs(rr).indexDOWNAll = rippleAllDOWN(rr,2);
        rippleAllDOWNs(rr).SWRdurDOWN = 1; %during DOWN
        rippleAllDOWNs(rr).SWRdurUP = 0; % not durng UP
        rippleAllDOWNs(rr).SWRpeak = ripplesAll.peaks(rippleAllDOWN(rr,1)); % timing peak ripple
        rippleAllDOWNs(rr).SWRstart = ripplesAll.timestamps(rippleAllDOWN(rr,1),1); % start ripple
        rippleAllDOWNs(rr).SWRdur = ripplesAll.timestamps(rippleAllDOWN(rr,1),2)-ripplesAll.timestamps(rippleAllDOWN(rr,1),1); % duration ripple
        rippleAllDOWNs(rr).SWRamp = ampRAll(rippleAllDOWN(rr,1)); % amp ripple FIX
        rippleAllDOWNs(rr).DOWNonset = SlowWaves_noneCutAll.ints.DOWN(rippleAllDOWN(rr,2),1); % DOWN onset
        rippleAllDOWNs(rr).DOWNoffset = SlowWaves_noneCutAll.ints.DOWN(rippleAllDOWN(rr,2),2); % DOWN offset
        rippleAllDOWNs(rr).DOWNdur = SlowWaves_noneCutAll.ints.DOWN(rippleAllDOWN(rr,2),2)-SlowWaves_noneCutAll.ints.DOWN(rippleAllDOWN(rr,2),1); % duration DOWN
        rippleAllDOWNs(rr).closestUPtoDOWNtrans = SlowWaves_noneCutAll.ints.DOWN(rippleAllDOWN(rr,2),1) - ripplesAll.peaks(rippleAllDOWN(rr,1)); % (- value, transition occurs first) distance swr peak from UP->DOWN transitoin
        rippleAllDOWNs(rr).closestDOWNtoUPtrans = SlowWaves_noneCutAll.ints.DOWN(rippleAllDOWN(rr,2),2) - ripplesAll.peaks(rippleAllDOWN(rr,1)); % (+ value, transition follows) distance swr peak from DOWN->UP transitoin
        if ischar(sharpWavesIn) % ripple features
            rippleAllDOWNs(rr).SWmag = SWAll.magnitude(rippleAllDOWN(rr,1)); % magnitude ripple
        end
    end
    % SWRs during UP
    for rr=1:length(rippleAllUP)
        rippleAllUPs(rr).indexRipplesAll = rippleAllUP(rr,1);
        rippleAllUPs(rr).indexUPAll = rippleAllUP(rr,2);
        rippleAllUPs(rr).SWRdurDOWN = 0; %during DOWN
        rippleAllUPs(rr).SWRdurUP = 1; % not durng UP
        rippleAllUPs(rr).SWRpeak = ripplesAll.peaks(rippleAllUP(rr,1)); % timing peak ripple
        rippleAllUPs(rr).SWRstart = ripplesAll.timestamps(rippleAllUP(rr,1),1); % start ripple
        rippleAllUPs(rr).SWRdur = ripplesAll.timestamps(rippleAllUP(rr,1),2)-ripplesAll.timestamps(rippleAllUP(rr,1),1); % duration ripple
        rippleAllUPs(rr).SWRamp = ampRAll(rippleAllUP(rr,1)); % amp ripple
        rippleAllUPs(rr).UPonset = SlowWaves_noneCutAll.ints.UP(rippleAllUP(rr,2),1); % UP onset
        rippleAllUPs(rr).UPoffset = SlowWaves_noneCutAll.ints.UP(rippleAllUP(rr,2),2); % UP offset
        rippleAllUPs(rr).UPdur = SlowWaves_noneCutAll.ints.UP(rippleAllUP(rr,2),2)-SlowWaves_noneCutAll.ints.UP(rippleAllUP(rr,2),1); % duration UP
        rippleAllUPs(rr).closestUPtoDOWNtrans = SlowWaves_noneCutAll.ints.UP(rippleAllUP(rr,2),2) - ripplesAll.peaks(rippleAllUP(rr,1)); % distance swr peak from closest UP->DOWN transitoin
        rippleAllUPs(rr).closestDOWNtoUPtrans = SlowWaves_noneCutAll.ints.UP(rippleAllUP(rr,2),1) - ripplesAll.peaks(rippleAllUP(rr,1)); % distance swr peak from DOWN->UP transitoin
        if ischar(sharpWavesIn) % ripple features
            rippleAllUPs(rr).SWmag = SWAll.magnitude(rippleAllUP(rr,1)); % magnitude ripple
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%
%%
% ii. Identify outline of DOWN state (very sloppy; must use up bool)
sizeTemp = size(CCGall);
% Invert boolean matrix so DOWN state = 1

traceDOWN = zeros(sizeTemp);
traceDOWN(CCGall_SlowWaveBool==0) = 1; %invert, DOWN state = 1
traceDOWN = traceDOWN(:,sortUse)'; %5976,617
% smooth
windowSize = 2;
kernel = ones(windowSize) / windowSize ^ 2;
blurryImage = conv2(single(traceDOWN), kernel, 'same');
traceDOWN = blurryImage >= .99; % Rethreshold
% Identify edge of object to be traced
r=100;
x=mean(traceDOWN(r:r+3,:),1);
y=diff(x);
[~,c] = max(y); % sample in time @edge
boundary = bwtraceboundary(logical(traceDOWN),[r c+1],'N'); % col 2 = x, col 1= y in imagesc
if isempty(boundary)
    r=110;
    x=mean(traceDOWN(r:r+3,:),1);
    y=diff(x);
    [~,c] = max(y); % sample in time @edge
    boundary = bwtraceboundary(logical(traceDOWN),[r c+1],'N'); % col 2 = x, col 1= y in imagesc
    if isempty(boundary)
        r=110;
        x=mean(traceDOWN(r:r+3,:),1);
        y=diff(x);
        [~,c] = max(y); % sample in time @edge
        boundary = bwtraceboundary(logical(traceDOWN),[r c-1],'N'); % col 2 = x, col 1= y in imagesc
        if isempty(boundary)
            r=115;
            x=mean(traceDOWN(r:r+3,:),1);
            y=diff(x);
            [~,c] = max(y); % sample in time @edge
            boundary = bwtraceboundary(logical(traceDOWN),[r c-1],'N'); % col 2 = x, col 1= y in imagesc
            if isempty(boundary)
                r=115;
                x=mean(traceDOWN(r:r+3,:),1);
                y=diff(x);
                [~,c] = max(y); % sample in time @edge
                boundary = bwtraceboundary(logical(traceDOWN),[r c+1],'N'); % col 2 = x, col 1= y in imagesc
            end
        end
    end
end



%% All together
figure(1) % FIGURE 1
clf
% MUA
subplot(4,3,[1 4 7 10])
% h1 = axes;
imagesc(squeeze(CCGall(:,sortUse))')
hold on
%caxis([.05e2 3e2])

% MUA + down outline
subplot(4,3,[2 5 8 11])
% h1 = axes;
imagesc(squeeze(CCGall(:,sortUse))')
hold on
%caxis([.05e2 3e2])
if ~isempty(boundary)
    patch(boundary(:,2),boundary(:,1),'r','faceAlpha',.2)
end
%imagesc(t_lag(:,1,1),1:length(sortUse),squeeze(CCGall{1,i}(:,sortUse))')
plot([sizeTemp(1)/2 sizeTemp(1)/2],get(gca,'ylim'),'--w')

% SWR bool
if ischar(ripplesIn)
    subplot(4,3,[3 6 9 12])
    %h1 = axes;
    spyMod((double(squeeze(CCGall_rippleBool(:,sortUse)))'),'k'); %t_lag(:,1,1),1:length(sortUse),
    hold on
    if ~isempty(boundary)
        patch(boundary(:,2),boundary(:,1),'r','faceAlpha',.2,'edgeColor','r')
    end
    plot([sizeTemp(1)/2 sizeTemp(1)/2],get(gca,'ylim'),'--w')
    %             h2 = axes;
    %             set(h2, 'Ydir', 'reverse')
    %             set(h2, 'YAxisLocation', 'Right')
    %             set(h2, 'Ytick', [])
end
xlabel('samples')

subtitle(basePath)

%% Figure - plot HPC and CTX MUA 
avCTX = mean(CCGall,2);
avCTX = bz_NormToRange(avCTX,[0 1]);

% figure(41)
% clf
% % CTX MUA
% subplot 322
% imagesc(squeeze(CCGall(:,sortUse))')
% hold on
% caxis([.05e2 3.5e2])
% if ~isempty(boundary)
%     patch(boundary(:,2),boundary(:,1),'r','faceAlpha',.2)
% end
% %imagesc(t_lag(:,1,1),1:length(sortUse),squeeze(CCGall{1,i}(:,sortUse))')
% plot([sizeTemp(1)/2 sizeTemp(1)/2],get(gca,'ylim'),'--w')
% % Overlay MUA
% yyaxis right
% plot(1:length(avCTX),avCTX./4,'linewidth',3,'color','w')
% axis tight
% ylim([0 1])
% ylabel('MUA CTX (0-1 Norm)')
% 
% subplot 324
% imagesc(squeeze(CCGall_HPC(:,sortUse))')
% hold on
% caxis([1.5 2.8])
% if ~isempty(boundary)
%    % patch(boundary(:,2),boundary(:,1),'r','faceAlpha',.2)
% end
% %imagesc(t_lag(:,1,1),1:length(sortUse),squeeze(CCGall{1,i}(:,sortUse))')
% plot([sizeTemp(1)/2 sizeTemp(1)/2],get(gca,'ylim'),'--w')
% 
% % Plot MUA 
% subplot 326
% % plot HPC av MUA
% avHPC = mean(CCGall_HPC,2);
% avHPC = bz_NormToRange(avHPC,[0 1]);
% yyaxis right
% plot(t_lag,avHPC,'color',[0.8500, 0.3250, 0.0980],'linewidth',3)
% axis tight
% set(gca,'yTick',[])
% ylabel('MUA HPC (0-1 Norm)')
% hold on
% 
% % Plot CTX av MUA
% yyaxis left
% plot(t_lag,avCTX,'linewidth',3,'color',[17 103 177]./255)
% axis tight
% ylabel('MUA CTX (0-1 Norm)')
% plot([0 0],get(gca,'ylim'),'--r')
% 
% subplot(3,2,[1 3 5])
% imagesc(squeeze(CCGall_HPC(:,sortUse))')
% hold on
% caxis([1.5 2.8])
% spyMod((double(squeeze(CCGall_rippleBool(:,sortUse)))'),'k');
% if ~isempty(boundary)
%     patch(boundary(:,2),boundary(:,1),'r','faceAlpha',.2)
% end
% %imagesc(t_lag(:,1,1),1:length(sortUse),squeeze(CCGall{1,i}(:,sortUse))')
% plot([sizeTemp(1)/2 sizeTemp(1)/2],get(gca,'ylim'),'--w')
% 
% suptitle('HPC and CTX MUA')
    
%% ripples only figure
if ischar(ripplesIn)
    figure(5) % FIGURE 5
    subplot(5,1,1:4)
    spyMod((double(squeeze(CCGall_rippleBool(:,sortUse)))'),'k'); %t_lag(:,1,1),1:length(sortUse),
    hold on
    if ~isempty(boundary)
        patch(boundary(:,2),boundary(:,1),'r','faceAlpha',.2,'edgeColor','r')
    end
    plot([sizeTemp(1)/2 sizeTemp(1)/2],get(gca,'ylim'),'--w')
    %             h2 = axes;
    %             set(h2, 'Ydir', 'reverse')
    %             set(h2, 'YAxisLocation', 'Right')
    %             set(h2, 'Ytick', [])
    xlabel('samples')
    
    subplot(5,1,5)
    avCCG = sum(CCGall_rippleBool,2);
    avCCG = avCCG./sum(avCCG);
        plot(t_lag,avCCG,'k','linewidth',3)
        hold on
        plot([0 0],get(gca,'ylim'),'r')
        xlabel('time(s)')
        box off 
        ylabel('P(SWR)')
end

% if strcmp(continuousDataUse,'spindlePower')
%     figure(21)
%     imagesc(squeeze(CCGall(:,sortUse))')
%     hold on
%     caxis([.1 800])
%     if ~isempty(boundary)
%         patch(boundary(:,2),boundary(:,1),'r','faceAlpha',.2)
%     end
%     %imagesc(t_lag(:,1,1),1:length(sortUse),squeeze(CCGall{1,i}(:,sortUse))')
%     plot([sizeTemp(1)/2 sizeTemp(1)/2],get(gca,'ylim'),'--w')
% end


    
%% TAKES VERY LONG FIGURE OUT HOW TO MAKE FASTER
% Replace each ripple INDEX with the amplitude or SW QUINTILE to which it belongs

if ~isempty(SWRFeatureSelection)
    
    % Deal with ripples being in root basepath if multiple regions
    if ~isempty(altPathSessionGeneral)
        use = altPathSessionGeneral;
    else
        use = basePath;
    end
    
    % SW Magnitude
    if strcmp(SWRFeatureSelection,'SWamp')
        if isfile(fullfile(use,'CCGVar_amplitudeSW.mat'))
            load(fullfile(use,'CCGVar_amplitudeSW.mat'))
            CCGVar = CCGVarBool.data;
        else
            CCGVar = intervalSWR;
            tic
            for kk=1:length(ripples.peakNormedPower) %for each ripple
                if mod(kk,100)==0; disp(num2str(kk)); end
                CCGVar(intervalSWR == kk) = xtileIndsSW(kk);
            end
            CCGVar(isnan(CCGVar)) = 0;
            CCGVarBool.data = CCGVar;
            CCGVarBool.time = ccgTimeCont; % length ripple lfp
            save('CCGVar_amplitudeSW.mat','CCGVarBool','-v7.3')
            disp('done SW amplitude sort')
        end
        
        % Amplitude - sort by quintiles - took 21 MINUTES
    elseif strcmp(SWRFeatureSelection,'SWRamp')
        if isfile(fullfile(use,'CCGVar_amplitudeSWR.mat'))
            load(fullfile(use,'CCGVar_amplitudeSWR.mat'))
            CCGVar = CCGVarBool.data;
        else
            CCGVar = intervalSWR;
            tic
            for kk=1:length(ripples.peakNormedPower)
                if mod(kk,100)==0; disp(num2str(kk)); end
                CCGVar(intervalSWR == kk) = xtileIndsA(kk);
            end
            CCGVar(isnan(CCGVar)) = 0;
            CCGVarBool.data = CCGVar;
            CCGVarBool.time = ccgTimeCont; % length ripple lfp
            save('CCGVar_amplitudeSWR.mat','CCGVarBool','-v7.3')
            disp('done SWR amplitude sort')
        end
        
        % Duration
    elseif strcmp(SWRFeatureSelection,'SWRdur')
        if isfile(fullfile(use,'CCGVar_durationSWR.mat'))
            load(fullfile(use,'CCGVar_durationSWR.mat'))
            CCGVar = CCGVarBool.data;
        else
            tic
            for kk=1:length(ripples.peakNormedPower)
                if mod(kk,100)==0; disp(num2str(kk)); end
                CCGVar(intervalSWR == kk) = xtileIndsD(kk);
            end
            CCGVar(isnan(CCGVar)) = 0;
            CCGVarBool.data = CCGVar;
            CCGVarBool.time = ccgTimeCont; % length ripple lfp
            save('CCGVar_durationSWR.mat','CCGVarBool','-v7.3')
            disp('done SWR duration sort')
        end
    end
    
    
    %% Plot SWR amplitude across ALL SWRs ... single plot not separate by DOWN dur xtile
    [CCGmean_SWRXtile,time,~,CCGall_SWRXtile,~ ] = EventVsContinousCCG2( CCGVar,ccgTimeCont,triggers,timelag);
    
    %% Plot SWR-feature version of figure
    colorSWRFeature = {'.c','.g','.b','.k','.r'};
    colorSWRFeature2 = {'c','g','b','k','r'};
    %figure(10)
    for kk=1:xtile % for each SWR xtile
        figure(40)
        % Raster plot by swr xtile
        subplot(5,1,1:4)
        hold on
        toPlot = double(CCGall_SWRXtile(:,sortUse)');
        toPlot(toPlot~=kk)=0;
        toPlot(toPlot~=0)=1;
        spyMod(sparse(toPlot),colorSWRFeature{kk}) %'color',colorSWRFeature(kk,:));
        hold on
        if ~isempty(boundary) & kk == xtile
            patch(boundary(:,2),boundary(:,1),'r','faceAlpha',.2,'edgeColor','r')
        end
        if kk==xtile
            ylabel('DOWNs')
            plot([sizeTemp(1)/2 sizeTemp(1)/2],get(gca,'ylim'),'--r')
        end
        
        subplot(5,1,5)
        avCCGSWR{kk} = sum(toPlot,1);
        avCCGSWR{kk} = avCCGSWR{kk}./sum(sum(CCGall_SWRXtile));
        hold on
        plot(time,avCCGSWR{kk},'color',colorSWRFeature2{kk},'linewidth',3)
        if kk==xtile
            xlabel('time(s)')
            plot([0 0],get(gca,'ylim'),'--r')
            ylabel('P(SWR)')
        end
        axis tight
        
    end
    
    %%%% WORK HERE ADD IN WHITE LINE 
   
    %%%%%%
    
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% Each Xtile

if ~isempty(xtile)
    
    % Structure of output vars: time x xtile x mouse
    figure(11)
    hold on
    
    for ii=1:xtile
        triggersXtile = triggers(xtileInds==ii);
        durXtile = durations(xtileInds==ii);
        [~,sortXtile]=sort(durXtile,'descend');
        xtileDurs{ii} = durXtile(sortXtile);
        % Calculate
        % MUA
        ccgDataCont = MUAAll.data;
        ccgTimeCont = MUAAll.timestamps;
        [XCCGmean(:,ii),t_lag,~,XCCGall{ii},~ ] = EventVsContinousCCG2( ccgDataCont,ccgTimeCont,triggersXtile,timelag); %CCG all dims: time x slow waves (trials); each cell = basePath
        
        % UP bool
        ccgDataCont = SlowWaveBoolAll;
        ccgTimeCont = MUAAll.timestamps;
        [XCCGmean_SlowWaveBool(:,ii),t_lag,~,XCCGall_SlowWaveBool{ii},~ ] = EventVsContinousCCG2( ccgDataCont,ccgTimeCont,triggersXtile,timelag); %CCG all dims: time x slow waves (trials); each cell = basePath
        
        % SWR bool
        if ischar(ripplesIn)
            % Ripple bool
            ccgDataCont = rippleBoolAll;
            ccgTimeCont = MUAAll.timestamps;
            [XCCGmean_rippleBool(:,ii),t_lag,~,XCCGall_rippleBool{ii},~ ] = EventVsContinousCCG2( ccgDataCont,ccgTimeCont,triggersXtile,timelag); %CCG all dims: time x slow waves (trials); each cell = basePath
            
            % Ripple indices
            ccgDataCont = intervalSWRAll;
            ccgTimeCont = MUAAll.timestamps;
            [XCCGmean_rippleIndices(:,ii),t_lag,~,XCCGall_rippleIndices{ii},~ ] = EventVsContinousCCG2( ccgDataCont,ccgTimeCont,triggersXtile,timelag); %CCG all dims: time x slow waves (trials); each cell = basePath
            
        end
        
        % Store
        sortXtileAll{ii} = sortXtile;
        Xall_rippleBool{ii}=single((XCCGall_rippleBool{ii}(:,sortXtile)'));
        Xall_SlowWaveBool{ii}=single(XCCGall_SlowWaveBool{ii}(:,sortXtile)');
        Xall_rippleIndices{ii}=single(XCCGall_rippleIndices{ii}(:,sortXtile)');
        Xall{ii}=squeeze(XCCGall{ii}(:,sortXtile))';
        
        sizeTemp = size(XCCGall{ii});
        
        % Make patches for DOWN
        % Identify outline of DOWN state
        % Invert boolean matrix so DOWN state = 1
        
        traceDOWN = zeros(sizeTemp);
        traceDOWN(XCCGall_SlowWaveBool{ii}==0) = 1; %invert, DOWN state = 1
        traceDOWN = traceDOWN(:,sortXtile)'; %5976,617
        % smooth
        windowSize = 2;
        kernel = ones(windowSize) / windowSize ^ 2;
        blurryImage = conv2(single(traceDOWN), kernel, 'same');
        traceDOWN = blurryImage >= .99; % Rethreshold
        
        % Identify edge of object to be traced
        r=100;
        x=mean(traceDOWN(r:r+3,:),1);
        y=diff(x);
        [~,c] = max(y); % sample in time @edge
        boundary2 = bwtraceboundary(logical(traceDOWN),[r c+1],'N'); % col 2 = x, col 1= y in imagesc
        if isempty(boundary2)
            r=110;
            x=mean(traceDOWN(r:r+3,:),1);
            y=diff(x);
            [~,c] = max(y); % sample in time @edge
            boundary2 = bwtraceboundary(logical(traceDOWN),[r c+1],'N'); % col 2 = x, col 1= y in imagesc
            if isempty(boundary2)
                r=110;
                x=mean(traceDOWN(r:r+3,:),1);
                y=diff(x);
                [~,c] = max(y); % sample in time @edge
                boundary2 = bwtraceboundary(logical(traceDOWN),[r c-1],'N'); % col 2 = x, col 1= y in imagesc
                if isempty(boundary2)
                    r=100;
                    x=mean(traceDOWN(r:r+3,:),1);
                    y=diff(x);
                    [~,c] = max(y); % sample in time @edge
                    boundary2 = bwtraceboundary(logical(traceDOWN),[r c-1],'N'); % col 2 = x, col 1= y in imagesc
                end
            end
        end
        
        
        boundaryUPBool{ii} = boundary2;
        %             else
        %                 boundaryUPBool{ii}=[];
        
        
        %% SWR Features (within DOWN xtile, go through each SWR xtile)
        % If want to calc SWR Amp xtile x Dur DOWN xtile
        %SWamp,SWRamp,SWRdur;
        if SWRFeatureSelection
            
            [XCCGmean_SWRXtile(:,ii),t_lag_SWRXtile,~,XCCGall_SWRXtile{ii},~ ] = EventVsContinousCCG2( CCGVar,ccgTimeCont,triggersXtile,timelag);
            
            % Plot SWR-feature version of figure
            colorSWRFeature = {'.c','.g','.b','.k','.r'};
            colorSWRFeature2 = {'c','g','b','k','r'};
            %figure(10)
            for kk=1:xtile % for each SWR xtile
                figure(10)
                % Raster plot by swr xtile
                subplot(6,xtile,[ii ii+xtile ii+xtile*2])
                toPlot = double(XCCGall_SWRXtile{ii}(:,sortXtile)');
                toPlot(toPlot~=kk)=0;
                toPlot(toPlot~=0)=1;
                spyMod(sparse(toPlot),colorSWRFeature{kk}) %'color',colorSWRFeature(kk,:));
                hold on
                if ~isempty(boundaryUPBool{ii}) & kk == xtile
                    patch(boundaryUPBool{ii}(:,2),boundaryUPBool{ii}(:,1),'r','faceAlpha',.2,'edgeColor','r')
                end
                if ii==1
                    ylabel('DOWNs')
                end
                plot([sizeTemp(1)/2 sizeTemp(1)/2],get(gca,'ylim'),'--r')
                
                % Line for each
                subplot(6,xtile,ii + xtile*3)
                avCCG = sum(toPlot,1);
                avCCG = avCCG./sum(avCCG);
                % yyaxis left
                plot(t_lag,avCCG,colorSWRFeature2{kk},'linewidth',1.2)
                hold on
                box off; axis tight; yy=get(gca,'ylim');
                if ii==xtile & kk==xtile
                    ylimzmax=get(gca,'ylim');
                end
                xlabel('time(s)')
                if ii==1
                    ylabel('P(SWR)')
                end
                
                figure(14) %all together
                subplot(1,xtile,kk)
                hold on
                % Raster plot by swr xtile
                toPlot = double(XCCGall_SWRXtile{ii}(:,sortXtile)');
                toPlot(toPlot~=kk)=0;
                toPlot(toPlot~=0)=1;
                spyMod(sparse(toPlot),colorSWRFeature{kk}) %'color',colorSWRFeature(kk,:));
            end
            %hold on
            %plot([0 0],[yy(1) yy(2)],'--r')
        end
        
        % Adjust axes to be same for lines
        title(['xtile ' num2str(ii)])
        
        % Collection variable; ii = xtile, iii = which variable
        %avCCG = CCGmean_SWRXtile(:,ii);
        %avCCGAll{iii,ii} = avCCG; %save avCCG for later plotting
        
    
    %%
    % Plot each xtile
    
    %% MUA
    figure(11)
    subplot(4,xtile,[ii ii+xtile ii+xtile*2])
    imagesc(t_lag,1:length(sortXtile),squeeze(XCCGall{ii}(:,sortXtile))')
    %caxis([.05e2 3e2])
    hold on
    if ~isempty(boundaryUPBool{ii})
        patch(boundaryUPBool{ii}(:,2),boundaryUPBool{ii}(:,1),'r','faceAlpha',.2)
    end
    plot([sizeTemp(1)/2 sizeTemp(1)/2],get(gca,'ylim'),'--r')
    
    if ii~=1
        set(gca, 'YTickLabel', [])
    end
    set(gca, 'XTickLabel', [])
    title(['xtile ' num2str(ii)])
    
    % Collection variable; ii = xtile, iii = which variable
    avCCG = XCCGmean(:,ii);
    avCCGAll{ii} = avCCG; %save avCCG for later plotting
    
    subplot(4,xtile,ii + xtile*3)
    plot(t_lag,avCCG,'linewidth',1.2,'color','k')
    
    box off; axis tight; yy=get(gca,'ylim');
    xlabel('time(s)')
    if ii==1
        ylabel('P(SWR)')
    end
    hold on
    plot([0 0],[yy(1) yy(2)],'--r')
    
    %% UP bool
    figure(12)
    subplot(5,xtile,[ii ii+xtile ii+xtile*2])
    spyMod(sparse(double((XCCGall_SlowWaveBool{ii}(:,sortXtile)'))),'k');
    hold on
    if ~isempty(boundaryUPBool{ii})
        patch(boundaryUPBool{ii}(:,2),boundaryUPBool{ii}(:,1),'r','faceAlpha',.2,'edgeColor','r')
    end
    if ii==1
        ylabel('DOWNs')
    end
    plot([sizeTemp(1)/2 sizeTemp(1)/2],get(gca,'ylim'),'--r')
    
    if ii~=1
        set(gca, 'YTickLabel', [])
    end
    set(gca, 'XTickLabel', [])
    title(['xtile ' num2str(ii)])
    
    % Collection variable; ii = xtile, iii = which variable
    avCCG = XCCGmean_SlowWaveBool(:,ii);
    avCCGAll_SlowWaveBool{ii} = avCCG; %save avCCG for later plotting
    
    subplot(4,xtile,ii + xtile*3)
    plot(t_lag,avCCG,'linewidth',1.2,'color','k')
    
    box off; axis tight; yy=get(gca,'ylim');
    xlabel('time(s)')
    if ii==1
        ylabel('P(SWR)')
    end
    hold on
    plot([0 0],[yy(1) yy(2)],'--r')
    %% ripple bool
    figure(13)
    
    subplot(5,xtile,[ii ii+xtile ii+xtile*2])
    spyMod(sparse(double((XCCGall_rippleBool{ii}(:,sortXtile)'))),'k');
    hold on
    if ~isempty(boundaryUPBool{ii})
        patch(boundaryUPBool{ii}(:,2),boundaryUPBool{ii}(:,1),'r','faceAlpha',.2,'edgeColor','r')
    end
    if ii==1
        ylabel('DOWNs')
    end
    plot([sizeTemp(1)/2 sizeTemp(1)/2],get(gca,'ylim'),'--r')
    title(['Xtile' num2str(ii) '; ' num2str(size(XCCGall{ii},2)) ' DOWNs'])
    
    if ii~=1
        set(gca, 'YTickLabel', [])
    end
    set(gca, 'XTickLabel', [])
    title(['xtile ' num2str(ii)])
    
    % Collection variable; ii = xtile, iii = which variable
    avCCG = XCCGmean_rippleBool(:,ii);
    avCCGAll_rippleBool{ii} = avCCG; %save avCCG for later plotting
    
    subplot(4,xtile,ii + xtile*3)
    yyaxis left
    plot(t_lag,avCCG,'linewidth',1.2,'color','k')
    
    box off; axis tight; yy=get(gca,'ylim');
    xlabel('time(s)')
    if ii==1
        ylabel('P(SWR)')
    end
    hold on
    plot([0 0],[yy(1) yy(2)],'--r')
end
end



%% Done collecting data for each xtile; add remaining stuff

if SWRFeatureSelection
    
    % Plot outline on top of all together x amplitude swr
    for kk=1:xtile
    figure(14)
    subplot(1,xtile,kk)
    hold on
    if ~isempty(boundary)
        patch(boundary(:,2),boundary(:,1),'r','faceAlpha',.2,'edgeColor','r')
        plot([sizeTemp(1)/2 sizeTemp(1)/2],get(gca,'ylim'),'--w')
    end
    end
    %% Adjust ylim so is same across each DOWN
    figure(10)
    % Set ripple bounds to min/max of last
    for ii=1:xtile
        mm(ii)=max(avCCGAll_rippleBool{ii});
        nn(ii)=min(avCCGAll_rippleBool{ii});
    end
    
    for ii=1:xtile
        subplot(6,xtile,ii + xtile*3)
        ylim([ylimzmax(1) ylimzmax(2)])
        plot([0 0],[ylimzmax(1) ylimzmax(2)],'--r')
        if ii~=1
            set(gca, 'YTickLabel', [])
        end
        set(gca, 'XTickLabel', [])
        xlabel([])
        
        % Plot SWRs
        subplot(6,xtile,ii + xtile*4)
        yyaxis left
        plot(t_lag,avCCGAll_rippleBool{ii},'linewidth',1.4,'color','k')
        ylim([min(nn) max(mm)])
        if ii==1
            ylabel('P(SWR)')
        else
            set(gca,'yTick',[])
        end
        set(gca,'xTick',[])
        hold on
        plot([0 0],[min(nn) max(mm)],'--r')
        
        iii=2; %MUA
        % Plot MUA
        yyaxis right
        plot(t_lag,avCCGAll{ii},'color',[0.8500, 0.3250, 0.0980],'linewidth',1)
        x=min(avCCGAll{ii});
        y=max(avCCGAll{ii});
        ylim([x y])
        set(gca,'yTick',[])
        if ii==xtile
            ylabel('MUA(0-1 Norm)')
        end
        
        % P UP
        iii=1;
        subplot(6,xtile,ii + xtile*5)
        plot(t_lag,avCCGAll_SlowWaveBool{ii},'linewidth',1.4,'color','k')
        hold on
        plot([0 0],get(gca,'ylim'),'--r')
        if ii==1
            ylabel('P(UP)')
        else
            set(gca,'yTick',[])
        end
        xlabel('time(s)')
        
    end
    
    suptitle(basePath)
    
end

%% Add P(UP) and MUA to ripples figure
figure(13)
for ii=1:xtile
    subplot(5,xtile,ii+xtile*3)
    plot(t_lag,avCCGAll_SlowWaveBool{ii},'linewidth',1.4,'color','k')
    hold on
    plot([0 0],get(gca,'ylim'),'--r')
    if ii==1
        ylabel('P(UP)')
    else
        set(gca,'yTick',[])
    end
    set(gca,'xTick',[])
    axis tight
end

% set bounds to top xtile
for ii=1:xtile
    mm(ii)=max(avCCGAll_rippleBool{ii});
    nn(ii)=min(avCCGAll_rippleBool{ii});
end

for ii=1:xtile
    subplot(4,xtile,ii+xtile*3)
    hold on
    % Plot MUA
    yyaxis right
    plot(t_lag,avCCGAll{ii},'color',[0.8500, 0.3250, 0.0980],'linewidth',1)
    x=min(avCCGAll{ii});
    y=max(avCCGAll{ii});
    ylim([x y])
    set(gca,'yTick',[])
    if ii==xtile
        ylabel('MUA(0-1 Norm)')
    end
    % Plot SWRs
    yyaxis left
    plot(t_lag,avCCGAll_rippleBool{ii},'linewidth',1.4,'color','k')
    ylim([min(nn) max(mm)])
    if ii==1
        ylabel('P(SWR)')
    else
        set(gca,'yTick',[])
    end
    plot([0 0],[min(nn) max(mm)],'--r')
end

%% %%%%%%%%%%%%%%%%%%

%%
if ~isempty(savePathMouse)
    % Save figures 1-4
    namesave = trigger;
    mkdir(fullfile(savePathMouse,'DOWNCenteredAnalysis','figureSetFx'))
    savePathMouse = fullfile(savePathMouse,'DOWNCenteredAnalysis','figureSetFx');
    figure(1)
    % suptitle(['Mouse' num2str(mice(i)) ': ' num2str(size(CCGall{1,i},2)) namesave])
    NiceSave('AllTogether',savePathMouse,restrictIn)
    figure(12)
    NiceSave([namesave '_UpDownLogical'],savePathMouse,restrictIn)
    figure(11)
    NiceSave([namesave '_MUA'],savePathMouse,restrictIn)
    figure(13)
    NiceSave([namesave '_CumPlotSWR'],savePathMouse,restrictIn)
    figure(5)
    NiceSave([namesave '_SWRsOnly'],savePathMouse,restrictIn)
    figure(41)
    NiceSave([namesave '_CTXandHPCMUA'],savePathMouse,restrictIn)
    
    
    %     figure(30)
    %     NiceSave([namesave '_SynchronyPlots_DowntoUP'],savePathMouse,num2str(mouse))
    if ~isempty(SWRFeatureSelection)
        figure(10)
        NiceSave([namesave '_' SWRFeatureSelection],savePathMouse,restrictIn)
        figure(40)
        NiceSave([namesave '_' SWRFeatureSelection '_ampPlot'],savePathMouse,restrictIn)
    end
end

%% Output structures

% All CCGs
CCGAll.t_lag = t_lag;
CCGAll.MUA = single(squeeze(CCGall_rippleBool(:,sortUse)))';
CCGAll.rippleBool = single(squeeze(CCGall_rippleBool(:,sortUse)))';
CCGAll.rippleIndex = single(squeeze(CCGall_rippleIndex(:,sortUse)))';
CCGAll.slowWaveBool = single(squeeze(CCGall_SlowWaveBool(:,sortUse)))';
CCGAll.boundary = boundary;
CCGAll.durDOWNs = downSort; %duration of DOWNs in same order of CCG

% Xtile CCGs
if ~isempty(xtile)
    CCGAllXtile.t_lag = t_lag;
    CCGAllXtile.MUA = XCCGall;
    CCGAllXtile.rippleBool = Xall_rippleBool;
    CCGAllXtile.rippleIndices = Xall_rippleIndices;
    CCGAllXtile.slowWaveBool = XCCGall_SlowWaveBool;
    CCGAllXtile.averages.MUA = avCCGAll;
    CCGAllXtile.averages.rippleBool = avCCGAll_rippleBool;
    CCGAllXtile.averages.slowWaveBool = avCCGAll_SlowWaveBool;
    CCGAllXtile.boundary = boundaryUPBool;
    CCGAllXtile.durDOWNs = xtileDurs;
else
    CCGAllXtile = [];
end

% Stuff in synchInfo
%     synchInfo(:,1) = downSort; %duration down; sorted
%     synchInfo(:,2) = sum(windowCCG_50ms,2); %50 ms post UP
%     synchInfo(:,3) = sum(windowCCG_20ms,2); %20 ms post UP
%     synchInfo(:,4) = durUP(sortUse); %duration of following UP; sort by same thing as DOWN, works bc = UP following the DOWN
%     temp=sortUse-1;
%     temp=temp(temp>=1);
%     synchInfo(:,5) = [nan; durUP(temp)];%duration preceding UP;
%     %PSS value at DOWN->UP transition
%     %triggersPSS=interp1(timePSS,PSS,triggers);
%     %synchInfo(:,6) = triggersPSS(sortUse); %PSS value at each DOWN state (too short timescale to be relevant)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADD ADDITIONAL FIGURES using structs made
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Summary synchInfo figure - WORKING ON THIS ; need to update to new
% synchInfo

if strcmp(trigger,'DOWNtoUP')
    
    % Plot
    colors1 = autumn(xtile+3); colors1=colors1(3:end,:);
    colors2 = gray(xtile+3); colors2=colors2(2:end,:); %colors2=flip(colors2,1);
    figure(30)
    for ii=1:xtile
        subplot(4,3,1)
        iii=1; %P(UP)
        plot(CCGAll.t_lag,CCGAllXtile.averages.slowWaveBool{ii},'color',colors1(ii,:),'linewidth',1.4)
        hold on
        plot([0 0],get(gca,'ylim'),'--r')
        ylabel('P(UP)')
        box off
        
        
        subplot(4,3,2)
        hold on
        plot(CCGAll.t_lag,CCGAllXtile.averages.MUA{ii},'color',colors1(ii,:),'linewidth',1.4)
        axis tight
        ylabel('MUA')
        xlabel('time(s)')
        %legend(num2str(1:xtile))
        
        subplot(4,3,3)
        hold on
        plot(CCGAll.t_lag,CCGAllXtile.averages.rippleBool{ii},'color',colors1(ii,:),'linewidth',1.4)
        axis tight
        ylabel('P(SWR)')
        xlabel('time(s)')
        %legend(num2str(1:xtile))
        %
        %         subplot(4,3,4)
        %         title('DurDOWN x MUA(20msPost)')
        %         scatter(log10(synchInfo(:,1)),synchInfo(:,3),'.k')
        %         axis tight
        %         LogScale('x',10)
        %         xlabel('DurDown(s)')
        %         ylabel('Sum MUA (20 ms bin)')
        %
        %         subplot(4,3,5)
        %         title('DurDOWN x MUA(50msPost)')
        %         scatter(log10(synchInfo(:,1)),synchInfo(:,2),'.k')
        %         axis tight
        %         LogScale('x',10)
        %         xlabel('DurDown(s)')
        %         ylabel('Sum MUA (50 ms bin)')
        %
        %         subplot(4,3,6)
        %         title('DurDOWN x DurUP(Post)')
        %         scatter(log10(synchInfo(:,1)),log10(synchInfo(:,3)),'.k')
        %         axis tight
        %         LogScale('xy',10)
        %         xlabel('DurDown n (s)')
        %         ylabel('DurUP n+1 (s)')
        %
        %         subplot(4,3,7)
        %         title('DurDOWN x DurUP(Pre)')
        %         scatter(log10(synchInfo(:,1)),log10(synchInfo(:,5)),'.k')
        %         axis tight
        %         LogScale('xy',10)
        %         xlabel('DurDown n (s)')
        ylabel('DurUP n-1 (s)')
        
        % Add in below once figure out what state variable wanna use
        %         subplot(4,3,8)
        %         title('DurDOWN x PSS')
        %         scatter(log10(synchInfo(:,1)),synchInfo(:,6),'.k')
        %         axis tight
        %         LogScale('x',10)
        %         xlabel('DurDown (s)')
        %         ylabel('PSS (light->deep NREM)')
        %
        %         subplot(4,3,9)
        %         title('DurUP(Post) x PSS')
        %         scatter(log10(synchInfo(:,4)),synchInfo(:,6),'.k')
        %         axis tight
        %         LogScale('x',10)
        %         xlabel('DurUP (s)')
        %         ylabel('PSS (light->deep NREM)')
        
        subplot(4,3,10)
        title('Synchrony x P(SWR)')
        
        subplot(4,3,11)
        title('Synchrony x Mag(SWR)')
        
    end
    suptitle('Aligned to DOWN->UP transition')
end

%% Remove duplicates (ripples that fall during UP and DOWN bc right at transition)
% Strategy: remove from 'DOWN' ripples, add to UP

swrs = unique(CCGAll.rippleIndex);
swrs = swrs(swrs ~= 0);
rippleAllDOWNsDeleteTest = rippleAllDOWNs; % duplicate variable bc deleting things

% Make temp matrix w ripple index + category
tempUP=[];
tempUP=[rippleAllUPs(:).indexRipplesAll]';
tempDOWN=[];
tempDOWN=[rippleAllDOWNsDeleteTest(:).indexRipplesAll]';
rippleCat = vertcat(tempUP,tempDOWN);

% Look for duplicates
if length(rippleCat)>length(swrs)
    swrsUnique = unique(rippleCat(:,1));
    [binc,ind] = histc(rippleCat(:,1),swrsUnique);
    x=find(binc==2);
    dups=swrsUnique(x); % SWR number
    
    for i=1:length(dups)
        r=find([rippleAllDOWNsDeleteTest(:).indexRipplesAll]==dups(i));
        disp(r)
        rippleAllDOWNsDeleteTest(r) = [];
    end
    
    % Remake with no duplicates
    tempUP=[];
    tempUP=[rippleAllUPs(:).indexRipplesAll]';
    tempDOWN=[];
    tempDOWN=[rippleAllDOWNsDeleteTest(:).indexRipplesAll]';
    rippleCat = vertcat(tempUP,tempDOWN);
    
    if length(rippleCat) == length(swrs)
        rippleAllDOWNs = rippleAllDOWNsDeleteTest;
        disp(['All ripples accounted for, total swr number: ' num2str(length(rippleCat))])
    end
    
end

%% Add PSS and infra @ each ripple

if ischar(includeState)
    
    % Pull out SWR times (using SWR peak for now); add to SWR structs
    tempUP=[];
    tempUP=[rippleAllUPs(:).SWRpeak]';
    tempDOWN=[];
    tempDOWN=[rippleAllDOWNs(:).SWRpeak]';
    
    % PSS
    UPtriggersPSS=interp1(timePSS,PSS,tempUP);
    DOWNtriggersPSS=interp1(timePSS,PSS,tempDOWN);
    
    % Infra
    %UPtriggersInfra=interp1(timeInfra,phaseInfra,tempUP);
    %DOWNtriggersInfra=interp1(timeInfra,phaseInfra,tempDOWN);
    
    % Add to output vars
    for i=1:length(rippleAllUPs)
        rippleAllUPs(i).PSS = UPtriggersPSS(i);
        %rippleAllUPs(i).infra = UPtriggersInfra(i);
    end
    
    for i=1:length(rippleAllDOWNs)
        rippleAllDOWNs(i).PSS = DOWNtriggersPSS(i);
        %rippleAllDOWNs(i).infra = DOWNtriggersInfra(i);
    end
    
    
end

%% Combine UP and DOWN into one struct

% Combine (same format struct as others)
% indexRipplesAll = vertcat([rippleAllUPs(:).indexRipplesAll]',[rippleAllDOWNs(:).indexRipplesAll]');
% indexUPorDOWNAll = vertcat([rippleAllUPs(:).indexUPAll]',[rippleAllDOWNs(:).indexDOWNAll]');
% SWRdurDOWN = vertcat([rippleAllUPs(:).SWRdurDOWN]',[rippleAllDOWNs(:).SWRdurDOWN]');
% SWRdurUP = vertcat([rippleAllUPs(:).SWRdurUP]',[rippleAllDOWNs(:).SWRdurUP]');
% SWRpeak =vertcat([rippleAllUPs(:).SWRpeak]',[rippleAllDOWNs(:).SWRpeak]');
% SWRstart = vertcat([rippleAllUPs(:).SWRstart]',[rippleAllDOWNs(:).SWRstart]');
% SWRdur = vertcat([rippleAllUPs(:).SWRdur]',[rippleAllDOWNs(:).SWRdur]');
% SWRamp = vertcat([rippleAllUPs(:).SWRamp]',[rippleAllDOWNs(:).SWRamp]');
% UPorDOWNonset = vertcat([rippleAllUPs(:).UPonset]',[rippleAllDOWNs(:).DOWNonset]');
% UPorDOWNoffset = vertcat([rippleAllUPs(:).UPoffset]',[rippleAllDOWNs(:).DOWNoffset]');
% UPorDOWNdur = vertcat([rippleAllUPs(:).UPdur]',[rippleAllDOWNs(:).DOWNdur]');
% closestUPtoDOWNtrans = vertcat([rippleAllUPs(:).closestUPtoDOWNtrans]',[rippleAllDOWNs(:).closestUPtoDOWNtrans]');
% closestDOWNtoUPtrans = vertcat([rippleAllUPs(:).closestDOWNtoUPtrans]',[rippleAllDOWNs(:).closestDOWNtoUPtrans]');
% if ischar(sharpWavesIn) % ripple features
%     SWmag = vertcat([rippleAllUPs(:).SWmag]',[rippleAllDOWNs(:).SWmag]');
% end
%
% for rr=1:length(indexRipplesAll);
%     rippleAll(rr).indexRipplesAll = indexRipplesAll(rr);
%     rippleAll(rr).indexUPorDOWNAll = indexUPorDOWNAll(rr);
%     rippleAll(rr).SWRdurDOWN = SWRdurDOWN(rr);
%     rippleAll(rr).SWRdurUP = SWRdurUP(rr);
%     rippleAll(rr).SWRpeak = SWRpeak(rr);
%     rippleAll(rr).SWRstart = SWRstart(rr);
%     rippleAll(rr).SWRdur = SWRdur(rr);
%     rippleAll(rr).SWRamp = SWRamp(rr);
%
%     rippleAll(rr).UPorDOWNonset = UPorDOWNonset(rr);
%     rippleAll(rr).UPorDOWNoffset = UPorDOWNoffset(rr);
%     rippleAll(rr).UPorDOWNdur = UPorDOWNdur(rr);
%     rippleAll(rr).closestUPtoDOWNtrans = closestUPtoDOWNtrans(rr);
%     rippleAll(rr).closestDOWNtoUPtrans = closestDOWNtoUPtrans(rr);
%     if ischar(sharpWavesIn) % ripple features
%         rippleAll(rr).SWmag = SWmag(rr);
%     end
% end

% Combine
rippleAll.indexRipplesAll = vertcat([rippleAllUPs(:).indexRipplesAll]',[rippleAllDOWNs(:).indexRipplesAll]');
rippleAll.indexUPorDOWNAll = vertcat([rippleAllUPs(:).indexUPAll]',[rippleAllDOWNs(:).indexDOWNAll]');
rippleAll.SWRdurDOWN = vertcat([rippleAllUPs(:).SWRdurDOWN]',[rippleAllDOWNs(:).SWRdurDOWN]');
rippleAll.SWRdurUP = vertcat([rippleAllUPs(:).SWRdurUP]',[rippleAllDOWNs(:).SWRdurUP]');
rippleAll.SWRpeak =vertcat([rippleAllUPs(:).SWRpeak]',[rippleAllDOWNs(:).SWRpeak]');
rippleAll.SWRstart = vertcat([rippleAllUPs(:).SWRstart]',[rippleAllDOWNs(:).SWRstart]');
rippleAll.SWRdur = vertcat([rippleAllUPs(:).SWRdur]',[rippleAllDOWNs(:).SWRdur]');
rippleAll.SWRamp = vertcat([rippleAllUPs(:).SWRamp]',[rippleAllDOWNs(:).SWRamp]');
rippleAll.UPorDOWNonset = vertcat([rippleAllUPs(:).UPonset]',[rippleAllDOWNs(:).DOWNonset]');
rippleAll.UPorDOWNoffset = vertcat([rippleAllUPs(:).UPoffset]',[rippleAllDOWNs(:).DOWNoffset]');
rippleAll.UPorDOWNdur = vertcat([rippleAllUPs(:).UPdur]',[rippleAllDOWNs(:).DOWNdur]');
rippleAll.closestUPtoDOWNtrans = vertcat([rippleAllUPs(:).closestUPtoDOWNtrans]',[rippleAllDOWNs(:).closestUPtoDOWNtrans]');
rippleAll.closestDOWNtoUPtrans = vertcat([rippleAllUPs(:).closestDOWNtoUPtrans]',[rippleAllDOWNs(:).closestDOWNtoUPtrans]');
if ischar(sharpWavesIn) % ripple features
    rippleAll.SWmag = vertcat([rippleAllUPs(:).SWmag]',[rippleAllDOWNs(:).SWmag]');
end
if ischar(includeState)
    rippleAll.PSS = vertcat([rippleAllUPs(:).PSS]',[rippleAllDOWNs(:).PSS]');
end
