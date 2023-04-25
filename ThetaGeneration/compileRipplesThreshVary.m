% Compile data across all sessions
% Final ripple analysis 
% By IZ - 2/18/2022

% RippleMasterDetector MUST have been run already!

function compiledRipples = compileRipplesThreshVary(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'minEvents',50,@isnumeric); % each session should have atleast 50 ripples to be included in the stats
addParameter(p,'plotfig',false,@islogical);
parse(p,varargin{:});

expPath = p.Results.expPath;
minEvents = p.Results.minEvents;
plotfig = p.Results.plotfig;

%%
if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end
allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');

for ii = 1:8
    compiledRipples.rate{ii} = [];
    if ii <3
        compiledRipples.dur{ii} = [];
        compiledRipples.amp{ii} = [];
    end
end

%% Start collecting data
for ii = 1:size(allSess,1)

    cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
    [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);     

    if exist([sessionInfo.FileName '.ripples_restrict.events.mat'],'file') 
       load([sessionInfo.FileName '.ripples_restrict.events.mat']);
    else
       disp('First calculate .ripple file! Skipping');
       return
    end
    
    if exist([sessionInfo.FileName '.ripples_thresh.events.mat'],'file') 
       load([sessionInfo.FileName '.ripples_thresh.events.mat']);
    else
       disp('First calculate .ripple file! Skipping');
       return
    end
    
    % Load pulses
    disp('Getting analog-in inputs...');
    load([sessionInfo.FileName '.pulses.events.mat']);

    %Skip session if too few ripples
    if length(ripplesRestrict{2}.peaks)<minEvents
        continue
    end
    
    pulTr = (pulses.stimComb==2);
    
    %Only select the pulses that happened in the home cage,
    %i.e., which were 5 seconds long
    homeCagePulse = pulses.intsPeriods(2,:) - pulses.intsPeriods(1,:);
    homeCagePulseidx = homeCagePulse < 5.05 & homeCagePulse > 4.95;
    pulTr = pulTr & homeCagePulseidx;
    events = pulses.intsPeriods(:,pulTr)';

    ripple_pre = InIntervals(ripplesRestrict{1}.peaks,(events-5));
    ripple_post = InIntervals(ripplesRestrict{2}.peaks,events);

    rippleT1_pre = InIntervals(ripplesThresh{1}.peaks,(events-5));
    rippleT1_post = InIntervals(ripplesThresh{1}.peaks,events);
    
    rippleT2_pre = InIntervals(ripplesThresh{2}.peaks,(events-5));
    rippleT2_post = InIntervals(ripplesThresh{2}.peaks,events);
    
    rippleT3_pre = InIntervals(ripplesThresh{3}.peaks,(events-5));
    rippleT3_post = InIntervals(ripplesThresh{3}.peaks,events);    
    
    if sum(ripple_pre)<2
        continue
    end
    
    compiledRipples.rate{1} = [compiledRipples.rate{1} sum(ripple_pre)./(numel(events)*5)];
    compiledRipples.rate{2} = [compiledRipples.rate{2} sum(ripple_post)./(numel(events)*5)];
    compiledRipples.rate{3} = [compiledRipples.rate{3} sum(rippleT1_pre)./(numel(events)*5)];
    compiledRipples.rate{4} = [compiledRipples.rate{4} sum(rippleT1_post)./(numel(events)*5)];
    compiledRipples.rate{5} = [compiledRipples.rate{5} sum(rippleT2_pre)./(numel(events)*5)];
    compiledRipples.rate{6} = [compiledRipples.rate{6} sum(rippleT2_post)./(numel(events)*5)];
    compiledRipples.rate{7} = [compiledRipples.rate{7} sum(rippleT3_pre)./(numel(events)*5)];
    compiledRipples.rate{8} = [compiledRipples.rate{8} sum(rippleT3_post)./(numel(events)*5)];    

    dur = ripplesRestrict{1}.timestamps(:,2)-ripplesRestrict{1}.timestamps(:,1);  
    compiledRipples.dur{1} = [compiledRipples.dur{1}; dur(ripple_pre)];
    
    dur = ripplesRestrict{2}.timestamps(:,2)-ripplesRestrict{2}.timestamps(:,1);  
    compiledRipples.dur{2} = [compiledRipples.dur{2}; dur(ripple_post)];    
end

save([expPath '\Summ\' 'compiledRipplesThresh.mat'], 'compiledRipples');


