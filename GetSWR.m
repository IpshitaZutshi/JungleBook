function ripples = GetSWR(LFPs,opt)

%% default parameters
if ~exist('opt','var')
    opt.rippleband = [140,250];
    opt.threshold = [2,5];
    opt.duration = [20,200];
    opt.mergethr = 15;
    opt.nosiecorrthr = 0.9;
end

%% loa and filter LFPs and Get ripple band power
ripplelfp = bz_Filter(LFPs,'passband',opt.rippleband,'filter','butter');
zscoredpower = zscore(ripplelfp.amp',[],2);
[pmax,icmax] = max(zscoredpower,[],1);
pmax(1) = 0; pmax(end) = 0;

%% find events
epochs(:,1) = find(diff(pmax>opt.threshold(1)) == 1)+1;
epochs(:,2) = find(diff(pmax>opt.threshold(1)) == -1);
% merge close epochs
merie = find(epochs(2:end,1)-epochs(1:end-1,2) < opt.mergethr);
epochs(merie,2) = epochs(merie+1,2);
epochs(merie+1,:) = [];
% remove epochs that peak amplitude lower
rmlist = []; dt = 1/LFPs.samplingRate;
for ie = 1:size(epochs,1)
    [peakamp,peakidx] = max(pmax(epochs(ie,1):epochs(ie,2)));
    epochs(ie,3) = epochs(ie,1) + peakidx -1;
    epochs(ie,4) = icmax(epochs(ie,3));
    if peakamp < opt.threshold(2), rmlist = [rmlist,ie]; end
end
epochs(rmlist,:) = [];
% remove too short or too long epochs
epochs(diff(epochs(:,1:2),1,2)*dt*1000<opt.duration(1) | diff(epochs(:,1:2),1,2)*dt*1000>opt.duration(2),:) = [];

%% remove high correlation epochs
rmlist = [];
for ie = 4:size(epochs,1)
    rho = corr(double(LFPs.data(epochs(ie,1):epochs(ie,2),epochs(ie,4))),double(LFPs.data(epochs(ie,1):epochs(ie,2),:)));
    if mean(rho)>0.9, rmlist = [rmlist,ie]; end
end
epochs(rmlist,:) = [];
% Exclude ripples towards the beginning or end of the file
a = epochs(:,1)<5*1250;
epochs(a,:) = [];

a = epochs(:,1)>length(LFPs.timestamps)-5*1250;
epochs(a,:) = [];

ripples.timestamps = LFPs.timestamps(epochs(:,1:2));
ripples.peaks = LFPs.timestamps(epochs(:,3));
ripples.epochsidx = epochs;


%% manually select ripple
% init = LAMP_Initialization;
% basepath = LFPs.basepath; basepath(1:3) = init.root;
basename = bz_BasenameFromBasepath(pwd);
Inspector_CheckRipples(LFPs,ripples,pwd);

%% save events for neuroscope
eventtype = 'ripples'; %string for you that = event type you're saving
numeventtypes = 3; % you have 3 diff types of events, start, peak, stop

% Save as .evt for inspection
ripbasename = [basename '.' eventtype '.RIP.evt']; % you need the ROX bc neuroscope is buggy and uses this to parse files.
delete(fullfile(basepath,ripbasename));

% Below is for ripples specifically 
% Populate events.time field
lengthAll = numel(ripples.peaks)*numeventtypes;
events.time = zeros(1,lengthAll);
events.time(1:3:lengthAll) = ripples.timestamps(:,1);
events.time(2:3:lengthAll) = ripples.peaks;
events.time(3:3:lengthAll) = ripples.timestamps(:,2);

% Populate events.description field
events.description = cell(1,lengthAll);
events.description(1:3:lengthAll) = {'start'};
events.description(2:3:lengthAll) = {'peak'};
events.description(3:3:lengthAll) = {'stop'};

% Save .evt file for viewing in neuroscope - will save in your current directory
if exist(fullfile(basepath,ripbasename),'file'), delete(fullfile(basepath,ripbasename)); end
SaveEvents(fullfile(basepath,ripbasename),events) %Save and look at in neuroscope

end