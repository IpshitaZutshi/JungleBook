% Compile data across all sessions

function compiledRipples = compileRipples(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'analogEv',64,@isnumeric);
addParameter(p,'numAnalog',2,@isnumeric);
addParameter(p,'plotfig',false,@islogical);
parse(p,varargin{:});

expPath = p.Results.expPath;
analogEv = p.Results.analogEv;
numAnalog = p.Results.numAnalog;
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

compiledRipples.rate_pre = [];
compiledRipples.rate_prestim = [];
compiledRipples.rate_post = [];
compiledRipples.rate_poststim = [];
compiledRipples.rateAnalog = [];

%% Start collecting data
for ii = 1:size(allSess,1)

    cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
    [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);    

    if exist([sessionInfo.FileName '.hippocampalLayers.channelinfo.mat'],'file') 
        load([sessionInfo.FileName '.hippocampalLayers.channelinfo.mat']);
        pyrCh = hippocampalLayers.pyramidal;
    else
        disp('First calculate .region file to identify ripple channel! Skipping');
        return
    end
    
    load([sessionInfo.FileName '.session.mat']);
    if isfield(session.channelTags,'RippleNoise')
        noiseCh = session.channelTags.RippleNoise.channels-1;
    else
        noiseCh = [];
    end
    
    [ripples] = bz_FindRipples(pwd,pyrCh,'noise',noiseCh,'savemat',true);
    passband = [130 200];

    lfp = bz_GetLFP(pyrCh,'noPrompts', true);
    signal = bz_Filter(double(lfp.data),'filter','butter','passband',passband,'order', 3);
    timestamps = lfp.timestamps;
    [~,data] = bz_RippleStats(signal,timestamps,ripples);

    % Load pulses
    disp('Getting analog-in inputs...');
    [pulses] = bz_getAnalogPulsesSine('analogCh',analogCh);


    for i = 1:(numAnalog+1)
        if exist('pulses')
            if i<=numAnalog
                pulTr = (pulses.stimComb==i);
            else
                pulTr = (pulses.stimPerID'==1 & pulses.stimComb==i);
            end
        end

        events = pulses.intsPeriods(1,pulTr);
        events = events(((events + 5) <= max(ripples.peaks)) & ((events - 5) > 0));

        %Generate logicals for ripples in pre versus post
        ripple_pre = [];
        ripple_prestim = [];
        ripple_post = [];
        ripple_poststim = [];  
        psthAnalogEv= [];
        
        ripple_pre(1:length(ripples.peaks)) = 0;
        ripple_prestim(1:length(ripples.peaks)) = 0;
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
            elseif min(abs(tempDiff)) >=5 & min(abs(tempDiff)) < 10
                [~,idxmin] =  min(abs(tempDiff));
                if tempDiff(idxmin) > 0
                    ripple_poststim(pp) = 1;
               elseif tempDiff(idxmin) < 0
                   ripple_prestim(pp) = 1;
               end
            else
                continue
            end
        end

        ripple_pre = logical(ripple_pre);
        ripple_post = logical(ripple_post);
        ripple_poststim = logical(ripple_poststim);
        ripple_prestim = logical(ripple_prestim);
        psthAnalogEv(1:length(ripples.peaks)) = i;
        
        compiledRipples.ripple_pre = [compiledRipples.ripple_pre ripple_pre];
        compiledRipples.ripple_prestim = [compiledRipples.ripple_prestim ripple_prestim];
        compiledRipples.ripple_post = [compiledRipples.ripple_post ripple_post];
        compiledRipples.ripple_poststim = [compiledRipples.ripple_poststim ripple_poststim];
        compiledRipples.peakAmplitude = [compiledRipples.peakAmplitude data.peakAmplitude'];
        compiledRipples.duration = [compiledRipples.duration data.duration'];
        compiledRipples.peakFrequency = [compiledRipples.peakFrequency data.peakFrequency'];
        compiledRipples.psthAnalogEv = [compiledRipples.psthAnalogEv psthAnalogEv];
        compiledRipples.rate_pre = [compiledRipples.rate_pre sum(ripple_pre)./(numel(events)*5)];
        compiledRipples.rate_prestim = [compiledRipples.rate_prestim sum(ripple_prestim)./(numel(events)*5)];
        compiledRipples.rate_post = [compiledRipples.rate_post sum(ripple_post)./(numel(events)*5)];
        compiledRipples.rate_poststim = [compiledRipples.rate_poststim sum(ripple_poststim)./(numel(events)*5)];
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
