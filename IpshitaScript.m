% Using the ContinuousDataByDOWNXtile function

% 0. Things you must have detected:
% ripples, DOWN states, behavioral states

% 1. MUA OR population rate is required input for this function for one of the figures
%   (could make this optional, let me know if you want me to). You can
%   either detect MUA and save it as a variable in your basePath (option 1a, time
%   consuming), or you can create a population rate vector out of your
%   spikes.mat file and provide it as input to the function (option 1b, recommend). 

%   Option 1a: Detect MUA. Create MUA.lfp.mat file and save in basePath.
%   Change input MUAIn for ContinuousDataByDOWNXtile fx to a string that is
%   name or partial name of this file.
%   Note: this might be really big/time consuming if long recordings)
    
%     % Input
%     SWChan = %channel you want to use for extraction of MUA
%     basePath = %enter basePath
%     [ MUA ] = MUAfromDat( basePath,'channels',SWChan);

%   Option 1b: create population rate vector from your spikes.mat file and
%   provide as input to the ContinuousDataByDOWNXtile fx (input MUAIn can
%   be a string referring to the saved MUA file or a struct, for further
%   description see ContinuousDataByDOWNXtile fx). I include below how I
%   would make a popRate vector with spikes.mat file.
 expPath = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';
     basePaths = {'Z:\Homes\zutshi01\Recordings\CA1_silencing\IZ12\IZ12_288um_200121_sess2','Z:\Homes\zutshi01\Recordings\CA1_silencing\IZ12\IZ12_288um_200122_sess3','Z:\Homes\zutshi01\Recordings\CA1_silencing\IZ12\IZ12_288um_200127_sess4','Z:\Homes\zutshi01\Recordings\CA1_silencing\IZ12\IZ12_432um_200129_sess5','Z:\Homes\zutshi01\Recordings\CA1_silencing\IZ12\IZ12_576um_200131_sess6',...
    'Z:\Homes\zutshi01\Recordings\CA1_silencing\IZ13\IZ13_216um_200304_sess3','Z:\Homes\zutshi01\Recordings\CA1_silencing\IZ13\IZ13_360um_200306_sess4','Z:\Homes\zutshi01\Recordings\CA1_silencing\IZ13\IZ13_504um_200309_sess6','Z:\Homes\zutshi01\Recordings\CA1_silencing\IZ16\IZ16_144um_200616_sess1','Z:\Homes\zutshi01\Recordings\CA1_silencing\IZ18\IZ18_0um_200707_sess1',...
    'Z:\Homes\zutshi01\Recordings\CA1_silencing\IZ24\IZ24_0um_200903_sess1','Z:\Homes\zutshi01\Recordings\CA1_silencing\IZ24\IZ24_288um_200915_sess3','Z:\Homes\zutshi01\Recordings\CA1_silencing\IZ26\IZ26_0um_201003_sess1','Z:\Homes\zutshi01\Recordings\CA1_silencing\IZ26\IZ26_0um_201006_sess2',...
    'Z:\Homes\zutshi01\Recordings\CA1_silencing\IZ26\IZ26_0um_201015_sess3','Z:\Homes\zutshi01\Recordings\CA1_silencing\IZ27\IZ27_0um_201015_sess2','Z:\Homes\zutshi01\Recordings\CA1_silencing\IZ33\IZ33_0um_210222_sess1','Z:\Homes\zutshi01\Recordings\CA1_silencing\IZ34\IZ34_0um_210222_sess1'};
%   
%     for ii = 1:length(basePaths)
%         baseName = bz_BasenameFromBasepath(strcat(expPath,basePaths{ii}));
%         Input
%         load(fullfile(basePaths{ii},[baseName '.spikes.cellinfo.mat']))
%         binsize = .03;
%         overlap = 6;
%         lowThresh = .2;
% %         Fx
%         spikemat = bz_SpktToSpkmat(tmp_spikes, 'binsize', binsize,'overlap',overlap);
%         spikemat = bz_SpktToSpkmat(spikes, 'binsize', binsize,'overlap',overlap);
%         Store rate and moving average of rate for each region
%         MUA.data = sum(spikemat.data,2);
%         MUA.timestamps = spikemat.timestamps;
%         MUA.samplingRate = 30000;
%         savefile = fullfile(basePaths{ii},[baseName,'.MUA.lfp.mat']);
%         save(savefile,'MUA','-v7.3');
%    end
    
% 2. Make plots with ContinuousDataByDOWNXtile fx
    %Parms
     %'Z:\Homes\zutshi01\Recordings\CA1_silencing\IZ13\IZ13_216um_200304_sess3';%{'H:\mouse47\CumulativeVars'}; % If you include more than one, might run out of memory - try one for now
    MUAIn = 'MUA.lfp.mat'; % input could also be the string - 'MUA.lfp.mat' - if you've saved MUA variable in your basePath
    SlowWavesIn = 'SlowWaves.events.mat';%'SlowWavesMUA.events.mat'; % could be 'SlowWaves.events.mat'
    trigger = 'UPtoDOWN'; %trigger    'DOWNtoUP';
    timelag = [.75 .75]; %time surrounding trigger
    restrictIn = 'NREM';%'all'; %options: NREM, all, WAKE
    xtile = 5; % for making plots - only works w 5 right now
    lowerBoundCutoff = .08;%.08,.01,[]; %only take DOWNs this length and above, leave upperBoundCutoff empty
    upperBoundCutoff = [];%.08,[] %only take DOWNs this length and below, leave lowerBoundCutoff empty
    SWRFeatureSelection = 'SWRamp'; %SWamp
    %sharpWavesIn = 'SharpWaves.events.mat'; %can include as input if want to separte by SW magnitude
    ripplesIn = '.ripples.events.mat';
    includeState = 'PSS';

    % 2a. No separation by SWR feature 
    [CCGAll,CCGAllXtile,synchInfo,rippleAllDOWNs,rippleAllUPs,rippleAll,DOWNInfo] = ContinuousDataByDOWNXtile(basePaths,MUAIn,SlowWavesIn,trigger,timelag,'restrictIn',restrictIn,'xtile',xtile,'lowerBoundCutoff',lowerBoundCutoff,'upperBoundCutoff',upperBoundCutoff,'ripplesIn',ripplesIn);

%     % 2b. Separation by SWR feature 
%     %[CCGAll,CCGAllXtile,synchInfo,rippleAllDOWNs,rippleAllUPs,rippleAll,DOWNInfo] = ContinuousDataByDOWNXtile(basePaths,MUAIn,SlowWavesIn,trigger,timelag,'restrictIn',restrictIn,'xtile',xtile,'lowerBoundCutoff',lowerBoundCutoff,'upperBoundCutoff',upperBoundCutoff,'ripplesIn',ripplesIn,'SWRFeatureSelection',SWRFeatureSelection,'sharpWavesIn',sharpWavesIn);
