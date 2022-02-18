function calculateICAwrapper(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'limitTime',false,@islogical);
addParameter(p,'force',true,@islogical);
parse(p,varargin{:});

expPath = p.Results.expPath;
force = p.Results.force;
limitTime = p.Results.limitTime;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');

for ii = 2:size(allSess,1)
    fprintf(' ** Examining session %3.i of %3.i... \n',ii, size(allSess,1));
    cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
    [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
    load([sessionInfo.FileName '.MergePoints.events.mat']);
    file = dir(('*.hippocampalLayers.channelinfo.mat'));
    load(file.name);  
    pyrCh = hippocampalLayers.pyramidal; 
    layerInfo = hippocampalLayers.all;    

    %% Channels    
    for ch = 1:size(sessionInfo.AnatGrps,2)
        if ismember(pyrCh, sessionInfo.AnatGrps(ch).Channels)
            channelOrder = sessionInfo.AnatGrps(ch).Channels; 
        end
    end    
        %% Specify chanRange
    [~,layerInfoIdx] = ismember(layerInfo,channelOrder);
    startCh = layerInfoIdx(2)-2;
    endCh = min(length(channelOrder),layerInfoIdx(4)+2);
    chanRange = [startCh:endCh];
    channelOrder = channelOrder(chanRange);
    
    if limitTime
        intervals = MergePoints.timestamps([1 2 4],:);
        lfp = bz_GetLFP(channelOrder,'restrict',intervals);
        %concatenate lfp
        lfpCat.Filename = lfp(1).Filename;
        lfpCat.samplingRate = lfp(1).samplingRate;
        lfpCat.channels = lfp(1).channels;
        lfpCat.data = [];
        lfpCat.interval = [];
        lfpCat.timestamps = [];
        lfpCat.duration = 0;
        for int = 1:size(lfp,2)
            lfpCat.data = [lfpCat.data; lfp(int).data];
            lfpCat.timestamps = [lfpCat.timestamps; lfp(int).timestamps];
            lfpCat.interval = [lfpCat.interval;lfp(int).interval];
            lfpCat.duration = lfpCat.duration+lfp(int).duration;
        end
        clear lfp
        lfp = lfpCat;
        clear lfpCat
    else
        lfp = bz_GetLFP(channelOrder);
    end
    
    %% Calculate ICA
    bz_RunIca('lfp',lfp,'nICs',10,'force',force);
    close all
end
end