function lfp = bz_interpolateLFP(lfp, varargin)

% USAGE:
%     bz_interpolateLFP(lfp, badChannels, refChan)
%     This function interpolates bad channels by calculating the mean
%     between its flanking channels. It first looks for a user input vector
%     of channels, if not provided, looks within the basepath for a  .mat
%     'badChannels', ad if not present, skips the interpolation
%     It also re-references all channels by subtracting one channel(the 'reference' channel)
%     from other channels, if a reference channel is provided. 
%     It then creates a new .lfp file, with the tag '_corrected' in the
%     basepath.
% INPUT
%     lfp            a buzcode structure with fields lfp.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%
%     Note: this uses 0-based indexing

% TO DO: Find a solution for consecutive bad channels

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',[],@isfolder);
addParameter(p,'badChannels',[]);
addParameter(p,'refChan',[],@isnumeric);
parse(p,varargin{:});

basepath = p.Results.basepath;
badChannels= p.Results.badChannels;
refChan = p.Results.refChan;

if isempty(basepath)
    basepath = pwd;
end

if isempty(badChannels)
    if exist([basepath filesep 'badChannels.mat'],'file')
        load([basepath filesep 'badChannels.mat'])
    end
end

[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);

if ~isempty(badChannels)
    for ii = 1:length(badChannels)
        idxChan = [];
        chanPre = [];
        chanPost = [];
        curChan = badChannels(ii);
        for ch = 1:size(sessionInfo.AnatGrps,2)
            idxChan = find(sessionInfo.AnatGrps(ch).Channels==curChan);     
            if ~isempty(idxChan) && idxChan > 1 && idxChan <length(sessionInfo.AnatGrps(ch).Channels) % Find flanking channels
                chanPre = sessionInfo.AnatGrps(ch).Channels(idxChan-1);
                chanPost = sessionInfo.AnatGrps(ch).Channels(idxChan+1);
            end
        end
        if ~isempty(chanPre)
            idBadChan = find(lfp.channels==curChan);   
            idPreChan = find(lfp.channels==chanPre);  
            idPostChan = find(lfp.channels==chanPost);  
            if ~isempty(idBadChan)
                lfp.data(:,idBadChan) = (lfp.data(:,idPreChan)+ lfp.data(:,idPostChan))/2;     
            end
        end
    end
end

if ~isempty(refChan)
    idrefChan = find(lfp.channels==refChan);   
    lfp.data = lfp.data - repmat(lfp.data(:,idrefChan),1,size(lfp.data,2));
end


end