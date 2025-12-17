function [sessName, countStim] = estimateDeliberation_mPFC(varargin)
% Compare # deliberation events between stim==0 and stim==1 (and by type).
%
% Outputs:
%   sessName     [nSess x 1]  session labels
%   countStim    [nSess x 2 x 3] raw counts: (session, stim(0/1), type)

p = inputParser;
addParameter(p,'expPath','Z:\Homes\zutshi01\Recordings\Auditory_Task',@ischar);
addParameter(p,'sess',{},@(x) iscell(x) && ~isempty(x));
addParameter(p,'plotfig',true,@islogical);
parse(p,varargin{:});

expPath = p.Results.expPath;
sess    = p.Results.sess;
plotfig = p.Results.plotfig;

if isempty(sess)
    sess = {'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17', ...
            'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16', ...
            'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15', ...
            'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15'};
end

nSess = numel(sess);
sessName = string(sess(:));

countStim = zeros(nSess,2,3); % (session, stimIndex(1=0,2=1), type)

for ii = 1:nSess
    basepath = fullfile(expPath, sess{ii});
    if ~isfolder(basepath)
        warning('Missing folder: %s (skipping)', basepath);
        continue
    end

    Dec = findDecelerationPoints_mPFC('basepath', basepath, 'plotfig', false);

    if ~isfield(Dec,'decType') || isempty(Dec.decType) || ~isfield(Dec,'stim') || isempty(Dec.stim)
        continue
    end

    dType = Dec.decType(:); % should be 1/2/3
    stim  = Dec.stim(:);    % should be 0/1

    % Count by stim and type
    for t = 1:3
        countStim(ii,1,t) = sum(stim==0 & dType==t); % stim=0
        countStim(ii,2,t) = sum(stim==1 & dType==t); % stim=1
    end
end


if plotfig
    % 1) Per-session comparison: stim 0 vs stim 1 (total events)
    figure('Color','w');
    for decType = 1:3
        subplot(1,3,decType)
        stats{decType} = groupStats({countStim(:,1,decType),countStim(:,2,decType)},[],'inAxis',true,'plotType','BoxLinesSEM','repeatedMeasures',true);
    end
end

end
