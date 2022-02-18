
function [sessionPulses] = getSessionPulses(varargin)
% Compute and determine optogenetic stimulation during behavior corresponding to session arm choice over all session subfolders
%
% USAGE
%
%   [sessionPulses] = getSessionPulses(varargin)
%
% INPUTS
% basePath                      (default: pwd) basePath for the recording file, 
%                                    in buzcode format:
% task                          'alternation' and 'cudeSide'
% forceReload                   Force detection (boolean, default false)
% verbose                       Default false
% saveMat                       Default true
%
% OUTPUT
%       - sessionPulses.(subSessionFolder) output structure, with the fields:
% sessionPulses.timestamps          Pulse timestamps, per trial in seconds
% sessionPulses.stim                Stimulation logical, 0 is off, 1 is on
% sessionPulses.region              If multiple analog pulse inputs,
%                                       identify which, 0 is off, 1 is analog 64, 2 is 65, 3 is both
% sessionPulses.target              Zone of the maze targeted, 0 is off, 1
%                                       is stem, 2 is return
% sessionPulses.delay               If the behavior had a delay, 0, 5 or 10
% sessionPulses.performanceBlocks   Performance vector in blocks of 10, row 1 is whether the block was stim or not [1,0], row 2 is the performance
% sessionPulses.performance             Total performance, row 1 is stim or not [1,0], row 2 is the performance 

% Ipshita Zutshi 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'task','alternation',@ischar);
addParameter(p,'force',true,@islogical)
addParameter(p,'verbose',false,@islogical)
addParameter(p,'saveMat',true,@islogical)
parse(p,varargin{:});
task = p.Results.task;
forceReload = p.Results.force;
basepath = p.Results.basepath;
verbose = p.Results.verbose;
saveMat = p.Results.saveMat;


%% Deal with inputs
if ~isempty(dir([basepath filesep '*.SessionPulses.Events.mat'])) && ~forceReload
    disp('Session pulses already detected! Loading file.');
    file = dir([basepath filesep '*.SessionPulses.Events.mat']);
    load(file.name);
end

%% Find subfolder recordings
cd(basepath);
[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);
C = strsplit(sessionInfo.session.name,'_');
sess = dir(strcat(C{1},'_',C{2},'*')); % get session files
count = 0;
for ii = 1:size(sess,1)
    
    if sess(ii).isdir && ~isempty(dir([basepath filesep sess(ii).name filesep '*Basler*avi']))
        cd([basepath filesep sess(ii).name]);
        fprintf('Computing session pulses in %s folder \n',sess(ii).name);
        sessionPulses.(sess(ii).name)= getArmPulses('verbose',verbose,'task',task);
        count = count + 1;
    end
end
cd(basepath);

figure
flds = fieldnames(sessionPulses);
subplot(size(flds,1),1,1)
plot(sessionPulses.(flds{1}).performanceBlocks(1,:),'b')
hold on
plot(sessionPulses.(flds{1}).performanceBlocks(2,:),'k')
title(strcat('Region of stim = ',num2str(sessionPulses.(flds{1}).region), '  Zone targeted = ',num2str(sessionPulses.(flds{1}).target)));
if count > 1
    subplot(size(flds,1),1,2)
    plot(sessionPulses.(flds{2}).performanceBlocks(1,:),'b')
    hold on
    plot(sessionPulses.(flds{2}).performanceBlocks(2,:),'k')
    title(strcat('Region of stim = ',num2str(sessionPulses.(flds{2}).region), '  Zone targeted = ',num2str(sessionPulses.(flds{2}).target)));
end

if saveMat
    save([basepath filesep sessionInfo.FileName '.SessionPulses.Events.mat'],'sessionPulses');
end

end

