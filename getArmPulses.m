
function [armPulses] = getArmPulses(varargin)
% Compute T maze performance and according to TTL pulses 
%
% USAGE
%
%   [armPulses] = getArmPulses(varargin)
%
% INPUTS
% 
% armChoice                     structure which has computed the
%                               performance for the session
% pulses                        structure with analog input TTLs
%                               If not called, look for it.
%                               If not called, look for it.
% task                          'alternation' and 'cudeSide'
% force                         Force detection (boolean, default false)
% verbose                       default false
%
% OUTPUT
%       - armPulses output structure, with the fields:
% armPulses.timestamps          Pulse timestamps, per trial in seconds
% armPulses.stim                Stimulation logical, 0 is off, 1 is on
% armPulses.region              If multiple analog pulse inputs,
%                                       identify which, 0 is off, 1 is analog 64, 2 is 65, 3 is both
% armPulses.target              Zone of the maze targeted, 0 is off, 1
%                                       is stem, 2 is return
% armPulses.delay               If the behavior had a delay, 0, 5 or 10
% armPulses.performanceBlocks   Performance vector in blocks of 10, row 1 is whether the block was stim or not [1,0], row 2 is the performance
% armPulses.performance             Total performance, row 1 is stim or not [1,0], row 2 is the performance 

% Ipshita Zutshi 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Parms
p = inputParser;
addParameter(p,'armChoice',[]);
addParameter(p,'pulses',[]);
addParameter(p,'task','alternation');
addParameter(p,'force',true,@islogical)
addParameter(p,'verbose',false,@islogical)
parse(p,varargin{:});
armChoice = p.Results.armChoice;
pulses = p.Results.pulses;
task = p.Results.task;
force = p.Results.force;
verbose = p.Results.verbose;

if ~isempty(dir('*.ArmPulses.Events.mat')) && ~force 
    disp('Arm pulses already computed! Loading file.');
    file = dir('*.ArmPulses.Events.mat');
    load(file.name);
    return
end

%% Get TTL pulses
if isempty(pulses)
    disp('Loading pulses...');
    currdir = pwd;
    C = strsplit(pwd,'\');
    newdir = fullfile(C{1:(length(C)-1)});
    cd(newdir)
    [pulses] = bz_getAnalogPulses;
    if isempty(pulses)
        armPulses = [];
        return
    end
    cd(currdir)
end

if isempty(armChoice)
    disp('Loading SessionArmChoice...')
    currdir = pwd;
    C = strsplit(pwd,'\');
    newdir = fullfile(C{1:(length(C)-1)});
    cd(newdir)
    sessionArmChoice = getSessionArmChoice;
    armChoice = sessionArmChoice.(C{length(C)});
    if isempty(pulses)
        armPulses = [];
        return
    end
    cd(currdir)
end
%%
if strcmpi(task,'cueSide')
    
    % TO DO
    
elseif strcmpi(task,'alternation')
    
    armPulses.timestamps(1:length(armChoice.timestamps)) = 0;
    armPulses.stim(1:length(armChoice.timestamps)) = 0;
    armPulses.delay = armChoice.delay.dur;

    for ii = 1:length(pulses.timestamps)
        % first check if the pulses ON timestamp is aligned to the reward input
        tempdiffON = abs(armChoice.timestamps - pulses.timestamps(ii,1)); %armChoice.timestamps has the timestamps of the reward location
        tempdiffOFF = abs(armChoice.timestamps - pulses.timestamps(ii,2)); %armChoice.timestamps has the timestamps of the reward location
        if min(tempdiffON) < 1 % pulse is aligned, stim is return arm specific
            armPulses.timestamps(tempdiffON<1) = pulses.timestamps(ii,1);
            armPulses.stim(tempdiffON<1) =  1;
            armPulses.target =  2; % Return
            if sum(pulses.timestamps(:,1) == pulses.timestamps(ii,1)) > 1
                armPulses.region =  3;
            else
                armPulses.region =  pulses.eventID(ii);
            end

        elseif min(tempdiffOFF) < 1 % pulse is aligned to off, stim is stem specific
            armPulses.timestamps(find(tempdiffOFF<1)-1) = pulses.timestamps(ii,1);
            armPulses.stim(find(tempdiffOFF<1)-1) =  1;
            armPulses.target =  1; % Stem
            if sum(pulses.timestamps(:,2) == pulses.timestamps(ii,2)) > 1
                armPulses.region =  3;
            else
                armPulses.region =  pulses.eventID(ii);
            end
        end
    end

    numBlock = 1;
    for jj = (find(armPulses.stim==1,1)-10):10:length(armChoice.choice)
        if jj+10<= length(armChoice.choice) && jj >0
            armPulses.performanceBlocks(1,numBlock) = mean(armPulses.stim(jj:jj+9));
            if jj == 1
                armPulses.performanceBlocks(2,numBlock) = nansum(armChoice.choice(jj:jj+9))/9;
            else
                armPulses.performanceBlocks(2,numBlock) = nansum(armChoice.choice(jj:jj+9))/10;
            end
        elseif jj+8 <= length(armChoice.choice)
            armPulses.performanceBlocks(1,numBlock) = mean(armPulses.stim(jj:length(armChoice.choice)));
            armPulses.performanceBlocks(2,numBlock) = nansum(armChoice.choice(jj:length(armChoice.choice)))/(length(armChoice.choice)-jj);
        end
        numBlock = numBlock+1;
    end
    armPulses.performance(1,1:2) = [0 1];
    armPulses.performance(2,1) = nanmean(armPulses.performanceBlocks(2,(armPulses.performanceBlocks(1,:)<0.2)));
    armPulses.performance(2,2) = nanmean(armPulses.performanceBlocks(2,(armPulses.performanceBlocks(1,:)>=0.8)));
    
    C = strsplit(pwd,'\');
    save([C{end} '.ArmPulses.Events.mat'], 'armPulses');
end

end
