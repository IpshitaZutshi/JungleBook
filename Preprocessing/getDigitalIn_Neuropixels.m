
function [digitalIn] = getDigitalIn_Neuropixels(varargin)
% [pul, val, dur] = getPulses(d,varargin)
%
% Find digital In pulses
%
% INPUTS
% ch            Default all.
% <OPTIONALS>
% fs            Sampling frequency (in Hz), default 30000, or try to
%               recover for rhd 
% offset        Offset subtracted (in seconds), default 0.
% periodLag     How long a pulse has to be far from other pulses to be consider a different stimulation period (in seconds, default 5s)    
% filename      File to get pulses from. Default, digitalin.dat file with folder
%               name in current directory
%
%
% OUTPUTS
%               digitalIn - events struct with the following fields
% ints          C x 2  matrix with pulse times in seconds. First column of C 
%               are the beggining of the pulses, second column of C are the end of 
%               the pulses.
% dur           Duration of the pulses. Note that default fs is 30000.
% timestampsOn  Beggining of all ON pulses
% timestampsOff Beggining of all OFF pulses
% intsPeriods   Stimulation periods, as defined by perioLag


p = inputParser;
addParameter(p,'fs',30000,@isnumeric)
addParameter(p,'offset',0,@isnumeric)
addParameter(p,'filename',[],@isstring)
addParameter(p,'periodLag',0.4,@isnumeric)
addParameter(p,'force',false,@islogical)

parse(p, varargin{:});
fs = p.Results.fs;
offset = p.Results.offset;
filename = p.Results.filename;
lag = p.Results.periodLag;
force = p.Results.force;

if ~isempty(dir('*DigitalIn.events.mat')) && ~force
    disp('Pulses already detected! Loading file.');
    file = dir('*DigitalIn.events.mat');
    load(file.name);
    return
end

if isempty(filename)
    filename=dir('digitalIn.dat');
    filename = filename.name;
else
    disp('No digitalIn file found...');
end

disp('Loading digital channels...');
m = memmapfile(filename,'Format','uint16','writable',false);
data = double(m.Data);
clear m

%% Find pulses
data(data~=1) = 0;
test2 = diff(data);

pulses = find(test2 == 1);
pulses2 = find(test2 == -1);
digital_on = pulses;
digital_off = pulses2;


% take timestamp in seconds
digitalIn.timestampsOn = digital_on/fs;
digitalIn.timestampsOff = digital_off/fs;
        
% intervals
d = zeros(2,max([size(digitalIn.timestampsOn,1) size(digitalIn.timestampsOff,1)]));
d(1,1:size(digitalIn.timestampsOn,1)) = digitalIn.timestampsOn;
d(2,1:size(digitalIn.timestampsOff,1)) = digitalIn.timestampsOff;
if d(1,1) > d(2,1)
    d = flip(d,1);
end
if d(2,end) == 0; d(2,end) = nan; end
digitalIn.ints = d;
digitalIn.dur = digitalIn.ints(2,:) - digitalIn.ints(1,:); % duration
        
clear intsPeriods
intsPeriods(1,1) = d(1,1); % find stimulation intervals
intPeaks =find(diff(d(1,:))>lag);
for jj = 1:length(intPeaks)
    intsPeriods(jj,2) = d(2,intPeaks(jj));
    intsPeriods(jj+1,1) = d(1,intPeaks(jj)+1);
end
intsPeriods(end,2) = d(2,end);  
digitalIn.intsPeriods = intsPeriods;


if exist('digitalIn')==1
    xt = linspace(0,size(data,2)/fs,size(data,2));
    data = flip(data);
    data = data(1:size(digitalIn.intsPeriods,2),:);

    h=figure;
    imagesc(xt,1:size(data,2),data);
    xlabel('s'); ylabel('Channels'); colormap gray 
    mkdir('Pulses');
    saveas(h,'pulses\digitalIn.png');

    try save([sess.FileName '.DigitalIn.events.mat'],'digitalIn');
    catch
        save('digitalIn.events.mat','digitalIn');
    end
else
    digitalIn = [];
end
 
end