function [tracking] = toneTracking(varargin)

p = inputParser;
addParameter(p,'analogInputPos',2,@isnumeric);
addParameter(p,'analogInputTone',3,@isnumeric);
addParameter(p,'fs',150,@isnumeric);
addParameter(p,'trackLength',112,@isnumeric);
addParameter(p,'trackImgLength',410,@isnumeric);
addParameter(p,'freqRange',[1000 22000],@isnumeric);
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'forceReload',false,@islogical)

parse(p,varargin{:});
analogInputPos = p.Results.analogInputPos;
analogInputTone = p.Results.analogInputTone;
fs = p.Results.fs;
trackLength = p.Results.trackLength;
trackImgLength = p.Results.trackImgLength;
freqRange = p.Results.freqRange;
saveMat = p.Results.saveMat;
forceReload = p.Results.forceReload;

basepath = pwd;

%% Deal with inputs
if ~isempty(dir([basepath filesep '*Tracking.Behavior.mat'])) && ~forceReload
    disp('Trajectory already detected! Loading file.');
    file = dir([basepath filesep '*Tracking.Behavior.mat']);
    load(file.name);
    return
end

%% Get analog inputs
try [amplifier_channels, notes, aux_input_channels, spike_triggers,...
    board_dig_in_channels, supply_voltage_channels, frequency_parameters, board_adc_channels ] =...
    read_Intan_RHD2000_file_2;
catch
    disp('File ''info.rhd'' not found. (Type ''help <a href="matlab:help loadAnalog">loadAnalog</a>'' for details) ');
end

disp('Loading analogin.dat...');
filename = 'analogin.dat';  
num_channels = length(board_adc_channels); % ADC input info from header file
fileinfo = dir(filename);
num_samples = fileinfo.bytes/(num_channels * 2); % uint16 = 2 bytes
fid = fopen('analogin.dat', 'r');
v = fread(fid, [num_channels, num_samples], 'uint16');
fclose(fid);
posInfo = v(1:num_channels,:);
%posInfo = ((posInfo)./65535) * 3; % convert to volts
fsample = frequency_parameters.amplifier_sample_rate;
% Downsample to 1250
for ii = 1:num_channels
    posInfoDS(ii,:) = downsample(posInfo(ii,:),fsample/fs);
end

% Correct for position jumps greater than ~30 cm.
x = posInfoDS(analogInputPos,:);
tone = posInfoDS(analogInputTone,:);

[idxjumps] = find((diff(x)/3)>0.2);
x(idxjumps) = NaN;
%Rescale to cm
%x = (x/3)*trackLength;
xt = linspace(0,length(x)/fs,length(x));

%Do not rescale to frequency now
tone(idxjumps) = NaN;
tone(tone<0.08) = NaN;

[~,fbasename,~]=fileparts(pwd);

tracking.position.x = x;
tracking.position.tone = tone;
tracking.timestamps = xt;
tracking.description = 'Tone';
tracking.folder = fbasename;
tracking.samplingRate = fs;

if saveMat
    C = strsplit(pwd,'\');
    save([basepath filesep C{end} '.Tracking.Behavior.mat'],'tracking');
end

end