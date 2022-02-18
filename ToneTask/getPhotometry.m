function photometry = getPhotometry(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'forceReload',false,@islogical)
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'ds',1250,@isnumeric)

parse(p,varargin{:});
basepath = p.Results.basepath;
forceReload = p.Results.forceReload;
saveMat = p.Results.saveMat;
ds = p.Results.ds;

%% In case photometry already exists 
if ~isempty(dir([basepath filesep '*Photometry.Behavior.mat'])) && ~forceReload
    disp('Photometry already detected! Loading file.');
    file = dir([basepath filesep '*Photometry.Behavior.mat']);
    load(file.name);
    return
end

try [amplifier_channels, notes, aux_input_channels, spike_triggers,...
    board_dig_in_channels, supply_voltage_channels, frequency_parameters, board_adc_channels ] =...
    read_Intan_RHD2000_file_2;
catch
    disp('File ''info.rhd'' not found. (Type ''help <a href="matlab:help loadAnalog">loadAnalog</a>'' for details) ');
end

fsample = frequency_parameters.amplifier_sample_rate;
disp('Loading analogin.dat...');
filename = 'analogin.dat';  
num_channels = length(board_adc_channels); % ADC input info from header file
fileinfo = dir(filename);
num_samples = fileinfo.bytes/(num_channels * 2); % uint16 = 2 bytes
fid = fopen('analogin.dat', 'r');
v = fread(fid, [num_channels, num_samples], 'uint16');
fclose(fid);
v=v*0.000050354; %convert to volt

Fs = fsample;
Wn = 20/(Fs/2);
[b,a] = butter(3,Wn);

F_green=v(4,:);
F_lp = filtfilt(b,a,F_green);
F_lp = downsample(F_lp,Fs/ds);

t = 1:1:length(F_lp);
b = polyfit(1:1:length(F_lp),F_lp,5);
a = polyval(b,t);
dfByF = (F_lp-a)./a;

photometry.signal = F_lp;
photometry.corrSignal = dfByF;
photometry.timestamps = linspace(1,length(F_lp)/ds,length(F_lp));

if saveMat
    C = strsplit(pwd,'\');
    save([basepath filesep C{end} '.Photometry.Behavior.mat'],'photometry');
end

end