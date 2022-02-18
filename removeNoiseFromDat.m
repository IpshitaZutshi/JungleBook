function removeNoiseFromDat(basepath,varargin)
% Edit dat files with different options.
%
% USAGE
%   editDatFile(basepath,ints,varargin)
% 
% INPUT
% basepath      If not provided, takes pwd
% threshold     Intervals to be edited.
% method        'substractMedian'
%
% <optional>
% option        'remove' or 'zeroes' (default). 
%
% Manu Valero-BuzsakiLab 2020
% TO DO: Edit AnalogIn data!!! 
%
%% Defaults and Parms
p = inputParser;
addParameter(p,'ch','all');
addParameter(p,'method','substractMedian',@ischar);

parse(p,varargin{:});
ch = p.Results.ch;
method = p.Results.method;

% Get elements
prevPath = pwd;
cd(basepath);

xml = LoadParameters;
fileTargetAmplifier = dir('amplifier*.dat');
if isempty(fileTargetAmplifier)
    filename = split(pwd,filesep); filename = filename{end};
    fileTargetAmplifier = dir([filename '*.dat']);
end

if size(fileTargetAmplifier,1) == 0
    error('Dat file not found!!');
end

m = memmapfile(fullfile(basepath,fileTargetAmplifier.name),'Format','int16','Writable', true);
data=reshape(m.Data,xml.nChannels,[]);
if ischar(ch) && strcmpi(ch, 'all')
    ch = 1:size(data,1);
end
timestamps = linspace(0,size(data,2)/xml.rates.wideband,size(data,2));
if strcmpi('substractMedian',method)
    m_data = median(data(ch,:));
elseif strcmpi('substractMean',method)
    m_data = mean(data(ch,:));
end

for ii = 1:length(ch)
    data(ch(ii),:) = int16(double(data(ch(ii),:)) - double(m_data));
end
m.Data=data(:);

cd(prevPath);

end