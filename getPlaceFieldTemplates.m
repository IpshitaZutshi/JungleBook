function getPlaceFieldTemplates(varargin)


%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
parse(p,varargin{:});
basepath = p.Results.basepath;

%% Deal with inputs
if ~isempty(dir([basepath filesep '*.SessionPulses.Events.mat']))
    disp('Session pulses already detected! Loading file.');
    file = dir([basepath filesep '*.SessionPulses.Events.mat']);
    load(file.name);
end

%% Find subfolder recordings
flds = fieldnames(sessionPulses);

for ii = 1:size(flds,1)
    cd([basepath filesep flds{ii}]);
    fprintf('Computing place field templates in %s folder \n',flds{ii});
    bz_findPlaceFieldsTemplate_IZ
    close all
end

end