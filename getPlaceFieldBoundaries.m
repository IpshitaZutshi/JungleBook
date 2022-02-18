% Compile data across all sessions

function getPlaceFieldBoundaries(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'downsample',true,@islogical);
parse(p,varargin{:});

expPath = p.Results.expPath;
downsample = p.Results.downsample;

%%
if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end
allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');

%% Start collecting data
for ii = 1:size(allSess,1)

    cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
    [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);    
    file = dir(('*.SessionPulses.Events.mat'));
    load(file.name);
    
    efields = fieldnames(sessionPulses);    

    for jj = 1:length(efields)
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name,'\',efields{jj}));
        if ~downsample
            load(strcat(efields{jj},'\',allSess(ii).name,'.firingMapsAvg.cellinfo.mat'))
        else
            load(strcat(efields{jj},'\',allSess(ii).name,'_DS.firingMapsAvg.cellinfo.mat'))
        end
        bz_findPlaceFields1D_IZ('firingMaps',firingMaps,'downsample',downsample);
    end

end
