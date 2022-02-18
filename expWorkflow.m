function expWorkflow(varargin)

% Pipeline processing for dat + xml data:
%   1. Organizes folders data by sessions (session being all recording from same day).
%   2. Concatenate sessions data.
%   3. Spike sort sessions by kilosort.
%   4. Autocluster.
%   5. Makes a folder summary with spike-waveforms, autocorrelogram and spatial position; psth from analog-in inputs and
% psth around slow waves and ripples.

% INPUT
% No inputs, launch gui
%   <options>       optional list of property-value pairs (see table below)
%   expPath        - Basepath for experiment. It contains all session
%                       folders. If not provides, promt.
%   analogEv       - List of analog channels with pulses to be detected (it support Intan Buzsaki Edition).
%   forceSum       - Force make folder summary (overwrite, if necessary).
%                       Default false.
%   cleanArtifacts - Remove artifacts from dat file, by default, is there is
%                       analogEv, is true
%   analisysList   - Logical array to indicate summary analysis, according
%                       to the following list: 
%                           1. Spike-waveform, autocorrelogram and spatial position 
%                           2. Psth from analog-in inputs
%                   Example: [1 1] or 'all' make all. Default 'all'    
%   pullData       - Path for raw data. Look for not analized session to
%                       copy to the main folder basepath. To do...
%
%   MV-BuzsakiLab 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Parms
p = inputParser;
addParameter(p,'expPath',[],@isdir);
addParameter(p,'analogEv',[],@isnumeric);
addParameter(p,'forceSum',false,@islogical);
addParameter(p,'analisysList','all');
addParameter(p,'cleanArtifacts',false,@islogical);
% addParameter(p,'pullData',[],@isdir); To do... 
parse(p,varargin{:});

expPath = p.Results.expPath;
analogEv = p.Results.analogEv;
forceSum = p.Results.forceSum;
analisysList = p.Results.analisysList;
cleanArtifacts = p.Results.cleanArtifacts;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end
allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});

%% deals with xml.
disp('Check xml...');
if isempty(dir('global.xml')) 
    disp('No xml global file! Looking for it...');
    xmlFile = []; ii = 2;
    while isempty(xmlFile) && ii < size(allpath,2)
        disp(ii);
        cd(allpath{ii});
        xmlFile = dir('*.xml');
        ii = ii + 1;
    end
    if isempty(xmlFile)    
        [file, path] = uigetfile('*.xml','Select global xml file');
        copyfile(strcat(path,file),'global.xml');
    else
        copyfile(strcat(xmlFile.folder,filesep,xmlFile.name),strcat(allpath{1},filesep,'global.xml'));
    end
    cd(allpath{1});
end

%% deals with analysis list
if ischar(analisysList)
    if strcmpi(analisysList,'all')
        analisysList = [1 1 1 1];
        analisysList(2) = ~isempty(analogEv);
    else
        error('Analysis list format not recognized!');
    end
end

%% Build sessions
disp('Building session folders (It asumes session as all folder recordered same day)...');
allFolder = dir(pwd);
for ii = 1:length(allFolder)
    if strlength(allFolder(ii).name) > 12 && isfolder(allFolder(ii).name) % if looks like a data folder
        folderFiels = strsplit(allFolder(ii).name,'_');
        if isempty(dir(strcat('*',folderFiels{2},'_','sess*'))) % if there is no session folder yet
            mkdir(strcat(folderFiels{1},'_',folderFiels{2},'_','sess',num2str(size(dir('*sess*'),1) + 1))); % create folder
        end
        if ~contains(folderFiels{3},'sess') % if it is not a session folder
            targetfoder = dir(strcat(folderFiels{1},'_',folderFiels{2},'_','sess*'));
            movefile(strcat('*',folderFiels{2},'_',folderFiels{3}),targetfoder.name); % move to session folder
        end
    end
end

% If recordered with intan Buz Edit, change dat file name
allpath = strsplit(genpath(expPath),';'); % all folders again
for ii = 1:size(allpath,2)
    if strlength(allpath{ii}) > 12 && isfolder(allpath{ii}) % if looks like a data folder
        cd(allpath{ii});
        if ~isempty(dir('amplifier_analogin_auxiliary_int16.dat')) % if recordered with intan Buz Edit
            movefile('amplifier_analogin_auxiliary_int16.dat','amplifier.dat');
            try  movefile('amplifier_analogin_auxiliary_int16.xml','amplifier.xml');
                movefile('amplifier_analogin_auxiliary_int16.nrs','amplifier.nrs');
            end
        end
    end
end

%% Concatenate sessions
cd(allpath{1});
allSess = dir('*_sess*');
disp('Concatenate session folders...');
for ii = 1:size(allSess,1)
    fprintf(' ** Concatenating session %3.i of %3.i... \n',ii, size(allSess,1));
    cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
    delete(strcat(allSess(ii).name,'.xml'));% bring xml file
    copyfile(strcat(allpath{1},'\global.xml'),strcat(allSess(ii).name,'.xml'),'f');
    % [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
    bz_ConcatenateDats;
end

%% Kilosort all sessions
disp('Spike sort all sessions...');
for ii = 1:size(allSess,1)
    cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
    if  isempty(dir('*Kilosort*')) % if not kilosorted yet
    fprintf(' ** Kilosorting session %3.i of %3.i... \n',ii, size(allSess,1));
        KiloSortWrapper;
        kilosortFolder = dir('*Kilosort*');
        try PhyAutoClustering(strcat(kilosortFolder.folder,'\',kilosortFolder.name)); % autoclustering
        catch
            warning('PhyAutoClustering not possible!!');
        end
        if exist('phyLink') && ~isempty(phyLink) % move phy link to
            kilosort_path = dir('*Kilosort*');
            try copyfile(phyLink, strcat(kilosort_path.name,filesep,'LaunchPhy')); % copy pulTime to kilosort folder
            end
        end
    end
end

end