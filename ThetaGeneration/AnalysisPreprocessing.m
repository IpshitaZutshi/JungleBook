

%% KilosortBatch
% MV-BuzsakiLab 2019
% Edited by IZ-BuzsakiLab 2019

function AnalysisPreprocessing(varargin)

% AnalysisPreprocessing(varargin)

% Pipeline processing for dat + xml data:
%   1. Organizes folders data by sessions ( One session being all recordings at a certain depth from the same day).
%   2. Concatenate sessions data.
%   3. Spike sort sessions using kilosort.
%   4. Autocluster.
%   5. Makes a folder summary with LFP analysis (PSTH CSD, Theta amplitude, theta asymmetry, theta frequency) and 
%   Spike analysis (spike-waveforms, autocorrelogram and spatial position;
%   psth from analog-in inputs) & extract behavior (alternation)

% INPUT
% No inputs, launch gui
%   <options>       optional list of property-value pairs (see table below)
%   expPath        - Basepath for experiment. It contains all session
%                       folders. If not provides, promt.
%   analogEv       - List of analog channels with pulses to be detected (it supports Intan Buzsaki Edition) - usually 64 if input is going to ADC1.
%   numAnalog      - Number of analog channels in the xml file. 
%                       Default 1.
%   forceSum       - Force make folder summary (overwrite, if necessary).
%                       Default false.
%   forcesort      - Force to kilosort (overwrite, if necessary).
%                       Default false.
%   cleanArtifacts - Remove artifacts from dat file, by default, is there is
%                       analogEv, is true
%   analysisList   - Logical array to indicate summary analysis, according
%                       to the following list: 
%                    1. LFP power profile
%                    2. LFP peri-stimulus wavelet power
%                    3. LFP theta CSD
%                    4. LFP theta amplitude, frequency and asymmetry index
%                    5. Spike-waveform, autocorrelogram and spatial position 
%                    6. Psth from analog-in inputs
%                   Example: [1 1 1 1 1 1 1] or 'all' make all. Default 'all'    
%   pullData       - Path for raw data. Look for not analized session to
%                       copy to the main folder basepath. To do...
%
%   MV-BuzsakiLab 2019
%   Edited by IZ-BuzsakiLab 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Parms
p = inputParser;
addParameter(p,'expPath',[],@isdir);
addParameter(p,'analogEv',[],@isnumeric);
addParameter(p,'numAnalog',1,@isnumeric);
addParameter(p,'analysisPath',[]); % Local paht to run the anaysis, is em
addParameter(p,'forceSum',false,@islogical);
addParameter(p,'forcesort',false,@islogical);
addParameter(p,'medianSubstr',true,@islogical);
addParameter(p,'analysisList','all');
addParameter(p,'cleanArtifacts',false,@islogical);
% addParameter(p,'pullData',[],@isdir); To do... 
parse(p,varargin{:});

expPath = p.Results.expPath;
analogEv = p.Results.analogEv;
analysisPath = p.Results.analysisPath;
numAnalog = p.Results.numAnalog;
forceSum = p.Results.forceSum;
medianSubstr = p.Results.medianSubstr;
forcesort = p.Results.forcesort;
analysisList = p.Results.analysisList;
cleanArtifacts = p.Results.cleanArtifacts;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end
allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});

%waitMin(90)
%pause(90*60)
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
        copyfile(strcat(xmlFile(1).folder,filesep,xmlFile(1).name),strcat(allpath{1},filesep,'global.xml'));
    end
    cd(allpath{1});
end

%% deals with analysis list
if ischar(analysisList)
    if strcmpi(analysisList,'all')
        %analysisList = [1 1 1 1 1 0 1 1 1];
        analysisList = [1 1 1 0 0 0 0 0 1];
        if isempty(analogEv)
            analysisList([2 3 6]) = 0;
        end
    else 
        error('Analysis list format not recognized!');
    end
end

%% Build sessions
disp('Building session folders (It assumes session as all folder at a certain depth recordered on the same day)...');
allFolder = dir(pwd);

for ii = 1:length(allFolder)
    % Find max sessNumber already
    a = dir('*sess*');
    startSess = 1;
    for jj = 1:size(a)
        num = regexp(a(jj).name, 'sess(\d+)', 'tokens');
        num = str2double(num{1});
        if num>startSess
            startSess = num;
        end
    end
    if strlength(allFolder(ii).name) > 12 && isfolder(allFolder(ii).name) % if looks like a data folder
        folderFiels = strsplit(allFolder(ii).name,'_');
        if numel(folderFiels) == 4 % if the depth is included in the folder name
            if isempty(dir(strcat('*',folderFiels{2},'_',folderFiels{3},'_','sess*'))) % if there is no session folder yet
                mkdir(strcat(folderFiels{1},'_',folderFiels{2},'_',folderFiels{3},'_','sess',num2str(size(dir('*sess*'),1) + 1))); % create folder
            end
            if ~contains(folderFiels{4},'sess') % if it is not a session folder
                targetfoder = dir(strcat(folderFiels{1},'_',folderFiels{2},'_',folderFiels{3},'_','sess*'));
                movefile(strcat('*',folderFiels{2},'_',folderFiels{3},'_',folderFiels{4}),targetfoder.name); % move to session folder
            end
        elseif numel(folderFiels) == 3
            if isempty(dir(strcat('*',folderFiels{2},'_','sess*'))) % if there is no session folder yet
                mkdir(strcat(folderFiels{1},'_',folderFiels{2},'_','sess',num2str(startSess + 1))); % create folder
            end
            if ~contains(folderFiels{3},'sess') % if it is not a session folder
                targetfoder = dir(strcat(folderFiels{1},'_',folderFiels{2},'_','sess*'));
                movefile(strcat('*',folderFiels{2},'_',folderFiels{3}),targetfoder.name); % move to session folder
            end            
        end 
    end
end

% If recordered with intan Buz Edit, change dat file name
allpath = strsplit(genpath(expPath),';'); % all folders again
for ii = 1:size(allpath,2)
    if strlength(allpath{ii}) > 12 && isfolder(allpath{ii}) % if looks like a data folder
        cd(allpath{ii});
        if ~isempty(dir('amplifier_analogin_auxiliary_int16.dat')) % if recordered with intan Buz Edit
            movefile(strcat(allpath{ii},'\','amplifier_analogin_auxiliary_int16.dat'),strcat(allpath{ii},'\','amplifier.dat'));
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
    if ~exist(strcat(allSess(ii).name,'.xml'))
        delete(strcat(allSess(ii).name,'.xml'));% bring xml file
        copyfile(strcat(allpath{1},'\global.xml'),strcat(allSess(ii).name,'.xml'),'f');
    end
    % Concatenate sessions
    
    bz_ConcatenateDats(pwd,0,1);

end

%% Loading metadata
try
    session = sessionTemplate(pwd,'showGUI',false); % 
    session.channels = 1:session.extracellular.nChannels;    
    %session.channelTags.Bad.channels = [24:38 48:63];
    save([session.general.name,'.session.mat'],'session','-v7.3');
catch
    warning('it seems that CellExplorer is not on your path');
end

%% Get analog and digital pulses
for ii = 1:size(allSess,1)
    cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
    if  ~isempty(dir('analogin.dat'))
        try
            [pulses] = bz_getAnalogPulsesSine;%getAnalogPulses('analogChannelsList',analogChannelsList);
        catch
            warning('No analog pulses detected');
        end
    end
    if ~isempty(dir('*digitalin.dat'))
        digitalIn = getDigitalIn('all','fs',session.extracellular.sr); 
    end
    %% Make LFP
%     if isempty(dir('*.lfp'))
%         disp('Creating .lfp file. This could take a while...');
%         ResampleBinary(strcat(allSess(ii).name,'.dat'),strcat(allSess(ii).name,'.lfp'),...
%             session.extracellular.nChannels,1, session.extracellular.sr/session.extracellular.srLfp);
%     end
    
%     %% MEDIAN SUBS
%     if isempty(dir([session.general.name '_original.dat']))
%         if islogical(medianSubstr) && medianSubstr
%             medianSubtraction(pwd,'keepDat',true);
%         elseif medianSubstr
%             medianSubtraction(pwd,'ch',medianSubstr,'keepDat',true);
%         end
%     else
%         warning('Session was already median-subtracted. Spiking...');
%     end

end

%% LFP
[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true); sessionInfo.rates.lfp = 1250;  save(strcat(sessionInfo.session.name,'.sessionInfo.mat'),'sessionInfo');
if isempty(dir('*.lfp'))
    try 
        bz_LFPfromDat(pwd,'outFs',1250); % generating lfp
    catch
        disp('Problems with bz_LFPfromDat, resampling...');
        ResampleBinary(strcat(sessionInfo.session.name,'.dat'),strcat(sessionInfo.session.name,'.lfp'),...
            sessionInfo.nChannels,1,sessionInfo.rates.wideband/sessionInfo.rates.lfp);
    end
end

%% Kilosort all sessions
disp('Spike sort all sessions...');
for ii = size(allSess,1):-1:1
    cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
    if forcesort ||  isempty(dir('*Kilosort*')) % if not kilosorted yet
        fprintf(' ** Kilosorting session %3.i of %3.i... \n',ii, size(allSess,1));   
        KiloSortWrapper;
        kilosortFolder = dir('*Kilosort_*');
        try
            PhyAutoClustering(strcat(kilosortFolder.folder,'\',kilosortFolder.name)); % autoclustering
        catch err
            disp(err.message)
            warning('PhyAutoClustering not possible!!');
        end
        if exist('phyLink') && ~isempty(phyLink) % move phy link to
            kilosort_path = dir('*Kilosort*');
            try copyfile(phyLink, strcat(kilosort_path.name,filesep,'LaunchPhy')); % copy pulTime to kilosort folder
            end
        end
    end

end

% %% BatchAnalysis
% for ii = 1:size(allSess,1)
%     cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
%     getToneTracking
% %     fprintf(' ** Summary %3.i of %3.i... \n',ii, size(allSess,1));
% %     cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
% %     if forceSum || (~isempty(dir('*Kilosort*')) && isempty(dir('summ'))) % is kilosorted but no summ
% %       %  AnalysisBatchTheta;         
% %     end
% end

end