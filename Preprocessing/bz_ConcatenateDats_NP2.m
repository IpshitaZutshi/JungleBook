function concatenateNeuropixelsDats(varargin)
% bz_ConcatenateDats - Concatenate raw .dat files found in a session folder
% - for OpenEphys neuropixel2.0 recordings 
% 
% Ipshita Zutshi 2024

%% Handling inputs
% basic session name and and path
p = inputParser;
addParameter(p,'basepath',pwd,@isstr)

parse(p,varargin{:})
basepath = p.Results.basepath;

%% Get session info
sessionInfo = bz_getSessionInfo('basepath',basepath); 
basename = sessionInfo.session.name;

%% If the dats are already merged quit
if exist(fullfile(basepath,[basename,'.dat']),'file')
    disp('.dat already exists in session directory, not merging subdats')
    return
end

%% find all files and sort by time created
[fileNames, filePaths] = getAllFiles(pwd);

% Find all the files with the name continuous.dat
indices = find(strcmp(fileNames, 'continuous.dat'));

binFiles = checkFile('basepath',basepath,'fileType','.bin');
nFiles = length(binFiles);
startTimes = cell(nFiles,1);
fileDurations = zeros(nFiles,1); %seconds

for i = 1:nFiles
    tmpMeta = ReadMeta(binFiles(i).name,binFiles(i).folder);
    startTimes{i} = tmpMeta.fileCreateTime;
    fileDurations(i) = str2num(tmpMeta.fileTimeSecs);
end

[~,timeOrder] = sort(datetime(startTimes,'InputFormat','yyyy-MM-dd''T''HH:mm:ss'));
fileDurations = fileDurations(timeOrder);

%% write to new file
sr = session.extracellular.sr;
nChan = session.extracellular.nChannels;
chunkDuration = 60; %seconds

% concatenated neural signal file
newPath = [basepath, filesep, basename, '.dat'];

% concatenated input channel file
inpPath = [basepath, filesep, 'digitalin.dat'];

% concatenated dat of neural signals
fNew = fopen(newPath,'w');

% concatenated dat of input signal
fInp = fopen(inpPath,'w');

for fIdx = 1:nFiles
    
    tmpOldFile = binFiles(timeOrder(fIdx));
    disp(['writing ',tmpOldFile.name])
    
    % old dat in subdir
    fOld = fopen([tmpOldFile.folder,filesep,tmpOldFile.name],'r');
    
    % new dat of input signal in subdir
    fInSub = fopen([tmpOldFile.folder,filesep,'digitalin.dat'],'w');
    
    while true
        data = fread(fOld,[nChan sr*chunkDuration],'int16');
        
        if isempty(data)
            break;
        end
        
        fwrite(fNew,data(1:end-1,:),'int16');
        fwrite(fInp,data(end,:),'int16');
        fwrite(fInSub,data(end,:),'int16');
        
    end
    
    fclose(fOld);
    fclose(fInSub);
    
end

fclose(fNew);
fclose(fInp);

oldSize = sum([binFiles.bytes]);
newSize = dir(newPath).bytes;
inpSize = dir(inpPath).bytes;

if newSize+inpSize==oldSize
    disp('files concatenated and size checked')
    sizeCheck = true;
else
    error('New file not right size')
end

%% Compute merge points

%using fileDurations in seconds
transitionTimes_sec = zeros(nFiles,2);
transitionTimes_sec(:,2) = cumsum(fileDurations);
if nFiles>1
    transitionTimes_sec(2:end,1) = transitionTimes_sec(1:end-1,2)+1/sr;
end

transitionTimes_samp = transitionTimes_sec*sr;
firstLastSample = zeros(nFiles,2);
firstLastSample(:,2) = fileDurations*sr;

MergePoints.timestamps = transitionTimes_sec;
MergePoints.timestamps_samples = transitionTimes_samp;
MergePoints.firstlasttimepoints_samples = firstLastSample;
MergePoints.foldernames = {binFiles(timeOrder).folder};
MergePoints.sizecheck = sizeCheck;

eventsfilename = fullfile(basepath,[basename,'.MergePoints.events.mat']);
save(eventsfilename,'MergePoints');

%% Ammend session.mat file to reflect input channel cut out
session = removeChannelSess(session,nChan);
save([basename '.session.mat'],'session');

end