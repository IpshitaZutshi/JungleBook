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
sessionInfo = bz_getSessionInfo; 
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

if ~isempty(indices)
    for rr = 1:length(indices)
       datpaths.amplifier{rr} = filePaths{indices(rr)};
       t = dir(filePaths{indices(rr)});
       datsizes.amplifier(rr) = t.bytes;
       nameList = strsplit(filePaths{indices(rr)},'\');
       recordingnames{rr} = nameList{8}; %% Consistently number 8
       [folder, ~, ~] = fileparts(filePaths{indices(rr)});
       timestampsDir{rr} = strcat(folder,'\','timestamps.npy');
    end
else           
    disp(strcat('No .dats found in subfolders.  Skipping ',basename))
    return
end

%% Sort folders correctly
[dateTime,timeOrder] = sort(datetime(recordingnames,'InputFormat','yyyy-MM-dd_HH-mm-ss_SSS'));

%% Concatenate
nFiles = length(datpaths.amplifier);
startTimes = cell(nFiles,1);
fileDurations = zeros(nFiles,1); %seconds

%% Read start times from the timestamps.npy file
for i = timeOrder
    timeData = readNPY(timestampsDir{i});
    firsttimepoint(i) = timeData(1);
    lasttimepoint(i) = timeData(end);
    firstlasttimepoints(i,:) = [firsttimepoint(i) lasttimepoint(i)];
    if i==1
        transitiontimes_samp = firstlasttimepoints(i,:);
    else
        transitiontimes_samp(i,:) = firstlasttimepoints(i,:)+transitiontimes_samp(i-1,2)+1;
    end
    trueTime(i) = dateTime(i);
end

fileDurations = (firstlasttimepoints(timeOrder,2)-firstlasttimepoints(timeOrder,1));

%% write to new file
sr = sessionInfo.rates.wideband;
nChan = sessionInfo.nChannels;
chunkDuration = 60; %seconds

% concatenated neural signal file
newPath = [basepath, filesep, basename, '.dat'];

% concatenated input channel file
inpPath = [basepath, filesep, 'digitalin.dat'];

% concatenated dat of neural signals
fNew = fopen(newPath,'w');

% concatenated dat of input signal
fInp = fopen(inpPath,'w');

for fIdx = timeOrder
    
    tmpOldFile = datpaths.amplifier{fIdx};
    disp(['writing ',recordingnames{fIdx}])
    
    % old dat in subdir
    fOld = fopen(tmpOldFile,'r');
    
    % new dat of input signal in subdir
    [folder, ~, ~] = fileparts(tmpOldFile);
    fInSub = fopen([folder,filesep,'digitalin.dat'],'w');
    
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

oldSize = sum(datsizes.amplifier);
newSize = dir(newPath).bytes;
inpSize = dir(inpPath).bytes;

if newSize+inpSize==oldSize
    disp('files concatenated and size checked')
    sizeCheck = true;
else
    error('New file not right size')
end

%% Make the events.mat file that saves all the merge information

eventsfilename = fullfile(basepath,[basename,'.MergePoints.events.mat']);

MergePoints.timestamps = transitiontimes_samp;
MergePoints.timestamps_samples = transitiontimes_samp*sr;
MergePoints.firstlasttimpoints_samples = firstlasttimepoints;
MergePoints.foldernames = recordingnames;
MergePoints.filesmerged = datpaths;
MergePoints.filesizes = datsizes;
MergePoints.sizecheck = sizeCheck;
MergePoints.detectorinfo.detectorname = 'concatenateNeuropixelsDats';
MergePoints.detectorinfo.detectiondate = datestr(now,'yyyy-mm-dd');

%Saving 
save(eventsfilename,'MergePoints');

%% Ammend session.mat file to reflect input channel cut out
sessionInfo = removeChannelSess(sessionInfo,(nChan-1));%0 indexing in neuroscope
save([basename '.sessionInfo.mat'],'sessionInfo');

end