% Define source and destination folders
sourceFolder = 'Z:\Homes\zutshi01\Recordings\Cue task\T20\Final';
destFolder = 'C:\Data\DeepLabCut\Videos';

% Create the destination folder if it doesn't exist
if ~exist(destFolder, 'dir')
    mkdir(destFolder);
end

% List all folders in the source that start with 'T22'
folderList = dir(fullfile(sourceFolder, 'T20*'));
% Filter to keep only directories (ignore files if any)
folderList = folderList([folderList.isdir]);

% Loop through each folder
for i = 1:length(folderList)
    % Construct the full path of the current folder
    currentFolder = fullfile(sourceFolder, folderList(i).name);
    
    % Look for .avi files starting with 'Top' in the current folder
    videoFiles = dir(fullfile(currentFolder, 'Top*.avi'));
    
    if ~isempty(videoFiles)
        % Assuming one matching file per folder; use the first one if multiple found
        srcFile = fullfile(currentFolder, videoFiles(1).name);
        
        % Create the new filename using the folder name (e.g., T22_T22_241209_173821.avi)
        newFileName = ['T20_' folderList(i).name '.avi'];
        destFile = fullfile(destFolder, newFileName);
        
        % Copy the file to the destination folder with the new name
        copyfile(srcFile, destFile);
        
        fprintf('Copied %s to %s\n', srcFile, destFile);
    else
        fprintf('No video file starting with "Top" found in folder: %s\n', folderList(i).name);
    end
end
 