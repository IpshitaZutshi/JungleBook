
%Function that takes all behavioral files from one animal and returns a
%variable with all data
function allSessionsBehav = getAnimalBehav()

%1. Get patchBehavior folders
basepath = pwd;
allDirs = dir(pwd);
allDirs = allDirs([allDirs.isdir]); 

patchFolders = {};
for i = 1:length(allDirs)
    subfolder = fullfile(pwd, allDirs(i).name);
    % Check if there are any .avi files in this subfolder
    aviFiles = dir(fullfile(subfolder,'*', '*.avi'));
    % If there are any .avi files, add the subfolder name to the list
    if ~isempty(aviFiles)
        patchFolders{end+1} = aviFiles.folder; % Append subfolder name
    end
end
% If no folders with .avi files were found, display a message
if isempty(patchFolders)
    disp('No folders containing .avi files were found.');
end

%2. Read digitalIn file within pathFolders and extract behavior.
allSessionsBehav = {};  
for folder=1:length(patchFolders)
    try
        cd(patchFolders{folder});  % Change directory to current folder
        patchBehav = getPatchBehavior();  % Attempt to extract behavior data
        allSessionsBehav{folder} = patchBehav;  % Store the result
    catch ME
        disp(['Error in folder: ', patchFolders{folder}]);
        disp(['Error message: ', ME.message]);
        allSessionsBehav{folder} = [];  
        continue; 
    end
end
