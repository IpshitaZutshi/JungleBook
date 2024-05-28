% % Prompt the user to input the source folders
% sourceFolders = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ39\Final\,Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ40\Final\,Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43\Final\,Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ44\Final\,Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\,Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ48\Final\';
% 
% % Convert the comma-separated list of folders into a cell array
% sourceFolders = strsplit(sourceFolders, ',');
% 
% % Prompt the user to input the destination folder
% destinationFolder = 'C:\Data\DeepLabCut\Videos\';
% 
% % Loop through each source folder
% for i = 1:numel(sourceFolders)
%     % Get a list of all subfolders and files in the current source folder
%     files = dir(fullfile(sourceFolders{i}, '**', 'test*.avi'));
% 
%     % Loop through each file found
%     for j = 1:numel(files)
%         % Create the full path of the current file
%         sourceFile = fullfile(files(j).folder, files(j).name);
% 
%         % Create the destination path
%         [~, filename, ext] = fileparts(sourceFile);
%         destinationFile = fullfile(destinationFolder, [filename ext]);
% 
%         % Copy the file to the destination folder
%         copyfile(sourceFile, destinationFile);
%     end
% end
% 
% disp('Files copied successfully.');

% Specify the folder containing the CSV files
csvFolder = 'C:\Data\DeepLabCut\AuditoryTask-Ipshita-2024-03-15\videos';

% Specify the root folder containing the subdirectories with AVI files
rootFolder = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final';

% Get a list of all CSV files in the specified folder starting with 'test'
csvFiles = dir(fullfile(csvFolder, 'test*.csv'));

% Loop through each CSV file
for i = 1:numel(csvFiles)
    csvFile = fullfile(csvFiles(i).folder, csvFiles(i).name); % Full path of CSV file
    
    % Extract the base name of the CSV file (without extension)
    [~, baseName, ~] = fileparts(csvFile);
    
    baseName  = baseName(1:24);
    % Recursive function to search for AVI file in subdirectories
    destinationFolder = findMatchingFolder(rootFolder, baseName);
    
    if ~isempty(destinationFolder)
        % If an AVI file is found in a subdirectory, copy the CSV file to that subdirectory
        newCSVFileName = [baseName, '_DLC.csv']; % Append '_DLC' to the CSV file name
        destinationCSVFile = fullfile(destinationFolder, newCSVFileName);
        copyfile(csvFile, destinationCSVFile);
        disp(['Copied ', csvFiles(i).name, ' to ', destinationCSVFile]);
    else
        disp(['No matching AVI file found for ', csvFiles(i).name]);
    end
end

function destinationFolder = findMatchingFolder(currentFolder, baseName)
    % Recursive function to search for AVI file in subdirectories
    
    % Get a list of all subdirectories in the current folder
    subdirs = dir(currentFolder);
    
    % Initialize destination folder as empty
    destinationFolder = [];
    
    % Loop through each subdirectory
    for j = 1:numel(subdirs)
        if subdirs(j).isdir && ~strcmp(subdirs(j).name, '.') && ~strcmp(subdirs(j).name, '..')
            subdir = fullfile(currentFolder, subdirs(j).name); % Full path of current subdirectory
            
            % Search for an AVI file with the same base name in this subdirectory
            aviFiles = dir(fullfile(subdir, [baseName, '.avi']));
            
            if ~isempty(aviFiles)
                % If an AVI file is found, set the destination folder
                destinationFolder = subdir;
                return; % Exit the function once a matching AVI file is found
            else
                % If AVI file is not found, recursively call the function on this subdirectory
                destinationFolder = findMatchingFolder(subdir, baseName);
                
                if ~isempty(destinationFolder)
                    return; % Exit the function if a matching AVI file is found in any subdirectory
                end
            end
        end
    end
end


