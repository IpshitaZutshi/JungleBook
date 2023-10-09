% List of provided folders to search for "Final" subfolder
sourceFolders = {'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ40',...
    'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ43',...
    'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ44',...
    'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47',...
    'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ48',...
    'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ41',...
    'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ45',...
    'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ46'};

%'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ39'

% Destination folder where you want to copy the .mat files
destinationFolder = 'E:\';

% Loop through each source folder
for folderIdx = 1:length(sourceFolders)
    sourceFolder = sourceFolders{folderIdx};
    C = strsplit(sourceFolder,'\');
    
    % Check if the source folder exists
    if exist(sourceFolder, 'dir') ~= 7
        fprintf('Source folder "%s" does not exist.\n', sourceFolder);
        continue;
    end
    
    % Construct the path to the "Final" subfolder
    finalFolder = fullfile(sourceFolder, 'Final');
    
    % Check if the "Final" subfolder exists within the source folder
    if exist(finalFolder, 'dir') ~= 7
        fprintf('"%s" does not have a "Final" subfolder.\n', sourceFolder);
        continue;
    end
    
    % Create the corresponding destination folder structure
    destinationFinalFolder = fullfile(destinationFolder, C{end}, 'Final');
    
    % List all .mat files in the "Final" subfolder and its subdirectories
    matFiles = dir(fullfile(finalFolder, '**', '*.mat'));
    
    % Loop through each .mat file and copy it to the destination folder
    for fileIdx = 1:length(matFiles)
        D = strsplit(matFiles(fileIdx).folder,'Final\');
        matFileName = matFiles(fileIdx).name;
        sourceFilePath = fullfile(matFiles(fileIdx).folder, matFileName);
        destinationFilePath = fullfile(destinationFinalFolder,D{2},matFileName);
        
        % Create the destination folder if it does not exist
        destinationFileDir = fileparts(destinationFilePath);
        if ~exist(destinationFileDir, 'dir')
            mkdir(destinationFileDir);
        end
        
        % Copy the .mat file to the destination folder
        copyfile(sourceFilePath, destinationFilePath);
        
        fprintf('Copied file: %s\n', destinationFilePath);
    end
end

fprintf('Task completed successfully.\n');
