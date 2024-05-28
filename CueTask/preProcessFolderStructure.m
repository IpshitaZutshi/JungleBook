function preProcessFolderStructure

sourceFolder = pwd;
destBaseFolder = pwd;

% Get list of folders in the source directory
folders = dir(sourceFolder);
folders = folders([folders.isdir]);
folders = folders(~ismember({folders.name}, {'.', '..'}));

% Find folders with date pattern
datePattern = '\d{4}-\d{2}-\d{2}';
matchingFolders = regexp({folders.name}, datePattern, 'match');
matchingIndices = ~cellfun(@isempty, matchingFolders);

% Extract session numbers and find most recent allSessionNumbers
sessPattern = 'T\d+_\d{6}_sess\d+';
matchingSess = regexp({folders.name}, sessPattern, 'match');
matchingSess = [matchingSess{:}];  % Flatten cell array
sessNumbers = cellfun(@(x) str2double(x(strfind(x, 'sess')+4:end)), matchingSess);
maxSessionNumber = max(sessNumbers);
count=0;
% Loop through matching folders and move their contents

for i = find(matchingIndices)
    % Extract date from folder name
    folderName = folders(i).name;
    dateMatch = regexp(folderName, datePattern, 'match');
    if isempty(dateMatch)
        continue; % Skip if date not found
    end
    folderDate = dateMatch{1};
    
    % Determine matching sessions for the current folder
    sessPattern = sprintf('%s_\\d{2}-\\d{2}-\\d{2}_\\d{3}', folderDate);
    matchingSessions = regexp({folders(matchingIndices).name}, sessPattern, 'match');   
    matchingSessions = matchingSessions{i};
    if isempty(matchingSessions)
        continue; % Skip if no more sessions of that date
    end

    % Create destination folder name
    destDateFolder = strrep(folderDate, '-', '');
    count = count+1;
    % Set destination session folder
    destSessionFolder = sprintf('T17_%s_sess%d', destDateFolder(3:end), maxSessionNumber + 1 + count);
    destPath = fullfile(destBaseFolder, destSessionFolder,folderName);
    
    % Create destination folder if it doesn't exist
    if ~exist(destPath, 'dir')
        mkdir(destPath);
    end
    
    % Move contents of source folder to destination folder
    sourcePath = fullfile(sourceFolder, folderName);
    movefile(fullfile(sourcePath, '*'), destPath);

    % Remove the source folder
    rmdir(sourcePath);
end


end