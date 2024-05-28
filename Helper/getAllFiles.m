function [fileNames, filePaths] = getAllFiles(dirName)
    % Get the list of all files and folders in the given directory
    files = dir(dirName);
    fileNames = {};
    filePaths = {};

    % Loop through each item
    for i = 1:length(files)
        % Skip '.' and '..' folders
        if strcmp(files(i).name, '.') || strcmp(files(i).name, '..')
            continue;
        end
        
        % If it's a directory, recursively call getAllFiles
        if files(i).isdir
            subDir = fullfile(dirName, files(i).name);
            [subFileNames, subFilePaths] = getAllFiles(subDir);
            fileNames = [fileNames; subFileNames];
            filePaths = [filePaths; subFilePaths];
        else
            % If it's a file, add its name and path to the lists
            fileNames = [fileNames; files(i).name];
            filePaths = [filePaths; fullfile(dirName, files(i).name)];
        end
    end
end
