function preProcessNeuropixelData

%% First copy all folders within a day into one master folder
preProcessFolderStructure;

%% Now, look through all folders, rename files, and concatenate .continuous dat file
% Get list of folders in the source directory
folders = dir(pwd);
% Find folders with the correct pattern
sessPattern = 'T\d+_\d{6}_sess\d+';
matchingFolders = regexp({folders.name}, sessPattern, 'match');
matchingIndices = ~cellfun(@isempty, matchingFolders);

for ii = find(matchingIndices)
    foldername = folders(ii).name;
    cd(strcat(folders(ii).folder,'\',foldername))

    [fileNames, filePaths] = getAllFiles(pwd);
    
    % If there is no .xml file, find any and copy it. 
    if ~exist(strcat(foldername,'.xml'))
        indices = find(strcmp(fileNames, 'continuous.xml'));
        if isempty(indices)
            print(strcat('No xml file within folder ',foldername))
            continue
        else
            copyfile(filePaths{indices(1)},strcat(folders(ii).folder,'\',foldername,'\',foldername,'.xml'),'f');
        end
    end    

    % Make sessionInfo file
    if ~exist(strcat(foldername,'.sessionInfo.mat'))
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true,'saveMat',true);
    end
    
    %% Concatenate dats if not already concatenated
    if ~exist(strcat(foldername,'.dat'))
        concatenateNeuropixelsDats('basepath',pwd);
    else
        disp(strcat('Data already concatenated for session ',foldername))
    end
end

end

