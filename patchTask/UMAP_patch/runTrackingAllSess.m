%
% USAGE
%   runs tracking behavior for all folders
%
% INPUTS 
%    start in top directory of mouse
%
%    =========================================================================


direc = '\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\patchTask\N14';

folders = dir(direc);  % list all files/folders in direc

for i = 1:length(folders)
    if folders(i).isdir && ~ismember(folders(i).name, {'.', '..'}) % only folders, ignore '.' and '..'
        folderName = folders(i).name;
        if contains(folderName, 'sess')
            currFolder = fullfile(direc, folderName);
            fprintf('Processing folder: %s\n', currFolder);
            cd(currFolder);
            basepath = pwd;

            % get tracking
            if isempty (dir(fullfile(basepath, '*Tracking.Behavior.mat')))
                getPatchTracking() 
            end

            cd(direc);
        end
    end
end


