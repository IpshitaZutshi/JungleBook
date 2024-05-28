sess = {'IZ39\Final\IZ39_220714_sess18','IZ39\Final\IZ39_220707_sess17','IZ39\Final\IZ39_220702_sess14',...
    'IZ40\Final\IZ40_220707_sess16','IZ40\Final\IZ40_220708_sess17','IZ40\Final\IZ40_220714_sess18',...
    'IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220828_sess4','IZ43\Final\IZ43_220919_sess14',...
    'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220830_sess7','IZ44\Final\IZ44_220919_sess14',...
    'IZ47\Final\IZ47_230626_sess15','IZ47\Final\IZ47_230707_sess24','IZ47\Final\IZ47_230710_sess25','IZ47\Final\IZ47_230712_sess27',...
    'IZ48\Final\IZ48_230628_sess17','IZ48\Final\IZ48_230703_sess21','IZ48\Final\IZ48_230705_sess22','IZ48\Final\IZ48_230714_sess28'};

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';
umap_path = 'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions';

for ss = 1:length(sess)
    cd(strcat(expPath,sess{ss}))

    destPath = strcat(umap_path,'\',sess{ss}(12:end));

    %% If not exist, make the dest path folder
    if ~exist(destPath, 'dir')
        mkdir(destPath);
    end

    %% Copy each of the following files within the directory to the destPath folder
    files_to_copy = {'*.spikes.cellinfo.mat', '*.Tracking.Behavior.mat', '*TrialBehavior.Behavior.mat', ...
                     '*cell_metrics.cellinfo.mat', '*.DigitalIn.events.mat', '*.MergePoints.events.mat', ...
                     '*.session.mat', '*.sessionInfo.mat'};

    for f = 1:length(files_to_copy)
        file = dir(files_to_copy{f});
        for k = 1:length(file)
            copyfile(file(k).name, fullfile(destPath, file(k).name));
        end
    end

end
