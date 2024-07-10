    sess = {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
        'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
        'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   23
        'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...
        'IZ40\Final\IZ40_220708_sess17','IZ40\Final\IZ40_220714_sess18',...27
        'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
        'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
        'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
        'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    36
        'IZ44\Final\IZ44_220827_sess4', 'IZ44\Final\IZ44_220828_sess5',...
        'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...
        'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
        'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',... 45
        'IZ47\Final\IZ47_230626_sess15','IZ47\Final\IZ47_230707_sess24',...
        'IZ47\Final\IZ47_230710_sess25','IZ47\Final\IZ47_230712_sess27',...49
        'IZ48\Final\IZ48_230628_sess17','IZ48\Final\IZ48_230703_sess21',...
        'IZ48\Final\IZ48_230705_sess22','IZ48\Final\IZ48_230714_sess28'};  

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';
umap_path = 'E:\Data\';%'Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions';

for ss = 1:length(sess)
    cd(strcat(expPath,sess{ss}))

    destPath = strcat(umap_path,'\',sess{ss});

    %% If not exist, make the dest path folder
    if ~exist(destPath, 'dir')
        mkdir(destPath);
    end

    %% Copy each of the following files within the directory to the destPath folder
    % files_to_copy = {'*.spikes.cellinfo.mat', '*.Tracking.Behavior.mat', '*TrialBehavior.Behavior.mat', ...
    %                  '*cell_metrics.cellinfo.mat', '*.DigitalIn.events.mat', '*.MergePoints.events.mat', ...
    %                  '*.session.mat', '*.sessionInfo.mat','*rateMapsAvg.cellinfo.mat','*rateMapsAvg2D.cellinfo.mat',...
    %                  '*rateMapsAvgError.cellinfo.mat','*rateMapsAvgLickLocProbe.cellinfo.mat','*rateMapsAvgnotLog.cellinfo.mat',...
    %                  '*rateMapsTrial.cellinfo.mat','*.spikeData.cellinfo',};
    files_to_copy = {'*.thetaLFP.mat',};

    for f = 1:length(files_to_copy)
        file = dir(files_to_copy{f});
        for k = 1:length(file)
            copyfile(file(k).name, fullfile(destPath, file(k).name));
        end
    end

    % % Also copy video files
    % load(strcat(sess{ss}(12:end),'.MergePoints.events.mat'));
    % for ii = 1:size(MergePoints.foldernames,2)
    %      if ~isempty(dir([MergePoints.foldernames{ii} filesep 'test*']))   %if it has a video file
    %          list = dir([MergePoints.foldernames{ii} filesep 'test*']);
    %          if ~exist(fullfile(destPath,MergePoints.foldernames{ii}), 'dir')
    %             mkdir(fullfile(destPath,MergePoints.foldernames{ii}));
    %          end
    %          for k = 1:length(list)
    %             copyfile(fullfile(list(k).folder,list(k).name),fullfile(destPath,MergePoints.foldernames{ii}, list(k).name));
    %          end
    %          list = dir([MergePoints.foldernames{ii} filesep 'front*']);
    %          for k = 1:length(list)
    %             copyfile(fullfile(list(k).folder,list(k).name),fullfile(destPath,MergePoints.foldernames{ii}, list(k).name));
    %          end
    %      end
    % end
    % 
    % %Also copy maps folders
    % destPath1 = strcat(umap_path,'\',sess{ss},'\Maps');
    % if ~exist(destPath1, 'dir')
    %     mkdir(destPath1);
    % end
    % copyfile('Maps\', destPath1, 'f');

end
