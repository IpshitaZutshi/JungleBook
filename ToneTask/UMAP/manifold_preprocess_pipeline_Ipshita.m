
% need to run
basepath = 'C:\Data\Kilosort\T17_240421_sess6\'

%% get subfolder names in the path
% pathFolder = 'G:\Data\Ipshita\NatureRevisions\'
% d = dir(pathFolder);
% isub = [d(:).isdir]; %# returns logical vector
% nameFolds = {d(isub).name}';
%%
behavior_epoch_index =[2,0,1,3,3,3,3,2];
speed_lims = 0;
smooths = 5;
SPIKEbins_beh = 0.1;



%% preprocess behavior

nameFolds = {'T17_240421_sess6'}
pathFolder = 'C:\Data\Kilosort\T17_240421_sess6\'

for bp = 1:length(nameFolds)
    nameFold = nameFolds{bp}
    basepath =[pathFolder,nameFold] 
    % Add masks
    %ipshita_add_masks(basepath,behavior_epoch_index(bp));

    for sb = 1:length(SPIKEbins_beh)
        SPIKEbin_beh = SPIKEbins_beh(sb);
        for sm =1:length(smooths)
            smooth_win = smooths(sm);
            for sl = 1:length(speed_lims)
                speed_lim = speed_lims(sl);
                manifold_preprocessing_ipshita('basepath',basepath,'dataType','behavior',...
                    'behavior_epoch_index',[behavior_epoch_index(bp)],...
                    'maze_type',{'SoundFreq'},...
                     'speed_lim',speed_lim,'smooth_win',smooth_win,...
                     'SPIKEbin_beh',SPIKEbin_beh);
                close all;
            end
        end
    end
end

%%
% %% plot behavior one by one
% 
% for TRIAL_TYPE = 0:7
% 
%  manifold_load_plot_Umap_ipshita('load_name','behavior','save_name','behavior',...
%     'umap_name','behavior','epoch_type','beh','dim1',1,'dim2',2,'dim3',3,'A',63,'E',-6 ,'TRIAL_TYPE',TRIAL_TYPE)
% end
% %%
% save_name = 'all'
% TRIAL_TYPE = [0,1,2,3,4,5,6,7]
% 
%  manifold_load_plot_Umap_ipshita('load_name','behavior','save_name','tone_lin_Trials',...
%     'umap_name','behavior','epoch_type','beh','dim1',1,'dim2',2,'dim3',3,'A',-7,'E',40 ,'TRIAL_TYPE',TRIAL_TYPE,...
%     'save_name',save_name)
%  %%
% save_name = 'noReturn'
% TRIAL_TYPE = [0,1,2,3,4,5,6]
% 
%  manifold_load_plot_Umap_ipshita('load_name','behavior','save_name','tone_lin_Trials',...
%     'umap_name','behavior','epoch_type','beh','dim1',1,'dim2',2,'dim3',3,'A',-7,'E',40 ,'TRIAL_TYPE',TRIAL_TYPE,...
%     'save_name',save_name)
% %%
% TRIAL_TYPE = [7]
% 
%  manifold_load_plot_Umap_ipshita('load_name','behavior','save_name','tone_Trials',...
%     'umap_name','behavior','epoch_type','beh','dim1',1,'dim2',2,'dim3',3,'A',-7,'E',40 ,'TRIAL_TYPE',TRIAL_TYPE)
% 
%  %% plot behavior with shuffle
%  manifold_load_plot_Umap('load_name','behavior_shuffle_c','save_name','behavior_shuffle_c',...
%     'umap_name','behavior_shuffle_c','epoch_type','beh','dim1',1,'dim2',2,'dim3',3 )
% 
% %% plot certain trial
% trial_number = 35
% save_name = ['trial_',num2str(trial_number)]
% TRIAL_TYPE = [0,1,2,3,4,5,6,7]
% 
%  manifold_load_plot_Umap_ipshita('load_name','behavior','save_name','tone_lin_Trials',...
%     'umap_name','behavior','epoch_type','beh','dim1',1,'dim2',2,'dim3',3,'A',-7,'E',40 ,'TRIAL_TYPE',TRIAL_TYPE,...
%     'save_name',save_name,'trial_number',trial_number)