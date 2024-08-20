function manifold_preprocessing_ipshita(varargin)

p = inputParser;
addParameter(p,'basepath',pwd);
addParameter(p,'behavior_epoch_index',[]);
addParameter(p,'maze_type',{'Tmaze','Tmaze'});
addParameter(p,'speed_lim',0);
addParameter(p,'dataType','behavior'); % or 'behavior'
addParameter(p,'shuffle_type',{});
addParameter(p,'num_shuffle',10);
addParameter(p,'sub_cell','all');
addParameter(p,'downsample_cell',false);
addParameter(p,'down_cell_num',60);
addParameter(p,'save_name',[]);
addParameter(p,'theta_source','REM_cycles');
addParameter(p,'brain_region',{});
addParameter(p,'exclude_tuned_cell',0);
addParameter(p,'tuned_cell',0);
addParameter(p,'seed',0); % 
addParameter(p,'smooth_win',5); % 
addParameter(p,'SPIKEbin_beh',0.1); % 
addParameter(p,'correct_only',0); % 
addParameter(p,'error_only',0); % 


parse(p,varargin{:});
error_only = p.Results.error_only;
correct_only = p.Results.correct_only;
SPIKEbin_beh = p.Results.SPIKEbin_beh;
basepath = p.Results.basepath;
behavior_epoch_index = p.Results.behavior_epoch_index;
maze_type = p.Results.maze_type;
speed_lim = p.Results.speed_lim;
dataType = p.Results.dataType;
shuffle_type = p.Results.shuffle_type;
num_shuffle = p.Results.num_shuffle;
sub_cell = p.Results.sub_cell;
downsample_cell = p.Results.downsample_cell;
down_cell_num = p.Results.down_cell_num;
save_name = p.Results.save_name;
theta_source = p.Results.theta_source;
brain_region = p.Results.brain_region;
exclude_tuned_cell = p.Results.exclude_tuned_cell;
tuned_cell = p.Results.tuned_cell;
seed = p.Results.seed;
smooth_win = p.Results.smooth_win;
% behavior_epoch_index = [2];
% maze_type  ={'Tmaze'};
% speed_lim = 0;
% dataType = 'RSE';
% shuffle_type ={'irc'};
% num_shuffle = 10;


%%
basename = basenameFromBasepath(basepath);
if isempty(save_name)
%     save_name = ['behavior_speed_',num2str(speed_lim),'_smooth_',num2str(smooth_win),'_Tbin_',num2str(SPIKEbin_beh)];
    if correct_only
        save_name = ['behavior_speed_',num2str(speed_lim),'_smooth_',num2str(smooth_win),...
            '_bin_',num2str(SPIKEbin_beh),'_correct_only'];
    elseif error_only
        save_name = ['behavior_speed_',num2str(speed_lim),'_smooth_',num2str(smooth_win),...
            '_bin_',num2str(SPIKEbin_beh),'_error_only'];    
    else
        save_name = ['behavior_speed_',num2str(speed_lim),'_smooth_',num2str(smooth_win),...
            '_bin_',num2str(SPIKEbin_beh)];
    end

end

%% behavior only 

switch dataType
    case 'behavior'

        manifold_saveData_zscore_ipshita('basepath',basepath,...
            'behavior_epoch_index',behavior_epoch_index,...
            'save_name',save_name, 'maze_type',maze_type,...
           'save_spikemat',true,'speed_lim',speed_lim,...
           'sub_cell',sub_cell,'brain_region',brain_region,...
           'down_cell_num',down_cell_num,'downsample_cell',downsample_cell,'exclude_tuned_cell',exclude_tuned_cell,...
            'tuned_cell',tuned_cell,...
            'seed',seed,'smooth_win',smooth_win,'SPIKEbin_beh',SPIKEbin_beh,...
            'correct_only',correct_only,'error_only',error_only);

        % shuffle
        if ~isempty(shuffle_type)
            for sf = 1:length(shuffle_type)
                st = shuffle_type{sf};
    
                save_name_shuffle = [save_name,'_shuffle_', st];
                manifold_saveData_zscore_ipshita('behavior_epoch_index',behavior_epoch_index,...
                    'save_name', save_name_shuffle, 'maze_type',maze_type,...
                'save_spikemat',true,'speed_lim',speed_lim,...
                'shuffle_beh',true,'shuffle_type',st,'num_shuffle',num_shuffle,...
                'sub_cell',sub_cell,'brain_region',brain_region,'exclude_tuned_cell',exclude_tuned_cell,...
                'tuned_cell',tuned_cell,'smooth_win',smooth_win);
    
            end
        end




end







