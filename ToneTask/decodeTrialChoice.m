function decodeTrialChoice(varargin)

p = inputParser;
addParameter(p,'plotfig',true,@islogical);
addParameter(p,'shuff_num',10,@isnumeric);
addParameter(p,'eventType',0,@isnumeric); %if 0, choice, 1 toneGain, 2, home lick
addParameter(p,'bin_width',250,@isnumeric);
addParameter(p,'step_size',100,@isnumeric);

parse(p,varargin{:});
plotfig = p.Results.plotfig;
shuff_num = p.Results.shuff_num;
eventType = p.Results.eventType;
bin_width = p.Results.bin_width;
step_size = p.Results.step_size;

%% First generate the raster data
file = dir(['*.spikes.cellinfo.mat']);
load(file.name);

file = dir(['*.TrialBehavior.Behavior.mat']);
load(file.name);

sessname  = strsplit(file.folder,'\');
if ~exist(strcat('raster_data_ndt_',num2str(eventType)),'dir')
    mkdir(strcat('raster_data_ndt_',num2str(eventType)))
end

win = [-5000 5000];% in ms
%Generate raster windows 5 seconds before and after trial end

if eventType>0
    shuff_num = 0;
end

if isfield(behavTrials,'probe')
    if eventType==0
        st = behavTrials.timestamps(behavTrials.linTrial==0 & behavTrials.probe==0,2);
        trialType = behavTrials.lickLoc(behavTrials.linTrial==0 & behavTrials.probe==0);
    elseif eventType==1
        st = behavTrials.timestamps(behavTrials.linTrial==0 & behavTrials.probe==0,2);
        trialType = behavTrials.toneGain(behavTrials.linTrial==0 & behavTrials.probe==0);
    elseif eventType==2
        st = behavTrials.timestamps(behavTrials.linTrial==0 & behavTrials.probe==0,1);
        trialType = behavTrials.lickLoc(behavTrials.linTrial==0 & behavTrials.probe==0);
    end
else
    if eventType==0
        st = behavTrials.timestamps(behavTrials.linTrial==0,2);
        trialType = behavTrials.lickLoc(behavTrials.linTrial==0);
    elseif eventType==1
        st = behavTrials.timestamps(behavTrials.linTrial==0,2);
        trialType = behavTrials.toneGain(behavTrials.linTrial==0);
    elseif eventType==2
        st = behavTrials.timestamps(behavTrials.linTrial==0,1);
        trialType = behavTrials.lickLoc(behavTrials.linTrial==0);
    end
end

%% Make a separate raster matrix for each cell, which then can be processed by the toolbox
for ii = 1:size(spikes.UID,2)
    % Initiate the matrix, number of trials by time
    raster_data = zeros(length(st),(win(2)-win(1))+1);
    for kk = 1:length(st)
        temp_rast = (spikes.times{ii} - st(kk))*1000;% Convert to ms
        temp_rast = temp_rast(temp_rast>win(1) & temp_rast<win(2));
        temp_rast = round(temp_rast)+abs(win(1)); % Shift everything relative to the start of win(1)
        temp_rast(temp_rast==0) = 1; % If there's a random stray timepoint at timepoint 0, fix it so code does not break
        temp_rast(temp_rast>(win(2)-win(1))+1) = (win(2)-win(1))+1; % for the end of the timepoint
        raster_data(kk,temp_rast) = 1;
        
    end
    raster_site_info.sessionID = sessname{end};
    raster_site_info.unit = ii;
    raster_labels.trialType = trialType;
    save(strcat('raster_data_ndt_',num2str(eventType),'\',sessname{end},'_',num2str(ii),'_','raster_data.mat'),'raster_data','raster_labels','raster_site_info')
end
 
%% Everything below is just directly following the toolbox guidelines and code
%% Bin data
raster_file_directory_name = strcat('raster_data_ndt_',num2str(eventType),'\');
save_prefix_name = strcat('Binned_',sessname{end},'_',num2str(eventType));

create_binned_data_from_raster_data(raster_file_directory_name, save_prefix_name, bin_width, step_size);

%% Create a datasource object
binned_format_file_name = strcat('Binned_',sessname{end},'_',num2str(eventType),'_',num2str(bin_width),'ms_bins_',num2str(step_size),'ms_sampled.mat');

% will decode the identity of where the mouse licked
specific_label_name_to_use = 'trialType';
num_cv_splits = 3;

% create the CL object
the_classifier = libsvm_CL;%max_correlation_coefficient_CL;%
% the_classifier.kernel  = 'gaussian';
% the_classifier.gaussian_gamma  = 0.125;

% create a feature preprocessor that z-score normalizes each neuron
the_feature_preprocessors{1} = zscore_normalize_FP;
for num = 0:shuff_num
    if num == 0
        ds = basic_DS(binned_format_file_name, specific_label_name_to_use, num_cv_splits);
        ds.create_simultaneously_recorded_populations = 1;
        % create the CV object
        the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors);
        the_cross_validator.num_resample_runs = 20;

        % save the results
        save_file_name = strcat('Decoding\decoding_results_',num2str(eventType));

    else 
        %'Currently running shuffled label decoding results (data for the null distribution)'   
        ds = basic_DS(binned_format_file_name, specific_label_name_to_use, num_cv_splits);
        ds.create_simultaneously_recorded_populations = 1;        
        %ds.randomly_shuffle_labels_before_running = 1;% randomly shuffled the labels before running
        % Shuffling isn't working for some reason, so manually shuffle the
        % labels
        trialType1 = trialType(randperm(length(trialType)));
        for kk =1:length(ds.the_labels)
            ds.the_labels{kk} = trialType1;
        end
            
        % create the cross validator as before        
        the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors); 
        the_cross_validator.num_resample_runs = 10;  % only do 10 resample runs to save time

        % don't show progress to reduce visual clutter while the code is running
        the_cross_validator.display_progress.zero_one_loss = 0;  
        the_cross_validator.display_progress.resample_run_time = 0;

        % save the results with the appropriate name to the shuff_results/ directory
        save_file_name = ['Decoding\shuff_results\results_shuff_run_' num2str(num, '%03d')];  
    end

    % we will also greatly speed up the run-time of the analysis by not creating a full TC
    the_cross_validator.test_only_at_training_times = 0;

    % run the decoding analysis
    DECODING_RESULTS = the_cross_validator.run_cv_decoding;
    
    if ~exist('Decoding','dir')
        mkdir('Decoding')
        if ~exist('Decoding\shuff_results','dir')
            mkdir('Decoding\shuff_results')
        end
    end
    % save the results
    
    save(save_file_name, 'DECODING_RESULTS'); 
end

%% Plot the results
if plotfig
    % which results should be plotted (only have one result to plot here)
    % result_names{1} = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230710_sess25\Decoding\decoding_results';
    result_names{1} = strcat('Decoding\decoding_results_',num2str(eventType));
    % result_names{2} = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ48\Final\IZ48_230705_sess22\Decoding\decoding_results';
    % result_names{1} = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230710_sess25\Decoding\shuff_results\results_shuff_run_003';
    % result_names{2} = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ48\Final\IZ48_230705_sess22\Decoding\shuff_results\results_shuff_run_003';
    plot_obj = plot_standard_results_object(result_names);

    % % create the names of directories that contain the shuffled data for creating null distributions
    % pval_dir_name{1} = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ47\Final\IZ47_230710_sess25\Decoding\shuff_results\';
    % pval_dir_name{2} = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\IZ48\Final\IZ48_230705_sess22\Decoding\shuff_results\';
    pval_dir_name{1} = 'Decoding\shuff_results\';
    plot_obj.p_values = pval_dir_name;
    % 
    % % use data from all time bins when creating the null distribution
    plot_obj.collapse_all_times_when_estimating_pvals = 1;
    plot_obj.significant_event_times = 5001;
    %plot_obj.p_value_alpha_level = 0.001;
    plot_obj.plot_results;   % actually plot the results
end

end