function meanAccuracy = compileDecodedTrialChoice(varargin)

force = 0;
% 
% sess= {'IZ39\Final\IZ39_220714_sess18',...
%     'IZ40\Final\IZ40_220714_sess18',...
%     'IZ43\Final\IZ43_220919_sess14',...%'IZ43\Final\IZ43_220828_sess4',...
%     'IZ44\Final\IZ44_220830_sess7',...
%     'IZ47\Final\IZ47_230710_sess25',...
%     'IZ48\Final\IZ48_230705_sess22',...
%     }; 
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

p = inputParser;
addParameter(p,'plotfig',false,@islogical);
addParameter(p,'forFigure',true,@islogical);
addParameter(p,'bin_width',250,@isnumeric);
addParameter(p,'step_size',100,@isnumeric);

parse(p,varargin{:});
plotfig = p.Results.plotfig;
forFigure = p.Results.forFigure;
bin_width = p.Results.bin_width;
step_size = p.Results.step_size;

for ss = 1:length(sess)
    cd(strcat(expPath,sess{ss}))
    
    if ~exist('Decoding','dir') || force
        try
            decodeTrialChoice('plotfig',false,'eventType',0,'bin_width',bin_width,'step_size',step_size);
        catch
        end
        if ~forFigure
            decodeTrialChoice('plotfig',false,'eventType',1,'bin_width',bin_width,'step_size',step_size);
            decodeTrialChoice('plotfig',false,'eventType',2,'bin_width',bin_width,'step_size',step_size);
        end
    end
    
    result_names{ss} = strcat(expPath,sess{ss},'\Decoding\decoding_results_0');
       
    if exist(strcat(result_names{ss},'.mat'))
        load(strcat(result_names{ss},'.mat'))       
        mid_range = round(length(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results)/2);
        before = round(mid_range/10);
        meanAccuracy(ss) = max(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results(mid_range-before:mid_range));
    else
        meanAccuracy(ss) = nan;
    end
    
    if ~forFigure    
        result_names1{ss} = strcat(expPath,sess{ss},'\Decoding\decoding_results_1');
        result_names2{ss} = strcat(expPath,sess{ss},'\Decoding\decoding_results_2');
        pval_dir_name{ss} = strcat(expPath,sess{ss},'\Decoding\shuff_results\');

        conf0(:,:,ss) = nanmean(DECODING_RESULTS.NORMALIZED_RANK_RESULTS.confusion_matrix_results.rank_confusion_matrix,3);
        load(strcat(result_names1{ss},'.mat'))
        conf1(:,:,ss) = nanmean(DECODING_RESULTS.NORMALIZED_RANK_RESULTS.confusion_matrix_results.rank_confusion_matrix,3);

        load(strcat(result_names2{ss},'.mat'))
        conf2(:,:,ss) = nanmean(DECODING_RESULTS.NORMALIZED_RANK_RESULTS.confusion_matrix_results.rank_confusion_matrix,3);        
    end
    
end


if plotfig
    figure
    subplot(3,2,1)
    plot_obj = plot_standard_results_object(result_names);
    plot_obj.collapse_all_times_when_estimating_pvals = 1;
    plot_obj.significant_event_times = 5001;
    plot_obj.p_values = pval_dir_name;
    plot_obj.plot_results; 

    %Plot the confusion matrix
    subplot(3,2,2)
    confAvg = nanmean(conf0,3);
    imagesc(confAvg)
    colorbar

    subplot(3,2,3)
    plot_obj = plot_standard_results_object(result_names1);
    plot_obj.collapse_all_times_when_estimating_pvals = 1;
    plot_obj.significant_event_times = 5001;
    plot_obj.p_values = pval_dir_name;
    plot_obj.plot_results;

    %Plot the confusion matrix
    subplot(3,2,4)
    confAvg = nanmean(conf1,3);
    imagesc(confAvg)
    colorbar

    subplot(3,2,5)
    plot_obj = plot_standard_results_object(result_names2);
    plot_obj.collapse_all_times_when_estimating_pvals = 1;
    plot_obj.significant_event_times = 5001;
    plot_obj.p_values = pval_dir_name;
    plot_obj.plot_results; 

    %Plot the confusion matrix
    subplot(3,2,6)
    confAvg = nanmean(conf2,3);
    imagesc(confAvg)
    colorbar
end
end