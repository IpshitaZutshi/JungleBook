% This function calculates the percentage of cells for the chosen session
% that are tuned to each of the variables that we have chosen to test with
% the PGAM. 
% The function also identifies the cells that are tuned to licks and plots
% their PSTH. 

function [resultsPGAM] = postprocessPGAM_final(varargin)

close all 

p = inputParser;
addParameter(p,'basepath','Z:\Homes\zutshi01\Recordings\Auditory_Task',@isstr);
addParameter(p,'fits','Z:\Buzsakilabspace\LabShare\AthinaApostolelli\PGAM\results_IZ\new_vars2',@isstr);
addParameter(p,'data','C:\Data\PGAMAnalysis\data_analyzed',@isstr);
addParameter(p,'postprocess','Z:\Homes\zutshi01\Recordings\Auditory_Task\GAM fits\postprocess',@isstr);

parse(p,varargin{:});
basepath = p.Results.basepath;
fits = p.Results.fits;
data = p.Results.data;
postprocess = p.Results.postprocess;

sessions = {'IZ39\Final\IZ39_220622_sess8', 'IZ39\Final\IZ39_220624_sess10', 'IZ39\Final\IZ39_220629_sess12',...
    'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
    'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...  
    'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...
    'IZ40\Final\IZ40_220708_sess17','IZ40\Final\IZ40_220714_sess18',...
    'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
    'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    
    'IZ44\Final\IZ44_220827_sess4','IZ44\Final\IZ44_220828_sess5',...
    'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...
 'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
 'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15'};

for s = 1:length(sessions)
    session = extractAfter(sessions{s}, 'Final\');
    baseDir = fullfile(basepath, sessions{s});
    dataDir = fullfile(data);
    resultsDir = fullfile(fits, session);
    postprocessDir = fullfile(postprocess, session);

    if ~isfolder(postprocessDir)
        mkdir(postprocessDir)
    end
    
    % Load spiking data
    if ~isempty(dir([baseDir filesep '*.spikes.cellinfo.mat']))
        disp('Spikes already detected! Loading file.');
        file = dir([baseDir filesep '*.spikes.cellinfo.mat']);
        load(fullfile(baseDir, file.name));
    end
    
    % Load behavior data
    if ~isempty(dir([baseDir filesep '*TrialBehavior.Behavior.mat'])) 
        disp('Behavior already detected! Loading file.');
        file = dir([baseDir filesep '*TrialBehavior.Behavior.mat']);
        load(fullfile(baseDir, file.name));
    end
    
    % Load event data 
    if ~isempty(dir([dataDir filesep session '.sessionDataPGAM.mat']))
        disp('Variables already detected! Loading file.');
        file = dir([dataDir filesep session '.sessionDataPGAM.mat']);
        load(fullfile(dataDir, file.name));
    end
    
    % Load cell metrics
    if ~isempty(dir([baseDir filesep '*.cell_metrics.cellinfo.mat']))
        disp('Cell metrics already detected! Loading file.');
        file = dir([baseDir filesep '*.cell_metrics.cellinfo.mat']);
        load(fullfile(baseDir, file.name));
    end
    
    % Load digitalIn
    if ~isempty(dir([baseDir filesep '*.DigitalIn.events.mat']))
        disp('Digital input already detected! Loading file.');
        file = dir([baseDir filesep '*.DigitalIn.events.mat']);
        load(fullfile(baseDir, file.name));
    end
    
    % NOTE: numNeurons is different to number of files, because cells with
    % firing rate < 0.5 have been excluded from the PGAM fits.

    % Only select pyramidal neurons 
    resultsPGAM = {};
    
    numNeurons = spikes.numcells;
    pyramidalCells = find(contains(cell_metrics.putativeCellType, 'Pyramidal Cell'));
    
    resultsPGAM.stats.numNeurons = numNeurons;
    resultsPGAM.stats.pyramidalCells = pyramidalCells';
    
    % Load PGAM results
    filelist = dir(fullfile(resultsDir, '*k-fold.csv'));
    filename = {filelist.name}';
       
    % Select neurons that have been fit
    neurons = [];
    for i = 1:numNeurons
        if exist(fullfile(resultsDir, strcat('spatial_neuron_', ...
            string(i), '_fit_k-fold.csv')))

            neurons = [neurons, i];
        else
            continue
        end
    end

    % Select pyramidal cells that have been fit
    pyramidal = [];
    for i = 1:length(pyramidalCells)
        if exist(fullfile(resultsDir, strcat('spatial_neuron_', ...
            string(pyramidalCells(i)), '_fit_k-fold.csv')))

            pyramidal = [pyramidal, pyramidalCells(i)];
        else
            continue
        end
    end

    resultsPGAM.stats.neuronsFit = neurons;
    resultsPGAM.stats.pyramidalFit = pyramidal;

    %% Find the cells tuned to each variable 

    % Pyramidal cells
    tuned_y = zeros(length(pyramidal),1);
    tuned_ylin = zeros(length(pyramidal),1);
    tuned_yrev = zeros(length(pyramidal),1);
    tuned_relDistStop = zeros(length(pyramidal),1);
    tuned_licks = zeros(length(pyramidal),1);
    
    for i = 1:length(pyramidal)           
        results = readtable(fullfile(resultsDir, strcat('spatial_neuron_', ...
            string(pyramidal(i)), '_fit_k-fold.csv')));
        
        % Variables licks is in the last row and the pval is in column 10
        if table2array(results(1, 10)) < 0.001
            tuned_y(i, 1) = 1; 
        end
        if table2array(results(2, 10)) < 0.001
            tuned_ylin(i, 1) = 1; 
        end
        if table2array(results(3, 10)) < 0.001
            tuned_yrev(i, 1) = 1; 
        end        
        if table2array(results(4, 10)) < 0.001
            tuned_relDistStop(i, 1) = 1; 
        end
        if table2array(results(10, 10)) < 0.001
            tuned_licks(i, 1) = 1; 
        end
    end

    resultsPGAM.tuned_y = pyramidal(tuned_y == 1);
    resultsPGAM.tuned_ylin = pyramidal(tuned_ylin == 1);
    resultsPGAM.tuned_yrev = pyramidal(tuned_yrev == 1);
    resultsPGAM.tuned_relDistStop = pyramidal(tuned_relDistStop == 1);
    resultsPGAM.tuned_licks = pyramidal(tuned_licks == 1);
    
   
    % All neurons
    tuned_y_all = zeros(length(neurons),1);
    tuned_ylin_all = zeros(length(neurons),1);
    tuned_yrev_all = zeros(length(neurons),1);
    tuned_relDistStop_all = zeros(length(neurons),1);
    tuned_licks_all = zeros(length(neurons),1);
    
    for i = 1:length(neurons)
        results = readtable(fullfile(resultsDir, strcat('spatial_neuron_', ...
            string(neurons(i)), '_fit_k-fold.csv')));
    
        if table2array(results(1, 10)) < 0.001
            tuned_y_all(i, 1) = 1; 
        end
        if table2array(results(2, 10)) < 0.001
            tuned_ylin_all(i, 1) = 1; 
        end
        if table2array(results(3, 10)) < 0.001
            tuned_yrev_all(i, 1) = 1; 
        end        
        if table2array(results(4, 10)) < 0.001
            tuned_relDistStop_all(i, 1) = 1; 
        end
        if table2array(results(10, 10)) < 0.001
            tuned_licks_all(i, 1) = 1; 
        end
    end
    
    resultsPGAM.tuned_y_all = neurons(tuned_y_all == 1);
    resultsPGAM.tuned_ylin_all = neurons(tuned_ylin_all == 1);
    resultsPGAM.tuned_yrev_all = neurons(tuned_yrev_all == 1);
    resultsPGAM.tuned_relDistStop_all = neurons(tuned_relDistStop_all == 1);
    resultsPGAM.tuned_licks_all = neurons(tuned_licks_all == 1);

    %% Calculate proportion of cells tuned to each variable. 

    % Pyramidal cells
    yProp = length(find(tuned_y == 1)) / length(pyramidalCells);
    ylinProp = length(find(tuned_ylin == 1)) / length(pyramidalCells);
    yRevProp  = length(find(tuned_yrev == 1)) / length(pyramidalCells);
    relDistStopProp = length(find(tuned_relDistStop == 1)) / length(pyramidalCells);
    licksProp = length(find(tuned_licks == 1)) / length(pyramidalCells);
    
    resultsPGAM.proportions.yProp = yProp;
    resultsPGAM.proportions.ylinProp = ylinProp;
    resultsPGAM.proportions.yRevProp = yRevProp;
    resultsPGAM.proportions.relDistStopProp = relDistStopProp;
    resultsPGAM.proportions.licksProp = licksProp;
    
    % All neurons
    yProp_all = length(find(tuned_y_all == 1)) / length(numNeurons);
    ylinProp_all = length(find(tuned_ylin_all == 1)) / length(numNeurons);
    yrevProp_all = length(find(tuned_yrev_all == 1)) / length(numNeurons);
    relDistStopProp_all = length(find(tuned_relDistStop_all == 1)) / length(numNeurons);
    licksProp_all = length(find(tuned_licks_all == 1)) / length(numNeurons);

    resultsPGAM.proportions.yProp_all = yProp_all;
    resultsPGAM.proportions.ylinProp_all = ylinProp_all;
    resultsPGAM.proportions.yrevProp_all = yrevProp_all;
    resultsPGAM.proportions.relDistStopProp_all = relDistStopProp_all;
    resultsPGAM.proportions.licksProp_all = licksProp_all;

    %% Combine variables

    % Pyramidal cells - tunings
    tuned_y_relDistStop = find(tuned_y == 1 & tuned_relDistStop == 1);
    tuned_licks_relDistStop = find(tuned_relDistStop == 1 & tuned_licks == 1);
    tuned_y_licks = find(tuned_licks == 1 & tuned_y == 1);
    tuned_y_licks_relDistStop = find(tuned_y == 1 & tuned_licks == 1 & tuned_relDistStop == 1);

    resultsPGAM.y_relDistStop = pyramidal(tuned_y_relDistStop);
    resultsPGAM.licks_relDistStop = pyramidal(tuned_licks_relDistStop);
    resultsPGAM.y_licks = pyramidal(tuned_y_licks);
    resultsPGAM.y_licks_relDistStop = pyramidal(tuned_y_licks_relDistStop);
    
    % Pyramidal cells - proportions
    y_relDistStop_prop = length(tuned_y_relDistStop) / length(pyramidalCells);
    licks_relDistStop_prop = length(tuned_licks_relDistStop) / length(pyramidalCells);
    y_licks_prop = length(tuned_y_licks) / length(pyramidalCells);
    y_licks_relDistStop_prop = length(tuned_y_licks_relDistStop) / length(pyramidalCells);
    
    resultsPGAM.proportions.y_relDistStop_prop = y_relDistStop_prop;
    resultsPGAM.proportions.licks_relDistStop_prop = licks_relDistStop_prop;
    resultsPGAM.proportions.y_licks_prop = y_licks_prop;
    resultsPGAM.proportions.y_licks_relDistStop_prop = y_licks_relDistStop_prop;

    % All neurons - tunings
    tuned_y_relDistStop_all = find(tuned_y_all == 1 & tuned_relDistStop_all == 1);
    tuned_licks_relDistStop_all = find(tuned_relDistStop_all == 1 & tuned_licks_all == 1);
    tuned_y_licks_all = find(tuned_licks_all == 1 & tuned_y_all == 1);
    tuned_y_licks_relDistStop_all = find(tuned_y_all == 1 & tuned_licks_all == 1 & tuned_relDistStop_all == 1);

    resultsPGAM.y_relDistStop_all = neurons(tuned_y_relDistStop_all);
    resultsPGAM.licks_relDistStop_all = neurons(tuned_licks_relDistStop_all);
    resultsPGAM.y_licks_all = neurons(tuned_y_licks_all);
    resultsPGAM.y_licks_relDistStop_all = neurons(tuned_y_licks_relDistStop_all);

    % All neurons - proportions
    y_relDistStop_prop_all = length(tuned_y_relDistStop_all) / length(numNeurons);
    licks_relDistStop_prop_all = length(tuned_licks_relDistStop_all) / length(numNeurons);
    y_licks_prop_all = length(tuned_y_licks_all) / length(numNeurons);
    y_licks_relDistStop_prop_all = length(tuned_y_licks_relDistStop_all) / length(numNeurons);
    
    resultsPGAM.proportions.y_relDistStop_prop_all = y_relDistStop_prop_all;
    resultsPGAM.proportions.licks_relDistStop_prop_all = licks_relDistStop_prop_all;
    resultsPGAM.proportions.y_licks_prop_all = y_licks_prop_all;
    resultsPGAM.proportions.y_licks_relDistStop_prop_all = y_licks_relDistStop_prop_all;

    
    %% TODO: Calculate the PSTH    
   
    %% Select cells tuned to licks with positive kernel (pyramidal only).
    tunedLicks = resultsPGAM.tuned_licks;

    kernelStrength_licks = zeros(length(tunedLicks),1);
    
    for i = 1:length(tunedLicks)
        results = readtable(fullfile(resultsDir, strcat('spatial_neuron_', ...
            string(tunedLicks(i)), '_fit_k-fold.csv')));
           
        if table2array(results(9, 18)) > 0
            kernelStrength_licks(i, 1) = 1; 
        end
    end
    
    resultsPGAM.tuned_licksP = tunedLicks(kernelStrength_licks==1)';
    
    positiveLicksProp = length(tunedLicks(kernelStrength_licks==1)) / length(pyramidalCells);
    resultsPGAM.proportions.positiveLicksProp = positiveLicksProp;

%     %% Select cells tuned to y and ylin with fields in the forward run 
%     % (pyrmidal only). To determine that, we look at the raw firing rate
%     % used for training the PGAM. 
% 
%     tunedY = resultsPGAM.tuned_y;
%     tunedYlin = resultsPGAM.tuned_ylin;
% 
%     forward_y = zeros(length(tunedY),1);
%     forward_ylin = zeros(length(tunedYlin),1);
% 
%     for i = 1:length(tunedY)
%         results = readtable(fullfile(resultsDir, strcat('spatial_neuron_', ...
%             string(tunedY(i)), '_fit_k-fold.csv')));
% 
%         rates = table2array(results(1,14));
%         if ~iscell(rates)
%             continue
%         else
%             rates = rates{1};
%             stringData = strrep(rates, '[', '');
%             stringData = strrep(stringData, ']', '');
%             stringData = strrep(stringData, '↵', '');
%             raw_rate = regexp(stringData, '\s+', 'split');        
%             raw_rate = str2double(raw_rate);
%         end
% 
%         % The 240cm forward and reverse virtual linear track is split into
%         % 15 bins for the PGAM fits. To only consider the forward
%         % direction, we will look at the first half of the bins (8 bins -
%         % 128 cm). 
%         max_rate = max(raw_rate);
%         if intersect(max_rate, raw_rate(1:8))
%             forward_y(i) = 1;
%         end     
%     end
% 
%     for i = 1:length(tunedYlin)
%         results = readtable(fullfile(resultsDir, strcat('spatial_neuron_', ...
%             string(tunedYlin(i)), '_fit_k-fold.csv')));
%         
%         rates = table2array(results(2,14));
%         if ~iscell(rates)
%             continue
%         else
%             rates = rates{1};
%             stringData = strrep(rates, '[', '');
%             stringData = strrep(stringData, ']', '');
%             stringData = strrep(stringData, '↵', '');
%             raw_rate = regexp(stringData, '\s+', 'split');        
%             raw_rate = str2double(raw_rate);
%         end
% 
%         max_rate = max(raw_rate);
%         if intersect(max_rate, raw_rate(1:8))
%             forward_ylin(i) = 1;
%         end     
%     end
% 
%     resultsPGAM.tuned_YFor = tunedY(forward_y==1)';
%     resultsPGAM.tuned_YForLin = tunedYlin(forward_ylin==1)';
% 
%     forwardYProp = length(tunedY(forward_y==1)) / length(pyramidal);
%     forwardYLinProp = length(tunedYlin(forward_ylin==1)) / length(pyramidal);
% 
%     resultsPGAM.proportions.YForProp = forwardYProp;
%     resultsPGAM.proportions.YForLinProp = forwardYLinProp;

    %% Calculate new tunings and proportions 
    tuned_yFor_relDistStop = intersect(tunedY(forward_y==1), pyramidal(tuned_relDistStop==1))';
    tuned_licksP_relDistStop = intersect(tunedLicks(kernelStrength_licks == 1), pyramidal(tuned_relDistStop==1))';
    tuned_yFor_licksP = intersect(tunedLicks(kernelStrength_licks == 1), tunedY(forward_y==1))';
    tuned_y_licksP = intersect(tunedLicks(kernelStrength_licks == 1), tunedY)';
    tuned_yFor_licksP_relDistStop = intersect(tunedLicks(kernelStrength_licks == 1), ...
        intersect(tunedY(forward_y==1), pyramidal(tuned_relDistStop==1)));

    resultsPGAM.yFor_relDistStop = tuned_yFor_relDistStop;
    resultsPGAM.licksP_relDistStop = tuned_licksP_relDistStop;
    resultsPGAM.yFor_licksP = tuned_yFor_licksP;
    resultsPGAM.yFor_licksp_relDistStop = tuned_yFor_licksP_relDistStop;
    
    yFor_relDistStop_prop = length(tuned_yFor_relDistStop) / length(pyramidal);
    licksP_relDistStop_prop = length(tuned_licksP_relDistStop) / length(pyramidal);
    yFor_licksP_prop = length(tuned_yFor_licksP) / length(pyramidal);
    yFor_licksp_relDistStop_prop = length(tuned_yFor_licksP_relDistStop) / length(pyramidal);
    y_licksP_prop = length(tuned_y_licksP) / length(pyramidal);
    
    resultsPGAM.proportions.yFor_relDistStop_prop = yFor_relDistStop_prop;
    resultsPGAM.proportions.licksP_relDistStop_prop = licksP_relDistStop_prop;
    resultsPGAM.proportions.yFor_licksp_prop = yFor_licksP_prop;
    resultsPGAM.proportions.y_licksP_prop = y_licksP_prop;
    resultsPGAM.proportions.yFor_licksp_relDistStop_prop = yFor_licksp_relDistStop_prop;

    %% Select cells tuned to relDistStop based on the mutual information. 
    % Similarly, select cells tuned to y based on the mutual information 
    % (pyramidal only).

    % y vs relDistStop 
    cells_yFor_relDistStop = intersect(tunedY(forward_y==1), pyramidal(tuned_relDistStop==1));
    
    mi_yFor = zeros(length(cells_yFor_relDistStop), 1);
    mi_relDistStop = zeros(length(cells_yFor_relDistStop), 1);
    
    tuned_relDistStop_mi = zeros(length(cells_yFor_relDistStop), 1);
    tuned_yFor_mi = zeros(length(cells_yFor_relDistStop), 1);
    
    for i = 1:length(cells_yFor_relDistStop)
        results = readtable(fullfile(resultsDir, strcat('spatial_neuron_', ...
        string(cells_yFor_relDistStop(i)), '_fit_k-fold.csv')));

        mi_yFor(i, 1) = table2array(results(2, 11));
        mi_relDistStop(i, 1) = table2array(results(4, 11));
   
    end
    
    % Normalize the MI
    normMI_relDistStop = normalize(mi_relDistStop, "range");
    normMI_yFor = normalize(mi_yFor, "range");
    
    for i = 1:length(cells_yFor_relDistStop)
        if normMI_yFor(i, 1) < normMI_relDistStop(i, 1)
            tuned_relDistStop_mi(i, 1) = 1;
        else
            tuned_yFor_mi(i, 1) = 1;
        end 
    end
    
    resultsPGAM.MI.tuned_relDistStop_mi = cells_yFor_relDistStop(tuned_relDistStop_mi==1)';
    resultsPGAM.MI.tuned_yFor_mi = cells_yFor_relDistStop(tuned_yFor_mi==1)';

    resultsPGAM.MI.tuned_relDistStop_mi_licks = intersect(resultsPGAM.MI.tuned_relDistStop_mi, ...
        tunedLicks(kernelStrength_licks==1));
    resultsPGAM.MI.tuned_relDistStop_mi_nolicks = intersect(resultsPGAM.MI.tuned_relDistStop_mi, ...
        tunedLicks(kernelStrength_licks==1));

    relDistStop_mi_prop = length(resultsPGAM.MI.tuned_relDistStop_mi) / length(cells_yFor_relDistStop);
    yFor_mi_prop = length(resultsPGAM.MI.tuned_yFor_mi) / length(cells_yFor_relDistStop);
    relDistStop_licks_mi_prop = length(resultsPGAM.MI.tuned_relDistStop_mi_licks) / length(cells_yFor_relDistStop);

    resultsPGAM.MI.relDistStop_mi_prop = relDistStop_mi_prop;
    resultsPGAM.MI.yFor_mi_prop = yFor_mi_prop;
    resultsPGAM.MI.relDistStop_licks_mi_prop = relDistStop_licks_mi_prop;


    %% Plot the proportions
    
    % Plot A
    % Tuning to licks with positive kernel only 
    % Tuning to y and ylin in the forward direction only 

    variable_names = ["y", "ylin", "relDistStop", "licks", "y relDistStop", ...
            "licks relDistStop", "y licks", "y licks relDistStop"];

    variables = [yProp, ylinProp, relDistStopProp, positiveLicksProp, ...
        yFor_relDistStop_prop, licksP_relDistStop_prop, y_licksP_prop, ...
        yFor_licksp_relDistStop_prop];
        
    figure('Position', [100 100 1300 900]);
    subplot(1,2,1)
    hold on
    for v = 1:length(variables)
        h = bar(v, variables(v), 'EdgeColor', [0.5 0.5 0.5]);
        if v == 5
            h.FaceColor = [0.3 0.3 0.3];
        else
            h.FaceColor = [0.7 0.7 0.7];
        end
    end
    xlim([0, length(variables)+1])
    ylim([0,1])
    set(gca,'XTick', 1:length(variables), 'XTickLabel', variable_names);
    title('Percentage of pyramidal cells tuned to each variable and their combinations.')
    
    % Plot B 
    % Cells tuned to BOTH y and relDistStop 

    % Tuning to licks with positive kernel only 
    % Tuning to y in the forward direction only and MI > relDistStop
    % Tuning to relDistStop and MI > y
    variable_names = ["y > relDistStop", "relDistStop > y", "relDistStop > y & licks"];
    variables = [yFor_mi_prop, relDistStop_mi_prop, relDistStop_licks_mi_prop];
        
    subplot(1,2,2)
    hold on
    for v = 1:length(variables)
        bar(v, variables(v), 'EdgeColor', [0.5 0.5 0.5], 'FaceColor', [0.3 0.3 0.3])
    end
    xlim([0, length(variables)+1])
    ylim([0,1])
    set(gca,'XTick', 1:length(variables), 'XTickLabel', variable_names);
    title('Proportion of cells "preferentially" tuned to both y and relDistStop.')
    
    saveas(gcf, fullfile(postprocessDir, 'percVariableTuning'), 'fig');
    saveas(gcf, fullfile(postprocessDir, 'percVariableTuning'), 'png');

    close(figure(1));

    % All neurons
%     variables_all = [yProp_all, ylinProp_all, relDistStopProp_all, licksProp_all, ...
%         y_relDistStop_prop_all, licks_relDistStop_prop_all, y_licks_prop_all, ...
%         y_licks_relDistStop_prop_all];
%     
%     figure;
%     hold on
%     for v = 1:length(variables_all)
%         bar(v, variables_all(v), 'EdgeColor', [0.4 0.4 0.4], 'FaceColor', [0.7 0.7 0.7])
%     end
%     xlim([0, length(variables_all)+1])
%     ylim([0,1])
%     set(gca,'XTick', 1:length(variables_all), 'XTickLabel', variable_names);
%     title('Percentage of neurons tuned to each variable and their combinations.')
%     
%     saveas(gcf, fullfile(postprocessDir, 'percVariableTuning_neurons'), 'fig');
%     saveas(gcf, fullfile(postprocessDir, 'percVariableTuning_neurons'), 'png');
% 
%     close(figure(1));

    %% Save results
    if exist(fullfile(postprocessDir, strcat(session, '.postprocessPGAM.mat')))
        delete(fullfile(postprocessDir, strcat(session, '.postprocessPGAM.mat')));
    end
    save(fullfile(postprocessDir, strcat(session, '.postprocessPGAM.mat')), 'resultsPGAM'); 
   
end
end