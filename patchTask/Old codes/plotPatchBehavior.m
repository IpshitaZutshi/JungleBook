function patch_data = PatchPlot(makePlot)
    % makePlot = true;
    fileInfo = dir('PatchBehav*');
    
    if isempty(fileInfo)
            error('No file starting with "PatchBehav" found in the current directory.');
    elseif numel(fileInfo) > 1
            error('Multiple files starting with "PatchBehav" found. Please specify the file.');
    end
        
    filename = fileInfo.name;
    fid = fopen(filename, 'r');
    
    % Initialize variables
    trial_count = 0;
    trial_numbers = [];
    rewarded_trials = [];
    licked_ports = [];
    reward_probabilities = [];
    patch_numbers = [];
    highProb_licks_percent = [];
    timestamps_licks = [];
    previous_rewarded_count = NaN; % (NaN for first trial)
    

    while ~feof(fid)
        % Read the first line (last trials)
        trial_count = trial_count + 1; % Increment trial count
        lastTrials = fgetl(fid); 
        if startsWith(lastTrials, 'Updated lastTrials: ')
            % Remove the prefix 'Updated lastTrials: ' and keep only the numbers
            lastTrials = strrep(lastTrials, 'Updated lastTrials: ', '');
        end

        values = str2num(lastTrials); 
        total_num = numel(values);
        num_ones = sum(values == 1);
        highProb_licks_percent(trial_count) = (num_ones / total_num) * 100;
  
        % Skip the empty line
        empty_line = fgetl(fid);
        if ~isempty(strtrim(empty_line))
            % If not empty, warn and skip
            disp('Expected an empty line, skipping unexpected content...');
        end
    
        % Read the second line (trial data)
        trial_line = fgetl(fid);
       
    
        % Process second line
        parts = split(trial_line, '_');
        timestamp = parts{1};
        timestamps_licks(trial_count) = str2double(timestamp);
        trialCount = parts{2};
        rewardedTrialCount = parts{3};
        portLick = parts{4};
    
        trial_numbers(trial_count) = str2double(trialCount); % Trial number
        current_rewarded_count = str2double(rewardedTrialCount); % Current rewarded trial count
        licked_ports(trial_count) = str2double(portLick); % Port licked
    
        % Get reward probabilities
        prob_str = parts{5};
        reward_probs = str2double(split(prob_str, ','));
        reward_probabilities(trial_count, :) = reward_probs;
    
        % Determine if trial was rewarded based on previous trial
        if trial_count == 1
            rewarded_trials(trial_count) = str2double(rewardedTrialCount) == 1;
        else
            rewarded_trials(trial_count) = (current_rewarded_count == previous_rewarded_count + 1);
        end
        previous_rewarded_count = current_rewarded_count;
    
        % Read the third line (patch data)
        patch_line = fgetl(fid);

        % Process third line
        patch_parts = split(patch_line, '_');
        patch_trial(trial_count) = str2double(patch_parts{1});
        patch_type(trial_count) = str2double(patch_parts{2});
    end

    port_4_probs = zeros(length(reward_probabilities), 1);
    reward_probabilities_all = [reward_probabilities(:, 1:4-1), port_4_probs, reward_probabilities(:, 4:end)];
    reward_probabilities_all = reward_probabilities_all';

    fclose(fid);

    % create struct
    patch_data = struct;
    patch_data.trial_numbers = trial_numbers;            
    patch_data.rewarded_trials = rewarded_trials;  % whether trial was rewarded (1 = rewarded, 0 = not)
    patch_data.licked_ports = licked_ports + 1;             
    patch_data.reward_probabilities = reward_probabilities_all;  
    patch_data.patch_trials = patch_trial;    
    patch_data.patch_type = patch_type;
    patch_data.timestamps = timestamps_licks;
    patch_data.percentage = highProb_licks_percent;
    patch_diff = diff(patch_data.patch_type);
    
    %get patch numbers
    count = 1;
    trial_count = 1;
    for trial = 1:length(patch_data.patch_trials)
        if trial_count > 1 && patch_data.patch_type(trial_count) ~= patch_data.patch_type(trial_count - 1)
            count = count + 1;
        end
        patch_data.patch_number(trial_count) = count;
        trial_count = trial_count + 1;
    end



    if makePlot == 1  
        timestamps = patch_data.timestamps; 
        licked_ports = patch_data.licked_ports;
      
        
        %LICKS AND PROBABILITIES
        figure;
        colormap(flipud(gray)); 
        hold on;
        customGreen = [152, 194, 9] / 255;
        customRed = [238, 75, 43] / 255;
        num_trials = length(reward_probabilities);
        % x_limits = [min(timestamps), max(timestamps)];
        y_limits = [min(licked_ports), max(licked_ports)]; 
        timestamps_minutes = timestamps / 60000;
        x_limits = [min(timestamps_minutes), max(timestamps_minutes)];
        
        
        plot(timestamps_minutes, licked_ports, 'Color', [0.5, 0.5, 0.5]);  
        hold on;
        
        rewarded_indices = find(rewarded_trials == 1);
        not_rewarded_indices = find(rewarded_trials == 0);
        scatter(timestamps_minutes(rewarded_indices), licked_ports(rewarded_indices), 36, customGreen, 'filled');
        scatter(timestamps_minutes(not_rewarded_indices), licked_ports(not_rewarded_indices), 36, customRed, 'filled');

        start_idx = 1;
        while start_idx <= length(patch_data.patch_type)
            patch_value = patch_data.patch_type(start_idx);
            end_idx = start_idx;
            
            while end_idx < length(patch_data.patch_type) && patch_data.patch_type(end_idx + 1) == patch_value
                end_idx = end_idx + 1;
            end
            
            start_time = timestamps_minutes(start_idx);
            end_time = timestamps_minutes(end_idx);
 
            if patch_value == 1
                y_bottom = 5; 
                y_top = 7;
            else
                y_bottom = 1; 
                y_top = 3;
            end
            
            rectangle('Position', [start_time, y_bottom, end_time - start_time, y_top - y_bottom], ...
                      'FaceColor', [0.5, 0.5, 0.5], ... 
                      'EdgeColor', 'none', ...
                      'FaceAlpha', 0.2);  
            start_idx = end_idx + 1;
        end
        dummy_rewarded = scatter(NaN, NaN, 36, customGreen, 'filled');
        dummy_not_rewarded = scatter(NaN, NaN, 36, customRed, 'filled');
        dummy_patch_region = patch(NaN, NaN, [0.5, 0.5, 0.5], 'FaceAlpha', 0.2);        
        legend([dummy_rewarded, dummy_not_rewarded, dummy_patch_region], ...
        {'Rewarded', 'Not Rewarded','High Probability Patch'}, 'Location', 'best');
        xlabel('Time (min)');
        ylabel('Licked Ports');
    
      
        %LICKS OVER TIME
        licks_matrix = [licked_ports', timestamps_licks'];
        timeBin = 0.5*60*1000; %min x sec x ms
        
        max_time = max(timestamps);  
        num_bins = ceil(max_time / timeBin);
    
        bin_edges = 0:timeBin:(num_bins * timeBin);
    
        unique_ports = unique(licked_ports);
        num_ports = length(unique_ports);
        licks_per_bin = zeros(num_bins, num_ports);  % Matrix bins x ports
    
        for i = 1:num_ports
            port = unique_ports(i);
            port_indices = licked_ports == port;
            port_timestamps = timestamps(port_indices);
            counts_per_bin = histcounts(port_timestamps, bin_edges);
            licks_per_bin(:, i) = counts_per_bin';
        end
    
        licksBin = licks_per_bin';
        licksBinFlipped = flipud(licksBin); % Flip rows of the data matrix
        
        % Create the heatmap
        figure;
        h = heatmap(licksBinFlipped, 'Colormap', parula);
        h.CellLabelColor = 'none';
        xlabel('Time bins');
        originalLabels = 1:size(licksBin, 1);
        flippedLabels = flip(originalLabels);
        h.YDisplayLabels = string(flippedLabels);
        ylabel('Ports');
           
        
    end
        
end