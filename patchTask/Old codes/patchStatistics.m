
patch_data = plotPatchBehavior(true);


% LICK PROBABILITIES VS REWARD PROBABILITIES - PORTS
% Initialize an empty structure array
results_struct = struct('Port', {}, 'Patch', {}, 'PortRewardProbability', {}, 'PortLickProbability', {});

for port = unique(patch_data.licked_ports)

    %create patch mask
    patch_type_mask = zeros(length(patch_data.patch_type), 1);
    count = 1;
    patch_type_mask(1) = count;
    
    for i = 2:length(patch_data.patch_type)
        if patch_data.patch_type(i) ~= patch_data.patch_type(i - 1);
            count = count + 1;
        end
        patch_type_mask(i) = count;
    end
    
    %calcualte probability to lick in port X while port probability is Y
    port_probabilities = patch_data.reward_probabilities(port, :); %port's probability throughout session
    patch_type_mask = patch_type_mask';
    for patch = unique(patch_type_mask);
        licked_ports_count = nnz((patch_type_mask == patch) & (patch_data.licked_ports == port));
        total_patch_counts = nnz((patch_type_mask == patch));
        port_lick_probability = licked_ports_count / total_patch_counts;
        port_probability_patch = port_probabilities(find(patch_type_mask == patch, 1));
        
        result_entry = struct('Port', port, 'Patch', patch, ...
                              'PortRewardProbability', port_probability_patch, ...
                              'PortLickProbability', port_lick_probability);
                          
        results_struct(end+1) = result_entry;
    end

end

% Create a multiplot figure organized in a single row
unique_ports = unique([results_struct.Port]);
num_ports = length(unique_ports);
num_cols = num_ports;
num_rows = 1; 

figure;
for idx = 1:num_ports
    port = unique_ports(idx);
    
    subplot(num_rows, num_cols, idx);
    
    logicalIndex = [results_struct.Port] == port;
    filtered_values = results_struct(logicalIndex);
    bar([filtered_values.Patch], [filtered_values.PortLickProbability], 'FaceColor', [.7 .7 .7], 'EdgeColor', [.7 .7 .7]); 

   % plot([filtered_values.Patch], [filtered_values.PortLickProbability], 'k-', 'LineWidth', 2); 
    hold on;
    plot([filtered_values.Patch], [filtered_values.PortRewardProbability], 'r-', 'LineWidth', 2);

    xlabel('Patch');
    ylabel('Probability');
    title(['Port: ', num2str(port)]);
    axis tight;      
end

for idx = 1:num_ports
    subplot(num_rows, num_cols, idx);
    ylim([0 1]); 
end
%legend({'Lick Probability', 'Reward Probability'}, 'Location', 'BestOutside', 'Orientation', 'vertical');

sgtitle('Lick and Reward Probabilities Across Ports');

%PREDOMINANT PATCH OVER TIME
lickedPatch = patch_data.licked_ports;
lickedPatch(ismember(lickedPatch, [1, 2, 3])) = -1; 
lickedPatch(lickedPatch == 4) = 0;                  
lickedPatch(ismember(lickedPatch, [5, 6, 7])) = 1; 

%perform moving average
window_size = 15;
smoothed_ports = movmean(lickedPatch, window_size);
patchChangePoint = double(diff(patch_type_mask) ~= 0);

figure;
hold on;
%plot(licked_ports, 'k', 'DisplayName', 'Original Data');   % Original data in black
plot(smoothed_ports, 'b', 'DisplayName', 'Smoothed Data'); % Smoothed data in red
change_indices = find(patchChangePoint);

for i = 1:length(change_indices)
    xline(change_indices(i), 'r--', 'DisplayName', 'Patch Change'); 
end

xlabel('Index');
ylabel('Value');
legend;
grid on;
hold off;


%PROBABILITY OF CHANGING PATCH AS A FUNCTION OF CONSECUTIVE ERRORS
decisionOutcome = patch_data.rewarded_trials;

trialsNumber = length(lickedPatch);
nonRewardedTrials = zeros(1, trialsNumber); % Number of consecutive trials without reward
patchChange = zeros(1, trialsNumber - 1); 

% Calculate accumulated trials without reward
for i = 2:trialsNumber
    if decisionOutcome(i-1) == 0
        nonRewardedTrials(i) = nonRewardedTrials(i-1) + 1;
    else
        nonRewardedTrials(i) = 0;
    end
    
    % Detect decision changes
    if lickedPatch(i) ~= lickedPatch(i-1)
        patchChange(i-1) = 1; 
    end
end

% Group trials by the number of trials without reward
maxNonRewarded = max(nonRewardedTrials);
changeProbability = zeros(1, length(maxNonRewarded) + 1);

for n = 0:maxNonRewarded
    trials_n = (nonRewardedTrials == n);
    
    % Count how many of these trials have a decision change
    changes = sum(patchChange(trials_n(2:end)));
    totalTrials = sum(trials_n(2:end));
    if totalTrials > 0
        changeProbability(n + 1) = changes / totalTrials;
    else
        changeProbability(n + 1) = NaN;
    end
end

figure;
bar(0:maxNonRewarded, changeProbability, 'FaceColor', [0.2 0.6 0.8]);
xlabel('Trials without reward');
ylabel('Change probability');
title('Change probability based on trials without reward');
grid on;



% %LICKS VS PROBABILITIES IN EACH PORT
% unique_patch = unique([results_struct.Patch]);
% num_patch = length(unique_patch);
% num_cols = num_patch;
% num_rows = 1; 
% 
% figure;
% for idx = 1:num_patch
%     patch = unique_patch(idx);
% 
%     subplot(num_rows, num_cols, idx);
% 
%     % Filter the data for the current patch
%     logicalIndex = [results_struct.Patch] == patch;
%     filtered_values = results_struct(logicalIndex);
% 
%     yyaxis left; % Left y-axis for lick probabilities
%     bar([filtered_values.Port], [filtered_values.PortLickProbability], ...
%         'FaceColor', [.7 .7 .7], 'EdgeColor', [.7 .7 .7]); 
%     if idx == 1
%         yyaxis left; 
%         ax = gca; 
%         ax.YColor = [0.7, 0.7, 0.7]; 
%         ylabel('Lick Probability');
%     end
%     yyaxis right; 
%     plot([filtered_values.Port], [filtered_values.PortRewardProbability], ...
%         'r-', 'LineWidth', 2);
%     if idx ==num_patch
%         ylabel('Reward Probability');
%     end
% 
%     xlabel('Port');
%     title(['Patch: ', num2str(patch)]);
%     axis tight;      
% end
% 
% sgtitle('Lick and Reward Probabilities Across Ports');
% 

%LICKS HISTOGRAM OVER SESSION
% leftHighProb = (patch_data.patch_type == 1);
% rightHighProb = (patch_data.patch_type == 0);
% lickedPatch = patch_data.licked_ports;
% lickedPatch(ismember(lickedPatch, [0, 1, 2])) = -1;  %right
% lickedPatch(lickedPatch == 3) = 0;                   %middle
% lickedPatch(ismember(lickedPatch, [4, 5, 6])) = 1;   %left
% 
% bin_size = 10;
% num_bins = ceil(length(lickedPatch) / bin_size);
% bin_indices = ceil((1:length(lickedPatch)) / bin_size);
% 
% % Initialize
% sum_left_licks = zeros(1, num_bins);
% sum_middle_licks = zeros(1, num_bins);
% sum_right_licks = zeros(1, num_bins);
% left_high = zeros(1, num_bins);
% right_high = zeros(1, num_bins);
% 
% for bin = 1:num_bins
%     bin_data = lickedPatch(bin_indices == bin);
%     left_data = leftHighProb(bin_indices == bin);
%     right_data = rightHighProb(bin_indices == bin);
%     sum_left_licks(bin) = sum(bin_data == 1);
%     sum_middle_licks(bin) = sum(bin_data == 0);
%     sum_right_licks(bin) = sum(bin_data == -1);
%     left_high(bin) = round(mean(left_data), 0);
%     right_high(bin) = round(mean(right_data), 0);
% end


% figure;
% 
% % Subplot 1: Stacked bar plot
% subplot(2, 1, 1);
% hold on;
% bar_data = [left_high(:), right_high(:)];
% b = bar(bar_data, 'stacked', 'BarWidth', 1);
% b(1).FaceColor = [1, 0.4, 0]; b(2).FaceColor = [0, 0.6, 1]; 
% b(1).EdgeColor = 'none'; 
% b(2).EdgeColor = 'none'; 
% axis off; 
% grid on;
% pos = get(gca, 'Position'); 
% pos(4) = pos(4) * 0.1; 
% set(gca, 'Position', pos);
% x_limits = xlim;
% 
% % Subplot 2: Line plots
% subplot(2, 1, 2);
% hold on;
% plot(sum_left_licks, 'Color', [1, 0.4, 0], 'DisplayName', 'Left', 'LineWidth', 1); 
% plot(sum_middle_licks, 'Color', [0.7, 0.7, 0.7], 'DisplayName', 'Middle', 'LineWidth', 1); 
% plot(sum_right_licks, 'Color', [0, 0.6, 1], 'DisplayName', 'Right', 'LineWidth', 1); 
% xlabel('Time');
% ylabel('Number of licks');
% legend;
% grid on;
% xlim(x_limits);
% pos = get(gca, 'Position'); 
% pos(4) = pos(4) * 1.4; 
% set(gca, 'Position', pos);
% hold off;


%PERCENTAGE OF RIGHT PATCH OVER TIME
trialPercentage = patch_data.percentage;
leftHighProb = (patch_data.patch_type == 1);
rightHighProb = (patch_data.patch_type == 0);


figure;
subplot(2, 1, 1); 
hold on;
bar_data = [leftHighProb(:), rightHighProb(:)];
b = bar(bar_data, 'stacked', 'BarWidth', 1);
b(1).FaceColor = [1, 0.4, 0]; % Left bars (orange)
b(2).FaceColor = [0, 0.6, 1]; % Right bars (blue)
b(1).EdgeColor = 'none'; 
b(2).EdgeColor = 'none';
axis off; 
grid on;

subplot(2, 1, 1);
pos = get(gca, 'Position'); 
pos(4) = pos(4) * 0.1;
set(gca, 'Position', pos);
x_limits = xlim;

subplot(2, 1, 2);
hold on;
% plot(trialPercentage, 'Color', [0.7, 0.7, 0.7], 'DisplayName', 'Licks in High Probability Patch', 'LineWidth', 1);
scatter(1:length(trialPercentage), trialPercentage, 10, 'MarkerFaceColor', [0.7, 0.7, 0.7], 'MarkerEdgeColor', 'none', 'DisplayName', 'Licks in High Probability Patch');
xlabel('Trials');
ylabel('Percentage');
legend;
grid on;
subplot(2, 1, 2);
pos = get(gca, 'Position'); 
pos(4) = pos(4) * 1.4;
set(gca, 'Position', pos);
xlim(x_limits); 

%NUMBER OF TRIALS X PATCH
patches = struct();
for patch_n = unique(patch_data.patch_number)
    currentPatch = patch_data.patch_number == patch_n;
    number_trials = sum(currentPatch);
    patchType = mean(patch_data.patch_type(currentPatch == 1));
    if patchType == 0
        meanProb = mean(patch_data.reward_probabilities(1:3, currentPatch), 'all');
    else
        meanProb = mean(patch_data.reward_probabilities(5:7, currentPatch), 'all');
    end
    patches.probs(patch_n) = meanProb;
    patches.number(patch_n) = patch_n;
    patches.trials_num(patch_n) = number_trials;
    patches.type(patch_n) = patchType;
end

% Plot histogram
figure;
hold on;
for i = 1:length(patches.number)
    if patches.type(i) == 1
        % Orange for left patch
        bar(patches.number(i), patches.trials_num(i), 'FaceColor', [1, 0.4, 0]);
    else
        % Blue for right patch
        bar(patches.number(i), patches.trials_num(i), 'FaceColor', [0, 0.6, 1]);
    end
end
xlabel('Patch Number');
ylabel('Number of Trials');
title('Histogram of Trials per Patch with Left/Right Colors');
hold off;

% % EFFECT OF AVG PROBS OF HIGH PROB PATCH
% figure;
% scatter(patches.probs, patches.trials_num, 'filled'); % Scatter plot
% hold on;
% % Linear regression
% p = polyfit(patches.probs, patches.trials_num, 1);
% yfit = polyval(p, patches.probs); 
% plot(patches.probs, yfit, '-r', 'LineWidth', 2); 
% 
% % % Correlation
% % correlation_value = corr(patches.probs(:), patches.trials_num(:)); % Get correlation
% % text(min(patches.probs), max(patches.trials_num), ...
% %     sprintf('Correlation: %.2f', correlation_value), 'FontSize', 12, 'Color', 'blue');
% 
% xlabel('Average Probabilities'); 
% ylabel('Number of Trials');
% grid on;
% legend({'Data Points', 'Regression Line'}, 'Location', 'best'); % Add legend
% hold off;


% EFFECT OF REWARD IN LOW PROB PATCH
patch_lowProb = struct();
for patch_n = unique(patch_data.patch_number)
    currentPatch = patch_data.patch_number == patch_n;
    patchType = mean(patch_data.patch_type(currentPatch == 1));
    if patchType == 1
        ports = [1, 2, 3]; 
    else
        ports = [5, 6, 7];
    end
    number_trials = sum(currentPatch);
    lowProbLicked = currentPatch & ismember(patch_data.licked_ports, ports);
    rewardedLowPatch = sum(patch_data.rewarded_trials(lowProbLicked) == 1);
    patch_lowProb.number(patch_n) = patch_n;
    patch_lowProb.trials_num(patch_n) = number_trials;
    patch_lowProb.type(patch_n) = patchType;
    patch_lowProb.rewLowProb(patch_n) = rewardedLowPatch;
end

%Plot correlation
figure;
scatter(patch_lowProb.rewLowProb, patch_lowProb.trials_num, 'filled'); % Scatter plot
hold on;
% Linear regression
p = polyfit(patch_lowProb.rewLowProb, patch_lowProb.trials_num, 1);
yfit = polyval(p, patch_lowProb.rewLowProb); 
plot(patch_lowProb.rewLowProb, yfit, '-r', 'LineWidth', 2); 

xlabel('Rewarded trials in Low probability patch'); 
ylabel('Number of Trials');
grid on;
legend({'Data Points', 'Regression Line'}, 'Location', 'best'); % Add legend
hold off;
