function plotPatchBehaviorSummary

col1 = [238/255 67/255 69/255;...
    241/255 114/255 42/255;...
    247/255 149/255 33/255;...
    0.5 0.5 0.5;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];

basepath = pwd;
behaviorFile = dir(fullfile(basepath, '*.TrialBehavior.mat'));
trackFile = dir(fullfile(basepath, '*Tracking.Behavior.mat'));
load(behaviorFile.name);
load(trackFile.name);

%% First plot one example session
figure
set(gcf,'Color','w')
set(gcf,'Renderer','painters')
subplot(3,4,1:4)
plot(tracking.timestamps, tracking.position.y,'Color',[0.3, 0.3 0.3],'LineWidth',1)
hold on

% for kk = 1:7
%     idx = find(behavTrials.port==kk);
%     for ii= 1:length(idx)
%         lickTS = behavTrials.timestamps(idx(ii));
%         [~,ixTS] = min(abs(tracking.timestamps-lickTS));    
%         if tracking.position.vy(ixTS)<0
%             scatter(tracking.timestamps(ixTS),tracking.position.y(ixTS),40,col1(kk,:),'filled')        
%         else
%             scatter(tracking.timestamps(ixTS),tracking.position.y(ixTS),40,col1(kk,:))        
%         end
%     end
% end
xlabel('Time')
ylabel('Position')
ylim([0 125])
xlim([tracking.timestamps(1) tracking.timestamps(end)])
box off

%% Plot licks - rewarded/unrewarded based on changing patch probability
subplot(3,4,5:8)
box off
hold on
col = [0/255 151/255 150/255; 128/255 0/255 128/255];
for ii= 1:length(behavTrials.timestamps)
   % current high patch
   patch1 = mean(behavTrials.ports_probability(ii,1:3));
   patch2 = mean(behavTrials.ports_probability(ii,5:7));
   if patch1>patch2
       highPatch = 1;
   else
       highPatch = 2;
   end
   if behavTrials.reward_outcome(ii) == 1
        scatter(behavTrials.timestamps(ii),behavTrials.port(ii),45,col(highPatch,:),"filled") 
   else
        scatter(behavTrials.timestamps(ii),behavTrials.port(ii),45,col(highPatch,:)) 
   end
end
%xlim([0 524])
%xlim([tracking.timestamps(1) tracking.timestamps(end)])
xlabel('Time')
ylabel('Port')

%% Running percentage of licks in a patch
% Define the number of trials for the running window
windowSize = 30;

% Identify high-probability patch at each trial
for ii = 1:behavTrials.num_trials
    highProbPatch(ii) = mean(behavTrials.ports_probability(ii,1:3)) > mean(behavTrials.ports_probability(ii,5:7));
end

% Determine whether each lick was in the high-probability patch
licksInHighProbPatch = ismember(behavTrials.port, 1:3) & highProbPatch' | ...
                        ismember(behavTrials.port, 5:7) & ~highProbPatch';

% Compute running percentage over a sliding window
runningPercentage = zeros(size(licksInHighProbPatch));

for t = 1:length(licksInHighProbPatch)
    if t >= windowSize
        windowLicks = licksInHighProbPatch(t-windowSize+1:t);
    else
        windowLicks = licksInHighProbPatch(1:t);
    end
    runningPercentage(t) = mean(windowLicks) * 100;
end

% Plot the running percentage
subplot(3,4,9:12)
plot(behavTrials.timestamps, runningPercentage, 'b', 'LineWidth', 1);
hold on
%line([tracking.timestamps(1) tracking.timestamps(end)],[79 79],'Color','k','LineWidth',1)
plot(behavTrials.timestamps,highProbPatch*100,'Color','r','LineWidth',1)
%xlim([tracking.timestamps(1) tracking.timestamps(end)])
xlabel('Trial Number');
ylabel('Running % of Licks in High-Probability Patch');


%% Across sessions: Distribution of average trials in a patch
expPath = 'Z:\Buzsakilabspace\LabShare\ZutshiI\patchTask\';
sess = {'N7\N7_241205_sess17','N7\N7_241206_sess18','N7\N7_241210_sess19','N7\N7_241211_sess20',...
    'N7\N7_241212_sess21','N7\N7_241213_sess22','N7\N7_241216_sess23','N7\N7_241217_sess24',...
    'N9\N9_241206_sess13','N9\N9_241210_sess12','N9\N9_241211_sess14','N9\N9_241212_sess15',...
    'N9\N9_241213_sess16'};

patch1_N7 = []; 
patch2_N7 = [];
patch1_N9 = []; 
patch2_N9 = [];
patch_N7 = []; % Total trials in a patch (combined Patch 1 & 2)
patch_N9 = [];
switches_N7 = []; % Store patch switches per session for N7
switches_N9 = []; % Store patch switches per session for N9

for ss = 1:length(sess)

    % Load behavior file
    cd(strcat(expPath,sess{ss}));
    behaviorFile = dir(fullfile(pwd, '*.TrialBehavior.mat'));
    load(behaviorFile.name);

    % Identify high-probability patch switches
    highProbPatch = mean(behavTrials.ports_probability(:,1:3), 2) < mean(behavTrials.ports_probability(:,5:7), 2);
    
    % Find consecutive trials in the same patch
    patch_switches = [1; find(diff(highProbPatch) ~= 0); behavTrials.num_trials]; 
    patch_durations = diff(patch_switches); % Compute duration of each patch

    % Separate durations by patch type and mouse
    if contains(sess{ss}, 'N7')
        patch1_N7 = [patch1_N7; patch_durations(highProbPatch(patch_switches(1:end-1)) == 1)];
        patch2_N7 = [patch2_N7; patch_durations(highProbPatch(patch_switches(1:end-1)) == 0)];
        patch_N7 = [patch_N7; patch_durations]; % Store total trials in patch
        switches_N7 = [switches_N7; length(patch_switches)-2]; % Store number of switches per session
    else
        patch1_N9 = [patch1_N9; patch_durations(highProbPatch(patch_switches(1:end-1)) == 1)];
        patch2_N9 = [patch2_N9; patch_durations(highProbPatch(patch_switches(1:end-1)) == 0)];
        patch_N9 = [patch_N9; patch_durations]; % Store total trials in patch
        switches_N9 = [switches_N9; length(patch_switches)-2]; % Store number of switches per session
    end
end

%% **Plot Box Plots**
col = [0/255 151/255 150/255;128/255 0/255 128/255];
figure
set(gcf, 'Color', 'w', 'Renderer', 'painters')

% Trials in Patch 1 vs Patch 2 (All Sessions)
subplot(1,4,1)
hold on;
data = [patch1_N7; patch1_N9; patch2_N7; patch2_N9];
groups = [repmat(1, length([patch1_N7; patch1_N9]), 1); ...
          repmat(2, length([patch2_N7; patch2_N9]), 1)];

stats.patchDiff = groupStats(data,groups,'inAxis',true,'plotType','boxplot','color',col,'labelSummary',false);
for ii = 1:2
    hold on
    scatter((ones*ii)-0.5, data(groups == ii),10,'MarkerFaceColor',col(ii,:),'MarkerEdgeColor',col(ii,:),"XJitter","rand")    
end
xticklabels(["Patch 1", "Patch 2"])
ylabel('Trials in Patch');
hold off;

% Box Plot: Trials in Patch 1 vs Patch 2 (N7 vs N9)
col = [0/255 151/255 150/255; 128/255 0/255 128/255;130/255 196/255 195/255;150/255 112/255 159/255];
subplot(1,4,2)
hold on;
data = [patch1_N7; patch2_N7; patch1_N9; patch2_N9];
groups = [repmat(1, length(patch1_N7), 1); ...
          repmat(2, length(patch2_N7), 1); ...
          repmat(3, length(patch1_N9), 1); ...
          repmat(4, length(patch2_N9), 1)];

stats.patchDiffMouse = groupStats(data,groups,'inAxis',true,'plotType','boxplot','color',col,'labelSummary',false);
for ii = 1:4
    hold on
    scatter((ones*ii)-0.5, data(groups == ii),10,'MarkerFaceColor',col(ii,:),'MarkerEdgeColor',col(ii,:),"XJitter","rand")    
end
ylabel('Trials in Patch');
hold off;


% Box Plot: Total Trials in Patch (N7 vs N9)
col = [0.8 0.8 0.8; 0.3 0.3 0.3];
subplot(1,4,3)
hold on;
data = [patch_N7; patch_N9];
groups = [repmat(1, length(patch_N7), 1); repmat(2, length(patch_N9), 1)];

stats.mouse = groupStats(data,groups,'inAxis',true,'plotType','boxplot','color',col,'labelSummary',false);
for ii = 1:2
    hold on
    scatter((ones*ii)-0.5, data(groups == ii),10,'MarkerFaceColor',col(ii,:),'MarkerEdgeColor',col(ii,:),"XJitter","rand")    
end
ylabel('Trials in Patch');
xticklabels(["N7", "N9"])
hold off;

% Box Plot: Average Patch Switches per Session (N7 vs N9)
subplot(1,4,4)
hold on;
data = [switches_N7; switches_N9];
groups = [repmat(1, length(switches_N7), 1); repmat(2, length(switches_N9), 1)];

stats.mouseSwitches = groupStats(data,groups,'inAxis',true,'plotType','boxplot','color',col,'labelSummary',false);
for ii = 1:2
    hold on
    scatter((ones*ii)-0.5, data(groups == ii),10,'MarkerFaceColor',col(ii,:),'MarkerEdgeColor',col(ii,:),"XJitter","rand")    
end
ylabel('Average switches in a session');
hold off;