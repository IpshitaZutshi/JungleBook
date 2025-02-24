function [behavTrials] = getCueBehavior_newSetup(varargin)

%% Updated on 3/15/2024 to match new digital inputs.

p = inputParser;
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'plotfig',true,@islogical)
addParameter(p,'forceRun',true,@islogical)
addParameter(p,'updatedIntan',true,@islogical)

parse(p,varargin{:});
saveMat = p.Results.saveMat;
plotfig = p.Results.plotfig;
forceRun = p.Results.forceRun;
updatedIntan = p.Results.updatedIntan;

basepath = pwd;

%% Deal with inputs
if ~isempty(dir([basepath filesep '*.TrialBehavior.Events.mat'])) && ~forceRun
    disp('Trial behavior already detected! Loading file.');
    file = dir([basepath filesep '*.TrialBehavior.Events.mat']);
    load(file.name);
    return
end

%% Get digital inputs
if exist('settings.xml')
    delete 'settings.xml'
end
disp('Loading digital In...');
digitalIn = bz_getDigitalIn; % Figure out why periodLag is crashing
if isempty(digitalIn)
    toneBehav = [];
    return
end

% Start of trial coincides with cue onset - take both left and right
a = digitalIn.dur{4}>0.001;
b = digitalIn.dur{5}>0.001;
timestamps = [digitalIn.timestampsOn{4}(a); digitalIn.timestampsOn{5}(b)];
[~,idxTS] = sort(timestamps);

if updatedIntan % base solenoid should have this information
    b = digitalIn.dur{7}>0.001;
    timestampsEnd = digitalIn.timestampsOn{7}(b);
    %end of trial coincides with baseIR
    if length(timestampsEnd) == length(idxTS)
        behavTrials.timestamps(:,1) = timestamps(idxTS);
        behavTrials.timestamps(:,2) = timestampsEnd;    
    elseif length(timestampsEnd) < length(idxTS)
        a = timestamps(idxTS);
        behavTrials.timestamps(:,1) = a(1:(end-1));
        behavTrials.timestamps(:,2) = timestampsEnd;     
    else
        a = timestampsEnd;  
        behavTrials.timestamps(:,1) = timestamps(idxTS);
        behavTrials.timestamps(:,2) = a(1:(end-1));   
    end    
    solChoice = digitalIn.ints{6}(1,:);
else
    %see which solenoid trials coincide with base
    solTS = digitalIn.ints{6}(1,:);
    solBase = InIntervals(solTS,[digitalIn.intsPeriods{7}(:,1)-0.5 digitalIn.intsPeriods{7}(:,1)]); % base pulses must be within 0.5 seconds of base IR
    %end of trial coincides with baseIR
    if sum(solBase) == length(idxTS)
        behavTrials.timestamps(:,1) = timestamps(idxTS);
        behavTrials.timestamps(:,2) = solTS(solBase);    
    elseif sum(solBase) < length(idxTS)
        a = timestamps(idxTS);
        behavTrials.timestamps(:,1) = a(1:(end-1));
        behavTrials.timestamps(:,2) = solTS(solBase);    
    else
        a = solTS(solBase); 
        behavTrials.timestamps(:,1) = timestamps(idxTS);
        behavTrials.timestamps(:,2) = a(1:(end-1));   
    end
    solChoice = digitalIn.ints{6}(1,~solBase);
end


a = digitalIn.dur{4}>0.001;
b = digitalIn.dur{5}>0.001;
arm = [zeros(sum(a),1);ones(sum(b),1)];
cueDur = [digitalIn.dur{4}(a)'; digitalIn.dur{5}(b)'];

behavTrials.cue = arm(idxTS);
behavTrials.cueDur = cueDur(idxTS);
durtemp  = digitalIn.dur{11}>0.001;
behavTrials.delayDur = digitalIn.dur{11}(durtemp)'-behavTrials.cueDur;

behavTrials.cue = behavTrials.cue(1:size(behavTrials.timestamps,1));
behavTrials.cueDur = behavTrials.cueDur(1:size(behavTrials.timestamps,1));
behavTrials.delayDur = behavTrials.delayDur(1:size(behavTrials.timestamps,1));

if updatedIntan
    a = digitalIn.dur{9}>0.001;
    b = digitalIn.dur{10}>0.001;
    choiceTS = [digitalIn.timestampsOn{9}(a); digitalIn.timestampsOn{10}(b)];
    %choiceTS = choiceTS(2:end);
    choice = [zeros(sum(a),1);ones(sum(b),1)];
    [~,idxTS] = sort(choiceTS);
    behavTrials.choiceTS = choiceTS(idxTS);
    behavTrials.choice = choice(idxTS);
    behavTrials.choiceTS = behavTrials.choiceTS(1:length(behavTrials.cue));
    behavTrials.choice = behavTrials.choice(1:length(behavTrials.cue));
    behavTrials.correct = ~xor(behavTrials.choice,behavTrials.cue); % 1 if they match, 0 if not
end

for ii = 1:size(behavTrials.timestamps,1)    
    curTrialTime = behavTrials.timestamps(ii,:);
    if ~updatedIntan 
        %correct solenoid Times are in solChoice        
        % check if choice solenoid was triggered during this trial   
        idx = solChoice>curTrialTime(1) & solChoice<curTrialTime(2);
        if sum(idx)>0
           behavTrials.correct(ii) = 1; 
           behavTrials.choice(ii) = behavTrials.cue(ii); 
           behavTrials.choiceTS(ii,1) = solChoice(idx);
        else
           behavTrials.correct(ii) = 0; 
           behavTrials.choice(ii) = ~behavTrials.cue(ii); 
           behavTrials.choiceTS(ii,1) = nan;
        end
    end
    % Check if it was a stim trial, and if so, stimdur and when was
    % delivered % to add - when was it delivered
    digiChan = 14;
    if ~isempty(digitalIn.ints{digiChan})
        stimTS = digitalIn.ints{digiChan}(1,:);
        idx = stimTS>curTrialTime(1) & stimTS<curTrialTime(2);
        %idx = 0;
        if sum(idx)>0
            behavTrials.stim(ii) = 1;
            a = digitalIn.dur{digiChan}(idx);
            behavTrials.stimDur(ii) = a(1);
        else
            behavTrials.stim(ii) = 0;
            behavTrials.stimDur(ii) = nan;
        end
    else
        behavTrials.stim(ii) = 0;
        behavTrials.stimDur(ii) = nan;        
    end
end

numTrials = size(behavTrials.timestamps,1);
performance = [];
for kk = 1:10:numTrials
    if (kk+10)<numTrials
        performance = [performance sum(behavTrials.correct(kk:kk+9))/10];
    else 
        performance = [performance sum(behavTrials.correct(kk:end))/(numTrials-kk)];
    end
end
behavTrials.performance10 = performance;    
behavTrials.performance = [sum(behavTrials.correct(behavTrials.stim==0))./sum(behavTrials.stim==0) sum(behavTrials.correct(behavTrials.stim==1))./sum(behavTrials.stim==1)];

if plotfig
    %plot the data
    figure
    set(gcf,'Renderer','painters')
    set(gcf,'Color','w')    
    set(gcf,'Position',[100 72 1662 888])
    
    subplot(2,2,[1 2])
    hold on
    plot(behavTrials.choiceTS, behavTrials.choice,'color',[.7 .7 .7]);
    scatter(behavTrials.choiceTS(find(behavTrials.correct == 1)),...
        behavTrials.choice(find(behavTrials.correct== 1)),100,[.6 .9 .7],'filled');
    scatter(behavTrials.choiceTS(find(behavTrials.correct== 0)),...
        behavTrials.choice(find(behavTrials.correct == 0)),100,[.9 .6 .7],'filled');
    scatter(behavTrials.choiceTS(find(behavTrials.stim== 1)),...
        ones(1,sum(behavTrials.stim== 1))+0.1,100,[8/243 133/243 161/243],'filled');    
    ylim([-0.1 1.2])
    title('Stim: bilateral hippocampus, delay+choice')
    
    subplot(2,2,3)
    plot(behavTrials.performance10,'LineWidth',1.5)
    ylabel('Performance')
    xlabel('Trials (blocks of 10)')
    ylim([0 1])
    xlim([0.5 13.5])
    hold on
    line([0.5 10.5],[0.5 0.5],'Color',[0.5 0.5 0.5],'LineStyle','--')
    title(strcat('Average performance:',num2str(sum(behavTrials.correct)./length(behavTrials.correct))))
    
    subplot(2,2,4)    
    plot(behavTrials.performance,'LineWidth',1.5)
    ylim([0 1])
    xlim([0 3])
    hold on
    line([0.5 10.5],[0.5 0.5],'Color',[0.5 0.5 0.5],'LineStyle','--')
    title(strcat('Mean stimdur:',num2str(nanmean(behavTrials.stimDur)),'ms'))
    
%     subplot(2,2,4)
%     curDurTypes = [1 2 5 10];
%     for kk = 1:4
%         behavCue(kk) = sum(behavTrials.correct(floor(behavTrials.cueDur*10)==curDurTypes(kk)))./length(behavTrials.correct(floor(behavTrials.cueDur*10)==curDurTypes(kk)));
%         numCue(kk) = length(behavTrials.correct(floor(behavTrials.cueDur*10)==curDurTypes(kk)));
%     end
%     plot(behavCue)
%     title(num2str(numCue))
%     xlim([0 4])
        
    mkdir('Behavior');
    saveas(gcf,'Behavior\cueBehavior.png');
end

% Always skip the first trial
behavTrials.timestamps = behavTrials.timestamps(2:end,:);
behavTrials.cue = behavTrials.cue(2:end);
behavTrials.cueDur = behavTrials.cueDur(2:end);
behavTrials.delayDur = behavTrials.delayDur(2:end);
behavTrials.choiceTS = behavTrials.choiceTS(2:end);
behavTrials.choice = behavTrials.choice(2:end);
behavTrials.correct = behavTrials.correct(2:end);
behavTrials.stim = behavTrials.stim(2:end);
behavTrials.stimDur = behavTrials.stimDur(2:end);

if saveMat
    C = strsplit(pwd,'\');
    save([basepath filesep C{end} '.TrialBehavior.Events.mat'],'behavTrials');
end

end