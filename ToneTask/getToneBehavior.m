function [behavTrials] = getToneBehavior(varargin)

p = inputParser;
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'plotfig',true,@islogical)
addParameter(p,'forceRun',true,@islogical)

parse(p,varargin{:});
saveMat = p.Results.saveMat;
plotfig = p.Results.plotfig;
forceRun = p.Results.forceRun;

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
digitalIn = bz_getDigitalIn;
if isempty(digitalIn)
    toneBehav = [];
    return
end


%% Read text file
disp('Reading behavior text file')
fileinfo = dir('ToneBehav*');
fid = fopen(fileinfo.name, 'r');
behavInfo = textscan(fid,'%s','delimiter','\n');
behavTS = cellfun(@(x)strsplit(x,'_'),behavInfo{1},'UniformOutput',false);
behavTS = vertcat(behavTS{:});
behavTS = str2double(behavTS);
behavTS(:,1) = behavTS(:,1)-behavTS(1,1);

% Start of each trial coincides with the solenoid valve 1, i.e., digital input 2
% End of trial coincides with digital input 10
if length(digitalIn.timestampsOn{2})==length(digitalIn.timestampsOn{10})
    behavTrials.timestamps = [digitalIn.timestampsOn{2}(1:(end-1)) digitalIn.timestampsOn{10}(2:end)];
elseif length(digitalIn.timestampsOn{2})< length(digitalIn.timestampsOn{10})
    behavTrials.timestamps = [digitalIn.timestampsOn{2}(1:end) digitalIn.timestampsOn{10}(2:end)];
end
% 
% %Assign gain to every trial
% posGain = posInfoDS(2,:)./posInfoDS(1,:);
% posGain(posGain==Inf) = nan;

% Initialize
behavTrials.linTrial = zeros(size(behavTrials.timestamps,1),1);
behavTrials.toneTrial = zeros(size(behavTrials.timestamps,1),1);
behavTrials.toneGain = zeros(size(behavTrials.timestamps,1),1);
behavTrials.numLicks = zeros(size(behavTrials.timestamps,1),4);
behavTrials.correct = zeros(size(behavTrials.timestamps,1),1);
behavTrials.lickLoc = zeros(size(behavTrials.timestamps,1),1);

[~,idx] = unique(behavTS(:,2));
numTrials = size(behavTrials.timestamps,1);
%there is some lag between trial update an gain update, add a buffer of 10
%idxes
behavTrials.toneGain = behavTS(idx(2:(numTrials+1))+10,3);
if size(behavTS,2)>5
    behavTrials.linTrial =behavTS(idx(2:(numTrials+1))+10,5);
    behavTrials.toneTrial =behavTS(idx(2:(numTrials+1))+10,6);
end

for ii = 1:size(behavTrials.timestamps,1)
    
    trialWin = behavTrials.timestamps(ii,:)*1250;
    %if the solenoid was triggered, it was a correct trial
    if sum(InIntervals(digitalIn.timestampsOn{3}*1250,[trialWin(1) trialWin(2)+1250]))|| sum(InIntervals(digitalIn.timestampsOn{4}*1250,[trialWin(1) trialWin(2)+1250])) ||...
            sum(InIntervals(digitalIn.timestampsOn{5}*1250,[trialWin(1) trialWin(2)+1250]))
        if behavTrials.linTrial(ii) == 0
            behavTrials.correct(ii) = 1;
        else
            behavTrials.correct(ii) = 0;
        end
        behavTrials.lickLoc(ii) = behavTrials.toneGain(ii);
    end  
    % if incorrect, where did they lick?
    for lickSens = 6:9
        behavTrials.numLicks(ii,lickSens-5) = sum(InIntervals(digitalIn.timestampsOn{lickSens}*1250,[trialWin(1) trialWin(2)+1250]));
    end   
    if behavTrials.correct(ii) == 0
        if sum(behavTrials.numLicks(ii,2:4)) == 0 
            behavTrials.lickLoc(ii) = nan;
        else
            behavTrials.lickLoc(ii) = find(behavTrials.numLicks(ii,2:4),1)-1;
        end
    end
end



if plotfig
    %plot the data
    figure
    set(gcf,'Renderer','painters')
    subplot(2,3,[1 4])
    colMap = cbrewer('seq','BuPu',20);
    imagesc(behavTrials.numLicks(:,2:4))
    colormap(colMap)
    colorbar
    xlim([0.5 4.5])
    ylabel(strcat('Trials (',num2str(1),'-',num2str(size(behavTrials.timestamps,1)),')'))
    xticklabels({'Lick1','Lick2','Lick3',''})
    h=gca; h.XAxis.TickLength = [0 0];
    hold on
    scatter(ones(1,sum(behavTrials.correct))*3.5,find(behavTrials.correct),25,[70/243 148/243 73/243],'filled');
    scatter(ones(1,sum(behavTrials.correct==0))*3.75,find(behavTrials.correct==0),25,[175/243 54/243 60/243],'filled');
    scatter(ones(1,sum(behavTrials.linTrial==1))*3.9,find(behavTrials.linTrial==1),25,[50/243 50/243 50/243],'filled');
    scatter(ones(1,sum(behavTrials.toneTrial>0))*3.9,find(behavTrials.toneTrial>0),25,[150/243 150/243 150/243],'filled');
    scatter(ones(1,sum(behavTrials.toneGain==0))*4.1,find(behavTrials.toneGain==0),25,[103/243 189/243 170/243],'filled');
    scatter(ones(1,sum(behavTrials.toneGain==1))*4.1,find(behavTrials.toneGain==1),25,[8/243 133/243 161/243],'filled');
    scatter(ones(1,sum(behavTrials.toneGain==2))*4.1,find(behavTrials.toneGain==2),25,[56/243 60/243 150/243],'filled');

    % For incorrect trials, where does the lick happen
    subplot(2,3,2)
    idxIncorrect = find(behavTrials.correct==0 & behavTrials.toneTrial==0 & behavTrials.linTrial==0);
    numLicks_inc = behavTrials.numLicks(idxIncorrect,2:4);
    distLick = [sum(numLicks_inc(:,1)>0)/length(idxIncorrect) sum(numLicks_inc(:,2)>0)/length(idxIncorrect) sum(numLicks_inc(:,3)>0)/length(idxIncorrect)];
    bar(distLick,'FaceColor',[175/243 54/243 60/243])
    title('Distribution of licks for incorrect trials')
    xticklabels({'Lick1','Lick2','Lick3'})
    ylim([0 1])

    subplot(2,3,3)
    idx = find(behavTrials.correct==1);
    numLicks_inc = behavTrials.numLicks(idx,2:4);
    distLick = [sum(numLicks_inc(:,1)>0)/length(idx) sum(numLicks_inc(:,2)>0)/length(idx) sum(numLicks_inc(:,3)>0)/length(idx)];
    bar(distLick,'FaceColor',[70/243 148/243 73/243])
    title('Distribution of licks for correct trials')
    xticklabels({'Lick1','Lick2','Lick3'})
    ylim([0 1])

    subplot(2,3,5)
    percentCorrect = [sum(behavTrials.correct(behavTrials.toneGain==0 & behavTrials.toneTrial==0 & behavTrials.linTrial==0))./sum((behavTrials.toneGain==0 & behavTrials.toneTrial==0 & behavTrials.linTrial==0)) ...
        sum(behavTrials.correct(behavTrials.toneGain==1 & behavTrials.toneTrial==0 & behavTrials.linTrial==0))./sum((behavTrials.toneGain==1 & behavTrials.toneTrial==0 & behavTrials.linTrial==0)) ...
        sum(behavTrials.correct(behavTrials.toneGain==2 & behavTrials.toneTrial==0 & behavTrials.linTrial==0))./sum((behavTrials.toneGain==2 & behavTrials.toneTrial==0 & behavTrials.linTrial==0))];
    bar(percentCorrect,'FaceColor',[70/243 148/243 73/243])
    ylim([0 1])
    title(strcat('% Correct:',num2str(sum(behavTrials.correct==1 & behavTrials.toneTrial==0)/sum(behavTrials.toneTrial==0)),',',num2str(sum(behavTrials.toneGain==0 & behavTrials.toneTrial==0)),',',...
        num2str(sum(behavTrials.toneGain==1 & behavTrials.toneTrial==0)),',',num2str(sum(behavTrials.toneGain==2 & behavTrials.toneTrial==0))))

    subplot(2,3,6)
    performance = [];
    for kk = 1:10:numTrials
        if (kk+10)<numTrials
            performance = [performance sum(behavTrials.correct(kk:kk+10))/10];
        else 
            performance = [performance sum(behavTrials.correct(kk:end))/(numTrials-kk)];
        end
    end
    plot(performance,'LineWidth',1.5)
    ylabel('Performance')
    xlabel('Trials (blocks of 10)')
    ylim([0 1])
    xlim([0.5 10.5])
    hold on
    line([0.5 10.5],[0.33 0.33],'Color',[0.5 0.5 0.5],'LineStyle','--')
        
    mkdir('Behavior');
    saveas(gcf,'Behavior\toneBehavior.png');

end

if saveMat
    C = strsplit(pwd,'\');
    save([basepath filesep C{end} '.TrialBehavior.Events.mat'],'behavTrials');
end

end