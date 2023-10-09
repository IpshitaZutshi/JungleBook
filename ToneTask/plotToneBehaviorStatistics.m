function Summary = plotToneBehaviorStatistics

sesstoAnalyze = {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
    'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
    'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   
    'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...
    'IZ40\Final\IZ40_220708_sess17','IZ40\Final\IZ40_220714_sess18',...
    'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
    'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
    'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
    'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    
    'IZ44\Final\IZ44_220827_sess4', 'IZ44\Final\IZ44_220828_sess5',...
    'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...
    'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
    'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',... 
    'IZ47\Final\IZ47_230626_sess15','IZ47\Final\IZ47_230707_sess24',...
    'IZ47\Final\IZ47_230710_sess25','IZ47\Final\IZ47_230712_sess27',...
    'IZ48\Final\IZ48_230628_sess17','IZ48\Final\IZ48_230703_sess21',...
    'IZ48\Final\IZ48_230705_sess22','IZ48\Final\IZ48_230714_sess28',... 
    }; 

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task';

Summary.PortDiffDist = [];

for ii = 1:length(sesstoAnalyze)
    cd(strcat(expPath,'\',sesstoAnalyze{ii}))
    file = dir('*.TrialBehavior.Behavior.mat');
    load(file.name);
    
    if ~isfield(behavTrials,'probe')
        behavTrials.probe(1:length(behavTrials.correct),1) = 0;
    end
    
    if strcmp(sesstoAnalyze{ii}(1:4),'IZ39')==1
        Summary.mouseID(ii) = 1;
    elseif strcmp(sesstoAnalyze{ii}(1:4),'IZ40')==1
        Summary.mouseID(ii) = 2;
    elseif strcmp(sesstoAnalyze{ii}(1:4),'IZ43')==1
        Summary.mouseID(ii) = 3;
    elseif strcmp(sesstoAnalyze{ii}(1:4),'IZ44')==1
        Summary.mouseID(ii) = 4;
    elseif strcmp(sesstoAnalyze{ii}(1:4),'IZ47')==1
        Summary.mouseID(ii) = 5;
    elseif strcmp(sesstoAnalyze{ii}(1:4),'IZ48')==1
        Summary.mouseID(ii) = 6;
    end
    
    
    numToneTrials = sum(behavTrials.linTrial == 0 & behavTrials.stim == 0 & behavTrials.probe == 0);
    
    for kk = 1:6
        % Collect distribution of tone gains
        Summary.trialType(ii,kk)  = sum(behavTrials.toneGain == (kk-1) & behavTrials.linTrial == 0 & behavTrials.stim == 0 & behavTrials.probe == 0)./numToneTrials; 
        
        % Collect distribution of lick locations
        Summary.lickChoice(ii,kk)  = sum(behavTrials.lickLoc == (kk-1) & behavTrials.linTrial == 0 & behavTrials.stim == 0 & behavTrials.probe == 0)./numToneTrials; 
        
        % Collect distribution of correct by port (hit rate)
        Summary.portCorrect(ii,kk) = sum(behavTrials.toneGain == (kk-1) & behavTrials.correct == 1 & behavTrials.linTrial == 0 & behavTrials.stim == 0 & behavTrials.probe == 0)./...
            sum(behavTrials.toneGain == (kk-1) & behavTrials.linTrial == 0 & behavTrials.stim == 0 & behavTrials.probe == 0);      
        
        % Collect false alarm rate
        Summary.falseAlarm(ii,kk) = sum(behavTrials.correct == 0 & behavTrials.lickLoc == (kk-1) & behavTrials.linTrial == 0 & behavTrials.stim == 0 & behavTrials.probe == 0)./...
            sum(behavTrials.toneGain ~= (kk-1) & behavTrials.linTrial == 0 & behavTrials.stim == 0 & behavTrials.probe == 0);   
    end
    
    Summary.PortDiffDist = [Summary.PortDiffDist; (behavTrials.toneGain(behavTrials.linTrial == 0 & behavTrials.stim == 0 & behavTrials.probe == 0)-behavTrials.lickLoc(behavTrials.linTrial == 0 & behavTrials.stim == 0 & behavTrials.probe == 0))];
    
    % Collect performance 
    Summary.performance(ii) = sum(behavTrials.correct==1 & behavTrials.linTrial == 0 & behavTrials.stim == 0 & behavTrials.probe == 0)./numToneTrials;
    
    
end

% colMap = cbrewer('seq','Blues',18);
% col = [colMap(5,:);colMap(8,:);colMap(10,:);colMap(13,:);colMap(16,:);colMap(18,:)];
% 
% figure
% set(gcf,'Color','w')
% set(gcf,'Renderer','painters')
% set(gcf,'Position', [1 41 1920 970])
% 
% subplot(2,6,1)
% stats.trialType = groupStats([{trialType(:,1)},{trialType(:,2)},{trialType(:,3)},{trialType(:,4)},{trialType(:,5)},{trialType(:,6)}],[],'repeatedMeasures',true,'inAxis',true,'plotType','boxplot','color',col);
% title('Trial probability')
% for ii = 1:6
%     hold on
%     scatter((ones*ii)-0.4, trialType(:,ii),[],'MarkerFaceColor',col(ii,:),'MarkerEdgeColor',col(ii,:))    
% end
% 
% subplot(2,6,2)
% stats.licks = groupStats([{lickChoice(:,1)},{lickChoice(:,2)},{lickChoice(:,3)},{lickChoice(:,4)},{lickChoice(:,5)},{lickChoice(:,6)}],[],'repeatedMeasures',true,'inAxis',true,'plotType','boxplot','color',col);
% title('Lick probability')
% for ii = 1:6
%     hold on
%     scatter((ones*ii)-0.4, lickChoice(:,ii),[],'MarkerFaceColor',col(ii,:),'MarkerEdgeColor',col(ii,:))    
% end
% 
% subplot(2,6,3)
% stats.performance = groupStats([{portCorrect(:,1)},{portCorrect(:,2)},{portCorrect(:,3)},{portCorrect(:,4)},{portCorrect(:,5)},{portCorrect(:,6)}],[],'repeatedMeasures',true,'inAxis',true,'plotType','boxplot','color',col);
% for ii = 1:6
%     hold on
%     scatter((ones*ii)-0.4, portCorrect(:,ii),[],'MarkerFaceColor',col(ii,:),'MarkerEdgeColor',col(ii,:))    
% end
% ylim([0.4 1.2])
% title('Performance by trial type')
% 
% subplot(2,6,4)
% histogram(-PortDiffDist,'Normalization','probability')
% title('Distance from target (incorrect)')
% ylabel('Probability')
% xlim([-5 5])
% 
% subplot(2,6,5)
% for ii = 1:6
%     data = performance(mouseID==ii);
%     scatter(ones*ii, data,[],'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.2)
%     hold on
%     line([ii-0.3 ii+0.3],[mean(data) mean(data)], 'Color','b','LineWidth',4)
% end
% xlim([0 7])
% ylim([0.4 1])
% 
% for ii = 1:4
%     subplot(2,6,6+ii)
%     data = portCorrect(mouseID==ii,:);
%     groupStats([{data(:,1) data(:,2) data(:,3) data(:,4) data(:,5) data(:,6)}],[],'repeatedMeasures',true,'inAxis',true,'plotType','BoxLinesSEM','color',col)
%     ylim([0.4 1.2])
% end
% 
% saveas(gcf,strcat(expPath,'\Compiled\behaviorSummary.png'));
% saveas(gcf,strcat(expPath,'\Compiled\behaviorSummary.eps'),'epsc');
% saveas(gcf,strcat(expPath,'\Compiled\behaviorSummary.fig'));
% save(strcat(expPath,'\Compiled\behaviorSummary.mat'),'stats'); 

end

