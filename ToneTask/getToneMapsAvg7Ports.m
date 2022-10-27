function getToneMapsAvg7Ports(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'plotfig',true,@islogical);
addParameter(p,'toneMap',true,@islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
toneMap = p.Results.toneMap;
plotfig = p.Results.plotfig;


%% Deal with inputse
if ~isempty(dir([basepath filesep '*.Tracking.Behavior.mat'])) 
    disp('Loading tracking');
    file = dir([basepath filesep '*.Tracking.Behavior.mat']);
    load(file(1).name);
end

if ~isempty(dir([basepath filesep '*TrialBehavior.Behavior.mat'])) 
    disp('Behavior already detected! Loading file.');
    file = dir([basepath filesep '*TrialBehavior.Behavior.mat']);
    load(file(1).name);
end

if ~isempty(dir([basepath filesep '*.spikes.cellinfo.mat']))
    disp('Spikes already detected! Loading file.');
    file = dir([basepath filesep '*.spikes.cellinfo.mat']);
    load(file.name);
end

[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);

fprintf('Computing place fields\n');

% Compute three sets of maps. Forward mapped to space, forward, mapped to
% tone, and reverse. 
%Get the index for different conditions
idx{1} = behavTrials.linTrial ==1;
idx{2} = behavTrials.toneGain ==0 & behavTrials.correct==1 & behavTrials.linTrial ==0;
idx{3} = behavTrials.toneGain ==1 & behavTrials.correct==1 & behavTrials.linTrial ==0;
idx{4} = behavTrials.toneGain ==2 & behavTrials.correct==1 & behavTrials.linTrial ==0;
idx{5} = behavTrials.toneGain ==3 & behavTrials.correct==1 & behavTrials.linTrial ==0;
idx{6} = behavTrials.toneGain ==4 & behavTrials.correct==1 & behavTrials.linTrial ==0;
idx{7} = behavTrials.toneGain ==5 & behavTrials.correct==1 & behavTrials.linTrial ==0;

%gain =[420/55, 420/130, 420/210, 420/290, 420/370, 420/420];
gain = [13, 4.2, 2.2068, 1.5205, 1.1698, 1.0000];
freqExp = log10(22000/1000);

for ii = 1:length(idx)    
    [idxPos] = InIntervals(tracking.timestamps,behavTrials.timestamps(idx{ii}(1:(end-1)),:));
    y = tracking.position.y(idxPos);
    y = y-7;    
    positions.forward{ii} = [tracking.timestamps(idxPos) y];
    if ~isempty(positions.forward{ii})
        positions.forward{ii} = [positions.forward{ii};[positions.forward{ii}(end,1)+0.000001 112]];% Add a  fake 112
    end
    if ii==1
        positions.tone{ii} = [tracking.timestamps(idxPos) y*nan];
        kk = 0;
    else
        kk = ii-1;
    end
    
    if kk >0
        y = tracking.position.y(idxPos);
        y = y-7;    
        tonepos = [];
        for jj = 1:length(y)
            freq = (y(jj)*gain(kk))/111;
            tonepos(jj) = (1000*(10.^(freqExp*freq)));%/22000)*111;
        end
        tonepos(tonepos>25000) = nan;        
        positions.tone{ii} = [tracking.timestamps(idxPos) tonepos'];    
    end  
end

firingMaps.forward = bz_firingMapAvg_IZ(positions.forward,spikes,'minTime',0.0001,'plotFig',false,'saveMat',false);
if toneMap
    firingMaps.tone = bz_firingMapAvg_IZ(positions.tone,spikes,'minTime',0.0001,'plotFig',false,'saveMat',false);
end

firingMaps.linTrial = behavTrials.linTrial(1:(end-1));
firingMaps.toneTrial = behavTrials.toneTrial(1:(end-1));
firingMaps.toneGain = behavTrials.toneGain(1:(end-1));
firingMaps.correct = behavTrials.correct(1:(end-1));
firingMaps.numLicks = behavTrials.numLicks(1:(end-1),:);

save([sessionInfo.FileName '.firingMapsAvg.cellinfo.mat'],'firingMaps'); 

labels = {'forward','tone'};
colmat{1}  = [175/243 54/243 60/243; 160/243 160/243 160/243; 85/243 85/243 85/243;52/243 52/243 52/243;...
    160/243 160/243 160/243; 85/243 85/243 85/243;52/243 52/243 52/243];
colmat{2} = [0 0 0; 103/243 189/243 170/243; 8/243 133/243 161/243;56/243 61/243 150/243;...
    103/243 189/243 170/243; 8/243 133/243 161/243;56/243 61/243 150/243];

if plotfig
    
    if ~isfolder('FiringMap')
        mkdir('FiringMap')
    end
    
    figure(1)
    set(gcf,'Renderer','painters')
    set(gcf,'Position',[1921 41 1920 970])
    figure(2)
    set(gcf,'Renderer','painters')
    set(gcf,'Position',[1921 41 1920 970])
    corrMaps = [];
    for pf = 1:length(firingMaps.forward.rateMaps)
%          figure(1)
%          subplot(7,ceil(length(firingMaps.forward.rateMaps)/7),pf);                 
%          hold on
%          figure(2)
%          subplot(7,ceil(length(firingMaps.forward.rateMaps)/7),pf);                 
%          hold on         
%          for ll = 1:length(labels)
%              figure(ll)
%              for ii = 1:7
%                 datamat = firingMaps.(labels{ll}).rateMaps{pf}{ii};                
%                 if ii ==2 && ll == 1
%                     datamat(:,12:end) = nan;
%                 elseif ii == 3 && ll == 1
%                     datamat(:,21:end) = nan;
%                 elseif ii ==4 && ll == 1
%                     datamat(:,32:end) = nan;
%                 elseif ii == 5 && ll == 1
%                     datamat(:,43:end) = nan;
%                 elseif ii ==6 && ll == 1
%                     datamat(:,52:end) = nan;                 
%                 end
%                 plot(datamat, 'color',colmat{ll}(ii,:),'LineWidth',1.5);
%              end                       
%          end               
         a = corrcoef(firingMaps.forward.rateMaps{pf}{5},firingMaps.forward.rateMaps{pf}{7},'Rows','pairwise');
         corrMaps(pf,1) = a(1,2);
         a = corrcoef(firingMaps.tone.rateMaps{pf}{5},firingMaps.tone.rateMaps{pf}{7},'Rows','pairwise');
         corrMaps(pf,2) = a(1,2);
         
         
         a = corrcoef(firingMaps.forward.rateMaps{pf}{1},firingMaps.forward.rateMaps{pf}{7},'Rows','pairwise');
         corrMaps2(pf,1) = a(1,2);
         a = corrcoef(firingMaps.tone.rateMaps{pf}{1},firingMaps.tone.rateMaps{pf}{7},'Rows','pairwise');
         corrMaps2(pf,2) = a(1,2);
    end
         figure(1)
         saveas(gcf,['FiringMap',filesep ,'SpaceAvg.png'],'png');
         saveas(gcf,['FiringMap',filesep ,'SpaceAvg.fig'],'fig');
         
         figure(2)
         saveas(gcf,['FiringMap',filesep ,'ToneAvg.png'],'png');
         saveas(gcf,['FiringMap',filesep ,'ToneAvg.fig'],'fig');                  

%         saveas(gcf,['FiringMap',filesep ,'cell_' num2str(pf) '.png'],'png');
%         saveas(gcf,['FiringMap',filesep ,'cell_' num2str(pf) '.fig'],'fig');
%         close all; 
    
    figure
    subplot(1,2,1)    
    scatter(corrMaps(:,1),corrMaps(:,2),'o','filled')
    refline(1)
    xlim([-0.6 1])
    ylim([-0.6 1])
    ylabel('Tone correlation')
    xlabel('Space correlation')
    subplot(1,2,2)

    scatter(corrMaps2(:,1),corrMaps2(:,2),'o','filled')
    refline(1)
    xlim([-0.6 1])
    ylim([-0.6 1])
    ylabel('Tone correlation')
    xlabel('Space correlation')
    saveas(gcf,['FiringMap',filesep ,'MapCorr.png'],'png');
    saveas(gcf,['FiringMap',filesep ,'MapCorr.fig'],'fig');    
end
end