function getToneMapsAvgPop(varargin)

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
behavTrials.correct(24) = 0;
idx{1} = behavTrials.linTrial ==1;
idx{2} = behavTrials.toneGain ==0 & behavTrials.toneTrial ==1 & behavTrials.linTrial ==0;
idx{3} = behavTrials.toneGain ==1 & behavTrials.toneTrial ==1 & behavTrials.linTrial ==0;
idx{4} = behavTrials.toneGain ==2 & behavTrials.toneTrial ==1 & behavTrials.linTrial ==0;
idx{5} = behavTrials.toneGain ==0 & behavTrials.toneTrial ==0 & behavTrials.correct==1 & behavTrials.linTrial ==0;
idx{6} = behavTrials.toneGain ==1 & behavTrials.toneTrial ==0 & behavTrials.correct==1 & behavTrials.linTrial ==0;
idx{7} = behavTrials.toneGain ==2 & behavTrials.toneTrial ==0 & behavTrials.correct==1 & behavTrials.linTrial ==0;
idx{8} = behavTrials.toneGain ==0 & behavTrials.toneTrial ==0 & behavTrials.correct==0 & behavTrials.linTrial ==0;
idx{9} = behavTrials.toneGain ==1 & behavTrials.toneTrial ==0 & behavTrials.correct==0 & behavTrials.linTrial ==0;
idx{10} = behavTrials.toneGain ==2 & behavTrials.toneTrial ==0 & behavTrials.correct==0 & behavTrials.linTrial ==0;

for ii = 1:length(idx)    
    [idxPos] = InIntervals(tracking.timestamps,behavTrials.timestamps(idx{ii}(1:(end-1)),:));
    positions.forward{ii} = [tracking.timestamps(idxPos) tracking.position.y(idxPos)];
    if ~isempty(positions.forward{ii})
        positions.forward{ii} = [positions.forward{ii};[positions.forward{ii}(end,1)+0.000001 112]];% Add a  fake 112
    end
    if ii==1
        positions.tone{ii} = [tracking.timestamps(idxPos) tracking.position.y(idxPos)*nan];
    elseif ii == 2 || ii == 5 || ii == 8
        tonepos = tracking.position.y(idxPos)*(400/220);
        tonepos(tonepos>112) = nan;
        positions.tone{ii} = [tracking.timestamps(idxPos) tonepos];
    elseif ii == 3 || ii == 6 || ii == 9
        tonepos = tracking.position.y(idxPos)*(400/320);
        tonepos(tonepos>112) = nan;
        positions.tone{ii} = [tracking.timestamps(idxPos) tonepos];
    elseif ii == 4 || ii == 7 || ii == 10
        tonepos = tracking.position.y(idxPos)*(400/400);
        tonepos(tonepos>112) = nan;
        positions.tone{ii} = [tracking.timestamps(idxPos) tonepos];       
    end  
end

firingMaps.forward = bz_firingMapAvg_IZ(positions.forward,spikes,'minTime',0.1,'plotFig',false,'saveMat',false);
if toneMap
    firingMaps.tone = bz_firingMapAvg_IZ(positions.tone,spikes,'minTime',0.1,'plotFig',false,'saveMat',false);
end

firingMaps.linTrial = behavTrials.linTrial(1:(end-1));
firingMaps.toneTrial = behavTrials.toneTrial(1:(end-1));
firingMaps.toneGain = behavTrials.toneGain(1:(end-1));
firingMaps.correct = behavTrials.correct(1:(end-1));
firingMaps.numLicks = behavTrials.numLicks(1:(end-1),:);

save([sessionInfo.FileName '.firingMapsAvg.cellinfo.mat'],'firingMaps'); 

idx_imp = [1 5 6 7 8 9 10];
linTrack = linspace(0,112,55);
YlGnBu=cbrewer('seq', 'YlGnBu', 11);
if plotfig   
    figure
    set(gcf,'Renderer','painters')
    for ii = 1:length(idx_imp)
        subplot(2,length(idx_imp),ii)
        dataMat{ii} = [];
        for pf = 1:length(firingMaps.forward.rateMaps)
            dataMat{ii} = [dataMat{ii};firingMaps.forward.rateMaps{pf}{idx_imp(ii)}];
        end
        if ii == 1
            [~,idx] = max(dataMat{ii},[],2);
            [~,sortidx] = sort(idx);
        end
        
        if ~isempty(dataMat{ii})            
            if ii == 2 || ii == 5
               data = dataMat{ii}(sortidx,1:33);
            elseif ii == 3 || ii == 6
               data = dataMat{ii}(sortidx,1:47);
            else
               data = dataMat{ii}(sortidx,:);
            end
            %data(isnan(data)) = 0;
            %imagesc(data)
            a = zscore(data,[],2);
            if ii == 2 || ii == 5
               a(:,34:size(dataMat{1},2)) = nan;
            elseif ii == 3 || ii == 6
               a(:,48:size(dataMat{1},2)) = nan;
            end                       
            imagesc(linTrack, 1:1:size(a), a)
            colormap(YlGnBu)
            colorbar
            caxis([-1.5 3])
            ylabel('Cell ID')
            xlabel('Position on track (cm)') 
            
            subplot(2,length(idx_imp),ii+length(idx_imp))
            plot(nanmean(data))  
            xlim([0 size(dataMat{1},2)])
            ylim([0 5.5])
            
            dataAll{ii} = a;
        end
    end
end
figure
set(gcf,'Color','w')
hold on
a  = corr(dataAll{2},dataAll{4},'Rows','pairwise','Type','Spearman');
plot(linTrack,diag(a),'color',[56/243 61/243 150/243],'LineWidth',2)
a  = corr(dataAll{3},dataAll{4},'Rows','pairwise','Type','Spearman');
plot(linTrack,diag(a),'color',[8/243 133/243 161/243],'LineWidth',2)
a  = corr(dataAll{1},dataAll{4},'Rows','pairwise','Type','Spearman');
plot(linTrack,diag(a),'color',[85/243 85/243 85/243],'LineWidth',2)
xlabel('Position on track (cm)');
ylabel('PV correlation');
line([62 62],[-0.1 0.9],'Color','r','LineWidth',1.5)
line([87 87],[-0.1 0.9],'Color','r','LineWidth',1.5)
line([110 110],[-0.1 0.9],'Color','r','LineWidth',1.5)
end