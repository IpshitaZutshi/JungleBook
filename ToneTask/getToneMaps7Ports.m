function getToneMaps7Ports(varargin)

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

%Assign spike position to each spike
if ~isempty(dir([basepath filesep '*.spikeData.cellinfo.mat']))
    file = dir([basepath filesep '*.spikeData.cellinfo.mat']);
    load(file.name);
else
    for unit = 1:length(spikes.UID)
        [idx] = InIntervals(spikes.times{unit},[tracking.timestamps(1) tracking.timestamps(end)]); 
        tsBehav = spikes.times{unit}(idx);
        for tt = 1:length(tsBehav)
            [~,closestIndex] = min(abs(tracking.timestamps-tsBehav(tt)));
            spikeData.posIdx{unit}(tt) = closestIndex;
        end
        spikeData.pos{unit} = tracking.position.y(spikeData.posIdx{unit});
    end
    save([sessionInfo.FileName '.spikeData.cellinfo.mat'],'spikeData'); 
end

gain =[420/55, 420/130, 420/210, 420/290, 420/370, 420/420];
freqExp = log10(22000/1000);

for pf = 1:(size(behavTrials.timestamps,1)-1)    
    [idx] = InIntervals(tracking.timestamps,behavTrials.timestamps(pf,:));
    positions.forward{pf} = [tracking.timestamps(idx) tracking.position.y(idx)];
    positions.forward{pf} = [positions.forward{pf};[behavTrials.timestamps(pf,2) 118]];% Add a  fake 118
    if toneMap
        if behavTrials.linTrial(pf)==1
            positions.tone{pf} = [tracking.timestamps(idx) tracking.position.y(idx)*nan];
        else
            y = tracking.position.y(idx);
            tonepos = [];
            for ii = 1:length(y)
                freq = (y(ii)*gain(behavTrials.toneGain(pf)+1))/118;
                tonepos(ii) = 1000*(10.^(freqExp*freq));
            end
            tonepos(tonepos>25000) = nan;
            positions.tone{pf} = [tracking.timestamps(idx) tonepos'];
        end
    end    
    [idx] = InIntervals(tracking.timestamps,[behavTrials.timestamps(pf,2) behavTrials.timestamps(pf+1,1)]);    
    positions.reverse{pf} = [tracking.timestamps(idx) tracking.position.y(idx)];
end

firingMaps.forward = bz_firingMapAvg_IZ(positions.forward,spikes,'minTime',0.0001,'plotFig',false,'saveMat',false);
if toneMap
    firingMaps.tone = bz_firingMapAvg_IZ(positions.tone,spikes,'minTime',0.0001,'plotFig',false,'saveMat',false);
end
firingMaps.reverse = bz_firingMapAvg_IZ(positions.reverse,spikes,'minTime',0.0001,'plotFig',false,'saveMat',false);

firingMaps.linTrial = behavTrials.linTrial(1:(end-1));
firingMaps.toneTrial = behavTrials.toneTrial(1:(end-1));
firingMaps.toneGain = behavTrials.toneGain(1:(end-1));
firingMaps.correct = behavTrials.correct(1:(end-1));
firingMaps.numLicks = behavTrials.numLicks(1:(end-1),:);

save([sessionInfo.FileName '.firingMapsTrial.cellinfo.mat'],'firingMaps'); 

labels = {'forward','tone','reverse'};

clear idx
if plotfig
    
    if ~isfolder('FiringMap')
        mkdir('FiringMap')
    end
    
    %Get the index for different conditions
    idx{1} = firingMaps.linTrial ==1;
    idx{2} = firingMaps.toneGain ==0 & firingMaps.correct==1 & firingMaps.linTrial ==0;
    idx{3} = firingMaps.toneGain ==1 & firingMaps.correct==1 & firingMaps.linTrial ==0;
    idx{4} = firingMaps.toneGain ==2 & firingMaps.correct==1 & firingMaps.linTrial ==0;
    idx{5} = firingMaps.toneGain ==3 & firingMaps.correct==1 & firingMaps.linTrial ==0;
    idx{6} = firingMaps.toneGain ==4 & firingMaps.correct==1 & firingMaps.linTrial ==0;
    idx{7} = firingMaps.toneGain ==5 & firingMaps.correct==1 & firingMaps.linTrial ==0;
    idx{8} = firingMaps.toneGain ==0 & firingMaps.correct==0 & firingMaps.linTrial ==0;
    idx{9} = firingMaps.toneGain ==1 & firingMaps.correct==0 & firingMaps.linTrial ==0;
    idx{10} = firingMaps.toneGain ==2 & firingMaps.correct==0 & firingMaps.linTrial ==0;
    idx{11} = firingMaps.toneGain ==3 & firingMaps.correct==0 & firingMaps.linTrial ==0;
    idx{12} = firingMaps.toneGain ==4 & firingMaps.correct==0 & firingMaps.linTrial ==0;
    idx{13} = firingMaps.toneGain ==5 & firingMaps.correct==0 & firingMaps.linTrial ==0;    
    
    for pf = 1:length(firingMaps.forward.rateMaps)
         %  First get lin track
         figure
         set(gcf,'Renderer','painters')
         set(gcf,'Position',[2200 200 1185 712])
         for ll = 1:length(labels)
             for ii = 1:length(idx)
                datamat = [];
                numtrials = find(idx{ii}==1);
                for kk = 1:length(numtrials)
                    datamat = [datamat;firingMaps.(labels{ll}).rateMaps{pf}{numtrials(kk)}];
                end
                subplot(3,13,13*(ll-1)+ii) % linear trials
                
                if ii ==2 && ll == 1
                    datamat(:,12:end) = nan;
                elseif ii == 3 && ll == 1
                    datamat(:,21:end) = nan;
                elseif ii ==4 && ll == 1
                    datamat(:,32:end) = nan;
                elseif ii == 5 && ll == 1
                    datamat(:,43:end) = nan;
                elseif ii ==6 && ll == 1
                    datamat(:,52:end) = nan;                 
                end
                    
                if ~isempty(datamat)
                    h = imagesc(datamat);
                    %caxis([0 cmax])
                    set(h, 'AlphaData', ~isnan(datamat))
                    
                else
                    axis off
                end
             end                       
         end                  

        saveas(gcf,['FiringMap',filesep ,'cell_' num2str(pf) '.png'],'png');
        saveas(gcf,['FiringMap',filesep ,'cell_' num2str(pf) '.fig'],'fig');
        close all;
    end    
end
end